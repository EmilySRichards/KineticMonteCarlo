# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Julia (6 threads) 1.8.2
#     language: julia
#     name: julia-_6-threads_-1.8
# ---

# ### Set up the geometry

@everywhere function BathSetup(Δ, Nx, Sx, W)
    Bh = [[] for d in 1:Δ.maxdim]
    Bc = [[] for d in 1:Δ.maxdim]    
    # BATH SETUP
    
    # decide which vertices are in which region
    for (i, Δi) in enumerate(Δ.cells[1])
        # cold region
        if Δi.x[1] < Sx*W
            push!(Bc[1], i)
        end

        # hot region
        if Δi.x[1] >= Sx*(Nx-W)
            push!(Bh[1], i)
        end
    end
    
    # iteratively define the regions for higher-dim cells by whether they're in the coboundary
    for d in 1:Δ.maxdim-1
        for (σ, Δσ) in enumerate(Δ.cells[d+1])
            if any(in(Bh[d]), Δσ.∂)
                push!(Bh[d+1], σ)
            elseif any(in(Bc[d]), Δσ.∂)
                push!(Bc[d+1], σ)
            end
        end
    end
    
    
    
    # STRIP SETUP
    # set up strips for averaging over to find T(x)
    strips = [[[], []] for n in 1:Nx]
    # strip "positions" are at x=(n-1)*Lx/Nx
       
    #vertices
    for (i, Δi) in enumerate(Δ.cells[1])
        x = 1 + Δi.x[1]/Sx
        n = floor(Int, x)
        push!(strips[n][1], i)
    end
    
    # edges
    for (e, Δe) in enumerate(Δ.cells[2])
        xEdge = 0.5 .* (Δ.cells[1][Δe.∂[1]].x[1] + Δ.cells[1][Δe.∂[2]].x[1])
        x = 1 + xEdge/Sx
        n = floor(Int, x)
        push!(strips[n][2], e)
    end
    
    return Bh, Bc, strips
end

# ### Canonical Thermalisation Routine

@everywhere function CanonBath(S, D, Δ, therm_runtime, bathEdges, T)
    for t in 1:therm_runtime
        for _ in eachindex(bathEdges)
            e = rand(bathEdges) # pick a random edge e
            
            ΔE, _ = EnergyChange(S, Δ, [e], 0, D)
            
            if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/T) # if favourable, do the spin flip
                S.vals[e] = -S.vals[e]
            end
        end
    end
end

# ### Demon dynamics routine 

@everywhere function DemonBath(S, D, Δ, runtime, Th, Tc, hotEdges, coldEdges)
    
    ΔEh =  zeros(runtime)
    ΔEc =  zeros(runtime)
    
    Demons = zeros(length(Δ.cells[2]), runtime+1)
    Energy = zeros(length(Δ.cells[2]), runtime+1) # no need to set initial energy b/c we just need the variance
    for e in eachindex(Δ.cells[2]) # set initial demon energy
        Demons[e,1] = D.vals[e]
    end
    
    for t in 1:runtime
        Demons[:,t+1] = Demons[:,t]
        Energy[:,t+1] = Energy[:,t]
        for _ in Δ.cells[2]
            e = rand(eachindex(Δ.cells[2]))
            ΔE, _ = EnergyChange(S, Δ, [e], 0, D)
            
            if e in hotEdges # if edge lies in the hot bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Th)
                    S.vals[e] = -S.vals[e]
                    
                    ΔEh[t] += ΔE
                    Energy[e,t+1] += ΔE
                end
                
            elseif e in coldEdges # if edge lies in the cold bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Tc)
                    S.vals[e] = -S.vals[e]
                    
                    ΔEc[t] += ΔE
                    Energy[e,t+1] += ΔE
                end
                
            else # otherwise...
                if D.vals[e] >= ΔE
                    S.vals[e] = -S.vals[e]
                    D.vals[e] -= ΔE
                    
                    Demons[e,t+1] -= ΔE
                    Energy[e,t+1] += ΔE
                end
            end
        end
    end
    
    return ΔEh, ΔEc, Demons[:,2:end], Energy[:,2:end] # cut off start point for consistency
end



# ### Single Simulation Run

@everywhere function BathSingle(S, D, Δ, Length, Area, Tc, Th, hotEdges, coldEdges, strips, therm_runtime, runtime, t_therm, t_autocorr, N_blocks)
    
    # -- -1. Define Observables --
    Δx = Length/length(strips)
    
    Dfun = (T) -> δE/(exp(δE/T)-1)
    Tfun = (D) -> δE/log(1.0 + δE/mean(D))
    
    Jfun = (ΔEc, ΔEh) -> mean((ΔEc-ΔEh)/2/Area) # dividing by Area relies on baths having same number of edges!
    κfun = (ΔEc, ΔEh, Dl, Dr) -> -2*Δx * Jfun(ΔEc, ΔEh) / (Tfun(Dr) - Tfun(Dl))
    
    CDfun = (N, D) -> (δE/Tfun(D))^2 * exp(δE/Tfun(D))/(exp(δE/Tfun(D))-1)^2
    C0fun = (N, D, E) -> Var(E) / Tfun(D)^2 / N
    Cfun = (N, D, E) -> 1/(1/C0fun(N,D,E) - 1/CDfun(N,D))
    
    # -- 0. Run Simulation --
    
    # thermalise hot & cold baths to right temperature    
    CanonBath(S, D, Δ, therm_runtime, hotEdges, Th)
    CanonBath(S, D, Δ, therm_runtime, coldEdges, Tc)
    
    ΔEh, ΔEc, Demons, Energy = DemonBath(S, D, Δ, runtime, Th, Tc, hotEdges, coldEdges)
    
    # cut out thermalisation time
    ΔEh = ΔEh[t_therm+1:end]
    ΔEc = ΔEc[t_therm+1:end]
    Demons = Demons[:,t_therm+1:end]
    Energy = Energy[:,t_therm+1:end]
    
    tmax = runtime - t_therm
    
    # Calculate strip energies
    avgD = zeros(Float64, (length(strips), tmax))
    totE = zeros(Float64, (length(strips), tmax))
    
    NumSpins = zeros(Float64, (length(strips)))
    for x in eachindex(strips)
        NumSpins[x] = length(strips[x][2])
        
        tot_D_x = zeros(size(Demons, 2))
        tot_E_x = zeros(size(Energy, 2))
        for e in strips[x][2]
            if e in coldEdges
                avgD[x,:] .+= Dfun(Tc)
            elseif e in hotEdges
                avgD[x,:] .+= Dfun(Th)
            else
                avgD[x,:] += Demons[e,:]
            end
            
            totE[x,:] += Energy[e,:] 
        end
    end
    avgD ./= NumSpins
    
    # Functions
    T_μ = zeros(length(strips))
    T_σ = zeros(length(strips))
    C_μ = zeros(length(strips))
    C_σ = zeros(length(strips))
    κ_μ = zeros(length(strips))
    κ_σ = zeros(length(strips))
    D_μ = zeros(length(strips))
    D_σ = zeros(length(strips))
    
    for x in 2:length(strips)-1
        Cfunx = (D, E) -> Cfun(NumSpins[x], D, E)
        Difffunx = (ΔEc, ΔEh, Dl, Dr, D, E) -> κfun(ΔEc, ΔEh, Dl, Dr) ./ Cfunx(D, E)
        
        T_μ[x], T_σ[x] = Estimator(Bootstrap, [avgD[x,:]], Tfun, t_autocorr, N_blocks)
        C_μ[x], C_σ[x] = Estimator(Bootstrap, [avgD[x,:], totE[x,:]], Cfunx,  t_autocorr, N_blocks)
        κ_μ[x], κ_σ[x] = Estimator(Bootstrap, [ΔEc, ΔEh, avgD[x-1,:], avgD[x+1,:]], κfun, t_autocorr, N_blocks)
        D_μ[x], D_σ[x] = Estimator(Bootstrap, [ΔEc, ΔEh, avgD[x-1,:], avgD[x+1,:], avgD[x,:], totE[x,:]], Difffunx, t_autocorr, N_blocks)
    end
    
    result = zeros(2, 4, length(strips))
    result[:,1,:] = hcat(T_μ, T_σ.^2)'
    result[:,2,:] = hcat(κ_μ, κ_σ.^2)'
    result[:,3,:] = hcat(C_μ, C_σ.^2)'
    result[:,4,:] = hcat(D_μ, D_σ.^2)'
    
    return result[:,:,2:end-1] # cut off ends where κ ill-defined 
end

# ### Overall simulation routine

function BathSimulation(L, PBC, Basis, W, Tc, Th, num_histories, therm_runtime, runtime, t_therm, t_autocorr, N_blocks)
    
    # set up graph and demarcate baths and strips
    Δ, Scale = LatticeGrid(L, PBC, Basis);
    Length = L[1]*Scale[1] # length of sample
    Area = prod(L[2:end])*prod(Scale[2:end]) # cross-sectional area of sample
    
    hotCells, coldCells, strips = BathSetup(Δ, L[1], Scale[1], W)
    
    # initialise spins & demons in ground state
    S = CreateField(Δ, 1)
    D = CreateField(Δ, 1)
    GroundState!(S, Δ)
    
    ks = range(1,2*num_histories)
    Hs = [num_histories for k=ks]
    args = [[deepcopy(S), deepcopy(D), Δ, Length, Area, Tc, Th, hotCells[2], coldCells[2], strips, therm_runtime, runtime, t_therm, t_autocorr, N_blocks] for k=ks]
    
    function hfun(k, H, args)
        n = div(k-1,H) + 1 # unif/rand index
        
        if n==2 # if random initial state
            S = (-1).^rand(Bool, length(Δ.cells[2]))
        end
        
        return BathSingle(args...)
    end
    
    if multiProcess
        results = pmap(hfun, ks, Hs, args)
    else
        results = Array{Any}(undef, length(ks))
        for k in ks
            results[k] = hfun(k, Hs[k], args[k])
        end
    end 
        
    tmp = zeros(2, 4, 2, length(strips)-2, num_histories) # estimates for T,κ,C,D
    for k in ks
        n,h = divrem(k-1,num_histories) .+ (1,1)
        
        tmp[:,:,n,:,h] = results[k]
    end
    tmp = sum(tmp, dims=5)
    
    # average over observables for all histories - okay b/c iid random variables
    tmp[2,:,:,:] = sqrt.(tmp[2,:,:,:])
    tmp ./= num_histories
        
    return tmp[1,1,:,:], tmp[1,2,:,:], tmp[1,3,:,:], tmp[1,4,:,:], tmp[2,1,:,:], tmp[2,2,:,:], tmp[2,3,:,:], tmp[2,4,:,:]
end
