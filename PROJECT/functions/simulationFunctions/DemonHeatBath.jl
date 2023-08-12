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

@everywhere function BathSetup(vertices, edges, Nx, Sx, W)
    Bh_j = []
    Bc_j = []
    
    Bh_α = []
    Bc_α = []
    
    # BATH SETUP
    for j in eachindex(vertices)
        # cold region
        if vertices[j].x[1] < Sx*W
            push!(Bc_j, j)
        end

        # hot region
        if vertices[j].x[1] >= Sx*(Nx-W)
            push!(Bh_j, j)
        end
    end
    
    # rough BCs on both sides for more symmetric definition of baths
    for α in eachindex(edges)
        if any(in(Bh_j), edges[α].∂)
            push!(Bh_α, α)
        elseif any(in(Bc_j), edges[α].∂)
            push!(Bc_α, α)
        end
    end
    
    # STRIP SETUP
    # set up strips for averaging over to find T(x)
    strips = [[[], []] for n in 1:Nx]
    # strip "positions" are at x=(n-1)*Lx/Nx
       
    #vertices
    for j in eachindex(vertices)
        x = 1 + vertices[j].x[1]/Sx
        n = floor(Int, x)
        push!(strips[n][1], j)
    end
    
    # edges
    for α in eachindex(edges)
        x = 1 + edges[α].x[1]/Sx
        n = floor(Int, x)
        push!(strips[n][2], α)
    end
    
    return Bh_j, Bc_j, Bh_α, Bc_α, strips
end

# ### Canonical Thermalisation Routine

@everywhere function CanonBath(vertices, edges, therm_runtime, B_α, T, 𝒽)
    for t in 1:therm_runtime
        for n in eachindex(B_α)
            β = B_α[rand(eachindex(B_α))]
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)

            if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/T)
                edges[β].σ = !edges[β].σ
            end
        end
    end
end

# ### Demon dynamics routine 

@everywhere function DemonBath(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α, 𝒽)
    
    ΔEh =  zeros(runtime)
    ΔEc =  zeros(runtime)
    
    D = zeros(length(edges), runtime+1)
    E = zeros(length(edges), runtime+1) # no need to set initial energy b/c we just need the variance
    for α in eachindex(edges) # set initial demon energy
        D[α,1] = edges[α].D
    end
    
    for t in 1:runtime
        D[:,t+1] = D[:,t]
        E[:,t+1] = E[:,t]
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)
            
            if β in Bh_α # if edge lies in the hot bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Th)
                    edges[β].σ = !edges[β].σ
                    
                    ΔEh[t] += ΔE
                    E[β,t+1] += ΔE
                end
                
            elseif β in Bc_α # if edge lies in the cold bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Tc)
                    edges[β].σ = !edges[β].σ
                    
                    ΔEc[t] += ΔE
                    E[β,t+1] += ΔE
                end
                
            else # otherwise...
                if edges[β].D >= ΔE
                    edges[β].σ = !edges[β].σ
                    edges[β].D -= ΔE
                    
                    D[β,t+1] -= ΔE
                    E[β,t+1] += ΔE
                end
            end
        end
    end
    
    return ΔEh, ΔEc, D[:,2:end], E[:,2:end] # cut off start point for consistency
end

# ### Spin-swap demon dynamics routine 

@everywhere function DemonBath_2flip(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α, 𝒽)
    
    ΔEh =  zeros(runtime)
    ΔEc =  zeros(runtime)
    
    D = zeros(length(edges), runtime+1)
    for α in eachindex(edges) # set initial demon energy
        D[α,1] = edges[α].D
    end
    
    for t in 1:runtime
        D[:,t+1] = D[:,t]
        for _ in vertices
            # propose flips
            i = rand(eachindex(vertices)) # shared vertex
            𝜷 = sample(vertices[i].δ, 2; replace=true) # two nearest-neighbour spins to flip (in order)
            
            𝒊 = [edges[𝜷[n]].∂[findfirst(edges[𝜷[n]].∂ .!= i)] for n in 1:2] # outer vertices (but may still coincide)
            
            # calculate overall energy change and current density between the two unshared vertices
            ΔE = ΔE_2flip(vertices, edges, 𝜷, 𝒊, i, 𝒽)
            
            if β in Bh_α # if edge lies in the hot bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Th)
                    edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                    edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                    
                    ΔEh[t] += ΔE
                end

            elseif β in Bc_α # if edge lies in the cold bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Tc)
                    edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                    edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                    
                    ΔEc[t] += ΔE
                end

            else # otherwise...
                if edges[𝜷[1]].D >= ΔE/2 && edges[𝜷[2]].D >= ΔE/2 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ
                    edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                    edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                    
                    edges[𝜷[1]].D -= ΔE/2
                    edges[𝜷[2]].D -= ΔE/2
                    
                    D[𝜷[1],t+1] -= ΔE/2
                    D[𝜷[2],t+1] -= ΔE/2
                end
            end
        end
    end
    
    return ΔEh, ΔEc, D[:,2:end] # cut off start point for consistency
end

# ### Single Simulation Run

@everywhere function BathSingle(vertices, edges, Length, Area, Tc, Th, Bh_α, Bc_α, strips, therm_runtime, runtime, t_therm, t_autocorr, N_blocks, 𝒽)
    
    # -- -1. Define Observables --
    g = 2*𝒽 - δE*ceil(2*𝒽/δE)
    Δx = Length/length(strips)
    
    Dfun = (T) -> δE/(exp(δE/T)-1) - g/(exp(-g/T)+1)
    Tfun = (D) -> (𝒽==0) ? δE/log(1.0 + δE/mean(D)) : find_zero((T) -> sign(T)*Dfun(abs(T)) - mean(D), (-2*Th, 2*Th))
    
    Jfun = (ΔEc, ΔEh) -> mean((ΔEc-ΔEh)/2/Area) # dividing by Area relies on baths having same number of edges!
    κfun = (ΔEc, ΔEh, Dl, Dr) -> -2*Δx * Jfun(ΔEc, ΔEh) / (Tfun(Dr) - Tfun(Dl))
    
    CDfun = (N, D) -> ((δE/Tfun(D))^2 * exp(δE/Tfun(D))/(exp(δE/Tfun(D))-1)^2 + (g/Tfun(D))^2 * exp(g/Tfun(D))/(exp(g/Tfun(D))+1)^2)
    C0fun = (N, D, E) -> Var(E) / Tfun(D)^2 / N
    Cfun = (N, D, E) -> 1/(1/C0fun(N,D,E) - 1/CDfun(N,D))
    
    # -- 0. Run Simulation --
    
    # thermalise hot & cold baths to right temperature
    CanonBath(vertices, edges, therm_runtime, Bh_α, Th, 𝒽)
    CanonBath(vertices, edges, therm_runtime, Bc_α, Tc, 𝒽)
    
    ΔEh, ΔEc, D, E = DemonBath(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α, 𝒽)
    
    #𝒽 = twoFlip ? 1 : 0
    #
    ## thermalise hot & cold baths to right temperature
    #CanonBath(vertices, edges, therm_runtime, Bh_α, Th, 𝒽)
    #CanonBath(vertices, edges, therm_runtime, Bc_α, Tc, 𝒽)
    #
    ## run simulation for whole system
    #if twoFlip
    #    ΔEh, ΔEc, D = DemonBath_2flip(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α, 𝒽)
    #else
    #    ΔEh, ΔEc, D = DemonBath(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α)
    #end
    #
    #δE0 = twoFlip ? δE/2 : δE
    
    # cut out thermalisation time
    ΔEh = ΔEh[t_therm+1:end]
    ΔEc = ΔEc[t_therm+1:end]
    D = D[:,t_therm+1:end]
    E = E[:,t_therm+1:end]
    
    tmax = runtime - t_therm
    
    # Calculate strip energies
    avgD = zeros(Float64, (length(strips), tmax))
    totE = zeros(Float64, (length(strips), tmax))
    
    NumSpins = zeros(Float64, (length(strips)))
    for x in eachindex(strips)
        NumSpins[x] = length(strips[x][2])
        
        tot_D_x = zeros(size(D, 2))
        tot_E_x = zeros(size(E, 2))
        for α in strips[x][2]
            if α in Bc_α
                avgD[x,:] .+= Dfun(Tc)
            elseif α in Bh_α
                avgD[x,:] .+= Dfun(Th)
            else
                avgD[x,:] += D[α,:]
            end
            
            totE[x,:] += E[α,:] 
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
        Difffunx = (ΔEc, ΔEh, Dl, Dr, D0, E) -> κfun(ΔEc, ΔEh, Dl, Dr) ./ Cfunx(D, E)
        
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

function BathSimulation(L, PBC, Basis, W, Tc, Th, num_histories, therm_runtime, runtime, t_therm, t_autocorr, N_blocks, 𝒽)
    
    # set up graph and demarcate baths and strips
    cells, Scale = LatticeGrid(L, PBC, Basis);
    vertices = cells[1]
    edges = cells[2]
    Length = L[1]*Scale[1] # length of sample
    Area = prod(L[2:end])*prod(Scale[2:end]) # cross-sectional area of sample
    
    Bh_j, Bc_j, Bh_α, Bc_α, strips = BathSetup(vertices, edges, L[1], Scale[1], W)
    
    # initialise spins in ground state
    GroundState!(cells)
    
    ks = range(1,2*length(𝒽)*num_histories)
    Hs = [num_histories for k=ks]
    args = [[deepcopy(vertices), deepcopy(edges), Length, Area, Tc, Th, Bh_α, Bc_α, strips, therm_runtime, runtime, t_therm, t_autocorr, N_blocks, 𝒽[rem(div(k-1,num_histories), length(𝒽))+1]] for k=ks]
    
    function hfun(k, H, args)
        n = div(k-1,H) + 1 # unif/rand index
        
        if n==2 # if random initial state
            for edge in args[2]
                edge.σ = rand(Bool)
            end
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
        
    tmp = zeros(2, 4, 2, length(strips)-2, length(𝒽), num_histories) # estimates for T,κ,C,D
    for k in ks
        ni,h = divrem(k-1,num_histories) .+ (1,1)
        n,i = divrem(ni-1,length(𝒽)) .+ (1,1)
        
        tmp[:,:,n,:,i,h] = results[k]
    end
    tmp = sum(tmp, dims=6)
    
    # average over observables for all histories - okay b/c iid random variables
    tmp[2,:,:,:,:] = sqrt.(tmp[2,:,:,:,:])
    tmp ./= num_histories
        
    return tmp[1,1,:,:,:], tmp[1,2,:,:,:], tmp[1,3,:,:,:], tmp[1,4,:,:,:], tmp[2,1,:,:,:], tmp[2,2,:,:,:], tmp[2,3,:,:,:], tmp[2,4,:,:,:]
end
