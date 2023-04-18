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

@everywhere function BathSetup(vertices, edges, Lx, W)
    Bh_j = []
    Bc_j = []
    
    Bh_α = []
    Bc_α = []
    
    # BATH SETUP
    for j in eachindex(vertices)
        # cold region
        if vertices[j].x[1] <= W
            push!(Bc_j, j)
        end

        # hot region
        if vertices[j].x[1] > Lx-W
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
    strips = [[[], []] for n in 1:Lx]
    
    #vertices
    for j in eachindex(vertices)
        x = vertices[j].x[1] + 0.5 # + 0.5 needed so as not to include some bath edges in the rightmost strip
        n = floor(Int, x)
        push!(strips[n][1], j)
    end
    
    # edges
    for α in eachindex(edges)
        x = edges[α].x[1] + 0.5 # + 0.5 needed so as not to include some bath edges in the rightmost strip
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
    for α in eachindex(edges) # set initial demon energy
        D[α,1] = edges[α].D
    end
    
    for t in 1:runtime
        D[:,t+1] = D[:,t]
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)
            
            if β in Bh_α # if edge lies in the hot bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Th)
                    edges[β].σ = !edges[β].σ
                    
                    ΔEh[t] += ΔE
                end

            elseif β in Bc_α # if edge lies in the cold bath...
                if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/Tc)
                    edges[β].σ = !edges[β].σ
                    
                    ΔEc[t] += ΔE
                end

            else # otherwise...
                if edges[β].D >= ΔE
                    edges[β].σ = !edges[β].σ
                    edges[β].D -= ΔE
                    
                    D[β,t+1] -= ΔE
                end
            end
        end
    end
    
    return ΔEh, ΔEc, D[:,2:end] # cut off start point for consistency
end

# ### Single Simulation Run

@everywhere function BathSingle(vertices, edges, Area, Tc, Th, Bh_α, Bc_α, strips, 𝒽, therm_runtime, runtime, t_therm, t_autocorr, N_blocks)
    
    # -- 0. Run Simulation --
    
    # thermalise hot & cold baths to right temperature
    CanonBath(vertices, edges, therm_runtime, Bh_α, Th, 𝒽)
    CanonBath(vertices, edges, therm_runtime, Bc_α, Tc, 𝒽)

    # run simulation for whole system
    ΔEh, ΔEc, D = DemonBath(vertices, edges, runtime, Th, Tc, Bh_α, Bc_α, 𝒽)
    
    # cut out thermalisation time
    ΔEh = ΔEh[t_therm+1:end]
    ΔEc = ΔEc[t_therm+1:end]
    D = D[:,t_therm+1:end]
    
    tmax = runtime - t_therm
    
    # Calculate strip energies
    totD = zeros(Float64, (length(strips), tmax))
    NumSpins = zeros(Float64, (length(strips)))
    for x in eachindex(strips)
        NumSpins[x] = length(strips[x][2])
        
        tot_D_x = zeros(size(D, 2))
        for α in strips[x][2]
            if α in Bc_α
                totD[x,:] .+= δE/(exp(δE/Tc)-1) + 2*𝒽/(exp(2*𝒽/Tc)+1)
            elseif α in Bh_α
                totD[x,:] .+= δE/(exp(δE/Th)-1) + 2*𝒽/(exp(2*𝒽/Th)+1)
            else
                totD[x,:] += D[α,:]
            end
        end
    end
    avgD = totD ./ NumSpins
    
    
    # Functions
    Δx = 1 # FOR NOW FINE BUT WILL DEPEND HOW STRIPS ARE DEFINED IN GENERAL

    Jfun = (ΔEc, ΔEh) -> mean((ΔEc-ΔEh)/2/Area) # dividing by Area relies on baths having same number of edges!
    Dfun = (T) -> δE/(exp(δE/T)-1) + 2*𝒽/(exp(2*𝒽/T)+1)
    D2fun = (T) -> sign(T)*Dfun(abs(T)) # makes D(T) an odd function so bracketing works better
    Tfun = (D) -> 𝒽>0 ?  find_zero((T) -> D2fun(T)-mean(D), 5.0) : δE/log(1.0 + δE/mean(D))
    κfun = (ΔEc, ΔEh, Dl, Dr) -> -2*Δx * Jfun(ΔEc, ΔEh) / (Tfun(Dr) - Tfun(Dl))
    
    
    T_μ = zeros(length(strips))
    T_σ = zeros(length(strips))
    C_μ = zeros(length(strips))
    C_σ = zeros(length(strips))
    κ_μ = zeros(length(strips))
    κ_σ = zeros(length(strips))
    
    for x in 2:length(strips)-1
        CDfun = (D) -> NumSpins[x] * (δE/Tfun(D))^2 * exp(δE/Tfun(D))/(exp(δE/Tfun(D))-1)^2
        Cfun = (D, E) -> CDfun(D) * Var(E) /(CDfun(D)*Tfun(D)^2 - Var(E)) / NumSpins[x] # THIS ALSO NEEDS TO BE CHANGED FOR 
        
        T_μ[x], T_σ[x] = MyBootstrap([avgD[x,:]], Tfun, t_autocorr, N_blocks)
        C_μ[x], C_σ[x] = MyBootstrap([avgD[x,:], totD[x,:]], Cfun, t_autocorr, N_blocks)
        κ_μ[x], κ_σ[x] = MyBootstrap([ΔEc, ΔEh, avgD[x-1,:], avgD[x+1,:]], κfun, t_autocorr, N_blocks)
    end
    
    result = zeros(2, 3, length(strips))
    result[:,1,:] = hcat(T_μ, T_σ.^2)'
    result[:,2,:] = hcat(κ_μ, κ_σ.^2)'
    result[:,3,:] = hcat(C_μ, C_σ.^2)'
    
    return result[:,:,2:end-1] # cut off ends where κ ill-defined 
end

# ### Overall simulation routine

function BathSimulation(L, PBC, Basis W, Tc, Th, num_histories, therm_runtime, runtime, t_therm, t_autocorr, N_blocks, 𝒽)
    
    # set up graph and demarcate baths and strips
    vertices, edges = LatticeGrid(L, PBC, Basis);
    Lx = L[1] # length of sample
    Area = prod(L[2:end]) # cross-sectional area of sample
    
    Bh_j, Bc_j, Bh_α, Bc_α, strips = BathSetup(vertices, edges, L[1], W)
    
    # initialise spins in ground state
    if sixVertex
        for edge in edges # gives ~GS (exact for PBCs) for square lattice
            edge.σ = vertices[edge.∂[1]].x[1]-vertices[edge.∂[2]].x[1]==0
        end
    end
    
    
    ks = range(1,2*length(𝒽)*num_histories)
    Hs = [num_histories for k=ks]
    ℋs = [length(𝒽) for k=ks]
    args = [[deepcopy(vertices), deepcopy(edges), Area, Tc, Th, Bh_α, Bc_α, strips, 𝒽[rem(div(k-1,num_histories),length(𝒽))+1], therm_runtime, runtime, t_therm, t_autocorr, N_blocks] for k=ks]
    
    function hfun(k, H, ℋ, args)
        n = div(div(k-1,H), ℋ) + 1 # unif/rand index
        
        if n==2 # if random initial state
            for edge in args[2]
                edge.σ = rand(Bool)
            end
        end
        
        return BathSingle(args...)
    end
    
    if multiProcess
        results = pmap(hfun, ks, Hs, ℋs, args)
    else
        results = Array{Any}(undef, length(ks))
        for k in ks
            results[k] = hfun(k, Hs[k], ℋs[k], args[k])
        end
    end 
        
    tmp = zeros(2, 3, 2, length(strips)-2, length(𝒽), num_histories) # estimates for T,κ,C
    for k in ks
        ni,h = divrem(k-1,num_histories) .+ (1,1)
        n,i = divrem(ni-1,length(𝒽)) .+ (1,1)
        
        tmp[:,:,n,:,i,h] = results[k]
    end
    tmp = sum(tmp, dims=6)
    
    # average over observables for all histories - okay b/c iid random variables
    tmp[2,:,:,:,:] = sqrt.(tmp[2,:,:,:,:])
    tmp ./= num_histories
        
    return tmp[1,1,:,:,:], tmp[1,2,:,:,:], tmp[1,3,:,:,:], tmp[2,1,:,:,:], tmp[2,2,:,:,:], tmp[2,3,:,:,:]
end
