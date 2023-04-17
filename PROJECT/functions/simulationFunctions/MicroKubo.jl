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

@everywhere function MicroKuboSetup(vertices, edges, therm_runtime, T, 𝒽, isRandom)
    
    if isRandom # initialise entire system in random state
        for edge in edges
            edge.σ = rand(Bool)
        end
    else # initialise entire system in ground state
        for edge in edges
            if sixVertex
                edge.σ = vertices[edge.∂[1]].x[1]-vertices[edge.∂[2]].x[1]==0 # gives ~GS ONLY for PBCs on square lattice
            else
                edge.σ = false
            end
            edge.D = 0
        end
    end
    
    E = zeros(therm_runtime+1) # just set initial energy to zero since we only need the variance
    
    # thermalise entire system
    for t in 1:therm_runtime
        E[t+1] = E[t]
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)

            if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/T)
                edges[β].σ = !edges[β].σ
                E[t+1] += ΔE
            end
        end
    end
    
    return E[2:end] # cut out initial energy for consistency with other observables
end

# ### Single spin-flip dynamics routine 

@everywhere function MicroKubo(vertices, edges, runtime, 𝒽)
    J = zeros(Float64, (length(vertices[1].x), runtime))
    
    for t in 1:runtime
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)
            
            if ΔE == 0
                Δj_β = Δj_flip(vertices, edges, β)
                edges[β].σ = !edges[β].σ
            
                # update x-current
                r_β = vertices[edges[β].∂[1]].x - vertices[edges[β].∂[2]].x
                for d in 1:length(r_β) # if vector has any axis displacement > 1, normalise to handle PBCs
                    r_β[d] /= (abs(r_β[d])>1) ? -abs(r_β[d]) : 1
                end
                
                J[:,t] += r_β * Δj_β # note no factor of 1/2 b/c only sum each edge once
            end
        end
    end
    
    return J
end

# ### Double spin-flip dynamics routine

@everywhere function MicroKubo_2flip(vertices, edges, runtime, 𝒽)
    J = zeros(Float64, (length(vertices[1].x), runtime))
    
    for t in 1:runtime
        for _ in 1:floor(UInt32, length(edges)/2)
            
            # propose flips
            i = rand(eachindex(vertices)) # shared vertex
            𝜷 = sample(vertices[i].δ, 2; replace=true) # two nearest-neighbour spins to flip (in order)
            
            𝒊 = [edges[𝜷[n]].∂[findfirst(edges[𝜷[n]].∂ .!= i)] for n in 1:2] # outer vertices (but may still coincide)
            
            ΣA = A(edges, vertices[i]) + A(edges, vertices[𝒊[1]]) + A(edges, vertices[𝒊[2]])
            
            # calculate overall energy change and current density between the two unshared vertices
            ΔE = ΔE_2flip(vertices, edges, 𝜷, 𝒊, 𝒽)
            Δj = Δj_2flip(vertices, edges, 𝜷, 𝒊, 𝒽)
                
            # decide whether to accept and perform the move
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA>0 # energy AND magnetisation conserved AND NO pair diffusion moves (i.e. no particle at central site i)
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA<0 # energy AND magnetisation conserved AND ONLY pair diffusion moves (i.e. no particle at central site i)
            if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ # energy AND magnetisation conserved
            #if ΔE == 0 # energy conserved
                
                edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                
                # get path of current flow
                r_β1 = vertices[i].x - vertices[𝒊[1]].x
                for d in 1:length(r_β1) # if vector has any axis displacement > 1, normalise to handle PBCs
                    r_β1[d] /= (abs(r_β1[d])>1) ? -abs(r_β1[d]) : 1
                end
                
                r_β2 = vertices[𝒊[2]].x - vertices[i].x
                for d in 1:length(r_β2) # if vector has any axis displacement > 1, normalise to handle PBCs
                    r_β2[d] /= (abs(r_β2[d])>1) ? -abs(r_β2[d]) : 1
                end
                
                J[:,t] += (r_β1 + r_β2) * Δj # note no factor of 1/2 b/c only sum each pair of sites once
            end
        end
    end
    
    return J
end

# ### Single Simulation Run

@everywhere function MKuboSingle(vertices, edges, runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽)
    
    Cfun = (E) -> var(E) / T^2 / length(edges)
    κfun = (S) -> mean(S) / T^2 / length(edges)
    Dfun = (E,S) -> κfun(S) / Cfun(E)
    
    tmax = runtime-t_therm
    
    # -- 0. Run Simulation --
    E = MicroKuboSetup(vertices, edges, therm_runtime, T, 𝒽, false)
    
    M = 0
    for edge in edges
        M += (-1)^edge.σ
    end
    M /= length(edges)
    
    bondNumber = 0
    maxBondNumber = 0
    for vertex in vertices
        z = length(vertex.δ)
        z₋ = 0
        for α in vertex.δ
            z₋ += edges[α].σ ? 1 : 0
        end
        z₊ = z - z₋
        
        bondNumber += z₊*z₋/2
        maxBondNumber += z*(z-1)/2
    end
    ℙ = bondNumber/maxBondNumber
    
    if twoFlip
        J = MicroKubo_2flip(vertices, edges, runtime, 𝒽)
    else
        J = MicroKubo(vertices, edges, runtime, 𝒽)
    end
    
    # cut out thermalisation time
    J = J[:,t_therm+1:end]
    E = E[t_therm+1:end]
    
    # -- 1. Heat Capacity --
    C_μ, C_σ = MyBootstrap([E], Cfun, t_autocorr, N_blocks)
    
    # -- 2. Thermal Conductivity and Diffusivity--s
    κ_μ = 0
    κ_v = 0
    D_μ = 0
    D_v = 0
    for τ in 0:t_cutoff
        statistic = (τ==0 ? 0.5 : 1.0) .* J[1,:] .* circshift(J[1,:], -τ)
        
        tmp1, tmp2 = MyBootstrap([statistic[1:end-τ]], κfun, t_autocorr, N_blocks)
        κ_μ += tmp1
        κ_v += tmp2^2
        
        tmp1, tmp2 = MyBootstrap([E[1:end-τ], statistic[1:end-τ]], Dfun, t_autocorr, N_blocks)
        D_μ += tmp1
        D_v += tmp2^2
    end
    
    #push!(testing, [T, 𝒽, IntAutocorrTime([E, J[1,:], J[2,:]])])
    
    return [κ_μ C_μ D_μ abs.(M) ℙ; κ_v C_σ^2 D_v 0 0]
end

# ### Overall simulation routine

function MKuboSimulation(vertices, edges, num_histories, runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽)
    
    ks = range(1,length(T)*length(𝒽)*num_histories)
    args = [[deepcopy(vertices), deepcopy(edges), runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T[div(div(k-1,num_histories),length(𝒽))+1], 𝒽[rem(div(k-1,num_histories),length(𝒽))+1]] for k=ks]
    
    function hfun(args)
        return MKuboSingle(args...)
    end
    
    
    if multiProcess
        results = pmap(hfun, args)
    else
        results = Array{Any}(undef, length(ks))
        Threads.@threads for k in ks
            results[k] = hfun(args[k])
        end
    end 
    
    tmp = zeros(2,5,length(T),length(𝒽),num_histories) # rows for mean and stdv of κ,C
    for k in ks
        ni,h = divrem(k-1,num_histories) .+ (1,1)
        n,i = divrem(ni-1,length(𝒽)) .+ (1,1)
        
        tmp[:,:,n,i,h] = results[k]
    end
    tmp = sum(tmp, dims=5)
    
    # average over observables for all histories - okay b/c iid random variables
    tmp[2,:,:,:] = sqrt.(tmp[2,:,:,:])
    tmp ./= num_histories
        
    return tmp[1,1,:,:], tmp[1,2,:,:], tmp[1,3,:,:], tmp[1,4,:,:], tmp[1,5,:,:], tmp[2,1,:,:], tmp[2,2,:,:], tmp[2,3,:,:], tmp[2,4,:,:], tmp[2,5,:,:]
end
