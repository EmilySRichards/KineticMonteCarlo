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

@everywhere function MicroKuboSetup(Δ, therm_runtime, T, 𝒽, isRandom)
    
    S = CreateField(Δ, 1)
    
    if isRandom # initialise entire system in random state
        S.vals = (-1).^rand(Bool, length(Δ.cells[2]))
    else # initialise entire system in ground state
        GroundState!(S, Δ)
    end
    
    Energy = zeros(therm_runtime+1) # just set initial energy to zero since we only need the variance
    
    # thermalise entire system
    for t in 1:therm_runtime
        Energy[t+1] = Energy[t]
        for _ in Δ.cells[2]
            e = rand(eachindex(Δ.cells[2]))
            ΔE, _ = EnergyChange(S, Δ, [e], 𝒽)

            if ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/T)
                S.vals[e] = -S.vals[e]
                Energy[t+1] += ΔE
            end
        end
    end
    
    return S, Energy[2:end] # cut out initial energy for consistency with other observables
end

# ### Set up the geometry - with magnetisation-conserving stage

@everywhere function MicroKuboSetup_2flip(S, Δ, therm_runtime, T, 𝒽, isRandom)
    
    # initial thermalisation WITHOUT conserving magnetisation
    S, _ = MicroKuboSetup(S, Δ, therm_runtime, T, 𝒽, isRandom)
    
    # additional thermalisation step WITH FIXED MAGNETISATION!
    Energy = zeros(therm_runtime+1) # just set initial energy to zero since we only need the variance
    
    for t in 1:therm_runtime
        Energy[t+1] = Energy[t]
        for _ in Δ.cells[1]
            i = rand(eachindex(Δ.cells[1])) # random vertex
            edges = sample(Δ.cells[1][i].∂ᵀ, 2; replace=false) # pick two different spins in its coboundary
            
            ΔE, J = EnergyChange(S, Δ, edges, 𝒽) 

            if (ΔE <= 0 || rand(Uniform(0,1)) < exp(-ΔE/T)) && (S.vals[edges[1]] + S.vals[edges[2]] == 0)
                for e in edges
                    S.vals[e] = -S.vals[e]
                end
                
                Energy[t+1] += ΔE
            end
        end
    end
    
    return S, Energy[2:end] # cut out initial energy for consistency with other observables
end

# ### Single spin-flip dynamics routine 

@everywhere function MicroKubo(S, Δ, runtime, 𝒽)
    Current = zeros(Float64, length(Δ.cells[2][1].x), runtime)
    Polarisation = zeros(Float64, length(Δ.cells[2][1].x), runtime)
    
    for t in 1:runtime
        for _ in Δ.cells[2] 
            e = rand(eachindex(Δ.cells[2]))
            ΔE, J = EnergyChange(S, Δ, [e], 𝒽)
            
            if ΔE == 0
                S.vals[e] = -S.vals[e]
            
                Current[:,t] += Δ.cells[2][e].x * J[1]
            end
        end

        if t<runtime
            Polarisation[:,t+1] = copy(Polarisation[:,t]) + Current[:,t]
        end
    end
    
    return Current, Polarisation
end

# ### Double spin-flip dynamics routine

@everywhere function MicroKubo_2flip(S, Δ, runtime, 𝒽)
    Current = zeros(Float64, length(Δ.cells[2][1].x), runtime)
    Polarisation = zeros(Float64, length(Δ.cells[2][1].x), runtime)
    
    for t in 1:runtime
        for _ in Δ.cells[1]
            # propose flips
            i = rand(eachindex(Δ.cells[1])) # random vertex
            edges = sample(Δ.cells[1][i].∂ᵀ, 2; replace=false) # pick two different spins in its coboundary
            
            # calculate overall energy change and current density between the two unshared vertices
            ΔE, J = EnergyChange(S, Δ, edges, 𝒽) 
            
            # decide whether to accept and perform the move
            if ΔE == 0 && (S.vals[edges[1]] + S.vals[edges[2]] == 0) # energy AND magnetisation conserved
                for e in edges
                    S.vals[e] = -S.vals[e]
                    
                    Current[:,t] += Δ.cells[2][e].x * J[1] # note no factor of 1/2 b/c only sum each pair of sites once
                end
            end
            
            if t<runtime
                Polarisation[:,t+1] = copy(Polarisation[:,t]) + Current[:,t]
            end
        end
    end
    
    return Current, Polarisation
end

# ### Single Simulation Run

@everywhere function MKuboSingle(Δ, runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽, allComponents)
    
    dim = allComponents ? length(Δ.cells[2][1].x) : 1
    
    Cfun = (E) -> var(E) / T^2 / length(Δ.cells[2])
    κfun = (S) -> mean(S) / T^2 / length(Δ.cells[2])
    Dfun = (E,S) -> κfun(S) / Cfun(E)
    
    tmax = runtime-t_therm
    
    # -- 0. Run Simulation --
    #if twoFlip
    #    S, E = MicroKuboSetup_2flip(Δ, therm_runtime, T, 𝒽, false)
    #else
    S, Energy = MicroKuboSetup(Δ, therm_runtime, T, 𝒽, false)
    #end
        
    M = mean(S.vals)
    
    if twoFlip
        Current, Polarisation = MicroKubo_2flip(S, Δ, runtime, 𝒽)
    else
        Current, Polarisation = MicroKubo(S, Δ, runtime, 𝒽)
    end
    
    
    
    
    #print(mean(Current, dims=2), "\n\n")
    #print(Polarisation[1,:] - cumsum(Current[1,:]), "\n\n")
    #if mean(cumsum(Current, dims=2), dims=2) ≈ zeros(size(Current, 1))
    #    print(T, " ")
    #end
    
    
    
    
    
    # cut out thermalisation time
    Current = Current[:,t_therm+1:end]
    Polarisation = Polarisation[:,t_therm+1:end]
    Energy = Energy[t_therm+1:end]
    
    # -- 1. Heat Capacity --
    C_μ, C_s = Estimator(Bootstrap, [Energy], Cfun, t_autocorr, N_blocks)
    
    
    # -- 2. Thermal Conductivity and Diffusivity --
    result = zeros(dim, dim, 2, 3)
    
    if allComponents
        statistic = zeros(Float64, tmax, dim, dim)
        for t in 1:tmax
            for τ in 0:min(tmax-t, t_cutoff)
                # symmetric part: #statistic[t,:,:] += (τ==0 ? 0.5 : 1.0) * 0.5 * (J[:,t+τ] .* J[:,t]' + J[:,t] .* J[i,t+τ]') * tmax/(tmax-τ)
                statistic[t,:,:] += 0.5 .* Current[:,t+τ] .* Current[:,t]' .* tmax/(tmax-τ)
            end
            statistic[t,:,:] -= 0.5 .* Current[:,t] .* Polarisation[:,t]'
        end
        #statistic .*= prod(scale) # rescaling to correct for scaling of unit cells

        κ_μ = zeros(dim, dim)
        κ_s = zeros(dim, dim)
        D_μ = zeros(dim, dim)
        D_s = zeros(dim, dim)  
        for i in 1:dim
            for j in 1:dim
                κ_μ[i,j], κ_s[i,j] = Estimator(Bootstrap, [statistic[:,i,j]], κfun, t_autocorr, N_blocks)
                D_μ[i,j], D_s[i,j] = Estimator(Bootstrap, [Energy, statistic[:,i,j]], Dfun, t_autocorr, N_blocks)
            end
        end
    else
        statistic = zeros(Float64, tmax)
        for t in 1:tmax
            for τ in 0:min(tmax-t, t_cutoff)
                statistic[t] += (τ==0 ? 0.5 : 1.0) * Current[1,t+τ] * Current[1,t] * tmax/(tmax-τ)
            end
        end
        #statistic .*= prod(scale) # rescaling to correct for scaling of unit cells
        
        κ_μ, κ_s = Estimator(Bootstrap, [statistic], κfun, t_autocorr, N_blocks) # note rescaling b/c propto V
        D_μ, D_s = Estimator(Bootstrap, [Energy, statistic], Dfun, t_autocorr, N_blocks)
    end
    
    result[:,:,1,1] .= κ_μ
    result[:,:,2,1] .= κ_s.^2
    result[:,:,1,2] .= C_μ
    result[:,:,2,2] .= C_s.^2
    result[:,:,1,3] .= D_μ
    result[:,:,2,3] .= D_s.^2
    
    return result
end

# +
#κ_μ = 0
#κ_s = 0
#D_μ = 0
#D_s = 0
#for τ in 0:t_cutoff
#    statistic = (τ==0 ? 0.5 : 1.0) .* J[1,:] .* circshift(J[1,:], -τ)
#    statistic /= length(statistic)
#    
#    tmp1, tmp2 = Estimator(Bootstrap, [statistic[1:end-τ]], κfun, t_autocorr, N_blocks)
#    κ_μ += tmp1
#    κ_s += tmp2
#    
#    tmp1, tmp2 = Estimator(Bootstrap, [E[1:end-τ], statistic[1:end-τ]], Dfun, t_autocorr, N_blocks)
#    D_μ += tmp1
#    D_s += tmp2
#end
# -

# ### Overall simulation routine

function MKuboSimulation(L, PBC, Basis, num_histories, runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽, allComponents)
    
    Δ, scale = LatticeGrid(L, PBC, Basis)
    
    ks = range(1,length(T)*length(𝒽)*num_histories)
    args = [[Δ, runtime, therm_runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T[div(div(k-1,num_histories),length(𝒽))+1], 𝒽[rem(div(k-1,num_histories),length(𝒽))+1], allComponents] for k=ks]
    
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
    
    dim = allComponents ? length(Δ.cells[1][1].x) : 1
    tmp = zeros(dim, dim, 2, 3, length(T), length(𝒽), num_histories) # rows for mean and stdv of κ,C
    for k in ks
        ni,h = divrem(k-1,num_histories) .+ (1,1)
        n,i = divrem(ni-1,length(𝒽)) .+ (1,1)
        
        tmp[:,:,:,:,n,i,h] = results[k]
    end
    
    # average over observables for all histories - okay b/c iid random variables
    tmp = sum(tmp, dims=7)
    tmp[:,:,2,:,:,:] = sqrt.(tmp[:,:,2,:,:,:])
    tmp ./= num_histories
        
    return tmp[:,:,1,1,:,:], tmp[1,1,1,2,:,:], tmp[:,:,1,3,:,:], tmp[:,:,2,1,:,:], tmp[1,1,2,2,:,:], tmp[:,:,2,3,:,:]
end
