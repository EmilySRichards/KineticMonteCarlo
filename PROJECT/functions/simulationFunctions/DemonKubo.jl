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

@everywhere function DemonKuboSetup(cells, z, T, 𝒽)
    vertices = cells[1]
    edges = cells[2]
    
    g = 2*𝒽 - δE*ceil(2*𝒽/δE)
    
    Dfun = (T) -> δE/(exp(δE/T)-1) - g/(exp(-g/T)+1)
    
    # REinitialise entire system in ground state
    GroundState!(cells)

    # calculate total demon energy for given temperature T
    D_tot = 0
    for v in vertices
        D_tot -= ϵ(vertices, edges, v, 𝒽)
    end
    
    D_tot += length(edges) * Dfun(T)
    D_tot += length(vertices) * (-λ*Asv([T], 𝒽, z)[1] + ξ*Bsv([T], 𝒽, z)[1] - 𝒽*Magnetisation([T], 𝒽, z)[1])
    
    # randomly increment demon energies
    idxs = collect(eachindex(edges)) 
    while D_tot>0 # while loop        
        hterm = (𝒽==0 || length(idxs)==0) ? false : rand(Bool) # randomly pick increment unit
        ΔD = hterm ? 2*𝒽 : δE
        
        α = hterm ? rand(idxs) : rand(eachindex(edges)) # pick a random (valid) edge
        
        edges[α].D += ΔD # increment its demon energy
        D_tot -= ΔD # decrement the total energy left to distribute
        
        if hterm # update which edges are valid
            deleteat!(idxs, findfirst(idxs .== α))
        end
    end
end

# ### Demon dynamics routine 

@everywhere function DemonKubo(vertices, edges, runtime, 𝒽)
    
    J = zeros(Float64, (length(vertices[1].x), runtime))
    D = zeros(Float64, (runtime+1))
    E = zeros(Float64, (runtime+1)) # just set zero of energy to 0 since we'll only use the variance
    
    # set initial demon energies
    for edge in edges
        D[1] += edge.D
    end
    
    for t in 1:runtime
        E[t+1] = E[t]
        D[t+1] = D[t]
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)
            
            if edges[β].D >= ΔE
                Δj_β = Δj_flip(vertices, edges, β)
                
                edges[β].σ = !edges[β].σ
                edges[β].D -= ΔE
                
                E[t+1] += ΔE
                D[t+1] -= ΔE
                
                # update current
                J[:,t] += edges[β].x * Δj_β # note no factor of 1/2 b/c only sum each edge once
            end
        end
    end
    
    D ./= length(edges)
    
    return J, D[2:end], E[2:end]
end

# ### Single Simulation Run

@everywhere function DKuboSingle(cells, scale, z, runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽)
    vertices = cells[1]
    edges = cells[2]
    
    # -- -1. Define Observables --
    g = 2*𝒽 - δE*ceil(2*𝒽/δE)
    
    Dfun = (T) -> δE/(exp(δE/T)-1) - g/(exp(-g/T)+1)
    Tfun = (D) -> (𝒽==0) ? δE/log(1.0 + δE/mean(D)) : find_zero((T) -> sign(T)*Dfun(abs(T)) - mean(D), (-20, 20))
    
    #CDfun = (D) -> ((δE/Tfun(D))^2 * exp(δE/Tfun(D))/(exp(δE/Tfun(D))-1)^2 + (g/Tfun(D))^2 * exp(g/Tfun(D))/(exp(g/Tfun(D))+1)^2)
    #C0fun = (D,E) -> Var(E) / Tfun(D)^2 / length(edges)
    #Cfun = (D,E) -> 1/(1/C0fun(D,E) - 1/CDfun(D))
    CDfun = (D) -> length(edges) * ((δE/Tfun(D))^2 * exp(δE/Tfun(D))/(exp(δE/Tfun(D))-1)^2 + (g/Tfun(D))^2 * exp(g/Tfun(D))/(exp(g/Tfun(D))+1)^2)
    Cfun = (D,E) -> CDfun(D) * Var(E) /(CDfun(D)*Tfun(D)^2 - Var(E)) / length(edges)
    
    κfun = (D,S) -> sum(S) / Tfun(D)^2 / length(edges)
    Difffun = (D,E,S) -> κfun(D, S) / Cfun(D, E)
    
    tmax = runtime-t_therm
    
    # -- 0. Run Simulation --
    DemonKuboSetup(cells, z, T, 𝒽)
    J, D, E = DemonKubo(vertices, edges, runtime, 𝒽)

    # cut out thermalisation time
    J = J[:,t_therm+1:end]
    D = D[t_therm+1:end]
    E = E[t_therm+1:end]
    
    #t_autocorr = IntAutocorrTime([D, E, J[1,:]])
    
    # -- 1. Temperature --
    T_μ, T_s = Estimator(Bootstrap, [D], Tfun, t_autocorr, N_blocks)
    
    # -- 2. Heat Capacity --
    C_μ, C_s = Estimator(Bootstrap, [D, E], Cfun, t_autocorr, N_blocks)
    
    # -- 3. Thermal Conductivity and Diffusivity--
    statistic = zeros(Float64, tmax)
    for t in 1:tmax
        for τ in 0:min(tmax-t, t_cutoff)
            statistic[t] += (τ==0 ? 0.5 : 1.0) * J[1,t] * J[1,t+τ] / (tmax-τ)
        end
    end
    #statistic .*= prod(scale) # rescaling to correct for scaling of unit cells
    
    κ_μ, κ_s = Estimator(Bootstrap, [D, statistic], κfun, t_autocorr, N_blocks)
    D_μ, D_s = Estimator(Bootstrap, [D, E, statistic], Difffun, t_autocorr, N_blocks)
    
    return [T_μ κ_μ C_μ D_μ; T_s^2 κ_s^2 C_s^2 D_s^2]
end

# +
# Old conductivity calculation (ASSUMES EACH CORRELATION TERM INDEP BUT DOESN'T SUFFER FROM NAN SAMPLES AS MUCH!!)
#κ_μ = 0
#κ_s = 0
#D_μ = 0
#D_s = 0
#for τ in 0:t_cutoff
#    statistic = (τ==0 ? 0.5 : 1.0) .* J[1,:] .* circshift(J[1,:], -τ)
#    statistic /= length(statistic)
#    
#    tmp1, tmp2 = Estimator(Bootstrap, [D[1:end-τ], statistic[1:end-τ]], κfun, t_autocorr, N_blocks)
#    κ_μ += tmp1
#    κ_s += tmp2
#    
#    tmp1, tmp2 = Estimator(Bootstrap, [D[1:end-τ], E[1:end-τ], statistic[1:end-τ]], Dfun, t_autocorr, N_blocks)
#    D_μ += tmp1
#    D_s += tmp2
#end
# -

# ### Overall simulation routine

function DKuboSimulation(L, PBC, Basis, num_histories, runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T, 𝒽)
    
    cells, scale = LatticeGrid(L, PBC, Basis)
    vertices = cells[1]
    edges = cells[2]
    z = Coordination(Basis)
    
    ks = range(1,length(T)*length(𝒽)*num_histories)
    args = [[deepcopy(cells), scale, z, runtime, t_therm, t_autocorr, N_blocks, t_cutoff, T[div(div(k-1,num_histories),length(𝒽))+1], 𝒽[rem(div(k-1,num_histories),length(𝒽))+1]] for k=ks]
    
    function hfun(args)
        return DKuboSingle(args...)
    end
    
    
    if multiProcess
        results = pmap(hfun, args)
    else
        results = Array{Any}(undef, length(ks))
        Threads.@threads for k in ks
            results[k] = hfun(args[k])
        end
    end 
    
    
    tmp = zeros(2, 4, length(T), length(𝒽), num_histories) # rows for mean and stdv of T,κ,C
    for k in ks
        ni,h = divrem(k-1,num_histories) .+ (1,1)
        n,i = divrem(ni-1,length(𝒽)) .+ (1,1)
        
        tmp[:,:,n,i,h] = results[k]
    end
    tmp = sum(tmp, dims=5)
    
    # average over observables for all histories - okay b/c iid random variables
    tmp[2,:,:,:] = sqrt.(tmp[2,:,:,:])
    tmp ./= num_histories
        
    return tmp[1,1,:,:], tmp[1,2,:,:], tmp[1,3,:,:], tmp[1,4,:,:], tmp[2,1,:,:], tmp[2,2,:,:], tmp[2,3,:,:], tmp[2,4,:,:]
end
