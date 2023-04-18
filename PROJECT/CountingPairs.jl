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

dir = dirname(pwd()) * "/PROJECT"
include(dir * "/functions/Preamble.jl")
@everywhere dir = dirname(pwd()) * "/PROJECT"

@everywhere global const sixVertex::Bool = false
@everywhere global const twoFlip::Bool = true
@everywhere global const δE::Int = sixVertex ? 8 : 4

@everywhere include(dir * "/functions/DataStructure.jl")
@everywhere include(dir * "/functions/Plotting.jl")
@everywhere include(dir * "/functions/Statistics.jl")
@everywhere include(dir * "/functions/Simulation.jl")

@everywhere include(dir * "/functions/simulationFunctions/MicroKubo.jl")

# +
L = [10, 10]
PBC = [true, true]
therm_runtime = 1000
runtime = 1000
num_histories = 10
tau = 2:floor(Int64, 0.75*runtime)
𝒽 = range(0.0, 2.0, length=15)

T = range(0.01, 10.0, length=50);

# +
pairCount = zeros(length(𝒽), length(T), num_histories)

for i in 1:num_histories
    for j in eachindex(T)
        for k in eachindex(𝒽)
            vertices, edges = CubicGrid(L, PBC);
            MicroKuboSetup(vertices, edges, therm_runtime, T[j], 𝒽[k], false)

            # Find all particles
            prtclIndices = []
            for (i,v) in enumerate(vertices)
                Ai = A(edges, v)
                if Ai < 0
                    push!(prtclIndices, i)
                end
            end
            
            if length(prtclIndices) == 0
                pairCount[k,j,i] = NaN
                continue
            end

            # Find pairs
            for p in prtclIndices
                for q in prtclIndices
                    if p != q
                        Δ = vertices[p].x - vertices[q].x
                        D = 0  
                        # correct for PBCs
                        for d in 1:length(Δ)
                            Δ[d] /= (abs(Δ[d])==L[d]-1) ? -abs(Δ[d]) : 1
                            D += Δ[d]^2
                        end

                        if D == 1
                            pairCount[k,j,i] += 1
                        end
                    end
                end
            end

            pairCount[k,j,i] /= 2*length(prtclIndices) # *fraction* of particles that are in pairs
        end
    end
end

# +
pairCount = mean(pairCount, dims=3)

for k in eachindex(𝒽)
    plot(T, pairCount[k, :])
end
savefig("figs/FractionOfParticlesInPairs.png")
# -

# Obviously the above will fail when there are multiple **overlapping** pairs, but at low temperatures it should give us a decent estimate.
#
# In terms of an analytic expectation, I'm really not sure. Recalling that particles are generated from the GS by creating blue spin-down strings, we should seek to minimise these strings and therefore increase the number of particles which exist in pairs.
