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

# +
dir = dirname(pwd()) * "/PROJECT"

using Pkg
Pkg.activate(dir)

using Distributed

# +
global const multiProcess = (nworkers()>1) ? true : false # only use multiprocessing if run with -p, otherwise use -t threads by default

if multiProcess
    print(nworkers())
    
    @everywhere using Pkg
    @everywhere Pkg.activate(dir)
    @everywhere using Distributed, StatsBase, Statistics, Distributions, Roots, PyPlot, LsqFit, Dates # , ProfileVega # , Bootstrap #, ProgressMeter, ProfileVega, JLD
else
    print(Threads.nthreads())
    
    using StatsBase, Statistics, Distributions, Roots, PyPlot, LsqFit, Dates # , ProfileVega # , Bootstrap #, ProgressMeter, ProfileVega, JLD
end

using Colors, PlotUtils, Graphs
# -

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
𝒽 = range(0.0, 2.0, length=7)

T = range(0.01, 10.0, length=20);

# +
pairCount = zeros(Int64, length(T), length(𝒽), num_histories)

for i in 1:num_histories
    for j in eachindex(T)
        for k in eachindex(𝒽)
            vertices, edges = CubicGrid(L, PBC);
            MicroKuboSetup(vertices, edges, therm_runtime, T[j], 𝒽, false)

            # Find all particles
            prtclIndices = []
            for (i,v) in enumerate(vertices)
                Ai = A(edges, v)
                if Ai < 0
                    push!(prtclIndices, i)
                end
            end

            # Find pairs

            pairCount = 0
            for i in prtclIndices
                for j in prtclIndices
                    if i != j
                        Δ = vertices[i].x - vertices[j].x
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

            #pairCount[k,j,i] /= 2 # correct for 2x-counting
            pairCount[k,j,i] /= length(prtclIndices) # *fraction* of particles that are in pairs
        end
    end
end

# +
pairCount = mean(pairCount, dims=3)

for k in eachindex(𝒽)
    plot(T, pairCount[k, :])
end
# -

# Obviously the above will fail when there are multiple **overlapping** pairs, but at low temperatures it should give us a decent estimate.
#
# In terms of an analytic expectation, I'm really not sure. Recalling that particles are generated from the GS by creating blue spin-down strings, we should seek to minimise these strings and therefore increase the number of particles which exist in pairs.
