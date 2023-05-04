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

# ### Single-Flip Dynamics

# #### 8-Vertex

@everywhere function A(edges, vertex) # calculates A at given vertex for 8-vertex model
    A = 1
    for α in vertex.δ # product of all adjacent spins
        A *= (-1)^edges[α].σ
    end

    return A
end

# #### 6-Vertex

@everywhere function Q(edges, vertex) # calculates B at given vertex for 6-vertex model
    Q = 0
    for α in vertex.δ # sum of all adjacent spins
        Q += (-1)^edges[α].σ
    end
    
    return Q
end

@everywhere function B(edges, vertex) # calculates A at given vertex for 6-vertex model
    return Q(edges, vertex)^2
end

# #### Both

@everywhere function ϵ(vertices, edges, vertex, 𝒽)
    ϵ = -λ*A(edges, vertex) + ξ*B(edges, vertex)
    
    for α in vertex.δ
        ϵ += 0.5 * (edges[α].D - 𝒽*edges[α].σ)
    end

    return ϵ
end

@everywhere function ΔE_flip(vertices, edges, β, 𝒽)
    v1 = vertices[edges[β].∂[1]]
    v2 = vertices[edges[β].∂[2]]
    σ = (-1)^edges[β].σ

    return 2*λ*(A(edges, v1) + A(edges, v2)) - 4*ξ*(σ*(Q(edges, v1) + Q(edges, v2)) - 2) + 2*𝒽*σ
end

@everywhere function Δj_flip(vertices, edges, β)
    v1 = vertices[edges[β].∂[1]]
    v2 = vertices[edges[β].∂[2]]
    σ = (-1)^edges[β].σ
    
    return λ*(A(edges, v2) - A(edges, v1)) - 2*ξ*σ*(Q(edges, v2) - Q(edges, v1))
end

# ### Double-Flip Dynamics

@everywhere function ΔE_2flip(vertices, edges, 𝜷, 𝒊, i, 𝒽)
    if 𝜷[1] == 𝜷[2]
        return 0
    end
    𝐯 = [vertices[𝒊[1]], vertices[𝒊[2]]]
    v = vertices[i]
    𝛔 = [(-1)^edges[𝜷[1]].σ, (-1)^edges[𝜷[2]].σ]
    
    return 2*λ*(A(edges, 𝐯[1]) + A(edges, 𝐯[2])) - 4*ξ*(𝛔[1]*Q(edges, 𝐯[1]) + 𝛔[2]*Q(edges, 𝐯[2]) - 2) + (2*𝒽 - 4*ξ*Q(edges, v) - 8*ξ*𝛔[1])*(𝛔[1] + 𝛔[2]) 
end

@everywhere function Δj_2flip(vertices, edges, 𝜷, 𝒊, 𝒽) # current flow from 𝒊[1] to 𝒊[2] via 𝜷[1] then 𝜷[2]
    𝐯 = [vertices[𝒊[1]], vertices[𝒊[2]]]
    𝛔 = [(-1)^edges[𝜷[1]].σ, (-1)^edges[𝜷[2]].σ]
    
    if ξ == 0
        return λ*(A(edges, 𝐯[1]) - A(edges, 𝐯[2])) + 0.5*𝒽*((-1)^edges[𝜷[2]].σ - (-1)^edges[𝜷[1]].σ)
    end
    
    return λ*(A(edges, 𝐯[2]) - A(edges, 𝐯[1])) - 2*ξ*(𝛔[2]*Q(edges, 𝐯[2]) - 𝛔[1]*Q(edges, 𝐯[1])) + 0.5*𝒽*((-1)^edges[𝜷[2]].σ - (-1)^edges[𝜷[1]].σ)
end
