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

@everywhere function Atilde(edges, vertex) # calculates A at given vertex for 6-vertex model
    A = 0
    for α in vertex.δ # product of all adjacent spins
        A += (-1)^edges[α].σ
    end
    
    return A^2
end

@everywhere function B(edges, vertex, β) # calculates B at given vertex for 6-vertex model
    B = 0
    for α in vertex.δ # sum of all adjacent spins EXCEPT β
        B += (α!=β) ? (-1)^edges[α].σ : 0
    end
    
    return B
end

# #### Both

@everywhere function ΔE_flip(vertices, edges, β, 𝒽)
    if sixVertex
        return -4*(-1)^edges[β].σ*(B(edges, vertices[edges[β].∂[1]], β) + B(edges, vertices[edges[β].∂[2]], β)) + 2*𝒽*(-1)^edges[β].σ
    else
        return 2*(A(edges, vertices[edges[β].∂[1]]) + A(edges, vertices[edges[β].∂[2]])) + 2*𝒽*(-1)^edges[β].σ
    end
end

@everywhere function Δj_flip(vertices, edges, β)
    if sixVertex
        Bi = B(edges, vertices[edges[β].∂[1]], β)
        Bj = B(edges, vertices[edges[β].∂[2]], β)
        
        return -2*(-1)^edges[β].σ*(Bj-Bi)
    else
        Ai = A(edges, vertices[edges[β].∂[1]])
        Aj = A(edges, vertices[edges[β].∂[2]])
        
        return (Aj-Ai)
    end
end

# ### Double-Flip Dynamics

@everywhere function ΔE_2flip(vertices, edges, 𝜷, 𝒊, 𝒽)
    if 𝜷[1] == 𝜷[2]
        return 0
    end
    return 2*(A(edges, vertices[𝒊[1]]) + A(edges, vertices[𝒊[2]])) + 2*𝒽*((-1)^edges[𝜷[1]].σ + (-1)^edges[𝜷[2]].σ) 
end

@everywhere function Δj_2flip(vertices, edges, 𝜷, 𝒊, 𝒽) # current flow from 𝒊[1] to 𝒊[2] via 𝜷[1] then 𝜷[2]
    return (A(edges, vertices[𝒊[2]]) - A(edges, vertices[𝒊[1]])) + 0.5*𝒽*((-1)^edges[𝜷[2]].σ - (-1)^edges[𝜷[1]].σ)
end
