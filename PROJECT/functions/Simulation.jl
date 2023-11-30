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

# ### Local energy density

@everywhere function ϵ(S, D, Δ, i, 𝒽)
    ϵi = 0
    
    # energy without demons
    ϵi += (λ!=0) ? -λ*Star(S, Δ, i) : 0
    ϵi += (ξ!=0) ? ξ*(-Boundary(S, Δ, i))^2 : 0
    
    Δi = Δ.cells[1][i]
    for e in Δi.∂ᵀ
        if isSpinIce
            ϵi -= 0.5 * 𝒽 * (Δ.cells[2][e].x[1]) * GetCpt(S, e, true) # magnetic field is a vector in the plane of the spins, taken to lie along [10...]
        else
            ϵi -= 0.5 * 𝒽 * GetCpt(S, e, false) # magnetic field is a scalar in a direc perp to the lattice
        end
    end
    
    # add on the local demon energy
    if D != nothing
        Δi = Δ.cells[1][i]
        for e in Δi.∂ᵀ
            ϵi += 0.5 * GetCpt(D, e, false)
        end
    end
    
    # note that term-by-term, all the above (De, Se_TC and re*Se_SI are EVEN under e -> -e so the above decompositions work)
    
    return ϵi
end

@everywhere function ϵ(S, D, Δ, 𝒽)
    ϵs = CreateField(Δ, 0)
    
    for i in eachindex(Δ.cells[1])
        ϵs.vals[i] = ϵ(S, D, Δ, i, 𝒽)
    end
    
    return ϵs
end




# ### Energy change and current function 

@everywhere function EnergyChange(S, Δ, edges, 𝒽, D=nothing) # takes an array of edges to flip in this move

    S′ = deepcopy(S)
    for e in edges
        S′.vals[e] = -S′.vals[e]
    end
        
    # affected vertices - union of the boundary sets of edges to flip
    visited = []
    ΔE = 0
    J = zeros(length(edges))
    for (n, e) in enumerate(edges)
        for (i, k) in zip(Δ.cells[2][e].∂, Δ.cells[2][e].η)
            Δϵi = ϵ(S′, D, Δ, i, 𝒽) - ϵ(S, D, Δ, i, 𝒽)
            
            if !(i in visited) # (avoids repeats!!)    
                ΔE += Δϵi
                push!(visited)
            end
            
            J[n] += k*Δϵi # ORIENTATION OF CURRENT HANDLED MANUALLY HERE --> BAD!!!
        end
    end
    
    J ./= 2 # dividing by the size of the coboundary of each edge (enforces continuity and energy cons)
        
    return ΔE, J
end