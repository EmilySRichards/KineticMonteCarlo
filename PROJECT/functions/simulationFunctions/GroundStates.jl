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

# ### Ground State Function

@everywhere function GroundState!(S, Δ, isPyrochlore)
    
    if !isSpinIce # ANY toric code
        for (e, Δe) in enumerate(Δ.cells[2])
            S.vals[e] = 1
        end
        
        return
    end
    
    if isPyrochlore # pyrochlore spin ice
        for (e, Δe) in enumerate(Δ.cells[2])
            S.vals[e] = (-1)^(Δe.x[1] ≈ Δe.x[2])
        end
        
        return
    end
    
    for (e, Δe) in enumerate(Δ.cells[2]) # square or kagome spin ice
        S.vals[e] = (-1)^(Δe.x[2]==0)
    end
    
end

