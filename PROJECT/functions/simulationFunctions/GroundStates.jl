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

@everywhere function GroundState!(cells, isPyrochlore)
    
    if !isSpinIce # ANY toric code
        for e in cells[2]
            e.σ = false
        end
        
        return
    end
    
    if isPyrochlore # pyrochlore spin ice
        for e in cells[2]
            e.σ = (e.x[1] ≈ e.x[2])
        end
        
        return
    end
    
    for e in cells[2] # square or kagome spin ice
        e.σ = (e.x[2]==0)
    end
    
end

