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
            e_dir = cells[1][e.∂[1]].x-cells[1][e.∂[2]].x

            for d in 1:length(e_dir) # make sure to correct for PBCs
                e_dir[d] /= (abs(e_dir[d])>1) ? -abs(e_dir[d]) : 1
            end

            e.σ = (e_dir[1] ≈ e_dir[2])
        end
        
        return
    end
    
    for e in cells[2] # square or kagome spin ice
        e.σ = (cells[1][e.∂[1]].x[2]-cells[1][e.∂[2]].x[2]==0)
    end
    
end

