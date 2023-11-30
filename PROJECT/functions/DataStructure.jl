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


# NOTE: we DEMAND that the Links in basis all contain boundary elements with a positive orientation (user-defined)  


# ### Functions to generate an arbitary lattice with vertices (0-cells) and edges (1-cells)

@everywhere function LatticeGrid(L, PBC, Basis)
    
    Verts0, Links0, n, Scale = Basis
    
    @assert length(L)==length(PBC)
    dim = length(L)
    
    N = prod(L) # number of unit cells
    
    cells = []
    
    for c in 1:length(Links0)+1
        push!(cells, [Cell(zeros(dim), [], [], [], []) for _ in 1:n[c]*N])
    end
    
    # define indexing convention
    function X_to_I(X, L)
        I = X[1]
        for d in 2:length(L)
            I += prod(L[1:d-1]) * (X[d]-1)
        end
        return I
    end

    function I_to_X(I, L)
        X = zeros(Int, length(L))
        J = I
        for d in 1:length(L)
            J, X[d] = divrem(J-1,L[d]) .+ (1,1)
        end
        return X
    end

    
    # place down all the vertices without connecting any edges/plaquettes/etc
    x0 = ones(length(L)) # origin for 1-indexing positions as done in X (just for convenience because Julia)
    for I in 1:N
        X = I_to_X(I, L) - x0 # don't forget to subtract the origin!
        
        for i in 1:n[1]
            𝐢 = n[1]*(I-1)+i # absolute index
            cells[1][𝐢].x = Verts0[i].x .+ X
        end
    end
    
    # for later
    cellsToKill = [[] for _ in 1:length(cells)]
    cellShifts = [zeros(length(cells[c])) for c in 1:length(cells)]
    
    # go back through and link up the hyperedges in order of dimension
    for c in 2:length(cells)
        
        for I in 1:N # for each unit cell
            X = I_to_X(I, L)

            # attach hyperedges (in +ve quadrant)
            for α in 1:n[c]
                𝛂 = n[c]*(I-1) + α # absolute index of the relevant hyperedge

                # if hyperedge crosses an OBC then don't link it up at all
                ifLink = true
                for hypervertex in Links0[c-1][α]
                    dir = hypervertex[2]
                    for d in 1:dim 
                        if dir[d]<0 && X[d]==1 && !PBC[d]
                            ifLink = false
                        elseif dir[d]>0 && X[d]==L[d] && !PBC[d]
                            ifLink = false
                        end
                    end
                end

                if ifLink # if NOT crossing an OBC, link up this hyperedge
                    for hypervertex in Links0[c-1][α]
                        dir = hypervertex[2]
                        orientation = hypervertex[3]

                        Y = copy(X)
                        for d in 1:dim
                            Y[d] = (dir[d]<0 && X[d]==1) ? L[d] : ((dir[d]>0 && X[d]==L[d]) ? 1 : X[d]+dir[d])
                        end
                        J = X_to_I(Y, L) # cell index of the relevant hypervertex
                        𝐣 = hypervertex[1] + n[c-1]*(J-1) # absolute index of the relevant hypervertex 

                        # update the relevant boundary list
                        push!(cells[c][𝛂].∂, 𝐣)
                        push!(cells[c][𝛂].η, orientation)
                        
                        # also update the coboundary lists
                        push!(cells[c-1][𝐣].∂ᵀ, 𝛂)
                        push!(cells[c-1][𝐣].ηᵀ, orientation)
                    end
                else # if it IS crossing an OBC, mark the hyperedge for deletion after the rest of the complex is constructed
                    push!(cellsToKill[c], 𝛂)
                    cellShifts[c][𝛂:end] .-= 1
                end
            end
        end
    end
    
    # define the displacements for each edge running from the -ve boundary element to the +ve one
    for I in 1:N
        for e in 1:n[2]
            𝐞 = n[2]*(I-1) + e # absolute index of the edge
                
            V1 = Links0[1][e][1]
            V2 = Links0[1][e][2]
            
            cells[2][𝐞].x = (Verts0[V2[1]].x + V2[2]) - (Verts0[V1[1]].x + V1[2]) # displacement of edge V1 -> V2
            # note: WE WILL GO BACK AND FIX THESE ORIENTATIONS ONCE WE'VE DELETED THE UNNECESSARY EDGES
        end
    end
    
    # Note: the higher-dim x's are just zero for our purposes...
    
    
    # kill off hyperedges crossing OBCs and shift indices to compensate in order of dimension
    for c in 2:length(cells)

        deleteat!(cells[c], cellsToKill[c])

        # fix the hyperedge indices in the coboundary of each hypervertex
        for v in cells[c-1]
            for i in eachindex(v.∂ᵀ)
                v.∂ᵀ[i] += cellShifts[c][v.∂ᵀ[i]]
            end
        end
        
        # fix the hyperedge indices in the boundary of each hyperface
        if c < length(cells)
            for v in cells[c+1]
                for i in eachindex(v.∂ᵀ)
                    v.∂ᵀ[i] += cellShifts[c][v.∂ᵀ[i]]
                end
            end
        end
        
    end
    
    
    # fix edge displacement orientations
    for edge in cells[2]
       edge.x .*= edge.η[2] # displacement of edge V- -> V+  
    end
    
    
    # rescale all the cell positions/displacements
    for dcells in cells
        for c in dcells
            c.x .*= Scale
        end
    end
    
    # iteratively assign hyperedge positions from hypervertex positions in order of dimension
    #for c in 2:length(cells)
    #    for e in cells[c]
    #        e.x = zeros(length(L))
    #        for i in e.∂
    #            e.x += cells[c-1][i].x
    #        end
    #        e.x ./= length(e.∂)
    #    end
    #end
    
    # 
    
    return CellComplex(dim, cells), Scale
end


function Displacement(Δe, i, j) # gets the displacement of an edge e i -> j
    if !(i in Δe.∂ && j in Δe.∂)
        error("The given edge does not connect the given vertices")
    end
    
    kj = Δe.η[findfirst(Δe.∂ .== j)]
       
    return Δe.x .* kj
end