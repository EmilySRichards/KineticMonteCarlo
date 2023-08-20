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

# ### Data structure for an arbitrary cell

@everywhere mutable struct Cell
    σ::Bool # false = +1, true = -1
    D::Float64 # demon energy in units of δE
    x::Array{Float64} # coords
    ∂::Array{UInt32} # boundary
    δ::Array{UInt32} # coboundary
end

# ### Functions to generate an arbitary lattice with vertices (0-cells) and edges (1-cells)

@everywhere function LatticeGrid(L, PBC, Basis)
    
    Verts0, Links0, n, Scale = Basis
    
    @assert length(L)==length(PBC)
    dim = length(L)
    
    N = prod(L) # number of unit cells
    
    cells = []
    
    for c in 1:length(Links0)+1
        push!(cells, [Cell(false, 0, zeros(dim), [], []) for _ in 1:n[c]*N])
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

                        Y = copy(X)
                        for d in 1:dim
                            Y[d] = (dir[d]<0 && X[d]==1) ? L[d] : ((dir[d]>0 && X[d]==L[d]) ? 1 : X[d]+dir[d])
                        end
                        J = X_to_I(Y, L) # cell index of the relevant hypervertex
                        𝐣 = hypervertex[1] + n[c-1]*(J-1) # absolute index of the relevant hypervertex 

                        # update the relevant boundary and coboundary lists
                        push!(cells[c][𝛂].∂, 𝐣)
                        push!(cells[c-1][𝐣].δ, 𝛂)
                    end
                else # if it IS crossing an OBC, mark the hyperedge for deletion after whole complex is constructed
                    push!(cellsToKill[c], 𝛂)
                    cellShifts[c][𝛂:end] .-= 1
                end
            end
        end
    end    
    
    
    # define the displacements for each edge - by doing this we don't need to worry about PBCs, we're in a basis relative to the torus
    for I in 1:N
        for α in 1:n[2]
            𝛂 = n[2]*(I-1) + α # absolute index of the edge
            
            V1 = Links0[1][α][1]
            V2 = Links0[1][α][2]
            
            cells[2][𝛂].x = (Verts0[V2[1]].x + V2[2]) - (Verts0[V1[1]].x + V1[2]) # displacement of edge V1->V2
        end
    end
    
    # Note: the higher-dim x's are just zero for our purposes...
    
    
    # kill off hyperedges crossing OBCs and shift indices to compensate in order of dimension
    for c in 2:length(cells)

        deleteat!(cells[c], cellsToKill[c])

        # fix the hyperedge indices in the coboundary of each hypervertex
        for v in cells[c-1]
            for i in eachindex(v.δ)
                v.δ[i] += cellShifts[c][v.δ[i]]
            end
        end
        
        # fix the hyperedge indices in the boundary of each hyperface
        if c < length(cells)
            for v in cells[c+1]
                for i in eachindex(v.∂)
                    v.∂[i] += cellShifts[c][v.∂[i]]
                end
            end
        end
        
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
    
    return cells, Scale
end


# ### Useful Functions

@everywhere function NearestNeighbourEdge(vertices, edges, α)
    neighbours = []
    for i in edges[α].∂
        δi = vertices[i].δ # coboundary of vertex i
        append!(neighbours, δi) # append coboundary to nn list
    end
    
    unique!(neighbours) # remove repeats
    deleteat!(neighbours, findfirst(neighbours .== α)) # remove the edge α itself!
    
    return neighbours
end

function RemoveEdges(vertices, edges, αs)
    
    for α in αs # for each edge α
        toKill = []
        for i in edges[α].∂ # for each vertex connected to α...

            deleteat!(vertices[i].δ ,findfirst(vertices[i].δ .== α)) # remove α from its coboundary
            push!(toKill, i)
        end

        edges[α].∂ = [] # set boundary of α to 0
    end
end

function LineGraph(vertices, edges)
    
    Lvertices = deepcopy(edges)
    Ledges = []
    
    α = 1
    for v in vertices
        pairs = combinations(v.δ,2)
        for pair in pairs
            push!(Ledges, Cell(false, 0, [], [pair[1], pair[2]], []))
            append!(Lvertices[pair[1]].δ, α)
            append!(Lvertices[pair[2]].δ, α)
            α += 1
        end
    end
    
    for α in eachindex(edges)
        Lvertices[α].x = 0.5 .* (vertices[edges[α].∂[1]].x + vertices[edges[α].∂[2]].x)
    end
    
    for i in eachindex(Ledges)
        Ledges[i].x = 0.5 .* (Lvertices[Ledges[i].∂[1]].x + Lvertices[Ledges[i].∂[2]].x)
    end
    
    return Lvertices, Ledges
end

# ### Interface with Graphs.jl package

# using Graphs, MetaGraphs, Plots, GraphRecipes

#function LatticeToGraph(vertices, edges)
#    # converts my custom data structure of vertices and edges to a structure matching the Graphs package
#
#    elist = []
#    for edge in edges
#        push!(elist, Tuple(edge.∂))
#    end
#
#    G = SimpleGraph(Graphs.SimpleEdge.(elist));
#
#    #G = MetaGraph(G)
#    #for i in eachindex(vertices)
#    #    set_prop!(G, i, :x, vertices[i].x)
#    #end
#    #for edge in edges
#    #    set_prop!(G, Edge(edge.∂...), :x, edge.x)
#    #    set_prop!(G, Edge(edge.∂...), :σ, edge.σ)
#    #    set_prop!(G, Edge(edge.∂...), :D, edge.D)
#    #end
#
#    return G
#end


#function GraphToLattice(G)
#    # converts a Graphs package graph to our data structure
#    
#    vertices = [Cell(false, 0, [], [], []) for i in 1:nv(G)] # list of vertices
#    edges = [Cell(false, 0, [], [], []) for α in 1:ne(G)] # list of edges
#    
#    for i in eachindex(vertices)
#        vertices[i].δ = neighbors(G, i)
#    end
#    
#    Gedges = edges(G)
#    
#    α = 1
#    for i in eachindex(vertices)
#        for j < i
#            if has_edge(i, j)
#                append!(edges[α].∂, [i, j])
#                α += 1
#            end
#        end
#    end
#    
#    return vertices, edges
#end

# Can then freely use e.g. ``graphplot(G, curves=false)``
