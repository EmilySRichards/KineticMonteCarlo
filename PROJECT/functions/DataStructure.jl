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

@everywhere function CreateCell(Vs0, Es0, I, X)
    Vs = deepcopy(Vs0)
    Es = deepcopy(Es0)
    
    for v in Vs
        v.x += X
        v.δ .+= length(Es)*(I-1)
    end
    for e in Es
        e.∂ .+= length(Vs)*(I-1)
    end
    
    return Vs, Es
end

@everywhere function LatticeGrid(L, PBC, Basis)
    
    Vs0, Es0, Bonds, Scale = Basis
    nv = length(Vs0)
    ne = length(Es0)
    
    @assert length(L)==length(PBC)
    dim = length(L)
    
    N = prod(L) # number of unit cells
    
    vertices = [Cell(false, 0, [], [], []) for j in 1:nv*N] # list of vertices
    edges = [Cell(false, 0, [], [], []) for α in 1:ne*N] # list of edges
    
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

    
    # place down all unit cells without connecting dangling edges
    x0 = ones(length(L)) # origin for 1-indexing positions as done in X (just for convenience because Julia)
    for I in 1:N
        X = I_to_X(I, L) - x0 # don't forget to subtract the origin!
        Vs, Es = CreateCell(Vs0, Es0, I, X)
        
        vertices[nv*(I-1)+1:nv*I] .= Vs
        edges[ne*(I-1)+1:ne*I] .= Es
    end
        
    # go back through and link up dangling edges vertices and edges arrays
    for I in 1:N
        X = I_to_X(I, L)
        
        # attach edges (in +ve quadrant)
        for Bond in Bonds
            dir = Bond[2] # relative displacement of linked cell, each a vector ∈ {0,±1}
            ev = Bond[1] # tuple of the edge and the vertex to be linked (within a unit cell)
            
            ifLink = true
            for d in 1:dim # if at a boundary and OBCs then don't link up the edge
                if dir[d]<0 && X[d]==1 && !PBC[d]
                    ifLink = false
                elseif dir[d]>0 && X[d]==L[d] && !PBC[d]
                    ifLink = false
                end
            end
            
            if ifLink
                Y = copy(X)
                for d in 1:dim
                    Y[d] = (dir[d]<0 && X[d]==1) ? L[d] : ((dir[d]>0 && X[d]==L[d]) ? 1 : X[d]+dir[d])
                end
                J = X_to_I(Y, L)
                
                push!(edges[ev[1]+ne*(I-1)].∂, ev[2]+nv*(J-1))
                push!(vertices[ev[2]+nv*(J-1)].δ, ev[1]+ne*(I-1))
            end
        end
    end
    
    # if there are any remaining dangling edges, they must be on an OBC => kill them
    if !all(PBC)
        toKill = []
        shifts = zeros(length(edges))
        for α in eachindex(edges)
            if length(edges[α].∂) < 2
                push!(toKill, α)
                shifts[α:end] .-= 1
            end
        end
        
        deleteat!(edges, toKill)

        # fix the edges in the coboundary of each vertex
        for v in vertices
            toKillv = []
            for i in eachindex(v.δ)
                if v.δ[i] in toKill
                    push!(toKillv, i)
                end
                v.δ[i] += shifts[v.δ[i]]
            end
            deleteat!(v.δ, toKillv)
        end
    end
    
    # rescale vertex positions
    for v in vertices
        v.x = v.x .* Scale
    end
    
    # calculate edge positions from vertex positions
    for e in edges
        e.x = zeros(length(L))
        for i in e.∂
            e.x += vertices[i].x
        end
        e.x ./= length(e.∂)
    end
    
    return vertices, edges, Scale
end

# ### Different Unit Cells

@everywhere function CubicBasis(dim)
    Nv = 1
    Ne = dim
    
    Vs = [Cell(false, 0, [], [], []) for j in 1:Nv]
    Es = [Cell(false, 0, [], [], []) for α in 1:Ne]
    
    Vs[1].x = zeros(dim)
    
    Bev = []
    for d in 1:dim
        push!(Es[d].∂, 1)
        push!(Vs[1].δ, d)
        dir = zeros(dim)
        dir[d] = 1
        push!(Bev, [(d, 1), dir])
    end
    
    Scale = ones(dim) # scale of the unit cell dimensions
    
    return Vs, Es, Bev, Scale
end

function HexBasis()
    Nv = 4
    Ne = 6
    
    Vs = [Cell(false, 0, [], [], []) for j in 1:Nv]
    Es = [Cell(false, 0, [], [], []) for α in 1:Ne]
    
    Vs[1].x = [0  , 0  ]
    Vs[2].x = [1/6, 1/2]
    Vs[3].x = [1/2, 1/2]
    Vs[4].x = [2/3, 0  ]
    
    tmp = [(1, 1), (1, 2), (2, 2), (3, 2), (3, 3), (4, 3), (5, 3), (5, 4), (6, 4)] # list of links (edge then vertex)
    
    for t in tmp
        push!(Es[t[1]].∂, t[2])
        push!(Vs[t[2]].δ, t[1])
    end
    
    Bev = [[(6, 1), [1 0]], [(2, 1), [0 1]], [(4, 4), [0 1]]] # dangling edge -> vertex bonds
    
    Scale = [3, sqrt(3)] # scale of the unit cell dimensions (such that bond length = 1)
    
    return Vs, Es, Bev, Scale
end

@everywhere function DiamondBasis()
    Nv = 8
    Ne = 16
    
    Vs = [Cell(false, 0, [], [], []) for j in 1:Nv]
    Es = [Cell(false, 0, [], [], []) for α in 1:Ne]
    
    Vs[1].x = [0,   0,   0  ]
    Vs[2].x = [1/4, 1/4, 1/4]
    Vs[3].x = [1/2, 1/2, 0  ]
    Vs[4].x = [1/2, 0  , 1/2]
    Vs[5].x = [0,   1/2, 1/2]
    Vs[6].x = [3/4, 3/4, 1/4]
    Vs[7].x = [3/4, 1/4, 3/4]
    Vs[8].x = [1/4, 3/4, 3/4]
    
    tmp = [(1, 1), (1, 2), (2, 2), (2, 3), (3, 2), (3, 4), (4, 2), (4, 5), (5, 3), (5, 6), (6, 4), (6, 7), (7, 5), (7, 8)] # list of INTERNAL links (edge then vertex)
    append!(tmp, [(8, 6), (9, 6), (10, 6), (11, 7), (12, 7), (13, 7), (14, 8), (15, 8), (16, 8)]) # list of EXTERNAL links
    
    for t in tmp
        push!(Es[t[1]].∂, t[2])
        push!(Vs[t[2]].δ, t[1])
    end
    
    Bev = [[(8, 1), [1 1 0]], [(9, 5), [1 0 0]], [(10, 4), [0 1 0]], [(11, 5), [1 0 0]], [(12, 1), [1 0 1]], [(13, 3), [0 0 1]], [(14, 4), [0 1 0]], [(15, 3), [0 0 1]], [(16, 1), [0 1 1]]] # dangling edge -> matching vertex
    
    Scale = ones(3) .* 4/sqrt(3) # scale of the unit cell dimensions (such that bond lengths=1)
    
    return Vs, Es, Bev, Scale
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

# ### Handy Graph Tools

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
