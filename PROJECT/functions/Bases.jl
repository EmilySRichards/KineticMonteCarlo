# ### Coordination from basis
@everywhere function Coordination(Basis)
    Vs  = Basis[1]
    Es  = Basis[2][1]
    
    count = zeros(Int64, length(Vs))
    for edge in Es
        
        for vertex in edge
            count[vertex[1]] += 1
        end
    end
    
    @assert all(count .== count[1]) # check that coordination is same for all all sites otherwise ill-defined
    
    return count[1]
end

# ### Different Unit Cells

@everywhere function CubicBasis(dim)
    n = [binomial(dim, 1) for d in 1:2] # ONLY UP TO EDGES FOR NOW!!!!
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    Verts[1].x = zeros(dim)

    Links = [] # all >0-dim cells
    
    # edges
    Edges = []
    for d in 1:dim
        dir = zeros(dim)
        dir[d] = 1
        
        push!(Edges, [(1, zeros(dim)), (1, dir)])
    end
    push!(Links, Edges)
    
    Scale = ones(dim) # scale of the unit cell dimensions
    
    return Verts, Links, n, Scale
end



@everywhere function SquareBasis()
    n = [1, 2, 1]
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    Verts[1].x = zeros(2)

    Links = [] # all >0-dim cells
    
    # edges
    push!(Links, [[(1, [0 0]), (1, [1 0])],
                  [(1, [0 0]), (1, [0 1])]])
    
    # faces
    push!(Links, [[(1, [0 0]), (2, [1 0]), (1, [0 1]), (2, [0 0])]])
    
    
    Scale = ones(2) # scale of the unit cell dimensions
    
    return Verts, Links, n, Scale
end



function HexBasis()
    # number of vertices, edges, faces, ...
    n = [4, 6, 2]
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    
    Verts[1].x = [0  , 0  ]
    Verts[2].x = [1/6, 1/2]
    Verts[3].x = [1/2, 1/2]
    Verts[4].x = [2/3, 0  ]
    
    Links = [] # all >0-dim cells
    
    # edges
    push!(Links, [[(1, [0 0]), (2, [0 0])],
                  [(2, [0 0]), (3, [0 0])],
                  [(3, [0 0]), (4, [0 0])],
                  [(4, [0 0]), (1, [1 0])],
                  [(3, [0 0]), (4, [0 1])],
                  [(2, [0 0]), (1, [0 1])]])
    
    # faces
    push!(Links, [[(6, [0 0]), (2, [0 0]), (5, [0 0]), (3, [0 1]), (2, [0 1]), (1, [0 1])],
                  [(3, [0 0]), (4, [0 0]), (1, [1 0]), (6, [1 0]), (4, [0 1]), (5, [0 0])]])
    
    Scale = [3, sqrt(3)] # scale of the unit cell dimensions (such that bond length = 1)
    
    return Verts, Links, n, Scale
end




@everywhere function DiamondBasis()
    # number of vertices, edges, faces, ...
    n = [8, 16]
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    
    Verts[1].x = [0,   0,   0  ]
    Verts[2].x = [1/4, 1/4, 1/4]
    Verts[3].x = [1/2, 1/2, 0  ]
    Verts[4].x = [1/2, 0  , 1/2]
    Verts[5].x = [0,   1/2, 1/2]
    Verts[6].x = [3/4, 3/4, 1/4]
    Verts[7].x = [3/4, 1/4, 3/4]
    Verts[8].x = [1/4, 3/4, 3/4]
    
    Links = [] # all >0-dim cells
    
    # edges
    push!(Links, [[(1, [0 0 0]), (2, [0 0 0])],
                  [(2, [0 0 0]), (3, [0 0 0])],
                  [(2, [0 0 0]), (4, [0 0 0])],
                  [(2, [0 0 0]), (5, [0 0 0])],
                  [(3, [0 0 0]), (6, [0 0 0])],
                  [(4, [0 0 0]), (7, [0 0 0])],
                  [(5, [0 0 0]), (8, [0 0 0])],
                  
                  [(6, [0 0 0]), (1, [1 1 0])],
                  [(6, [0 0 0]), (5, [1 0 0])],
                  [(6, [0 0 0]), (4, [0 1 0])],
            
                  [(7, [0 0 0]), (5, [1 0 0])],
                  [(7, [0 0 0]), (1, [1 0 1])],
                  [(7, [0 0 0]), (3, [0 0 1])],
            
                  [(8, [0 0 0]), (4, [0 1 0])],
                  [(8, [0 0 0]), (3, [0 0 1])],
                  [(8, [0 0 0]), (1, [0 1 1])]])
    
    
    Scale = ones(3) .* 4/sqrt(3) # scale of the unit cell dimensions (such that bond lengths=1)
    
    return Verts, Links, n, Scale
end



function SemiTriangBasis()
    # number of vertices, edges, faces, ...
    n = [2, 5]
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    
    Verts[1].x = [0  , 0]
    Verts[2].x = [1/2, 0]
    
    Links = [] # all >0-dim cells
    
    # edges
    push!(Links, [[(1, [0 0]), (2, [0 0])],
                  [(1, [0 0]), (1, [0 1])],
                  [(1, [0 0]), (2, [0 1])],
                  [(2, [0 0]), (2, [0 1])],
                  [(2, [0 0]), (1, [1 0])]])

    Scale = [3, sqrt(3)] # scale of the unit cell dimensions (such that bond length = 1)
    
    return Verts, Links, n, Scale
end


function KagomeBasis()
    # number of vertices, edges, faces, ...
    n = [6, 12]
    
    Verts = [Cell(false, 0, [], [], []) for j in 1:n[1]]
    
    Verts[1].x = [0, 0]
    Verts[2].x = [0.5, 0]
    Verts[3].x = [0.25, 0.25]
    Verts[4].x = [0, 0.5]
    Verts[5].x = [0.5, 0.5]
    Verts[6].x = [0.75, 0.75]
    
    Links = [] # all >0-dim cells
    
    # edges
    push!(Links, [[(1, [0 0]), (2, [0 0])],
                  [(2, [0 0]), (1, [1 0])],
                  [(1, [0 0]), (3, [0 0])],
                  [(2, [0 0]), (3, [0 0])],
                  [(3, [0 0]), (4, [0 0])],
                  [(3, [0 0]), (5, [0 0])],
                  [(4, [0 0]), (5, [0 0])],
                  [(5, [0 0]), (4, [1 0])],
                  [(5, [0 0]), (6, [0 0])],
                  [(6, [0 0]), (2, [0 1])],
                  [(6, [0 0]), (1, [1 1])],
                  [(6, [0 0]), (4, [1 0])]])
    
    Scale = [2, 2*sqrt(3)] # scale of the unit cell dimensions (such that bond lengths=1)
    
    return Verts, Links, n, Scale
end