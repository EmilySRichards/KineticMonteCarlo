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




function KagomeBasis()
    Nv = 6
    Ne = 12
    
    Vs = [Cell(false, 0, [], [], []) for j in 1:Nv]
    Es = [Cell(false, 0, [], [], []) for α in 1:Ne]
    
    Vs[1].x = [0, 0]
    Vs[2].x = [0.5, 0]
    Vs[3].x = [0.25, 0.25]
    Vs[4].x = [0, 0.5]
    Vs[5].x = [0.5, 0.5]
    Vs[6].x = [0.75, 0.75]
    
    tmp = [(1, 1), (1, 2), (2, 2), (3, 1), (3, 3), (4, 2), (4, 3), (5, 3), (5, 4), (6, 3), (6, 5), (7, 4), (7, 5), (8, 5), (9, 5), (9, 6), (10, 6), (11, 6), (12, 6)] # list of links (edge then vertex)
    
    for t in tmp
        push!(Es[t[1]].∂, t[2])
        push!(Vs[t[2]].δ, t[1])
    end
    
    Bev = [[(2, 1), [1 0]], [(8, 4), [1 0]], [(10, 2), [0 1]], [(11, 1), [1 1]], [(12, 4), [1 0]]] # dangling edge -> matching vertex
    
    return Vs, Es, Bev
end