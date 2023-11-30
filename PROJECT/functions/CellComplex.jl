# ### A single cell

@everywhere mutable struct Cell
    x::Array{Float64} # coords
    ∂::Array{Int32} # boundary cells
    η::Array{Int8} # boundary cell orientations (+1 or -1 -- a coefficient of 0 corresponds to not being in this array!)
    ∂ᵀ::Array{Int32} # coboundary cells
    ηᵀ::Array{Int8} # coboundary cell orientations (+1 or -1 -- a coefficient of 0 corresponds to not being in this array!)
end

# ### A cell complex

@everywhere struct CellComplex
    maxdim::Int8 # max dimension
    #embDims::Array{Int} # dimensions of embedding vector spaces for Δσ.x things
    #numCells::Array{Int} # array of numbers of k-cells
    cells::Array{Array{Cell}} # cells
end

# ### A field

@everywhere mutable struct Field
    dim::Int8 # dimension
    #isAdd::Bool # determines whether -F[σ] is an add. or mult. inverse (for SI and TC operators resp)
    vals::Array{Int16} # coefficients
end


# ### Function to create a new k-field

@everywhere function CreateField(Δ, k) # group A = ℤ
    @assert k>=0 && k<=Δ.maxdim
    
    return Field(k, zeros(length(Δ.cells[k+1])))
end

# ### Function to access a component of a field

@everywhere function GetCpt(F, σ, isAdd)
    if σ > 0
        return F.vals[σ]
    else
        return isAdd ? -F.vals[-σ] : F.vals[-σ]
    end
    # return sign(σ)^isAdd * F.vals[abs(σ)] # = -F[-σ] for SI and 1/F[-σ]=F[-σ] for Z2 TC
    # can think of isAdd as representing whether F is parallel to the edges or perpendicular to them, and how that affects the component on inverting σ
end

# ### Function to set a component of a field

@everywhere function SetCpt(F, σ, a)
    F.vals[abs(σ)] = sign(σ)^isAdd * a
end


# ### Function to take the inner product of two fields

@everywhere function InnerProduct(F, G)
    @assert F.dim==G.dim
    
    FG = 0
    
    for (f, g) in zip(F, G)
        FG += f*g
    end
    
    return FG
end


# ### Function to compute the boundary of a k-field

@everywhere function Boundary(F, Δ, i) # group A = ℤ
    if F.dim==0
        return 0
    end
    
    ∂F = 0
    
    Δi = Δ.cells[(F.dim+1)-1][i]
    for (e, k) in zip(Δi.∂ᵀ, Δi.ηᵀ)
        ∂F += GetCpt(F, k*e, true)
    end
    
   return ∂F
end

@everywhere function Boundary(F, Δ) # group A = ℤ
    ∂F = CreateField(Δ, F.dim-1)
    
    for i in eachindex(Δ.cells[(F.dim+1)-1])
        ∂F.vals[i] = Boundary(F, Δ, i)
    end
    
   return ∂F
end

# ### Function to compute the coboundary of a k-field

@everywhere function Coboundary(F, Δ, p)
    if F.dim==Δ.maxdim
        return 0
    end
    
    ∂ᵀF = 0
    
    Δp = Δ.cells[(F.dim+1)+1][p]
    for (e, k) in zip(Δp.∂, Δp.η)
        ∂ᵀF += GetCpt(F, k*e, true)
    end
    
    return ∂ᵀF
end

@everywhere function Coboundary(F, Δ)
    ∂ᵀF = CreateField(Δ, F.dim+1)
    
    for p in eachindex(Δ.cells[(F.dim+1)+1])
        ∂ᵀF[p] = Coboundary(F, Δ, p)
    end
    
    return ∂ᵀF
end

# ### Function to compute the star operator of a k-field

@everywhere function Star(F, Δ, i) # group A = ℤ
    if F.dim==0
        return 0
    end
    
    AF = 1
    
    Δi = Δ.cells[(F.dim+1)-1][i]
    for (e, k) in zip(Δi.∂ᵀ, Δi.ηᵀ)
        AF *= GetCpt(F, k*e, false)
    end
    
   return AF
end

@everywhere function Star(F, Δ) # group A = ℤ
    AF = CreateField(Δ, F.dim-1)
    
    for i in eachindex(Δ.cells[(F.dim+1)-1])
        AF.vals[i] = Star(F, Δ, i)
    end
    
   return AF
end

# ### Function to compute the plaquette operator of a k-field

@everywhere function Plaquette(F, Δ, p)
    if F.dim==Δ.maxdim
        return 0
    end
    
    BF = 1
    
    Δp = Δ.cells[(F.dim+1)+1][p]
    for (e, k) in zip(Δp.∂, Δp.η)
        BF *= GetCpt(F, k*e, false)
    end
    
    return BF
end

@everywhere function Plaquette(F, Δ)
    BF = CreateField(Δ, F.dim+1)
    
    for p in eachindex(Δ.cells[(F.dim+1)+1])
        BF[p] = Plaquette(F, Δ, p)
    end
    
    return BF
end
