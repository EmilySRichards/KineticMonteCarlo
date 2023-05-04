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

# ### Partition Function

function PartitionFunction(T, 𝒽, z)
    
    Z = zeros(size(T))
    for n in 0:z
        Z += binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end
        
    #Z  = 6 .* exp( λ ./ T)
    #Z += 2 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh.(4 .*  𝒽 ./ T)
    #Z += 8 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh.(2 .*  𝒽 ./ T)
    
    return  Z
end


# ### <A> Single Vertex

function Asv(T, 𝒽, z)
    
    A = zeros(size(T))
    for n in 0:z
        A += (-1)^n .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end
    
    #A  = 6 .* exp.( λ ./ T)
    #A += 2 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh(4 .*  𝒽 ./ T)
    #A -= 8 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh(2 .*  𝒽 ./ T)
    
    A /= PartitionFunction(T, 𝒽, z)
    
    return  A
end


# ### <B> Single Vertex

function Bsv(T, 𝒽, z)
    
    B = zeros(size(T))
    for n in 0:z
        B += - (z-2*n).^2 .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end
    
    #B  = 32 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh.(4 .* 𝒽 ./ T)
    #B += 32 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh.(2 .* 𝒽 ./ T)
    
    B /= PartitionFunction(T, 𝒽, z)
    
    return  B
end


# ### <Q> Single Vertex

function Qsv(T, 𝒽, z)
    
    Q = zeros(size(T))
    for n in 0:z
        Q += (z-2*n) .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end
    
    #Q  =  8 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* sinh.(4 .*  𝒽 ./ T)
    #Q -= 16 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* sinh.(2 .*  𝒽 ./ T)
    
    Q /= PartitionFunction.(T, 𝒽, z)
    
    return  Q
end


# ### Min-Energy Excitation Denstity

function ExcitationDensity(T, 𝒽, z)
    if λ==0 # spin ice case
        q = (mod(z, 2)==0) ? 2 : 1 # lowes-energy excitation charge
        
        Nq = zeros(size(T))
        
        Nq = binomial(z, (z-q)/2) * exp.(- (λ ./ T) - q^2 .* (ξ ./ T)) .* 2 .* cosh.(q .* (𝒽 ./ T))
        
        #for n in ns
        #    Nq += binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
        #end
          
        Nq /= PartitionFunction.(T, 𝒽, z)
            
        return Nq
    end
    
    return  0.5 .* (1 .- Asv(T, 𝒽, z)) # toric code case - easy!
end


# ### Magnetisation

function Magnetisation(T, 𝒽, z)
    return  Qsv(T, 𝒽, z) ./ z
end


# ### Heat Capacity

function HeatCapacity(T, 𝒽, z)
    
    Zfun = (β) -> PartitionFunction(1/β, 𝒽, z)
    #Zfun = (β) -> 6*exp(λ*β) + 2*exp(λ*β)*exp.(-16*ξ*β) * cosh(4*𝒽*β) + 8*exp(-λ*β)*exp(-4*ξ*β)*cosh(2*𝒽*β)
    Z1fun = (β) -> ForwardDiff.derivative(Zfun, β)
    Z2fun = (β) -> ForwardDiff.derivative(Z1fun, β)
    
    C = zeros(length(T))
    for n in eachindex(T)
        C[n]= Z2fun(1/T[n]) / Zfun(1/T[n]) - (Z1fun(1/T[n]) / Zfun(1/T[n])) ^ 2
    end
    C ./= 2 * T.^2
    
    return  C
end