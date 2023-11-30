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

    return  Z
end


# ### <A> Single Vertex

function Asv(T, 𝒽, z)
    
    A = zeros(size(T))
    for n in 0:z
        A += (-1)^n .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end

    A ./= PartitionFunction(T, 𝒽, z)
    
    return  A
end


# ### <B> Single Vertex

function Q²sv(T, 𝒽, z)
    
    Q² = zeros(size(T))
    for n in 0:z
        Q² += (z-2*n).^2 .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end

    Q² ./= PartitionFunction(T, 𝒽, z)
    
    return Q²
end


# ### <Q> Single Vertex

function Qsv(T, 𝒽, z)
    
    Q = zeros(size(T))
    for n in 0:z
        Q += (z-2*n) .* binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
    end
    
    Q ./= PartitionFunction(T, 𝒽, z)
    
    return  Q
end


# ### Min-Energy Excitation Denstity

function ExcitationDensity(T, 𝒽, z)
    if isSpinIce # spin ice case
        q = (mod(z, 2)==0) ? 2 : 3 # lowest-energy excitation charge above GS
        
        Nq = zeros(size(T))
        for n in 0:z
            if abs(z-2*n) == q # exclude GS states
                Nq += binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
            end
        end
        
        #Nq = binomial(z, (z-q)÷2) * exp.(- (λ ./ T) - q^2 .* (ξ ./ T)) .* 2 .* cosh.(q .* (𝒽 ./ T))
        
        
        Nq ./= PartitionFunction(T, 𝒽, z)
            
        return Nq
    else
        return  0.5 .* (1 .- Asv(T, 𝒽, z)) # toric code case - easy!
    end
end


# ### All-Energy Excitation Denstity

function AllExcitationDensity(T, 𝒽, z)
    if isSpinIce # spin ice case
        q = (mod(z, 2)==0) ? 0 : 1 # GS charge
        
        Nq = zeros(size(T))
        for n in 0:z
            if abs(z-2*n) != q # exclude GS states
                Nq += binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
            end
        end
        
        Nq ./= PartitionFunction(T, 𝒽, z)
            
        return Nq
    else
        return  0.5 .* (1 .- Asv(T, 𝒽, z)) # toric code case - easy!
    end
end


# ### Magnetisation

function Magnetisation(T, 𝒽, z)
    return  Qsv(T, 𝒽, z) ./ z
end


# ### Heat Capacity

function HeatCapacity(T, 𝒽, z)
    
    Zfun = (β) -> PartitionFunction([1/β], 𝒽, z)[1]
    Z1fun = (β) -> ForwardDiff.derivative(Zfun, β)
    Z2fun = (β) -> ForwardDiff.derivative(Z1fun, β)
    
    C = zeros(length(T))
    for n in eachindex(T)
        C[n]= Z2fun(1/T[n]) / Zfun(1/T[n]) - (Z1fun(1/T[n]) / Zfun(1/T[n])) ^ 2
    end
    C ./= T.^2
    C .*= 2 / z # want capacity per SPIN, not per VERTEX - should have z/2 = |E|/|V|
    
    return  C
end



# ### Heat Capacity Test

function HeatCapacityTest(T, 𝒽, z)
    
    if isSpinIce
        function tmp(T, h, z)
            q = (mod(z, 2)==0) ? 2 : 3 # GS charge

            Nq = zeros(size(T))
            for n in 0:z
                if abs(z-2*n) <= q # exclude extreme excited states
                    Nq += binomial(z, n) .* exp.((-1)^n .* (λ ./ T) - (z-2*n)^2 .* (ξ ./ T) + (z-2*n) .* (𝒽 ./ T))
                end
            end

            return Nq
        end

        Zfun = (β) -> tmp([1/β], 𝒽, z)[1]
        Z1fun = (β) -> ForwardDiff.derivative(Zfun, β)
        Z2fun = (β) -> ForwardDiff.derivative(Z1fun, β)

        C = zeros(length(T))
        for n in eachindex(T)
            C[n]= Z2fun(1/T[n]) / Zfun(1/T[n]) - (Z1fun(1/T[n]) / Zfun(1/T[n])) ^ 2
        end
        C ./= T.^2
        C .*= 2 / z # want capacity per SPIN, not per VERTEX - should have z/2 = |E|/|V|
        
        return  C
        
    else
        return HeatCapacity(T, 𝒽, z)
    end
end