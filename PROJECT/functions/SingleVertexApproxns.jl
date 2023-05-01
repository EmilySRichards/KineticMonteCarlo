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

function PartitionFunction(T, 𝒽)
    Z  = 6 .* exp( λ ./ T)
    Z += 2 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh.(4 .*  𝒽 ./ T)
    Z += 8 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh.(2 .*  𝒽 ./ T)
    
    return  Z
end


# ### <A> Single Vertex

function Asv(T, 𝒽)
    A  = 6 .* exp.( λ ./ T)
    A += 2 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh(4 .*  𝒽 ./ T)
    A -= 8 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh(2 .*  𝒽 ./ T)
    
    A /= PartitionFunction(T, 𝒽)
    
    return  A
end


# ### <B> Single Vertex

function Bsv(T, 𝒽)
    B  = 32 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh.(4 .* 𝒽 ./ T)
    B -= 32 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh.(2 .* 𝒽 ./ T)
    
    B /= PartitionFunction(T, 𝒽)
    
    return  B
end


# ### A=-1 Excitation Denstity

function ExcitationDensity(T, 𝒽)
    A   = 6 .* exp.( λ ./ T)
    A .+= 2 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* cosh.(4 .*  𝒽 ./ T)
    A .-= 8 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* cosh.(2 .*  𝒽 ./ T)
    
    A ./= PartitionFunction.(T, 𝒽)
    
    return  0.5 .* (1 .- A)
end


# ### Magnetisation

function Magnetisation(T, 𝒽)
    M   =  8 .* exp.( λ ./ T) .* exp.(-16 .* ξ ./ T) .* sinh.(4 .*  𝒽 ./ T)
    M .-= 16 .* exp.(-λ ./ T) .* exp.(- 4 .* ξ ./ T) .* sinh.(2 .*  𝒽 ./ T)
    
    M ./= 4 .* PartitionFunction.(T, 𝒽)
    
    return  M
end


# ### Heat Capacity

function HeatCapacity(T, 𝒽)
    
    Zfun = (β) -> 6*exp(λ*β) + 2*exp(λ*β)*exp.(-16*ξ*β) * cosh(4*𝒽*β) + 8*exp(-λ*β)*exp(-4*ξ*β)*cosh(2*𝒽*β)
    Z1fun = (β) -> ForwardDiff.derivative(Zfun, β)
    Z2fun = (β) -> ForwardDiff.derivative(Z1fun, β)
    
    C = zeros(length(T))
    for n in eachindex(T)
        C[n]= Z2fun(1/T[n]) / Zfun(1/T[n]) - (Z1fun(1/T[n]) / Zfun(1/T[n])) ^ 2
    end
    C ./= 2 * T.^2
    
    return  C
end