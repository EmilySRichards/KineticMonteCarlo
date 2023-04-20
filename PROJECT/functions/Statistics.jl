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
#     name: julia-(6-threads)-1.8
# ---

# ### Enum type for different methods
@enum Method begin
   Bootstrap = 1
   Jackknife = 2
   None = 3
end

# ### Jackknife & Bootstrap Algorithms

# https://www.physik.uni-leipzig.de/~spitzner/publications/Spitzner_bootstrap.pdf
@everywhere function MyJackknife(x, fun, W = 1)
    P = size(x, 1) # number of arguments to fun
    
    N = length(x[1])
    M = ceil(Int64, N/W) # number of blocks (ceil accounts for shorter blocks if W doesn't factor N)
    
    blocks = [[x[p][1+(m-1)*W : min(m*W,N)] for m=1:M] for p=1:P] # M x P x W nested vectors    
    samples = [[vcat(blocks[p][1:end .!= n]...) for p=1:P] for n=1:M] # M x P x W*(M-1) nested vectors
    
    fmavg = [fun(samples[n]...) for n in 1:M] # run function on all resamplings
    
    𝐟𝐦𝐚𝐯𝐠 = filter(!isnan, fmavg) # **filters out bad samples where the estimator is undefined**
    𝐌 = length(fmavg) # adjust M to account for that
    
    𝐟𝐚𝐯𝐠1 = mean(𝐟𝐦𝐚𝐯𝐠) # biased jackknife estimator
    𝐟𝐚𝐯𝐠 = (N*fun(x...) - (N-𝐌)*𝐟𝐚𝐯𝐠1)/𝐌 # UNbiased jackknife estimator
    𝐟𝐞𝐫𝐫 = sqrt(sum((𝐟𝐦𝐚𝐯𝐠 .- 𝐟𝐚𝐯𝐠1).^2)*(𝐌-1)/𝐌) # error in estimator

    return 𝐟𝐚𝐯𝐠, 𝐟𝐞𝐫𝐫
end

# https://www.physik.uni-leipzig.de/~spitzner/publications/Spitzner_bootstrap.pdf
@everywhere function MyBootstrap(x, fun, W, Nbps = -1) 
    P = size(x, 1) # number of arguments to fun
    
    N = length(x[1])
    M = ceil(Int64, N/W) # number of blocks (ceil accounts for shorter blocks if W doesn't factor N)
    if Nbps <= 0
        Nbps = M
    end
    
    blocks = [[x[p][1+(m-1)*W : min(m*W,N)] for m=1:M] for p=1:P] # M x P x W nested vectors
    rands = [rand(1:M, M) for n=1:Nbps] # have to precompute these so they're the same for all variables p
    samples = [[vcat(blocks[p][rands[n]]...) for p=1:P] for n=1:Nbps] # Nbps x P x W*M nested vectors
    
    fmavg = [fun(samples[n]...) for n in 1:Nbps] # construct Nbps randomly-chosen length-M resamplings and average each one
    
    𝐟𝐦𝐚𝐯𝐠 = filter(!isnan, fmavg) # **filters out bad samples where the estimator is undefined**
    𝐍𝐛𝐩𝐬 = length(fmavg) # adjust Nbps to account for that
    
    𝐟𝐚𝐯𝐠 = mean(𝐟𝐦𝐚𝐯𝐠) # biased estimator (bias hard to estimate and usually small => typically ignored)
    𝐟𝐞𝐫𝐫 = sqrt(sum((𝐟𝐦𝐚𝐯𝐠 .- 𝐟𝐚𝐯𝐠).^2)/(𝐍𝐛𝐩𝐬-1)) # error in estimator

    return 𝐟𝐚𝐯𝐠, 𝐟𝐞𝐫𝐫
end

# ### Add capability to override Bootstrap because it fails for heat bath method!!!

@everywhere function Estimator(method, x, fun, W=1, Nbps=-1)
    if method==Bootstrap
        return MyBootstrap(x, fun, W, Nbps)
    elseif method==Jackknife
        return MyJackknife(x, fun, W)
    else
        return fun(x...), 0
    end
end


# ### Autocorrelation function

@everywhere function MyAutocor(y, normalise=true) # shamelessly stolen from the StatsBase package - had to b/c of package issues on the TCM network
    D = size(y, 1)
    T = size(y, 2)
    
    lags = range(0,T-1)
    
    r = zeros(D, length(lags))
    for d in 1:D
        for (t, lag) in enumerate(lags)  # for each lag value
            r[d,t] = sum(y[d,1:T-lag] .* y[d,1+lag:T])
            r[d,t] /= T - lag
        end
    end
    
    if normalise
        r ./= r[:,1]
    end
    
    return r
end

# ### Integrated autocorrelation time

@everywhere function IntAutocorrTime(x)
    y = MyAutocor(hcat(x...)')
    y = y[:,2:end] # cut out τ=0 term
    y .*= ones(size(y)) .- repeat(collect(range(1,size(y,2))), 1, size(y,1))'./size(y,2)
    return ceil(Int64, 1 + 2 * maximum(sum(y, dims=2))) # exclude τ=0 term
end

# ### Useful alternative functions

@everywhere Var(x) = (length(x)>1) ? var(x) : 0.0
