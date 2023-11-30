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

# ### Set up the geometry

@everywhere function MicroDiffnSetup(Δ, numToFlip)
    
    # initialise entire system in ground state
    S = CreateField(Δ, 1)
    GroundState!(S, Δ)
    
    # flip numEdges random spins
    valid_edges = collect(eachindex(Δ.cells[2]))
    flipped_edges = []
    
    for n in 1:numToFlip
        if valid_edges==[]
            break
        end
        
        e = rand(valid_edges)
        
        push!(flipped_edges, e) # add edge to list to flip
        
        # remove e AND other edges which share vertices with it
        deleteat!(valid_edges, findall(valid_edges.==e))
        for i in Δ.cells[2][e].∂
            for f in Δ.cells[1][i].∂ᵀ
                deleteat!(valid_edges, findall(valid_edges.==f))
            end                
        end
    end

    for e in flipped_edges
        S.vals[e] = -S.vals[e]
    end
    
    return S
end

# ### Single spin-flip dynamics routine 

@everywhere function MicroDiffn(S, Δ, runtime, 𝒽, allowBckgdMoves)
    
    dim = length(Δ.cells[1][1].x)
    
    # find all the excitations
    excitations = []
    for i in eachindex(Δ.cells[1])
        A = Star(S, Δ, i)
        Q = abs(-Boundary(S, Δ, i))
        
        if (isSpinIce ? (Q == 3 || Q == 2) : A == -1)
            push!(excitations, i)
        end
    end

    
    δs = zeros(dim, length(excitations), runtime)
    xs = zeros(dim, length(excitations), runtime+1)
    for (n, i) in enumerate(excitations)
        xs[:,n,1] = Δ.cells[1][i].x
    end
    
    # actual simulation
    for t in 1:runtime
        xs[:,:,t+1] = xs[:,:,t]
        
        for _ in Δ.cells[2]
            e = rand(eachindex(Δ.cells[2]))
            ΔE, _ = EnergyChange(S, Δ, [e], 𝒽)
            
            Δe = Δ.cells[2][e]
            i = Δe.∂[1]
            j = Δe.∂[2]

            if ΔE == 0 && (allowBckgdMoves || (i in excitations || j in excitations))
                
                S.vals[e] = -S.vals[e]
                
                # displacement of edge (fixed to account for PBCs) from i->j
                δ = Displacement(Δe, i, j)
                    
                # if a prtcl is linked to this edge, move it - note ΔE=/=0 if prtcl on both vertices => ignore this case
                n = findfirst(excitations.==i)
                m = findfirst(excitations.==j)
                if n!=nothing     # i = excitations[n] = excitation
                    excitations[n] = j
                    xs[:,n,t+1] += δ # = Δ.cells[1][j].x # 
                    δs[:,n,t] += δ
                elseif m!=nothing # j = excitations[m] = excitation
                    excitations[m] = i
                    xs[:,m,t+1] -= δ # = Δ.cells[1][i].x # 
                    δs[:,m,t] -= δ
                end
                
                # else no prtcls => nothing to track!
            end
        end
    end
    
    return xs, δs
end

# ### Single spin-flip dynamics routine WITH PLAQUETTE FLIPS

@everywhere function MicroDiffn_plaqs(S, Δ, runtime, 𝒽, allowBckgdMoves, N)
    
    dim = length(Δ.cells[1][1].x)
    
    # find all the excitations
    excitations = []
    for i in eachindex(Δ.cells[1])
        A = Star(S, Δ, i)
        Q = abs(-Boundary(S, Δ, i))
        
        if (isSpinIce ? (Q == 3 || Q == 2) : A == -1)
            push!(excitations, i)
        end
    end
    
    δs = zeros(dim, length(excitations), runtime)
    xs = zeros(dim, length(excitations), runtime+1)
    for (n, i) in enumerate(excitations)
        xs[:,n,1] = Δ.cells[1][i].x
    end
    
    
    # actual simulation
    for t in 1:runtime
        xs[:,:,t+1] = xs[:,:,t]
        
        for n in eachindex(Δ.cells[2])
            
            # try to flip 1 edge
            for n in 1:N
                e = rand(eachindex(Δ.cells[2])) # pick a random edge
                ΔE, _ = EnergyChange(S, Δ, [e], 𝒽)
                
                Δe = Δ.cells[2][e]
                i = Δe.∂[1]
                j = Δe.∂[2]
                
                if ΔE == 0 && (allowBckgdMoves || (i in excitations || j in excitations))
                
                    S.vals[e] = -S.vals[e]

                     # displacement of edge (fixed to account for PBCs) from i->j
                    δ = Displacement(Δe, i, j)

                    # if a prtcl is linked to this edge, move it - note ΔE=/=0 if prtcl on both vertices => ignore this case
                    n = findfirst(excitations.==i)
                    m = findfirst(excitations.==j)
                    if n!=nothing     # i = excitations[n] = excitation
                        excitations[n] = j
                        xs[:,n,t+1] += δ # = Δ.cells[1][j].x # 
                        δs[:,n,t] += δ
                    elseif m!=nothing # j = excitations[m] = excitation
                        excitations[m] = i
                        xs[:,m,t+1] -= δ # j = Δ.cells[1][i].x # 
                        δs[:,m,t] -= δ
                    end
                end
            end
            
            # try to flip a plaquette every N moves
            #if false
            #for _ in 1:1
            if mod(n, N) == 0
                p = rand(eachindex(Δ.cells[3])) # pick a random plaquette
                edges = Δ.cells[3][p].∂ # all boundary edges to try and flip
                ΔE, _ = EnergyChange(S, Δ, edges, 𝒽)
                
                if ΔE == 0
                    for e in edges
                        S.vals[e] = -S.vals[e]
                    end
                end
                # either way, no particles move => no need to update any of xs, excitations or δs
            end
        end
    end
    
    return xs, δs
end

# ### Double spin-flip dynamics routine 

@everywhere function MicroDiffn_2flip(S, Δ, runtime, 𝒽)
    
    dim = length(Δ.cells[1][1].x)
    
    # find all the excitations
    excitations = []
    for i in eachindex(Δ.cells[1])
        A = Star(S, Δ, i)
        Q = abs(-Boundary(S, Δ, i))
        
        if (isSpinIce ? (Q == 3 || Q == 2) : A == -1)
            push!(excitations, i)
        end
    end
    
    δs = zeros(dim, length(excitations), runtime)
    xs = zeros(dim, length(excitations), runtime+1)
    for (n, i) in enumerate(excitations)
        xs[:,n,1] = Δ.cells[1][i].x
    end
    
    
    # actual simulation
    for t in 1:runtime
        xs[:,:,t+1] = xs[:,:,t]
        
        for _ in Δ.cells[1]
            # propose flips
            j = rand(eachindex(Δ.cells[1])) # random vertex
            edges = sample(Δ.cells[1][j].∂ᵀ, 2; replace=false) # two neighbouring edges
            ΔE, _ = EnergyChange(S, Δ, edges, 𝒽)
            
            
            Δe = Δ.cells[2][edges[1]] # first edge
            Δf = Δ.cells[2][edges[2]] # second edge
            
            i = Δe.∂[findfirst(Δe.∂ .!= j)] # outer vertex on e
            k = Δf.∂[findfirst(Δf.∂ .!= j)] # outer vertex on f
            
            ΣA = Star(S, Δ, i) + Star(S, Δ, j) + Star(S, Δ, k) # corrects for overcounting of Ai
            
            # decide whether to accept and perform the move
            if ΔE == 0 && (S.vals[edges[1]] + S.vals[edges[2]] == 0) # energy AND magnetisation conserved
            # && ΣA>0   NO hole diffusion moves
            # && ΣA<0 ONLY hole diffusion moves
                
                for e in edges
                    S.vals[e] = -S.vals[e]
                end
                
                # ΔE=0 => an excitation is linked to this edge => move it
                # we choose to assume the central particle is fixed => valid way of tracking them
                
                # vector from site i to site k (via j)
                δik = Displacement(Δe, i, j) + Displacement(Δf, j, k)
                
                # if a prtcl is linked to this edge, move it - note ΔE=/=0 if prtcl on both vertices => ignore this case
                n = findfirst(excitations.==i)
                m = findfirst(excitations.==k)
                if n!=nothing     # j = excitations[n] = excitation
                    excitations[n] = j
                    xs[:,n,t+1] += δik # = Δ.cells[1][k].x # 
                    δs[:,n,t] += δik
                elseif m!=nothing # k = excitations[m] = excitation
                    excitations[m] = k
                    xs[:,m,t+1] -= δik # = Δ.cells[1][j].x # 
                    δs[:,m,t] -= δik
                end
            end
        end
    end
    
    return xs, δs
end

# ### Observables

# #### Mean-Squared displacement functions

# https://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft
@everywhere function Msd_fft(x)
    # splits up the MSD calculation to allow for fft optimisations - fft not used yet though...
    
    T = length(x) # number of timesteps
    
    # calculate S1
    D = x.^2
    S1 = zeros(Float64, T)
    Q = 2*sum(D)
    for t in 0:T-1
        D1 = (t==0) ? 0 : D[t]
        D2 = (t==0) ? 0 : D[T-t+1]
        Q -= D1 + D2
        S1[t+1] = Q / (T-t)
    end
    
    # calculate S2
    S2 = MyAutocor(reshape(x, (1, T)), false)
    S2 = dropdims(S2, dims=1)
    
    return S1 .- 2*S2
end

@everywhere function Msd_ez(x)
    lags = range(0, length(x)-1)
    msd = zeros(length(x))    
    
    for (i, lag) in enumerate(lags)
        diffs = x[1:end-lag] .- x[lag+1:end]
        msd[i] += mean(diffs.^2)
    end
    
    return msd
end

@everywhere function Msd(x)
    D = size(x, 1)
    P = size(x, 2)
    T = size(x, 3)
    
    msd = zeros(T)
    for d in 1:D
        for p in 1:P
            msd .+= Msd_fft(x[d,p,:])
        end
    end
    msd ./= P
    
    return msd
end

# #### Step direction correlation functions

@everywhere function DirrCorr(dx)
    D = size(dx, 1)
    P = size(dx, 2)
    T = size(dx, 3)
    
    corr = zeros(T)
    for d in 1:D # sum over dimensions
        corr .+= dropdims(mean(MyAutocor(dx[d,:,:], false), dims=1), dims=1) # average over particles
    end
    #corr ./= corr[1] # normalise
    
    return corr
end

# ### Single diffusion routine

@everywhere function DiffSimSingle(Δ, therm_runtime, runtime, useT, ℓorT, 𝒽)
    
    # thermalise to correct temperature OR correct number of particles
    if useT
        S, _ = MicroKuboSetup(Δ, therm_runtime, ℓorT, 𝒽, true)
    else
        S = MicroDiffnSetup(Δ, ℓorT)
        
        if twoFlip # allow particles to separate before we start tracking them!
            MicroDiffn_2flip(S, Δ, therm_runtime, 𝒽)
        else
            MicroDiffn(S, Δ, therm_runtime, 𝒽, true)
            #MicroDiffn_plaqs(S, Δ, therm_runtime, 𝒽, true, 1)
        end
    end
    
    M = mean(S.vals)
    
    bondNumber = 0
    maxBondNumber = 0
    for (i, Δi) in enumerate(Δ.cells[1])
        z = length(Δi.∂ᵀ)
        z₋ = 0
        for e in Δi.∂ᵀ
            z₋ += S.vals[e]<0 ? 1 : 0
        end
        z₊ = z - z₋
        
        bondNumber += z₊*z₋/2
        maxBondNumber += z*(z-1)/2
    end
    ℙ = bondNumber/maxBondNumber
    
    # track the paths of the resulting excitations (can't annihiliate b/c microcanonical!)
    if twoFlip
        x, δ = MicroDiffn_2flip(S, Δ, runtime, 𝒽)
    else
        x, δ = MicroDiffn(S, Δ, runtime, 𝒽, true)
        #x, δ = MicroDiffn_plaqs(S, Δ, runtime, 𝒽, true, 1)
    end
    
    return x, δ, M, ℙ
end

# ### Overall diffusion routine

@everywhere function mpfun1(args)
    return DiffSimSingle(args...)
end

@everywhere function DiffSim(L, PBC, Basis, therm_runtime, runtime, ℓ, T, 𝒽)
    
    # set up lattice
    Δ, scale = LatticeGrid(L, PBC, Basis);
    
    useT = length(T)>0
    if !useT
       @assert length(ℓ)>0
    end
    M = useT ? length(T) : length(ℓ)   

    ns = 1:num_histories*length(𝒽)*M
    
    if useT
        args = [[Δ, therm_runtime, runtime, useT, T[rem(n-1,M)+1], 𝒽[rem(div(n-1,M),length(𝒽))+1]] for n in ns]
    else
        args = [[Δ, therm_runtime, runtime, useT, ℓ[rem(n-1,M)+1], 𝒽[rem(div(n-1,M),length(𝒽))+1]] for n in ns]
    end

    if multiProcess
        results = pmap(mpfun1, args)
    else
        results = Array{Any}(undef, length(ns))
        Threads.@threads for n in ns
            results[n] = mpfun1(args[n])
        end
    end


    x = [[[Array{Float64, 3}(undef,0,0,0) for _ in 1:num_histories] for _ in 1:length(𝒽)] for _ in 1:M]
    δ = [[[Array{Float64, 3}(undef,0,0,0) for _ in 1:num_histories] for _ in 1:length(𝒽)] for _ in 1:M]
    Mag = zeros(Float64, M, length(𝒽), num_histories)
    Perc = zeros(Float64, M, length(𝒽), num_histories)
    p = zeros(Int64, M, length(𝒽), num_histories)

    for n in ns
        m,t = divrem(n-1,M) .+ (1,1)
        h,i = divrem(m-1,length(𝒽)) .+ (1,1)

        x[t][i][h] = results[n][1]
        δ[t][i][h] = results[n][2]
        Mag[t,i,h]   = results[n][3]
        Perc[t,i,h]   = results[n][4]
        p[t,i,h]   = size(x[t][i][h], 2)
    end
    
    return x, δ, Mag, Perc, p, length(Δ.cells[1])
end

# ### Single analysis routine

@everywhere function DiffAnalysisSingle(p, x, δ, tau)
    dim = size(x[1], 1)
    T = size(x[1], 3)
    
    t = range(0,T)
    xfit = log.(t[tau])
    
    valid_histories = findall(p .> 0) # those for which there are particles! - - equiv to nanmean...
    
    ts = []
    sq_disp = []
    step_corr = []
    for h in valid_histories
        append!(ts, t[tau])       
        append!(sq_disp, Msd(x[h])[tau])
        append!(step_corr, DirrCorr(δ[h])[tau])
    end
    
    if length(ts) == 0
        return [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN for _ in 1:T], [NaN for _ in 1:T-1]
    end
    
    #xfit1 = log.(ts)
    #yfit1 = log.(sq_disp)
    xfit1 = ts
    yfit1 = sq_disp
    
    #idx = abs.(step_corr) .> 0
    #xfit2 = log.(ts[idx])
    #yfit2 = log.(abs.(step_corr[idx]))
    xfit2 = ts
    yfit2 = step_corr ./ sign(step_corr[findmax(abs.(step_corr))[2]])
    
    # linear fit function
    funlin = (x, p) -> p[1] .+ x .* p[2]
    funpow = (x, p) -> p[1] .* x .^ p[2]
    
    
    # MSD fit
    D = [NaN, NaN]
    α = [NaN, NaN]
    try
        #p1 = [log.(2*dim*Dself), 1.0]
        #fit1 = curve_fit(funlin, xfit1, yfit1, p1);

        #Est1 = fit1.param
        #Cov1 = estimate_covar(fit1)

        #D = [exp(Est1[1]), exp(Est1[1])*sqrt(Cov1[1,1])]
        #α = [Est1[2], sqrt(Cov1[2,2])]
        
        p1 = [2*dim*Dself, 1.0]
        fit1 = curve_fit(funpow, xfit1, yfit1, p1);

        Est1 = fit1.param
        Cov1 = estimate_covar(fit1)

        D = [Est1[1], sqrt(Cov1[1,1])]
        α = [Est1[2], sqrt(Cov1[2,2])]
        
        D ./= 2 * dim
    catch e
        print(e, "\n")
    end
    
    # DirrCorr fit
    C = [NaN, NaN]
    γ = [NaN, NaN]
    try
        #p2 = [log.(dim*Dself), -1.0]
        #fit2 = curve_fit(funlin, xfit2, yfit2, p2);
        
        #Est2 = fit2.param
        #Cov2 = estimate_covar(fit2)

        #C = [exp(Est2[1]), exp(Est2[1])*sqrt(Cov2[1,1])]
        #γ = [Est2[2], Cov2[2,2]]
        
        p2 = [dim*Dself, -1.0]
        fit2 = curve_fit(funpow, xfit2, yfit2, p2);
        
        Est2 = fit2.param
        Cov2 = estimate_covar(fit2)

        C = [Est2[1], Cov2[1,1]]
        γ = [Est2[2], Cov2[2,2]]
        
        C ./= dim
    catch e
        print(e, "\n")
    end
    
    
    MSD = zeros(T)
    VACF = zeros(T-1)
    for h in valid_histories
        MSD += Msd(x[h])
        VACF += DirrCorr(δ[h])
    end
    MSD ./= length(valid_histories)
    VACF ./= length(valid_histories)
    
    
    return D, α, C, γ, MSD, VACF
end

# ### Overall analysis routine

@everywhere function mpfun2(args)
    return DiffAnalysisSingle(args...)
end

@everywhere function DiffAnalysis(x, δ, p, runtime, ℓ, T, 𝒽)
    
    useT = length(T)>0
    if !useT
       @assert length(ℓ)>0
    end
    M = useT ? length(T) : length(ℓ)

    ns = 1:length(𝒽)*M

    args = []
    for n in ns
        i,t = divrem(n-1,M) .+ (1,1)
        push!(args, [p[t,i,:], x[t][i], δ[t][i], tau])
    end


    if multiProcess
        results = pmap(mpfun2, args)
    else
        results = Array{Any}(undef, length(ns))
        Threads.@threads for n in ns
            results[n] = mpfun2(args[n])
        end
    end


    D = zeros(2, M, length(𝒽))
    α = zeros(2, M, length(𝒽))
    C = zeros(2, M, length(𝒽))
    γ = zeros(2, M, length(𝒽))
    MSD = zeros(runtime+1, M, length(𝒽))
    DirrCorr = zeros(runtime, M, length(𝒽))
    
    for n in ns
        i,t = divrem(n-1,M) .+ (1,1)

        D[:,t,i]        = results[n][1]
        α[:,t,i]        = results[n][2]
        C[:,t,i]        = results[n][3]
        γ[:,t,i]        = results[n][4]
        MSD[:,t,i]      = results[n][5]
        DirrCorr[:,t,i] = results[n][6]
    end
    
    return D, α, C, γ, MSD, DirrCorr
end
