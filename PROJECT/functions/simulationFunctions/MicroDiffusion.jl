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

@everywhere function MicroDiffnSetup(vertices, edges, numToFlip)
    # initialise entire system in ground state
    GroundState!(vertices, edges)
    
    # flip numEdges random spins
    valid_edges = collect(eachindex(edges))
    flipped_edges = []
    
    for n in 1:numToFlip
        if valid_edges==[]
            break
        end
        
        α = rand(valid_edges)
        
        push!(flipped_edges, α) # add edge to list to flip
        
        # remove α AND other edges which share vertices with it
        deleteat!(valid_edges, findall(valid_edges.==α))
        for i in edges[α].∂
            for β in vertices[i].δ
                deleteat!(valid_edges, findall(valid_edges.==β))
            end                
        end
    end

    for α in flipped_edges
        edges[α].σ = !edges[α].σ
    end
end

# ### Single spin-flip dynamics routine 

@everywhere function MicroDiffn(vertices, edges, runtime, 𝒽, allowBckgdMoves)
    
    dim = length(vertices[1].x)
    
    # find all the excitations
    js = []
    for j in eachindex(vertices)
        Aj = A(edges, vertices[j])
        Qj = abs(Q(edges, vertices[j]))
        
        if (isSpinIce ? (Qj == 3 || Qj == 2) : Aj == -1)
            push!(js, j)
        end
    end
    
    xs = zeros(dim, length(js), runtime+1)
    δs = zeros(dim, length(js), runtime)
    for n in eachindex(js)
        xs[:,n,1] = vertices[js[n]].x
    end
    
    # actual simulation
    for t in 1:runtime
        xs[:,:,t+1] = xs[:,:,t]
        
        for _ in edges
            β = rand(eachindex(edges))
            ΔE = ΔE_flip(vertices, edges, β, 𝒽)
            
            j1 = edges[β].∂[1]
            j2 = edges[β].∂[2]
            
            if ΔE == 0 && (allowBckgdMoves || (j1 in js || j2 in js))
                
                edges[β].σ = !edges[β].σ
                
                # displacement of edge (fixed to account for PBCs)
                Δ = vertices[j2].x - vertices[j1].x
                for d in 1:length(Δ)
                    Δ[d] /= (abs(Δ[d])>1) ? -abs(Δ[d]) : 1 # note MINUS abs to ensure orientation is right (i.e. Δ>0 if going from RHS to LHS)
                end
                
                # if a prtcl is linked to this edge, move it - note ΔE=/=0 if prtcl on both vertices => ignore this case
                n1 = findfirst(js.==j1)
                n2 = findfirst(js.==j2)
                if n1!=nothing     # j1 = js[n1] = excitation
                    js[n1] = j2
                    xs[:,n1,t+1] += Δ # = vertices[j2].x # 
                    δs[:,n1,t] += Δ
                elseif n2!=nothing # j2 = js[n2] = excitation
                    js[n2] = j1
                    xs[:,n2,t+1] -= Δ # = vertices[j1].x # 
                    δs[:,n2,t] -= Δ
                end
                
                # else no prtcls => nothing to track!
            end
        end
    end
    
    return xs, δs
end

# ### Double spin-flip dynamics routine 

@everywhere function MicroDiffn_2flip(vertices, edges, runtime, 𝒽)
    
    dim = length(vertices[1].x)
    
    # find all the excitations
    js = []
    for j in eachindex(vertices)
        Aj = A(edges, vertices[j])
        Qj = abs(Q(edges, vertices[j]))
        
        if (isSpinIce ? (Qj == 3 || Qj == 2) : Aj == -1)
            push!(js, j)
        end
    end
    
    xs = zeros(dim, length(js), runtime+1)
    δs = zeros(dim, length(js), runtime)
    for n in eachindex(js)
        xs[:,n,1] = vertices[js[n]].x
    end
    
    # actual simulation
    for t in 1:runtime
        xs[:,:,t+1] = xs[:,:,t]
        
        for _ in vertices
            
            # propose flips
            i = rand(eachindex(vertices)) # shared vertex
            𝜷 = sample(vertices[i].δ, 2; replace=true) # two nearest-neighbour spins to flip (in order)
            
            𝒊 = [edges[𝜷[n]].∂[findfirst(edges[𝜷[n]].∂ .!= i)] for n in 1:2] # outer vertices (but may still coincide)
            
            ΣA = A(edges, vertices[i]) + A(edges, vertices[𝒊[1]]) + A(edges, vertices[𝒊[2]])
            
            # calculate overall energy change
            ΔE = ΔE_2flip(vertices, edges, 𝜷, 𝒊, i, 𝒽)

            # decide whether to accept and perform the move
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA>0 # energy AND magnetisation conserved AND no pair diffusion moves (i.e. no particle at central site i)
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA<0 # energy AND magnetisation conserved AND ONLY pair diffusion moves (i.e. no particle at central site i)
            if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ # energy AND magnetisation conserved
                
                edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                
                # ΔE=0 => an excitation is linked to this edge => move it
                # we choose to assume the central particle is fixed => valid way of tracking them
                
                # displacement of edge (fixed to account for PBCs)
                Δ1 = vertices[i].x - vertices[𝒊[1]].x
                for d in 1:length(Δ1)
                    Δ1[d] /= (abs(Δ1[d])>1) ? -abs(Δ1[d]) : 1 # note MINUS abs to ensure orientation is right (i.e. Δ>0 if going from RHS to LHS)
                end
                
                Δ2 = vertices[𝒊[2]].x - vertices[i].x
                for d in 1:length(Δ2)
                    Δ2[d] /= (abs(Δ2[d])>1) ? -abs(Δ2[d]) : 1 # note MINUS abs to ensure orientation is right (i.e. Δ>0 if going from RHS to LHS)
                end
                
                Δ = Δ2 + Δ1
                
                n1 = findfirst(js.==𝒊[1])
                n2 = findfirst(js.==𝒊[2])
                
                if n1!=nothing     # j1 = js[n1] = excitation
                    js[n1] = 𝒊[2]
                    xs[:,n1,t+1] += Δ # = vertices[𝒊[2]].x
                    δs[:,n1,t] += Δ
                elseif n2!=nothing # j2 = js[n2] = excitation
                    js[n2] = 𝒊[1]
                    xs[:,n2,t+1] -= Δ # = vertices[𝒊[1]].x
                    δs[:,n2,t] -= Δ
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

@everywhere function DiffSimSingle(vertices, edges, therm_runtime, runtime, useT, ℓorT, 𝒽)
    
    # thermalise to correct temperature OR correct number of particles
    if useT
        MicroKuboSetup(vertices, edges, therm_runtime, ℓorT, 𝒽, false)
    else
        MicroDiffnSetup(vertices, edges, ℓorT)
        
        if twoFlip # allow particles to separate before we start tracking them!
            MicroDiffn_2flip(vertices, edges, therm_runtime, 𝒽)
        else
            MicroDiffn(vertices, edges, therm_runtime, 𝒽, true)
        end
    end
    
    M = 0
    for edge in edges
        M += (-1)^edge.σ
    end
    M /= length(edges)
    
    bondNumber = 0
    maxBondNumber = 0
    for vertex in vertices
        z = length(vertex.δ)
        z₋ = 0
        for α in vertex.δ
            z₋ += edges[α].σ ? 1 : 0
        end
        z₊ = z - z₋
        
        bondNumber += z₊*z₋/2
        maxBondNumber += z*(z-1)/2
    end
    ℙ = bondNumber/maxBondNumber
    
    # track the paths of the resulting excitations (can't annihiliate b/c microcanonical!)
    if twoFlip
        x, δ = MicroDiffn_2flip(vertices, edges, runtime, 𝒽)
    else
        x, δ = MicroDiffn(vertices, edges, runtime, 𝒽, true)
    end
    
    return x, δ, M, ℙ
end

# ### Overall diffusion routine

@everywhere function mpfun1(args)
    return DiffSimSingle(args...)
end

@everywhere function DiffSim(L, PBC, Basis, therm_runtime, runtime, ℓ, T, 𝒽)
    
    # set up lattice
    vertices, edges, scale = LatticeGrid(L, PBC, Basis);
    
    useT = length(T)>0
    if !useT
       @assert length(ℓ)>0
    end
    M = useT ? length(T) : length(ℓ)   

    ns = 1:num_histories*length(𝒽)*M
    
    if useT
        args = [[deepcopy(vertices), deepcopy(edges), therm_runtime, runtime, useT, T[rem(n-1,M)+1], 𝒽[rem(div(n-1,M),length(𝒽))+1]] for n in ns]
    else
        args = [[deepcopy(vertices), deepcopy(edges), therm_runtime, runtime, useT, ℓ[rem(n-1,M)+1], 𝒽[rem(div(n-1,M),length(𝒽))+1]] for n in ns]
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
    
    return x, δ, Mag, Perc, p, length(vertices)
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
