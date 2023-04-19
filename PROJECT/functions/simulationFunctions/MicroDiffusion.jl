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
    for edge in edges
        if sixVertex
            edge.σ = vertices[edge.∂[1]].x[1]-vertices[edge.∂[2]].x[1]==0 # gives ~GS ONLY for PBCs on square lattice
        else
            edge.σ = false
        end
        edge.D = 0
    end
    
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
        deleteat!(valid_edges, findall(x->x==α, valid_edges))
        for i in edges[α].∂
            for β in vertices[i].δ
                deleteat!(valid_edges, findall(x->x==β, valid_edges))
            end                
        end
    end

    for α in flipped_edges
        edges[α].σ = !edges[α].σ
    end
end

# ### Single spin-flip dynamics routine 

@everywhere function MicroDiffn(vertices, edges, runtime, 𝒽)
    
    dim = length(vertices[1].x)
    
    # find all the excitations
    js = []
    for j in eachindex(vertices)
        if (sixVertex ? Atilde(edges, vertices[j])==4 : A(edges, vertices[j])<0) # in 6-vertex case, A_e = 2^2 = 4, A_2e = 4^2 = 16
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

            if ΔE == 0 # note ΔE NEVER zero if edge links two vertices (or none at all) => can ignore this case
                edges[β].σ = !edges[β].σ
                
                # ΔE=0 => an excitation is linked to this edge => move it 
                j1 = edges[β].∂[1]
                j2 = edges[β].∂[2]
                
                # displacement of edge (fixed to account for PBCs)
                Δ = vertices[j2].x - vertices[j1].x
                for d in 1:length(Δ)
                    Δ[d] /= (abs(Δ[d])>1) ? -abs(Δ[d]) : 1 # note MINUS abs to ensure orientation is right (i.e. Δ>0 if going from RHS to LHS)
                end
                
                n1 = findfirst(js.==j1)
                n2 = findfirst(js.==j2)
                if n1!=nothing     # j1 = js[n1] = excitation
                    js[n1] = j2
                    xs[:,n1,t+1] += Δ
                    δs[:,n1,t] += Δ
                elseif n2!=nothing # j2 = js[n2] = excitation
                    js[n2] = j1
                    xs[:,n2,t+1] -= Δ
                    δs[:,n2,t] -= Δ
                end
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
        if (sixVertex ? Atilde(edges, vertices[j])==4 : A(edges, vertices[j])<0) # in 2D spin ice case, A_e = 2^2 = 4, A_2e = 4^2 = 16
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
        
        for _ in 1:floor(Int64, length(edges)/2)
            # propose flips
            i = rand(eachindex(vertices)) # shared vertex
            𝜷 = sample(vertices[i].δ, 2; replace=false) # two nearest-neighbour spins to flip (in order)
            
            𝒊 = [edges[𝜷[n]].∂[findfirst(edges[𝜷[n]].∂ .!= i)] for n in 1:2] # outer vertices (but may still coincide)
            
            ΣA = A(edges, vertices[i]) + A(edges, vertices[𝒊[1]]) + A(edges, vertices[𝒊[2]])
            
            # calculate overall energy change and current density between the two unshared vertices
            ΔE = ΔE_2flip(vertices, edges, 𝜷, 𝒊, 𝒽)
            Δj = Δj_2flip(vertices, edges, 𝜷, 𝒊, 𝒽)

            # decide whether to accept and perform the move
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA>0 # energy AND magnetisation conserved AND no pair diffusion moves (i.e. no particle at central site i)
            #if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ && ΣA<0 # energy AND magnetisation conserved AND ONLY pair diffusion moves (i.e. no particle at central site i)
            if ΔE == 0 && edges[𝜷[1]].σ!=edges[𝜷[2]].σ # energy AND magnetisation conserved
            #if ΔE == 0 # energy conserved
                
                edges[𝜷[1]].σ = !edges[𝜷[1]].σ
                edges[𝜷[2]].σ = !edges[𝜷[2]].σ
                
                # ΔE=0 => an excitation is linked to this edge => move it
                # we choose to assume the moving particle always starts at one of the edge vertices => valid way of tracking them if we don't allow repeat edges
                
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
                    xs[:,n1,t+1] += Δ
                    δs[:,n1,t] += Δ
                elseif n2!=nothing # j2 = js[n2] = excitation
                    js[n2] = 𝒊[1]
                    xs[:,n2,t+1] -= Δ
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
    corr ./= corr[1] # normalise
    
    return corr
end

# ### Single diffusion routine

@everywhere function DiffSimSingle(vertices, edges, therm_runtime, runtime, useT, ℓorT, 𝒽)

    # thermalise to correct temperature OR correct number of particles
    if useT
        MicroKuboSetup(vertices, edges, therm_runtime, ℓorT, 𝒽, false)
    else
        MicroDiffnSetup(vertices, edges, ℓorT)
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
        x, δ = MicroDiffn(vertices, edges, runtime, 𝒽)
    end
    
    return x, δ, M, ℙ
end

# ### Overall diffusion routine

@everywhere function DiffSim(L, PBC, Basis, therm_runtime, runtime, ℓ, T, 𝒽)
    
    # set up lattice
    vertices, edges = LatticeGrid(L, PBC, Basis);
    
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

    @everywhere function mpfun1(args)
        return DiffSimSingle(args...)
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
    num_histories = size(p)
    T = size(x[1], 3)
    
    t = range(0,T)
    xfit = log.(t[tau])
    
    valid_histories = findall(p .> 0) # those for which there are particles! - - equiv to nanmean...
    
    if length(valid_histories) == 0
        return [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN for _ in 1:T], [NaN for _ in 1:T-1]
    end
    
    sq_disp = zeros(T)
    step_corr = zeros(T-1)
    for h in valid_histories
        sq_disp += Msd(x[h])
        step_corr += DirrCorr(δ[h])
    end
    sq_disp ./= length(valid_histories)
    step_corr ./= length(valid_histories)
    
    if sq_disp == zeros(size(sq_disp))
        return [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN, NaN], [NaN for _ in 1:T], [NaN for _ in 1:T-1]
    end
    
    # linear fit function
    fun = (x, p) -> p[1] .+ x .* p[2]
    
    
    # MSD fit
    p1 = [0.0, 1.0]
    yfit1 = log.(sq_disp[tau])
    fit1 = curve_fit(fun, xfit, yfit1, p1);

    Est = fit1.param
    Cov = estimate_covar(fit1)

    D = [exp(Est[1]), exp(Est[1])*sqrt(Cov[1,1])] ./4 # div by 4 b/c in 2 dims, x^2~4Dt and both x and t are measured in units of a=δt=1
    α = [Est[2], sqrt(Cov[2,2])]
    
    # DirrCorr fit
    #p2 = [2.0, 1.0] # 2 b/c in 2D
    #yfit2 = log.(abs.(step_corr[tau]))
    #fit2 = curve_fit(fun, xfit, yfit2, p2);

    #Est = fit2.param
    #Cov = estimate_covar(fit2)

    #C = [-Est[1], Cov[1,1]]
    #γ = [Est[2], Cov[2,2]]
    C = [NaN, NaN]
    γ = [NaN, NaN]
    
    return D, α, C, γ, sq_disp, step_corr
end

# ### Overall analysis routine

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

    @everywhere function mpfun2(args)
        return DiffAnalysisSingle(args...)
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

        D[:,t,i]    = results[n][1]
        α[:,t,i]    = results[n][2]
        C[:,t,i]    = results[n][3]
        γ[:,t,i]    = results[n][4]
        MSD[:,t,i] = results[n][5]
        DirrCorr[:,t,i] = results[n][6]
    end

    return D, α, C, γ, MSD, DirrCorr
end
