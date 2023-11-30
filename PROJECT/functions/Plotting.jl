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

# ### Graph plotting

function PlotGraph(S, Δ, toPlot=nothing) # plot the graph - NOTE will only nicely plot graphs embeddable into R^2 (i.e. won't show PBCs wel!)
    dim = length(Δ.cells[1][1].x)
    
    if dim==2
        fig = figure()
        ax = fig.add_subplot(111)
    else
        fig = figure()
        ax = fig.add_subplot(111, projection="3d") 
    end
    
    for (e, Δe) in enumerate(Δ.cells[2])
        if toPlot != nothing && !toPlot[2][e]
            continue
        end
        
        r1 = Δ.cells[1][Δe.∂[findfirst(Δe.η.<0)]].x
        r2 = Δ.cells[1][Δe.∂[findfirst(Δe.η.>0)]].x
        
        r0 = (GetCpt(S, e, true)>0) ? r1 : r2
        δr = (GetCpt(S, e, true)>0) ? r2-r1 : r1-r2
        δr *= 0.5
        
        if dim==2
            ax.plot([r1[1]; r2[1]], [r1[2]; r2[2]], color=(GetCpt(S, e, false)>0 ? :red : :blue), zorder=1)
            
            ax.arrow(r0[1], r0[2], δr[1], δr[2], color=(GetCpt(S, e, false)>0 ? :red : :blue), 
                  length_includes_head=true, head_length=0.15, head_width=0.15)
        else
            ax.plot3D([r1[1]; r2[1]], [r1[2]; r2[2]], [r1[3]; r2[3]], color=(GetCpt(S, e, false)>0 ? :red : :blue), zorder=1)
            
            ax.quiver(r0[1], r0[2], r0[3], δr[1], δr[2], δr[3], color=(GetCpt(S, e, false)>0 ? :red : :blue))
        end
    end
    
    for (i, Δi) in enumerate(Δ.cells[1])
        if toPlot != nothing && !toPlot[1][i]
            continue
        end
        
        if dim==2
            ax.scatter(Δi.x[1], Δi.x[2], color=:black, zorder=2)
        else
            ax.scatter3D(Δi.x[1], Δi.x[2], Δi.x[3], color=:black, zorder=2)
        end
    end
    
    axis("equal")
    
    return fig, ax
end

# ### Graphs.jl plotting

# + active=""
# edgewidth_dict = ones(length(verticesH), length(verticesH))
# for n in eachindex(b)
#     i = b[n]
#     j = (n<length(b)) ? b[n+1] : b[1]
#     
#     edgewidth_dict[i, j] = 5
#     edgewidth_dict[j, i] = 5
# end
# -

# ### Plotting with errors

function plotWithError(Y, X, C, mk, label="", YStd = zeros(size(Y)), XStd = zeros(size(X)))    
    p = errorbar(X, Y, xerr=XStd, yerr=YStd, ls="None", marker=mk, color=C, label=label)
    # fill_between(X, Y-YStd, Y+YStd, alpha=0.5)
    return p
end

# ### Colourmaps

function jetmap(N)
    c = cgrad(:jet);
    return (N==1) ? [(red(c[1]), green(c[1]), blue(c[1]))] : [(red(c[i]), green(c[i]), blue(c[i])) for i in range(0,1,length=N)]
end
