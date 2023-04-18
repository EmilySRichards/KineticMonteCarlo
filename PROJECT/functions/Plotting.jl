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

function PlotGraph(vertices, edges) # plot the graph - NOTE will only nicely plot graphs embeddable into R^2 (i.e. won't show PBCs wel!)
    figure()
    for e in edges
        r1 = vertices[e.∂[1]].x
        r2 = vertices[e.∂[2]].x
        
        plot([r1[1]; r2[1]], [r1[2]; r2[2]], color=(e.σ ? :red : :blue), zorder=1) 
    end
    
    for v in vertices
       scatter(v.x[1], v.x[2], color=:black, zorder=2) # color=(A(edges,v)<0 ? :yellow : :black)  
    end
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

function plotWithError(Y, X, C, YStd = zeros(size(Y)), XStd = zeros(size(X)))    
    errorbar(X, Y, xerr=XStd, yerr=YStd, ls="None", marker=".", color=C)
    # fill_between(X, Y-YStd, Y+YStd, alpha=0.5)
end

# ### Colourmaps

function jetmap(N)
    c = cgrad(:jet);
    return (N==1) ? [(red(c[1]), green(c[1]), blue(c[1]))] : [(red(c[i]), green(c[i]), blue(c[i])) for i in range(0,1,length=N)]
end
