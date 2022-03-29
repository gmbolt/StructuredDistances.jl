using NetworkLayout, Roots
export MultigraphPlot, relabel, get_curves

function relabel(edgelist::Vector{Tuple{Int,Int}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{Int,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

function relabel(edgelist::Vector{Tuple{String,String}})
    edgelist_relab = Tuple{Int,Int}[]
    vmap = Dict{String,Int}()
    i = 1
    for edge in edgelist
        for v in edge
            if v ∉ keys(vmap)
                vmap[v] = i 
                i+=1
            end 
        end 
        push!(edgelist_relab, map(x->vmap[x], edge))
    end 
    rev_vmap = Dict(v=>k for (k,v) in vmap)
    return edgelist_relab, rev_vmap
end 

node_edge_intercept(x,w,h,r) = (4*h/(w^2) * x*(w-x))^2 + x^2 - r^2

function get_arc(
    src, dst;
    bendprop=0.1
    )
    # f(x) = 4h/w² x(w-x) (equation for arc from (0,0)->(w,0))
    # We then then translate and rotate this 
    w = euclidean(src,dst)
    h = bendprop*euclidean(src,dst)
    θ = angle(complex(dst...)-complex(src...))
    cosθ, sinθ = (cos(θ), sin(θ))
    C = 4*h/w^2
    rawx = 0:0.01:w
    srcx,srcy = src
    xvals = [x*cosθ - C*x*(w-x)*sinθ + srcx for x in rawx]
    yvals = [x*sinθ + C*x*(w-x)*cosθ + srcy for x in rawx]
    return xvals, yvals
end


function get_arc_shorten(
    src, dst;
    bendprop=0.1,
    shorten=0.2
    )

    w = euclidean(src,dst)
    h = bendprop*euclidean(src,dst)
    x_shorten = find_zero(x->node_edge_intercept(x,w,h,shorten), (0.0,shorten))
    θ = angle(complex(dst...)-complex(src...))
    cosθ, sinθ = (cos(θ), sin(θ))
    C = 4*h/w^2
    rawx = x_shorten:0.01:(w-x_shorten)
    srcx,srcy = src
    xvals = [x*cosθ - C*x*(w-x)*sinθ + srcx for x in rawx]
    yvals = [x*sinθ + C*x*(w-x)*cosθ + srcy for x in rawx]
    return xvals, yvals
end

function get_curves(
    edges::AbstractVector,
    locs::AbstractVector;
    bendprop=0.1,
    shorten=0.2
    )

    xout, yout = (Vector{Float64}[],Vector{Float64}[])

    for (i,j) in edges
        xtmp, ytmp = get_arc_shorten(locs[i],locs[j],bendprop=bendprop,shorten=shorten)
        push!(xout,xtmp)
        push!(yout,ytmp)
    end 

    return xout, yout

end 

@userplot MultigraphPlot
@recipe function f(
    h::MultigraphPlot; 
    relabelnodes=true,
    shorten=0.1,
    bendprop=0.1,
    seed=1
    )
    edgelist = h.args[1]
    
    do_relabel = (relabelnodes) || (eltype(edgelist)==Tuple{String,String})
    if do_relabel
        edgelist, rev_vmap = relabel(edgelist)
        V = length(keys(rev_vmap))
        adjmat = zeros(Int, V, V)
        for (i,j) in edgelist
            adjmat[i,j] += 1 
        end 
    end 
    aspect_ratio --> 1
    locs = NetworkLayout.spring(adjmat,seed=seed)
    
    @series begin 
        # seriestype := :line
        linecolor := :blue
        arrow --> true
        get_curves(edgelist, locs, shorten=shorten, bendprop=bendprop)
    end 

    @series begin 
        seriestype := :scatter
        markercolor := :blue
        locs
    end 
end 