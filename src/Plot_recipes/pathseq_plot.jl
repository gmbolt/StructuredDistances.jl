using RecipesBase
export PathSeqPlot, NoisePlot 
export get_entry_xlocs

function get_entry_xlocs(
    lens::Vector{Int},
    width,
    entrymargin,
    pathmargin;
    origin=0.0
    )
    out = Float64[]
    x_tmp = origin
    for n in lens
        for i in 1:n 
            push!(out, x_tmp)
            x_tmp += width + entrymargin
        end 
        x_tmp -= entrymargin # Take of last entry margin
        x_tmp += pathmargin 
    end 
    return out
end 

@userplot PathSeqPlot
@recipe function f(
    h::PathSeqPlot; 
    entrymargin=0.1, 
    entryfontsize=20,
    pathmargin=0.5,
    entrycolor=:green,
    align=:center
    )
    obs = h.args[1]
    w,h = (1,1) 
    # Get entry centers 
    y = fill(h/2, sum(length, obs))
    x = get_entry_xlocs(
        length.(obs), 
        w, 
        entrymargin, 
        pathmargin, 
        origin=w/2
    )
    # Shift x according to alignment
    if align==:center 
        x .-= x[end]/2
    elseif align==:right 
        x .-= (x[end] + w/2)
    else
        if align!=:left 
            error("Align specification not supported. Must be, :center, :left or :right")
        end 
    end 
    annotations := [(x,y,string(i)) for (x,y,i) in zip(x,y,vcat(obs...))]
    annotationfontsize := entryfontsize
    x_cords, y_cords = rectangle_corners(x,y, w, h; anchor=:center)
    showaxis --> false 
    axis --> nothing 
    aspect_ratio --> 1 
    legend --> false 
    if typeof(entrycolor)==Symbol
        fillcolor --> entrycolor
    elseif typeof(entrycolor)<:Vector{Vector{T}} where {T}
        fillcolor --> permutedims(vcat(entrycolor...))
    end 
    @series begin
        seriestype := :shape
        x_cords,y_cords
    end 
end 

@userplot NoisePlot
@recipe function f(
    h::NoisePlot; 
    entrymargin=0.1, 
    entryfontsize=20,
    stdcolor=:green, errcolor=:magenta,
    pathmargin=0.5,
    align=:center
    )
    obs, err = (h.args[1], h.args[2])
    w,h = (1,1) 
    # Get entry centers 
    y = fill(h/2, sum(length, obs))
    x = get_entry_xlocs(
        length.(obs), 
        w, 
        entrymargin, 
        pathmargin, 
        origin=w/2
    )
    # Shift x according to alignment
    if align==:center 
        x .-= x[end]/2
    elseif align==:right 
        x .-= (x[end] + w/2)
    else
        if align!=:left 
            error("Align specification not supported. Must be, :center, :left or :right")
        end 
    end 
    annotations := [(x,y,string(i)) for (x,y,i) in zip(x,y,vcat(obs...))]
    annotationfontsize := entryfontsize
    x_cords, y_cords = rectangle_corners(x,y, w, h; anchor=:center)
    showaxis --> false 
    axis --> nothing 
    aspect_ratio --> 1 
    legend --> false 
    fillcolor --> permutedims([i ? stdcolor : errcolor for i in vcat(err...)])
    @series begin
        seriestype := :shape
        x_cords,y_cords
    end 
end 

