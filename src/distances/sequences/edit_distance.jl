using StatsBase, Printf, Distances

export EditDistance, EditDist
export FastEditDistance, FastEditDist
export FixedPenaltyEditDistance, FixPenEditDist
export FastFixedPenaltyEditDistance, FastFixPenEditDist
export AvgSizeFpEditDistance
export print_info, get_info, get_info_deep

# EditDistance
# ---
struct EditDistance{T<:Metric} <: Metric
    ground_dist::T
end

const EditDist = EditDistance

function (d::EditDistance)(S1::Vector{T}, S2::Vector{T}) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d₀ = d.ground_dist
        prev_row = pushfirst!(cumsum([d₀(nothing, p) for p in S2]), 0.0)
        curr_row = zeros(Float64, length(S2) + 1)

        for i = 1:length(S1)
            curr_row[1] = prev_row[1] + length(S1[i])
            for j = 1:(length(S2))
                curr_row[j+1] = min(prev_row[j] + d₀(S1[i], S2[j]),
                    prev_row[j+1] + d₀(nothing, S1[i]),
                    curr_row[j] + d₀(nothing, S2[j]))
            end
            # @show curr_row
            copy!(prev_row, curr_row)
        end
        return curr_row[end]
    end
end

# EditDistance with Memory
# ---------------

struct FastEditDistance{T<:Metric} <: Metric
    ground_dist::T
    curr_row::Vector{Float64}
    prev_row::Vector{Float64}
    function FastEditDistance(ground_dist::S, K::Int) where {S<:Metric}
        new{S}(ground_dist, zeros(Float64, K), zeros(Float64, K))
    end
end

const FastEditDist = FastEditDistance

function Base.show(io::IO, d::FastEditDistance{T}) where {T<:SemiMetric}
    print(io, "FastEditDistance{$(T),$(length(d.curr_row))}")
end

function (d::FastEditDistance)(S1::Vector{T}, S2::Vector{T}) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d₀ = d.ground_dist
        prev_row = view(d.prev_row, 1:(length(S2)+1))
        curr_row = view(d.curr_row, 1:(length(S2)+1))
        # prev_row = d.prev_row
        # curr_row = d.curr_row

        prev_row[1] = 0.0
        for i in 1:length(S2)
            prev_row[i+1] = prev_row[i] + d₀(nothing, S2[i])
        end
        curr_row .= 0.0


        @views for i = 1:length(S1)
            curr_row[1] = prev_row[1] + length(S1[i])
            for j = 1:(length(S2))
                # @show i, j, prev_row[j], d.ground_dist(S1[i], S2[j])
                curr_row[j+1] = min(
                    prev_row[j] + d₀(S1[i], S2[j]),
                    prev_row[j+1] + d₀(nothing, S1[i]),
                    curr_row[j] + d₀(nothing, S2[j])
                )
            end
            # @show curr_row
            copy!(prev_row, curr_row)
        end
        return curr_row[length(S2)+1]
    end
end

function (d::Union{EditDist,FastEditDist})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return sum(p -> d.ground_dist(nothing, p), Y)
end
function (d::Union{EditDist,FastEditDist})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    return d(Y, X)
end
function (d::Union{EditDist,FastEditDist})(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end


function print_info(
    d::Union{EditDistance,FastEditDistance},
    S1::Vector{T}, S2::Vector{T}
) where {T}

    d₀ = d.ground_dist
    # First find the substitution matrix
    C = zeros(Float64, length(S1) + 1, length(S2) + 1)
    C[:, 1] = pushfirst!(cumsum([d₀(Λ, p) for p in S1]), 0.0)
    C[1, :] = pushfirst!(cumsum([d₀(Λ, p) for p in S2]), 0.0)

    for j in 1:length(S2)
        for i in 1:length(S1)
            C[i+1, j+1] = minimum([
                C[i, j] + d₀(S1[i], S2[j]),
                C[i, j+1] + d₀(nothing, S2[j]),
                C[i+1, j] + d₀(nothing, S1[i])
            ])
        end
    end
    # Now retrace steps to determine an optimal matching
    i, j = size(C)
    pairs = Tuple{Int,Int}[]
    while (i ≠ 1) & (j ≠ 1)
        if C[i, j] == (C[i, j-1] + d₀(nothing, S2[j-1]))
            pushfirst!(pairs, (0, j - 1))
            j = j - 1
        elseif C[i, j] == (C[i-1, j] + d₀(nothing, S1[i-1]))
            pushfirst!(pairs, (i - 1, 0))
            i = i - 1
        else
            pushfirst!(pairs, (i - 1, j - 1))
            i = i - 1
            j = j - 1
        end
    end
    for k in Iterators.reverse(1:(i-1))
        pushfirst!(pairs, (k, 0))
    end
    for k in Iterators.reverse(1:(j-1))
        pushfirst!(pairs, (0, k))
    end
    max_len = maximum(map(x -> length(@sprintf("%s", x)), S1))
    # @show outputs
    title = "\nOptimal Matching"
    println(title)
    println("-"^length(title), "\n")
    for (k, l) in pairs
        if k == 0
            tmp_S1 = "Λ"
            tmp_S2 = S2[l]
        elseif l == 0
            tmp_S2 = "Λ"
            tmp_S1 = S1[k]
        else
            tmp_S1, tmp_S2 = (S1[k], S2[l])
        end
        pad = max_len - length(@sprintf("%s", tmp_S1))
        println("$tmp_S1" * " "^pad, " → $tmp_S2")
    end

end

function get_info(
    d::Union{EditDistance,FastEditDistance},
    x::Vector{T}, y::Vector{T}
) where {T}

    d_g = d.ground_dist
    # First find the substitution matrix
    C = zeros(Float64, length(x) + 1, length(y) + 1)
    C[:, 1] = pushfirst!(cumsum([d_g(nothing, p) for p in x]), 0.0)
    C[1, :] = pushfirst!(cumsum([d_g(nothing, p) for p in y]), 0.0)

    for j in 1:length(y)
        for i in 1:length(x)
            C[i+1, j+1] = minimum([
                C[i, j] + d_g(x[i], y[j]),
                C[i, j+1] + d_g(nothing, x[i]),
                C[i+1, j] + d_g(nothing, y[j])
            ])
        end
    end
    # Now retrace steps to find indices of matched paths 
    i, j = size(C)
    indx, indy = (zeros(Bool, length(x)), zeros(Bool, length(y)))
    while (i ≠ 1) & (j ≠ 1)
        if C[i, j] == (C[i, j-1] + d_g(nothing, y[j-1]))
            j = j - 1
        elseif C[i, j] == (C[i-1, j] + d_g(nothing, x[i-1]))
            i = i - 1
        else
            i = i - 1
            j = j - 1
            indx[i], indy[j] = (true, true)
        end
    end
    return indx, indy
end

struct FixedPenaltyEditDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    ρ::Real
end

const FixPenEditDist = FixedPenaltyEditDistance

function (d::FixPenEditDist)(S1::Vector{T}, S2::Vector{T}) where {T}
    if length(S1) < length(S2)
        d(S2, S1)
    else
        prev_row = d.ρ * collect(0:length(S2))
        curr_row = zeros(Float64, length(S2) + 1)

        for i in eachindex(S1)
            curr_row[1] = i * d.ρ
            for j in eachindex(S2)
                curr_row[j+1] = min(
                    prev_row[j] + d.ground_dist(S1[i], S2[j]),
                    prev_row[j+1] + d.ρ,
                    curr_row[j] + d.ρ
                )
            end
            copy!(prev_row, curr_row)
        end
        return curr_row[end]
    end
end

# FixedPenaltyEditDistance with Memory
# ------------------------------------

struct FastFixedPenaltyEditDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    ρ::Real
    curr_row::Vector{Float64}
    prev_row::Vector{Float64}
    function FastFixedPenaltyEditDistance(ground_dist::S, ρ::Float64, K::Int) where {S<:SemiMetric}
        new{S}(ground_dist, ρ, zeros(Float64, K), zeros(Float64, K))
    end
end

const FastFixPenEditDist = FastFixedPenaltyEditDistance

function Base.show(io::IO, d::FastFixPenEditDist{T}) where {T<:SemiMetric}
    print(io, "FastFixPenEditDist{$(T),$(length(d.curr_row))}")
end

function (d::FastFixPenEditDist)(S1::Vector{T}, S2::Vector{T}) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d₀ = d.ground_dist
        prev_row = view(d.prev_row, 1:(length(S2)+1))
        curr_row = view(d.curr_row, 1:(length(S2)+1))

        prev_row[1] = 0.0
        for i in eachindex(S2)
            prev_row[i+1] = prev_row[i] + d.ρ
        end
        curr_row .= 0.0


        @views for i in eachindex(S1)
            curr_row[1] = prev_row[1] + length(S1[i])
            for j in eachindex(S2)
                curr_row[j+1] = min(
                    prev_row[j] + d₀(S1[i], S2[j]),
                    prev_row[j+1] + d.ρ,
                    curr_row[j] + d.ρ
                )
            end
            copy!(prev_row, curr_row)
        end
        return curr_row[length(S2)+1]
    end
end

# Distance to nothings 

function (d::Union{FixPenEditDist,FastFixPenEditDist})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return d.ρ * length(Y)
end
function (d::Union{FixPenEditDist,FastFixPenEditDist})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    d(Y, X)
end
function (d::Union{FixPenEditDist,FastFixPenEditDist})(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end



function get_info(
    d::Union{FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T}
) where {T}
    C = zeros(Float64, length(x) + 1, length(y) + 1)

    C[:, 1] = [d.ρ * i for i = 0:length(x)]
    C[1, :] = [d.ρ * i for i = 0:length(y)]

    for j = 1:length(y)
        for i = 1:length(x)
            C[i+1, j+1] = minimum([
                C[i, j] + d.ground_dist(x[i], y[j]),
                C[i, j+1] + d.ρ,
                C[i+1, j] + d.ρ
            ])
        end
    end
    # Now retrace steps to determine an optimal matching
    i, j = size(C)
    indx, indy = (zeros(Bool, length(x)), zeros(Bool, length(y)))
    while (i ≠ 1) | (j ≠ 1)
        if C[i, j] == (C[i-1, j] + d.ρ)
            i = i - 1
        elseif C[i, j] == (C[i, j-1] + d.ρ)
            j = j - 1
        else
            i = i - 1
            j = j - 1
            indx[i], indy[j] = (true, true)
        end
    end
    return indx, indy
end


function print_info(
    d::Union{EditDistance,FastEditDistance,FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T}
) where {T}

    indx, indy = get_info(d, x, y)
    max_len = maximum(map(z -> length(@sprintf("%s", z)), x))
    ix, iy = (0, 0)
    title = "\nOptimal Matching"
    println(title)
    println("-"^length(title), "\n")
    while true
        ixerr = isnothing(findnext(indx, ix + 1)) ? length(x) : findnext(indx, ix + 1) - 1
        iyerr = isnothing(findnext(indy, iy + 1)) ? length(y) : findnext(indy, iy + 1) - 1
        for j in (ix+1):ixerr
            tmp1 = x[j]
            tmp2 = "Null"
            pad = max_len - length(@sprintf("%s", tmp1))
            println(" "^pad * "$tmp1", " → $tmp2")
        end
        for j in (iy+1):iyerr
            tmp1 = "Null"
            tmp2 = y[j]
            pad = max_len - length(@sprintf("%s", tmp1))
            println(" "^pad * "$tmp1", " → $tmp2")
        end
        if isnothing(findnext(indx, ix + 1))
            break
        end
        ix = findnext(indx, ix + 1)
        iy = findnext(indy, iy + 1)
        tmp1, tmp2 = (x[ix], y[iy])
        pad = max_len - length(@sprintf("%s", tmp1))
        println(" "^pad * "$tmp1", " → $tmp2")
    end
end

function get_info_deep(
    d::Union{EditDistance,FastEditDistance,FixPenEditDist,FastFixPenEditDist},
    x::Vector{T}, y::Vector{T}
) where {T}

    d_g = d.ground_dist
    indx, indy = get_info(d, x, y)
    outx, outy = (
        [zeros(Bool, i) for i in length.(x)],
        [zeros(Bool, i) for i in length.(y)]
    )
    ix, iy = (0, 0)
    while true
        if isnothing(findnext(indx, ix + 1))
            break
        end
        ix = findnext(indx, ix + 1)
        iy = findnext(indy, iy + 1)
        outx[ix], outy[iy] = get_info(d_g, x[ix], y[iy])
    end
    return outx, outy
end


struct AvgSizeFpEditDistance{T<:Metric} <: Metric
    ground_dist::T
    ρ::Real
end

function (d::AvgSizeFpEditDistance)(S1::Vector{T}, S2::Vector{T}) where {T}

    d_ed = FixPenEditDist(d.ground_dist, d.ρ)(S1, S2)

    return d_ed + (mean(length.(S1)) - mean(length.(S2)))^2

end

