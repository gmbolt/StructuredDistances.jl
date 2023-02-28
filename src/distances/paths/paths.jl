using Distances, InvertedIndices

export LCS, FastLCS
export print_info, get_info
export LSP, FastLSP
## Interaction Distances

struct LCS <: Metric end
struct FastLCS <: Metric
    curr_row::Vector{Int}
    prev_row::Vector{Int}
    function FastLCS(K::Int)
        new(zeros(Int, K), zeros(Int, K))
    end
end

function Base.show(io::IO, d::FastLCS)
    print(io, "FastLCS{$(length(d.curr_row))}")
end

# LCS
function (dist::LCS)(X::Vector{T}, Y::Vector{T})::Float64 where {T}

    n = length(X)
    m = length(Y)


    prevRow = 0:m
    currRow = zeros(Int, m + 1)
    @inbounds begin
        for i = 1:n
            currRow[1] = i
            for j = 1:m
                if X[i] == Y[j]
                    currRow[j+1] = prevRow[j]
                else
                    currRow[j+1] = min(currRow[j], prevRow[j+1]) + 1
                end
            end
            prevRow = copy(currRow)
        end
    end
    return currRow[end]
end

function (dist::LCS)(X::Tuple, Y::Tuple)::Float64

    n = length(X)
    m = length(Y)

    prevRow = 0:m
    currRow = zeros(Int, m + 1)
    @inbounds begin
        for i = 1:n
            currRow[1] = i
            for j = 1:m
                if X[i] == Y[j]
                    currRow[j+1] = prevRow[j]
                else
                    currRow[j+1] = min(currRow[j], prevRow[j+1]) + 1
                end
            end
            prevRow = copy(currRow)
        end
    end
    return currRow[end]
end

# With storage 
function (d::FastLCS)(
    X::Vector{T}, Y::Vector{T}
)::Float64 where {T}

    n = length(X)
    m = length(Y)

    @assert (n > 0) & (m > 0) "both paths must be either of type Nothing or of nonzero length."

    prev_row = view(d.prev_row, 1:(m+1))
    curr_row = view(d.curr_row, 1:(m+1))

    copy!(prev_row, 0:m)
    curr_row .= 0


    @inbounds begin
        for i = 1:n
            curr_row[1] = i
            for j = 1:m
                if X[i] == Y[j]
                    curr_row[j+1] = prev_row[j]
                else
                    curr_row[j+1] = min(curr_row[j], prev_row[j+1]) + 1
                end
            end
            copy!(prev_row, curr_row)
        end
    end
    return curr_row[m+1]
end

# Distances to the null 
function (dist::Union{LCS,FastLCS})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return length(Y)
end
function (dist::Union{LCS,FastLCS})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    return length(X)
end
function (dist::Union{LCS,FastLCS})(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end



# Get locations of longest common subseq (for visuals)

function get_info(
    d::Union{LCS,FastLCS},
    X::Vector{T}, Y::Vector{T}
) where {T}

    C = zeros(Int, length(X) + 1, length(Y) + 1)
    indx, indy = (zeros(Bool, length(X)), zeros(Bool, length(Y)))
    C[:, 1] = [i for i = 0:length(X)]
    C[1, :] = [i for i = 0:length(Y)]

    for j = 1:length(Y)
        for i = 1:length(X)
            C[i+1, j+1] = minimum([
                C[i, j] + 2 * (X[i] != Y[j]),
                C[i, j+1] + 1,
                C[i+1, j] + 1
            ])
        end
    end

    i = length(X) + 1
    j = length(Y) + 1
    while (i ≠ 1) & (j ≠ 1)
        if C[i, j] == (C[i-1, j] + 1)
            i = i - 1
        elseif C[i, j] == (C[i, j-1] + 1)
            j = j - 1
        elseif C[i, j] == (C[i-1, j-1])
            indx[i-1] = true
            indy[j-1] = true
            i = i - 1
            j = j - 1
        end
    end

    return indx, indy
end


# Longest Common Subpath (LSP)

struct LSP <: Metric end

function (dist::LSP)(X::Vector{T}, Y::Vector{T})::Float64 where {T}

    # Here we consider a Dynamic programming approach
    n = length(X)
    m = length(Y)

    @assert (n > 0) & (m > 0) "both paths must be either of type Nothing or of nonzero length."

    prev_row = zeros(Float64, m + 1)
    curr_row = zeros(Float64, m + 1)
    z = 0.0

    @inbounds begin
        for i = 1:n
            for j = 1:m
                if X[i] == Y[j]
                    curr_row[j+1] = prev_row[j] + 1.0 # Subpath length increment
                    if curr_row[j+1] > z
                        z = curr_row[j+1]
                    end
                else
                    curr_row[j+1] = 0.0
                end
            end
            copy!(prev_row, curr_row)
        end
    end
    return n + m - 2 * z
end

struct FastLSP <: Metric
    curr_row::Vector{Float64}
    prev_row::Vector{Float64}
    function FastLSP(K::Int)
        new(zeros(Float64, K), zeros(Float64, K))
    end
end

function Base.show(io::IO, d::FastLSP)
    print(io, "FastLSP{$(length(d.curr_row))}")
end

function (dist::FastLSP)(X::Vector{T}, Y::Vector{T})::Float64 where {T}
    # Here we take a Dynamic programming approach, but use pre-allocated arrays for storage.
    n = length(X)
    m = length(Y)

    @assert (n > 0) & (m > 0) "both paths must be either of type Nothing or of nonzero length."

    prev_row = view(dist.prev_row, 1:(m+1))
    curr_row = view(dist.curr_row, 1:(m+1))

    prev_row .= 0.0
    curr_row .= 0.0
    z = 0.0

    @inbounds begin
        for i = 1:n
            for j = 1:m
                if X[i] == Y[j]
                    curr_row[j+1] = prev_row[j] + 1.0 # Subpath length increment
                    if curr_row[j+1] > z
                        z = curr_row[j+1]
                    end
                else
                    curr_row[j+1] = 0.0
                end
            end
            copy!(prev_row, curr_row)
        end
    end
    return n + m - 2 * z
end

function (dist::Union{LSP,FastLSP})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return length(Y)
end
function (dist::Union{LSP,FastLSP})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    return length(X)
end
function (dist::Union{LSP,FastLSP})(X::Nothing, Y::Nothing) where {T}
    return 0.0
end

function get_info(
    d::Union{LSP,FastLSP},
    X::Vector{T}, Y::Vector{T}
) where {T}

    n, m = (length(X), length(Y))
    curr_row = zeros(Int, m + 1)
    prev_row = zeros(Int, m + 1)
    z = 0

    ix, iy = (0, 0) # Denotes end of common subpath in X and Y resp.

    begin
        for i = 1:n
            for j = 1:m
                if X[i] == Y[j]
                    curr_row[j+1] = prev_row[j] + 1 # Subpath length increment
                    if curr_row[j+1] > z
                        z = curr_row[j+1] # Found new longest subpath
                        ix, iy = (i, j)
                    end
                else
                    curr_row[j+1] = 0
                end
            end
            copy!(prev_row, curr_row)
        end
    end

    indx = [(ix - z + 1 ≤ i ≤ ix) ? true : false for i in 1:n]
    indy = [(iy - z + 1 ≤ i ≤ iy) ? true : false for i in 1:m]
    return indx, indy
end

function print_info(
    d::Union{LSP,FastLSP,LCS,FastLCS},
    x::Vector{T}, y::Vector{T}
) where {T}

    indx, indy = get_info(d, x, y)
    ix, iy = (0, 0)
    while true
        ixerr = isnothing(findnext(indx, ix + 1)) ? length(x) : findnext(indx, ix + 1) - 1
        iyerr = isnothing(findnext(indy, iy + 1)) ? length(y) : findnext(indy, iy + 1) - 1
        for j in (ix+1):ixerr
            println("$(x[j]) (delete)")
        end
        for j in (iy+1):iyerr
            println("$(y[j]) (add)")
        end
        if isnothing(findnext(indx, ix + 1))
            break
        end
        ix = findnext(indx, ix + 1)
        iy = findnext(indy, iy + 1)
        println("$(x[ix]) (keep)")
    end
end
