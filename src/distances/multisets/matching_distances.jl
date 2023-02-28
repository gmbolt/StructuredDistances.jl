using StatsBase, Distances, Hungarian, Printf

export MatchingDistance, MatchDist
export FastMatchingDistance, FastMatchDist
export FixedPenaltyMatchingDistance, FixPenMatchDist
export MinDistMatchingDistance, MinDistMatchDist
export AvgSizeMatchingDistance, AvgSizeMatchDist
export get_cost_matrix_dynamic, get_cost_matrix_fixed


struct MatchingDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
end

const MatchDist = MatchingDistance

function get_cost_matrix_dynamic(
    d::MatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    if length(X) < length(Y)
        return get_cost_matrix_dynamic(d, Y, X)
    elseif length(X) == length(Y)
        return pairwise_inbounds(d.ground_dist, X, Y)
    else
        C = pairwise_inbounds(d.ground_dist, X, Y)
        null_dists = [d.ground_dist(nothing, p) for p in X]
        size_diff = length(X) - length(Y)
        return [C [x for x ∈ null_dists, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::MatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    if length(X) < length(Y)
        C = pairwise_inbounds(d.ground_dist, X, Y)
        pentalty_vec = [d.ground_dist(nothing, p) for p in Y]
        size_diff = length(Y) - length(X)
        return [C; [y for i = 1:size_diff, y ∈ pentalty_vec]]
    elseif length(X) == length(Y)
        return pairwise_inbounds(d.ground_dist, X, Y)
    else
        C = pairwise_inbounds(d.ground_dist, X, Y)
        penalty_vec = [d.ground_dist(nothing, p) for p in X]
        size_diff = length(X) - length(Y)
        return [C [x for x ∈ penalty_vec, j = 1:size_diff]]
    end

end


struct FastMatchingDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    C::Matrix{Float64}
    function FastMatchingDistance(ground_dist::S, K::Int) where {S<:SemiMetric}
        new{S}(ground_dist, zeros(K, K))
    end
end

const FastMatchDist = FastMatchingDistance

Base.show(io::IO, d::FastMatchingDistance{T}) where {T<:SemiMetric} = print(io, "FastMatchingDistance{$(T),$(size(d.C,1))}")

function get_cost_matrix_dynamic(
    d::FastMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
    else
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for i in 1:N
            C[i, M+1] = d.ground_dist(nothing, X[i])
        end
        for j in (M+2):N
            for i in 1:N
                C[i, j] = C[i, j-1]
            end
        end
    end
    return C
end

function get_cost_matrix_fixed(
    d::FastMatchingDistance,
    X::Vector{T}, Y::Vector{T}
) where {T}

    N, M = (length(X), length(Y))
    if N < M
        C = view(d.C, 1:M, 1:M)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for j in 1:M
            C[N+1, j] = d.ground_dist(nothing, Y[j])
        end
        for i in (N+2):N
            for j in 1:M
                C[i, j] = C[i-1, j]
            end
        end
    elseif N == M
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
    else
        C = view(d.C, 1:N, 1:N)
        pairwise_inbounds!(C, d.ground_dist, X, Y)
        for i in 1:N
            C[i, M+1] = d.ground_dist(nothing, X[i])
        end
        for j in (M+2):M
            for i in 1:N
                C[i, j] = C[i, j-1]
            end
        end
    end
    return C

end

function (d::Union{MatchDist,FastMatchDist})(X::Nothing, Y::Vector{T})::Float64 where {T}
    return sum(p -> d.ground_dist(nothing, p), Y)
end
function (d::Union{MatchDist,FastMatchDist})(X::Vector{T}, Y::Nothing)::Float64 where {T}
    return d(Y, X)
end
function (d::Union{MatchDist,FastMatchDist})(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end

struct FixedPenaltyMatchingDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    penalty::Float64
end

const FixPenMatchDist = FixedPenaltyMatchingDistance
# const FixPenMatchDist{T} = FixedPenaltyMatchingDistance{T}
# FixPenMatchDist(d_g::SemiMetric, penalty::Float64) = FixedPenaltyMatchingDistance(d_g, penalty)

Base.show(io::IO, d::FixPenMatchDist) = print(io, "$(typeof(d))(ρ=$(d.penalty))")

function get_cost_matrix_dynamic(
    d::FixPenMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        return pairwise_inbounds(d_g, X, Y)
    else
        C = pairwise_inbounds(d_g, X, Y)
        size_diff = (N - M)
        return [C [d.penalty for i ∈ eachindex(X), j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::FixPenMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    C = pairwise_inbounds(d_g, X, Y)
    if N < M
        size_diff = (M - N)
        return [C; [d.penalty for i in 1:size_diff, j ∈ eachindex(Y)]]
    elseif N > M
        size_diff = (N - M)
        return [C [d.penalty for i ∈ eachindex(X), j in 1:size_diff]]
    else
        return pairwise_inbounds(d_g, X, Y)
    end

end

# function (d::FixPenMatchDist)(
#     S1::Vector{T}, S2::Vector{T}
#     ) where {T}

#     length(S1) < length(S2) ? (N_min, N_max)=(length(S1),length(S2)) : (N_min,N_max)=(length(S2),length(S1))

#     C = fill(d.penalty, N_max, N_max)
#     @views Distances.pairwise!(C, d.ground_dist, S1, S2)

#     # If we have max_dist > 2*ρ we must extend C to allow for interactions in BOTH observations being un-matched.
#     if maximum(view(C, 1:length(S1), 1:length(S2))) > 2*d.penalty
#         C = [
#             C fill(d.penalty/2, N_max, N_min); 
#             fill(d.penalty/2, N_min, N_max) fill(0.0, N_min, N_min)
#             ]
#     end  
#     return hungarian(C)[2]
# end

function (d::FixPenMatchDist)(X::Nothing, Y::Vector{T})::Float64 where {T}
    return d.penalty * length(Y)
end
function (d::FixPenMatchDist)(X::Vector{T}, Y::Nothing)::Float64 where {T}
    d(Y, X)
end
function (d::FixPenMatchDist)(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end

struct AvgSizeMatchingDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    penalty::Float64
end

const AvgSizeMatchDist = AvgSizeMatchingDistance
# const AvgSizeMatchDist{T} = AvgSizeMatchingDistance{T}

Base.show(io::IO, d::AvgSizeMatchDist) = print(io, "$(typeof(d))(ρ=$(d.penalty))")


function get_cost_matrix_dynamic(
    d::AvgSizeMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        return pairwise_inbounds(d_g, X, Y)
    else
        C = pairwise_inbounds(d_g, X, Y)
        mean_size = mean(y -> d_g(y, nothing), Y) # Find mean size in smaller set
        pentalty_vec = [abs(d_g(x, nothing) - mean_size) + d.penalty for x in X]
        size_diff = (N - M)
        return [C [x for x ∈ pentalty_vec, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::AvgSizeMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    C = pairwise_inbounds(d_g, X, Y)
    if N < M
        mean_size = mean(x -> d_g(x, nothing), X)
        penalty_vec = [abs(d_g(y, nothing) - mean_size) + d.penalty for y in Y]
        size_diff = (M - N)
        return [C; [y for i in 1:size_diff, y ∈ penalty_vec]]
    elseif N > M
        mean_size = mean(y -> d_g(y, nothing), Y)
        penalty_vec = [abs(d_g(x, nothing) - mean_size) + d.penalty for x in X]
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, i in 1:size_diff]]
    else
        return pairwise_inbounds(d_g, X, Y)
    end

end

function (d::AvgSizeMatchDist)(X::Nothing, Y::Vector{T})::Float64 where {T}
    return (d.penalty * length(Y)) + sum(x -> dist.ground_dist(x, nothing), Y)
end
function (d::AvgSizeMatchDist)(X::Vector{T}, Y::Nothing)::Float64 where {T}
    d(Y, X)
end
function (d::AvgSizeMatchDist)(X::Nothing, Y::Nothing)::Float64 where {T}
    return 0.0
end


struct MinDistMatchingDistance{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    penalty::Float64
end

const MinDistMatchDist{T} = MinDistMatchingDistance where {T}

Base.show(io::IO, d::MinDistMatchDist) = print(io, "$(typeof(d))(ρ=$(d.penalty))")


function get_cost_matrix_dynamic(
    d::MinDistMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    if N < M
        return get_cost_matrix_dynamic(d, Y, X)
    elseif N == M
        return pairwise_inbounds(d_g, X, Y)
    else
        C = pairwise_inbounds(d_g, X, Y)
        penalty_vec = map(minimum, eachrow(C))
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, j = 1:size_diff]]
    end
end

function get_cost_matrix_fixed(
    d::MinDistMatchDist,
    X::Vector{T}, Y::Vector{T}
) where {T}

    d_g = d.ground_dist
    N, M = (length(X), length(Y))
    C = pairwise_inbounds(d_g, X, Y)
    if N < M
        penalty_vec = map(minimum, eachcol(C))
        size_diff = (M - N)
        return [C; [y for i in 1:size_diff, y ∈ penalty_vec]]
    elseif N > M
        penalty_vec = map(minimum, eachrow(C))
        size_diff = (N - M)
        return [C [x for x ∈ penalty_vec, i in 1:size_diff]]
    else
        return pairwise_inbounds(d_g, X, Y)
    end

end

# Define type union for all these complete matching distances 

const CompleteMatchingDistance = Union{MatchingDistance,FastMatchingDistance,FixedPenaltyMatchingDistance,AvgSizeMatchingDistance,MinDistMatchingDistance}
