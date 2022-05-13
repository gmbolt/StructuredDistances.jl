export General
export get_info, print_info, get_info_deep
# General Matching
# ----------------

# * All matching distances find optimal COMPLETE matching by default 
# * To instead find general optimal matching one can construct the generalised distance 

struct General{T<:CompleteMatchingDistance} 
    d::T
end 

function get_cost_matrix_fixed(
    d_gen::General{S}, 
    X::Vector{T}, Y::Vector{T}
    ) where {S<:MatchDist, T}

    d_g = d_gen.d.ground_dist
    N, M = (length(X), length(Y))
    K = N + M
    C = zeros(K, K)
    for j in 1:K
        for i in 1:K
            if (i ≤ N) & (j ≤ M)
                C[i,j] = d_g(X[i],Y[j])
            elseif (i ≤ N) & (j > M)
                C[i,j] = d_g(X[i],nothing) # Penalty for non-inclusion of X entry
            elseif (i > N) & (j ≤ M)
                C[i,j] = d_g(Y[j],nothing) # Penalty for non-inclusion of Y entry
            end 
        end 
    end 
    return C
end 

function get_cost_matrix_fixed(
    d_gen::General{S}, 
    X::Vector{T}, Y::Vector{T}
    ) where {S<:FixPenMatchDist, T}

    d_g = d_gen.d.ground_dist
    penalty = d_gen.d.penalty
    N, M = (length(X), length(Y))
    K = N + M
    C = fill(penalty, K, K)
    for j in 1:K
        for i in 1:K
            if (i ≤ N) & (j ≤ M)
                C[i,j] = d_g(X[i],Y[j])
            elseif (i > N) & (j > M)
                C[i,j] = 0.0
            end 
        end 
    end 
    return C
end 

function get_cost_matrix_fixed(
    d_gen::General{S}, 
    X::Vector{T}, Y::Vector{T}
    ) where {S<:AvgSizeMatchDist, T} 

    d_g = d_gen.d.ground_dist
    penalty = d_gen.d.penalty
    N, M = (length(X), length(Y))
    K = N + M
    C = zeros(K, K)
    x_mean_size, y_mean_size = (
        mean(x->d_g(x,nothing), X),
        mean(y->d_g(y,nothing), Y)
    )
    for j in 1:K
        for i in 1:K
            if (i ≤ N) & (j ≤ M)
                C[i,j] = d_g(X[i],Y[j])
            elseif (i ≤ N) & (j > M)
                C[i,j] = abs(d_g(X[i],nothing)-y_mean_size) + penalty # Penalty for non-inclusion of X entry
            elseif (i > N) & (j ≤ M)
                C[i,j] = abs(d_g(Y[j],nothing)-x_mean_size) + penalty # Penalty for non-inclusion of Y entry
            end 
        end 
    end 
    return C

end 

function get_cost_matrix_fixed(
    d_gen::General{S}, 
    X::Vector{T}, Y::Vector{T}
    ) where {S<:MinDistMatchDist, T} 

    d_g = d_gen.d.ground_dist
    penalty = d_gen.d.penalty
    N, M = (length(X), length(Y))
    # C_v = view(C, 1:N, 1:M)
    C = pairwise_inbounds(d_g, X, Y)
    x_min_vec, y_min_vec = (
        map(minimum, eachrow(C)),
        map(minimum, eachcol(C))
    )
    C = [
        C [x + penalty for x ∈ x_min_vec, i in 1:N] ;
        [y + penalty for i in 1:M, y ∈ y_min_vec] zeros(M, N)
    ]
    return C

end 
