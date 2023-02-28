using StatsBase, Distances, Printf

export DTW, FixedPenaltyDTW, FixPenDTW

struct DTW{T<:SemiMetric} <: SemiMetric
    ground_dist::T
end

function (d::DTW)(
    S1::Vector{T}, S2::Vector{T}
) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d_g = d.ground_dist
        prev_row = pushfirst!(fill(Inf, length(S2)), 0.0)
        curr_row = fill(Inf, length(S2) + 1)

        for i in eachindex(S1)
            for j in eachindex(S2)
                cost = d_g(S1[i], S2[j])
                curr_row[j+1] = cost + min(prev_row[j], prev_row[j+1], curr_row[j])
            end
            copy!(prev_row, curr_row)
        end
        return curr_row[end]
    end
end

function (d::DTW)(
    S1::Nothing, S2::Vector{T}
) where {T}

    d_g = d.ground_dist
    z = 0.0
    for x in S2
        z += d_g(x, nothing)
    end
    return z + length(S2)

end
(d::DTW)(S1::Vector{T}, S2::Nothing) where {T} = d(S2, S1)
(d::DTW)(S1::Nothing, S2::Nothing) where {T} = 0.0


function print_info(
    d::DTW,
    S1::Vector{T}, S2::Vector{T}
) where {T}

    d = d.ground_dist
    # First find the substitution matrix
    C = fill(Inf, length(S1) + 1, length(S2) + 1)
    C[1, 1] = 0.0

    for j in 1:length(S2)
        for i in 1:length(S1)
            cost = d(S1[i], S2[j])
            C[i+1, j+1] = cost + min(C[i+1, j], C[i, j+1], C[i, j])
        end
    end
    # Now retrace steps to determine an optimal matching
    i, j = size(C)
    pairs = Tuple{Int,Int}[]
    pushfirst!(pairs, (i - 1, j - 1))
    while (i ≠ 2) | (j ≠ 2)
        i_tmp, j_tmp = argmin(view(C, (i-1):i, (j-1):j)).I
        i = i - 2 + i_tmp
        j = j - 2 + j_tmp
        pushfirst!(pairs, (i - 1, j - 1))
    end
    # @show outputs
    title = "Optimal Coupling"
    println(title)
    println("-"^length(title), "\n")
    # for statement in outputs
    #     println(statement)
    # end
    i_tmp, j_tmp = (0, 0)
    for (i, j) in pairs
        if i == i_tmp
            println(" "^length(@sprintf("%s", S1[i])) * " ↘ $(S2[j])")
        elseif j == j_tmp
            println("$(S1[i]) ↗ " * " "^length(@sprintf("%s", S2[j])))
        else
            println("$(S1[i]) → $(S2[j])")
        end
        i_tmp, j_tmp = (i, j)
    end

end

# Penalised 
# ---------

struct FixedPenaltyDTW{T<:SemiMetric} <: SemiMetric
    ground_dist::T
    ρ::Real
end

const FixPenDTW = FixedPenaltyDTW

function (d::FixPenDTW)(
    S1::Vector{T}, S2::Vector{T}
) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d_g = d.ground_dist
        prev_row = pushfirst!(fill(Inf, length(S2)), 0.0)
        curr_row = fill(Inf, length(S2) + 1)

        for i in eachindex(S1)
            for j in eachindex(S2)
                cost = d_g(S1[i], S2[j])
                curr_row[j+1] = cost + min(
                    prev_row[j], # New pair (no warping)
                    prev_row[j+1] + d.ρ, # Warping 
                    curr_row[j] + d.ρ)  # Warping
            end
            copy!(prev_row, curr_row)
        end
        return curr_row[end]
    end
end

function (d::FixPenDTW)(
    S1::Nothing, S2::Vector{T}
) where {T}

    d_g = d.ground_dist
    z = 0.0
    for x in S2
        z += d_g(x, nothing)
    end
    return z + length(S2) * d.ρ

end
(d::FixPenDTW)(S1::Vector{T}, S2::Nothing) where {T} = d(S2, S1)
(d::FixPenDTW)(S1::Nothing, S2::Nothing) where {T} = 0.0#


function print_info(
    d::Union{FixPenDTW},
    S1::Vector{T}, S2::Vector{T}
) where {T}

    d_g = d.ground_dist
    # First find the substitution matrix
    C = fill(Inf, length(S1) + 1, length(S2) + 1)
    C[1, 1] = 0.0

    for j in 1:length(S2)
        for i in 1:length(S1)
            cost = d_g(S1[i], S2[j])
            C[i+1, j+1] = cost + min(C[i+1, j] + d.ρ, C[i, j+1] + d.ρ, C[i, j])
        end
    end
    # Now retrace steps to determine an optimal matching
    i, j = size(C)
    pairs = Tuple{Int,Int}[]
    warp_shift_mat = [0 d.ρ; d.ρ 0] # Extra bit here different from DTW
    pushfirst!(pairs, (i - 1, j - 1))
    while (i ≠ 2) | (j ≠ 2)
        i_tmp, j_tmp = argmin(view(C, (i-1):i, (j-1):j) + warp_shift_mat).I
        i = i - 2 + i_tmp
        j = j - 2 + j_tmp
        pushfirst!(pairs, (i - 1, j - 1))
    end
    # @show outputs
    title = "Fixed Penalty DTW Print-out with $d_g Ground Distance and ρ=$(d.ρ)"
    println(title)
    println("-"^length(title), "\n")
    println("The cheapest way to do the tranformation...\n")
    println(S1, "---->", S2)
    println("\n...is the following series of edits...\n")
    # for statement in outputs
    #     println(statement)
    # end
    i_tmp, j_tmp = (0, 0)
    for (i, j) in pairs
        if i == i_tmp
            println(" "^length(@sprintf("%s", S1[i])) * " ↘ $(S2[j])")
        elseif j == j_tmp
            println("$(S1[i]) ↗ " * " "^length(@sprintf("%s", S2[j])))
        else
            println("$(S1[i]) → $(S2[j])")
        end
        i_tmp, j_tmp = (i, j)
    end
end