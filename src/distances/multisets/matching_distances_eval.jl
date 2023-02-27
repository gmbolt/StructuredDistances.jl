using PythonOT

export get_info, print_info, get_info_deep

# Distance evaluation & summarisation
# -----------------------------------

# Note each matching-based distance must have two methods 
# 1. get_cost_matrix_dynamic() - this will swap the dimensions of the cost matrix depending on size of passed objects
# 2. get_cost_matrix_fixed() - this will keep the demension fixed, so that the dimensions of cost matrix will depend on size of passed object

Base.show(io::IO, d::CompleteMatchingDistance) = print(io, typeof(d))

function (d::T where {T<:CompleteMatchingDistance})(
    X::Vector{S}, Y::Vector{S}
) where {S}
    C = get_cost_matrix_dynamic(d, X, Y)
    x = ones(size(C, 1))
    out = PythonOT.emd2(
        x, x, C
    )
    return out
end

function Base.show(io::IO, d_gen::General{T}) where {T<:CompleteMatchingDistance}
    b = IOBuffer()
    show(b, d_gen.d)
    d_str = String(take!(b))
    print(io, "General{$(d_str)}")
end

function (d::General{T} where {T<:CompleteMatchingDistance})(
    X::Vector{S}, Y::Vector{S}
) where {S}
    C = get_cost_matrix_fixed(d, X, Y)
    return hungarian(C)[2]
end

const MatchingBasedDistance = Union{CompleteMatchingDistance,General{T}} where {T<:CompleteMatchingDistance}

function get_info(
    d::T,
    X::Vector{S}, Y::Vector{S}
) where {T<:MatchingBasedDistance,S}

    C = get_cost_matrix_fixed(d, X, Y)
    assignment, cost = hungarian(C)
    indx, indy = (Int[], Int[])
    for i in eachindex(X)
        j = assignment[i]
        if j ≤ length(Y)
            push!(indx, i)
            push!(indy, j)
        end
    end
    return indx, indy
end

function print_info(
    d::T,
    X::Vector{S}, Y::Vector{S}
) where {T<:MatchingBasedDistance,S}

    title = "Matching Summary"
    println(title)
    println("-"^length(title) * "\n")
    indx, indy = get_info(d, X, Y)
    max_len = maximum(map(z -> length(@sprintf("%s", z)), X[indx]))
    for (i, j) in zip(indx, indy)
        xtmp, ytmp = (X[i], Y[j])
        pad = max_len - length(@sprintf("%s", xtmp))
        println(" "^pad * "$xtmp", " → $ytmp")
    end
    println("\nUnmatched entries (first):")
    for i in eachindex(X)
        if i ∉ indx
            print("$(X[i]), ")
        end
    end
    println("\nUnmatched entries (second):")
    for i in eachindex(Y)
        if i ∉ indy
            print("$(Y[i]), ")
        end
    end
end

function get_info_deep(
    d::T,
    X::Vector{S}, Y::Vector{S}
) where {T<:MatchingBasedDistance,S}

    d_g = d.ground_dist

    indx, indy = get_info(d, X, Y)

    outx, outy = (
        [zeros(Bool, i) for i in length.(X)],
        [zeros(Bool, i) for i in length.(Y)]
    )

    for (i, j) in zip(indx, indy)
        outx[i], outy[j] = get_info(d_g, X[i], Y[j])
    end

    return outx, outy

end

