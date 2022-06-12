using Distances, StatsBase, ProgressMeter
export pairwise_inbounds, pairwise_inbounds!, progress_pairwise


"""
In-place distance matrix calculation between elements of Vectors.
"""
function Distances.pairwise!(
    A::AbstractArray,
    metric::Metric,
    a::Vector{T} where T,
    b::Vector{T} where T
    )

    for j in 1:length(b)
        for i in 1:length(a)
            A[i,j] = metric(a[i],b[j])
        end
    end
end

function pairwise_inbounds!(
    A::AbstractArray,
    metric::Metric,
    a::Vector{T},
    b::Vector{T} 
    ) where {T}
    @inbounds begin
        for j in eachindex(b)
            for i in eachindex(a)
                A[i,j] = metric(a[i],b[j])
            end
        end
    end 
end

"""
In-place distance matrix calculation between elements of Vectors. (For SubArray)
"""
function Distances.pairwise!(
    A::SubArray,
    metric::Metric,
    a::Vector{T} where T,
    b::Vector{T} where T
    )

    for j in 1:length(b)
        for i in 1:length(a)
            A[i,j] = metric(a[i],b[j])
        end
    end
end

# function pairwise_inbounds!(
#     A::SubArray,
#     metric::Metric,
#     a::Vector{T},
#     b::Vector{T} 
#     ) where {T}
#     @inbounds begin
#         for j in eachindex(b)
#             for i in eachindex(a)
#                 A[i,j] = metric(a[i],b[j])
#             end
#         end
#     end 
# end

"""
In-place distance matrix calculation between elements of Vectors. (For SubArray)
"""
function Distances.pairwise!(
    A::AbstractMatrix,
    metric::Metric,
    a::Vector{T} where T,
    b::Vector{T} where T
    )

    for j in 1:length(b)
        for i in 1:length(a)
            A[i,j] = metric(a[i],b[j])
        end
    end
end


"""
Distance matrix calculation between elements of Vectors. This is a custom extension
of the function in the Distances.jl package to allow vectors of general type. The
function in Distances.jl is designed for univariate/multivariate data and so takes
as input either vectors or matrices (data points as rows).
"""
function Distances.pairwise(
    metric::SemiMetric,
    a::Vector{T},
    b::Vector{T} 
    ) where {T}
    D = Array{Float64,2}(undef, length(a), length(b))
    for j in 1:length(b)
        for i in 1:length(a)
            D[i,j] = metric(a[i],b[j])
        end
    end
    return D
end

function pairwise_inbounds(
    metric::SemiMetric,
    a::Vector{T},
    b::Vector{T} 
    ) where {T}
    D = Array{Float64,2}(undef, length(a), length(b))
    @inbounds begin
        for j in eachindex(b)
            for i in eachindex(a)
                D[i,j] = metric(a[i],b[j])
            end
        end
    end 
    return D
end

"""
Distance matrix calculation between elements of Vectors with progress bar.
"""
function progress_pairwise(
    d::SemiMetric,
    a::Vector{T}
    ) where {T}

    D = zeros(length(a), length(a))
    iter = Progress(Int(length(a)*(length(a)-1)/2), 1)
    for j in 1:length(a)
        for i in 1:(j-1)
            D[i,j] = d(a[i],a[j])
            next!(iter)
        end
    end
    D += D'
    return D

end
