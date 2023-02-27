using Distances
export PrecomputedDistance, pre_compute

"""
Distance which just looks-up a cache of stored distances.

"""
struct PrecomputedDistance{T<:Any} <: SemiMetric
    cache::Dict{Tuple{T,T},Float64}
end

function (d::PrecomputedDistance{T})(
    args...
) where {T}
    return d.cache[args]
end

"""
Constructor method which takes distance + collection of values and makes a pre-computed distance.
"""
function pre_compute(d::SemiMetric, data::Vector{T}) where {T}

    cache = Dict{Tuple{T,T},Float64}()
    for (a, b) in Base.Iterators.product(data, data)
        cache[(a, b)] = d(a, b)
    end
    return PrecomputedDistance(cache)

end
