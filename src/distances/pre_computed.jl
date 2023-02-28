using Distances
export pre_compute, Precomputed

"""
Distance which just looks-up a cache of stored distances.

"""
struct Precomputed{T<:SemiMetric,S} <: SemiMetric
    d::T
    cache::Dict{Tuple{Union{Nothing,S},Union{Nothing,S}},Float64}
end

function Base.show(io::IO, d_pre::Precomputed{T,S}) where {T<:SemiMetric,S}
    b = IOBuffer()
    show(b, d_pre.d)
    d_str = String(take!(b))
    print(io, "Precomputed{$(d_str),$(S)}")
end

# Distance eval is just a look-up of the cache
function (d::Precomputed{T,S})(
    args...
) where {T<:SemiMetric,S}
    return d.cache[args]
end

function Precomputed(d::SemiMetric, data::Vector{S}) where {S}
    TypeVal = Union{S,Nothing} # We include distances to nothing.
    cache = Dict{Tuple{TypeVal,TypeVal},Float64}()
    for (a, b) in Base.Iterators.product([data; nothing], [data; nothing])
        cache[(a, b)] = d(a, b)
    end
    return Precomputed(d, cache)
end

"""
Constructor method which takes distance + collection of values and makes a pre-computed distance.
"""
function pre_compute(d::SemiMetric, data::Vector{T}) where {T}

    TypeVal = Union{T,Nothing}
    cache = Dict{Tuple{TypeVal,TypeVal},Float64}()
    for (a, b) in Base.Iterators.product([data; nothing], [data; nothing])
        cache[(a, b)] = d(a, b)
    end
    return PrecomputedDistance(cache)

end
