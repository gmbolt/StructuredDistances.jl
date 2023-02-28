using Distances

export Normalised
"""
Normalised distance obtained by applying Steinhaus transform to given distance. 

Note this will be a metric given `d` is a metric.
"""
struct Normalised{T<:SemiMetric} <: SemiMetric
    d::T
end


function Base.show(io::IO, d_n::Normalised{T}) where {T<:SemiMetric}
    b = IOBuffer()
    show(b, d_n.d)
    d_str = String(take!(b))
    print(io, "Normalised{$(d_str)}")
end

function (d_n::Normalised{T})(x, y) where {T<:SemiMetric}
    d = d_n.d
    d_tmp = d(x, y)
    return 2 * d_tmp / (d(x, nothing) + d(y, nothing) + d_tmp)
end

(d_n::Normalised{T} where {T<:SemiMetric})(x::Nothing, y::Nothing) = 0.0 # To avoid NaN