using Distances

export Normalised

struct Normalised{T<:SemiMetric} <: SemiMetric
    d::T
end

function (d_n::Normalised{T})(x, y) where {T<:SemiMetric}
    d = d_n.d
    d_tmp = d(x, y)
    return 2 * d_tmp / (d(x, nothing) + d(y, nothing) + d_tmp)
end

(d_n::Normalised{T} where {T<:SemiMetric})(x::Nothing, y::Nothing) = 0.0 # To avoid NaN