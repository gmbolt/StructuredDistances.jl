using StatsBase
export get_aggregate_adj_mat

function StatsBase.counts(data::Vector{T}, ref::Vector{T}) where {T}
    mapper = Dict(val=>i for (i,val) in enumerate(ref))
    out = zeros(Int, length(ref))
    for x in data 
        out[mapper[x]] += 1
    end 
    return out
end 

function get_aggregate_adj_mat(data::Vector{Vector{T}}, ref::Vector{T}) where {T}
    mapper = Dict(val=>i for (i,val) in enumerate(ref))
    A = zeros(Int, length(ref), length(ref))
    for path in data
        for i in Iterators.rest(eachindex(path),1)
            A[mapper[path[i-1]], mapper[path[i]]] += 1
        end 
    end 
    return A
end 
