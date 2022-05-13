

struct EditDistance{T<:Metric} <: Metric
    ground_dist::T
end


function (d::EditDistance)(S1::Vector{T}, S2::Vector{T}) where {T}
    if length(S1) < length(S2)  # This ensures first seq is longest
        d(S2, S1)
    else
        d₀ = d.ground_dist
        prev_row = pushfirst!(cumsum([d₀(nothing, p) for p in S2]), 0.0);
        curr_row = zeros(Float64, length(S2) + 1);

        for i = 1:length(S1)
            curr_row[1] = prev_row[1] + length(S1[i])
            for j = 1:(length(S2))
                # @show i, j, prev_row[j], d.ground_dist(S1[i], S2[j])
                curr_row[j+1] = min(prev_row[j] + d₀(S1[i], S2[j]),
                                        prev_row[j+1] + d₀(nothing, S1[i]),
                                        curr_row[j] + d₀(nothing, S2[j]))
            end
            # @show curr_row
            copy!(prev_row, curr_row)
        end
        return curr_row[end]
    end
end