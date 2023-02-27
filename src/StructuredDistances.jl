module StructuredDistances

using Distances

# Write your package code here.

include("utils.jl")
include("distances/helpers.jl")

# Path distances 
include("distances/paths/paths.jl")

# Graph distances 
include("distances/graphs/graph_distances.jl")

# Multiset distances
include("distances/multisets/matching_distances.jl")
include("distances/multisets/matching_distances_generalised.jl")
include("distances/multisets/matching_distances_eval.jl")
include("distances/multisets/earth_movers_distance.jl")

# Sequence distances
include("distances/sequences/edit_distance.jl")
include("distances/sequences/dtw_distance.jl")

# Pre-computed distances
include("distances/pre_computed.jl")

# Basic set distances (Hamming + Jaccard)
include("distances/set_distances.jl")

# Normalised distances (normalises any metric)
include("distances/normalised.jl")

# Plotting recipies 
include("plot_recipes/seq_plot.jl")
include("plot_recipes/pathseq_plot.jl")

end
