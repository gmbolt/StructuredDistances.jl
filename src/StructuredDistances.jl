module StructuredDistances

using Distances 

# Write your package code here.

include("utils.jl")
include("threaded_pairwise.jl")

include("distances/helpers.jl")
include("distances/paths.jl")
include("distances/multisets/matching_distances.jl")
include("distances/multisets/matching_distances_generalised.jl")
include("distances/multisets/matching_distances_eval.jl")
include("distances/multisets/earth_movers_distance.jl")
include("distances/obj_sequences.jl")
# include("Distances/obj_multisets.jl")
include("distances/set_distances.jl")
include("distances/graph_distances.jl")
include("distances/normalised.jl")

include("plot_recipes/seq_plot.jl")
include("plot_recipes/pathseq_plot.jl")
# include("plot_recipes/multigraph_plot.jl")

end
