module StructuredDistances

using Distances 

# Write your package code here.

include("utils.jl")
include("threaded_pairwise.jl")

include("Distances/helpers.jl")
include("Distances/paths.jl")
include("Distances/multisets/matching_distances.jl")
include("Distances/multisets/matching_distances_generalised.jl")
include("Distances/multisets/matching_distances_eval.jl")
include("Distances/multisets/earth_movers_distance.jl")
include("Distances/obj_sequences.jl")
# include("Distances/obj_multisets.jl")
include("Distances/set_distances.jl")
include("Distances/graph_distances.jl")
include("Distances/normalised.jl")

include("Plot_recipes/seq_plot.jl")
include("Plot_recipes/pathseq_plot.jl")
# include("Plot_recipes/multigraph_plot.jl")

end
