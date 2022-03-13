module StructuredDistances

using Distances 

# Write your package code here.

include("aliases.jl")
include("utils.jl")
include("threaded_pairwise.jl")

include("Distances/helpers.jl")
include("Distances/sequences.jl")
include("Distances/obj_sequences.jl")
include("Distances/obj_multisets.jl")
include("Distances/set_distances.jl")
include("Distances/graph_distances.jl")
include("Distances/normalised.jl")

include("Plot_recipes/sequences.jl")

end
