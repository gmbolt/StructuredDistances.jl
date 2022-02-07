module StructuredDistances

using Distances 

# Write your package code here.

include("aliases.jl")

include("Distances/helpers.jl")
include("Distances/interactions.jl")
include("Distances/sequences.jl")
include("Distances/multisets.jl")
include("Distances/graph_distances.jl")
include("Distances/normalised.jl")

end
