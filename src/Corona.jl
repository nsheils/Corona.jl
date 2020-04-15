__precompile__()
module Corona

using Dates
using TimeSeries
using DataStructures

include("utils.jl")
include("DA.jl")
include("dataio.jl")
include("plot.jl")

end # module Corona
