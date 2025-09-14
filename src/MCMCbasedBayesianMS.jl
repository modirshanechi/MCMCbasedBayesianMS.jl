module MCMCbasedBayesianMS

using Random, Statistics, Distributions
using ConcreteStructs
using DataFrames
using Plots

include("bms_structs.jl")
include("posterior_mcmc.jl")
include("bor_sampling.jl")
include("full_inference.jl")
include("bms_plots.jl")

end # module MCMCbasedBayesianMS
