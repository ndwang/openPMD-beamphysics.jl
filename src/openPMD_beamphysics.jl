module openPMD_beamphysics

# External dependencies
using HDF5
using Plots
using StatsBase
using LinearAlgebra
using Distributions
using AtomicAndPhysicalConstants
@APCdef

# Include all source files in dependency order
include("statistics.jl")
include("particles.jl")

# Version info
const VERSION = v"0.1.0"

end # module openPMD_beamphysics
