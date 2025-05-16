module openPMD_beamphysics

# External dependencies
using HDF5
using Plots
using Statistics
using LinearAlgebra
using Distributions
using AtomicAndPhysicalConstants

# Include all source files in dependency order
include("particles.jl")

# Version info
const VERSION = v"0.1.0"

end # module openPMD_beamphysics
