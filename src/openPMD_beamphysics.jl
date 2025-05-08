module openPMD_beamphysics

# External dependencies
using HDF5
using Plots
using Statistics
using LinearAlgebra
using Distributions
using Dates
using TimeZones
using DelimitedFiles
using MAT
using FFTW
using Interpolations
using Unitful
using DifferentialEquations
using Optim

# Include all source files in dependency order
include("units.jl")
include("labels.jl")
include("status.jl")
include("species.jl")
include("tools.jl")
include("statistics.jl")
include("particles.jl")
include("readers.jl")
include("writers.jl")
include("plot.jl")
include("fields/fieldmesh.jl")
include("fields/analysis.jl")
include("fields/corrector_modeling.jl")
include("fields/conversion.jl")
include("fields/expansion.jl")

# Version info
const VERSION = v"0.1.0"

end # module openPMD_beamphysics
