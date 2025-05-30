module openPMD_beamphysics

# External dependencies
using HDF5
using Plots
using StatsBase
using LinearAlgebra
using Distributions
using AtomicAndPhysicalConstants
using Unitful
using Unitful.DefaultSymbols
@APCdef

# Include all source files in dependency order
include("units.jl")
include("readers.jl")
include("statistics.jl")
include("particles.jl")
include("labels.jl")
include("plot.jl")

# Export key types and functions
export particle_paths, field_paths, is_constant_component, constant_component_value
export component_data, offset_component_name, particle_array, all_components, component_str
export nice_scale_prefix, nice_array, pg_units
export load_field_attrs
export component_data, is_constant_component, constant_component_value
export norm_emit_calc, twiss_calc, twiss_dispersion
export ParticleGroup, single_particle, twiss, in_z_coordinates, in_t_coordinates, average_current
export drift!, drift_to_t!, drift_to_z!, split_particles
export limits, plottable_array
export mathlabel, texlabel


# Version info
const VERSION = v"0.1.0"

end # module openPMD_beamphysics
