module openPMD_beamphysics

# External dependencies
using HDF5
using Plots
using StatsBase: weights, fit, Histogram
using LinearAlgebra: det
using AtomicAndPhysicalConstants: Species, massof, chargeof, C_LIGHT, E_CHARGE
using Unitful
using Unitful.DefaultSymbols
import Statistics
import Statistics: mean, std
import StatsBase
import StatsBase: cov

# Include all source files in dependency order
include("utils.jl")
include("units.jl")
include("readers.jl")
include("particles.jl")
include("statistics.jl")
include("writers.jl")
include("labels.jl")
include("plot.jl")

# Export key types and constructors
export ParticleGroup, single_particle

# Derived property functions
export nalive, ndead, mass, species_charge, charge
export energy, kinetic_energy, momentum
export xp, yp, r, theta, pr, ptheta, Lz
export gamma, beta, beta_x, beta_y, beta_z
export norm_emit_x, norm_emit_y, norm_emit_4d
export x_bar, px_bar, y_bar, py_bar, Jx, Jy
export set_charge!, average_current

# Statistics & beam physics
export norm_emit_calc, twiss_calc, twiss_dispersion, twiss
export A_mat_calc, A_inverse_mat_calc, amplitude_calc
export particle_amplitude, normalized_particle_coordinate
export slice_statistics, delta, ptp

# Coordinates & operations
export in_z_coordinates, in_t_coordinates
export drift!, drift_to_t!, drift_to_z!, split_particles

# I/O
export particle_paths, field_paths, load_field_attrs
export is_constant_component, constant_component_value
export component_data, offset_component_name, particle_array, all_components, component_str
export pmd_init, pmd_field_init, write_component_data, write_pmd_bunch, write_pmd_field
export particle_data_dict, write_particle_group

# Units & display
export nice_scale_prefix, nice_array, pg_units, limits, plottable_array
export mathlabel, texlabel

# Plotting
export slice_plot, density_plot, marginal_plot, density_and_slice_plot

end # module openPMD_beamphysics
