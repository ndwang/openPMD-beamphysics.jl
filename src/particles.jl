"""
    ParticleGroup

Represents a collection of particles with their properties and methods for analysis.

# Fields
- `data::Dict`: Internal data dictionary containing particle properties
- `_settable_array_keys::Vector{String}`: Keys for array properties that can be set
- `_settable_scalar_keys::Vector{String}`: Keys for scalar properties that can be set
"""
mutable struct ParticleGroup
    data::Dict
    _settable_array_keys::Vector{String}
    _settable_scalar_keys::Vector{String}

    function ParticleGroup(h5=nothing, data=nothing)
        if !isnothing(h5) && !isnothing(data)
            error("Cannot initialize with both h5 and data")
        end

        if !isnothing(h5)
            # Handle file path
            if isa(h5, String) || isa(h5, Path)
                fname = expandvars(h5)
                @assert isfile(fname) "File does not exist: $fname"
                
                h5open(fname, "r") do hh5
                    pp = particle_paths(hh5)
                    @assert length(pp) == 1 "Number of particle paths in $h5: $(length(pp))"
                    data = load_bunch_data(hh5[pp[1]])
                end
            else
                # Try dict
                data = load_bunch_data(h5)
            end
        else
            # Fill out data. Exclude species.
            data = full_data(data)
            species = unique(data["species"])

            # Allow for empty data (len=0). Otherwise, check species.
            if length(species) >= 1
                @assert length(species) == 1 "mixed species are not allowed: $species"
                data["species"] = species[1]
            end
        end

        settable_array_keys = ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
        # Optional data
        if "id" in keys(data)
            push!(settable_array_keys, "id")
        end

        settable_scalar_keys = ["species"]
        settable_keys = vcat(settable_array_keys, settable_scalar_keys)

        # Internal data. Only allow settable keys
        internal_data = Dict(k => data[k] for k in settable_keys)

        new(internal_data, settable_array_keys, settable_scalar_keys)
    end
end

# Property accessors
for prop in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight"]
    @eval begin
        function Base.getproperty(pg::ParticleGroup, s::Symbol)
            if s == $(QuoteNode(Symbol(prop)))
                return pg.data[$prop]
            end
            return getfield(pg, s)
        end

        function Base.setproperty!(pg::ParticleGroup, s::Symbol, val)
            if s == $(QuoteNode(Symbol(prop)))
                pg.data[$prop] = full_array(length(pg), val)
                return
            end
            setfield!(pg, s, val)
        end
    end
end

# Special property accessors
function Base.getproperty(pg::ParticleGroup, s::Symbol)
    if s == :id
        if !("id" in keys(pg.data))
            assign_id(pg)
        end
        return pg.data["id"]
    elseif s == :species
        return pg.data["species"]
    elseif s == :data
        return pg.data
    end
    return getfield(pg, s)
end

function Base.setproperty!(pg::ParticleGroup, s::Symbol, val)
    if s == :id
        pg.data["id"] = full_array(length(pg), val)
    elseif s == :species
        pg.data["species"] = val
    elseif s == :data
        pg.data = val
    else
        setfield!(pg, s, val)
    end
end

# Derived properties
function Base.getproperty(pg::ParticleGroup, s::Symbol)
    if s == :n_particle
        return length(pg)
    elseif s == :n_alive
        return count(pg.status .== 1)
    elseif s == :n_dead
        return pg.n_particle - pg.n_alive
    elseif s == :mass
        return mass_of(pg.species)
    elseif s == :species_charge
        return charge_of(pg.species)
    elseif s == :charge
        return sum(pg.weight)
    elseif s == :p
        return sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2)
    elseif s == :energy
        return sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2 .+ pg.mass^2)
    elseif s == :kinetic_energy
        return pg.energy .- pg.mass
    elseif s == :xp
        return pg.px ./ pg.pz
    elseif s == :yp
        return pg.py ./ pg.pz
    elseif s == :r
        return hypot.(pg.x, pg.y)
    elseif s == :theta
        return atan.(pg.y, pg.x)
    elseif s == :pr
        theta = pg.theta
        return pg.px .* cos.(theta) .+ pg.py .* sin.(theta)
    elseif s == :ptheta
        theta = pg.theta
        return -pg.px .* sin.(theta) .+ pg.py .* cos.(theta)
    elseif s == :Lz
        return pg.x .* pg.py .- pg.y .* pg.px
    elseif s == :gamma
        return pg.energy ./ pg.mass
    elseif s == :beta
        return pg.p ./ pg.energy
    elseif s == :beta_x
        return pg.px ./ pg.energy
    elseif s == :beta_y
        return pg.py ./ pg.energy
    elseif s == :beta_z
        return pg.pz ./ pg.energy
    end
    return getfield(pg, s)
end

# Methods
function assign_id(pg::ParticleGroup)
    if !("id" in pg._settable_array_keys)
        push!(pg._settable_array_keys, "id")
    end
    pg.id = collect(1:pg.n_particle)
end

function units(pg::ParticleGroup, key::String)
    return pg_units(key)
end

function charge!(pg::ParticleGroup, val::Real)
    @assert val > 0 "charge must be >0. This is used to weight the particles."
    pg.weight .*= val / pg.charge
end

function higher_order_energy_calc(pg::ParticleGroup, order=2)
    if std(pg.z) < 1e-12
        # must be at a screen. Use t
        t = pg.t
    else
        # All particles at the same time. Use z to calc t
        t = pg.z ./ c_light
    end
    energy = pg.energy

    coeffs = polyfit(t, energy, order)
    best_fit = polyval(t, coeffs)
    return energy .- best_fit
end

# Statistical methods
function min(pg::ParticleGroup, key::String)
    return minimum(pg[key])
end

function max(pg::ParticleGroup, key::String)
    return maximum(pg[key])
end

function ptp(pg::ParticleGroup, key::String)
    return maximum(pg[key]) - minimum(pg[key])
end

function avg(pg::ParticleGroup, key::String)
    dat = pg[key]
    if isscalar(dat)
        return dat
    end
    return mean(dat, weights=pg.weight)
end

function std(pg::ParticleGroup, key::String)
    dat = pg[key]
    if isscalar(dat)
        return 0
    end
    avg_dat = avg(pg, key)
    return sqrt(mean((dat .- avg_dat).^2, weights=pg.weight))
end

function cov(pg::ParticleGroup, keys::String...)
    dats = [pg[key] for key in keys]
    return cov(dats, weights=pg.weight)
end

function histogramdd(pg::ParticleGroup, keys::String...; bins=10, range=nothing)
    H, edges = histogramdd(
        [pg[k] for k in keys],
        weights=pg.weight,
        bins=bins,
        range=range
    )
    return H, edges
end

# Beam statistics
function norm_emit_x(pg::ParticleGroup)
    return norm_emit_calc(pg, ["x"])
end

function norm_emit_y(pg::ParticleGroup)
    return norm_emit_calc(pg, ["y"])
end

function norm_emit_4d(pg::ParticleGroup)
    return norm_emit_calc(pg, ["x", "y"])
end

function twiss(pg::ParticleGroup, plane="x", fraction=1, p0c=nothing)
    d = Dict{String,Any}()
    for p in plane
        merge!(d, particle_twiss_dispersion(pg, plane=p, fraction=fraction, p0c=p0c))
    end
    return d
end

function twiss_match(pg::ParticleGroup, beta=nothing, alpha=nothing, plane="x", p0c=nothing, inplace=false)
    return matched_particles(pg, beta=beta, alpha=alpha, plane=plane, p0c=p0c, inplace=inplace)
end

# Coordinate checks
function in_z_coordinates(pg::ParticleGroup)
    return length(unique(pg.z)) == 1
end

function in_t_coordinates(pg::ParticleGroup)
    return length(unique(pg.t)) == 1
end

function average_current(pg::ParticleGroup)
    dt = ptp(pg, "t")
    if dt == 0
        # must be in t coordinates. Calc with
        dt = ptp(pg, "z") / (avg(pg, "beta_z") * c_light)
    end
    return pg.charge / dt
end

# Drift methods
function drift!(pg::ParticleGroup, delta_t::Real)
    pg.x .+= pg.beta_x .* c_light .* delta_t
    pg.y .+= pg.beta_y .* c_light .* delta_t
    pg.z .+= pg.beta_z .* c_light .* delta_t
    pg.t .+= delta_t
end

function drift_to_z!(pg::ParticleGroup, z=nothing)
    if isnothing(z)
        z = avg(pg, "z")
    end
    dt = (z .- pg.z) ./ (pg.beta_z .* c_light)
    drift!(pg, dt)
    # Fix z to be exactly this value
    pg.z .= z
end

function drift_to_t!(pg::ParticleGroup, t=nothing)
    if isnothing(t)
        t = avg(pg, "t")
    end
    dt = t .- pg.t
    drift!(pg, dt)
    # Fix t to be exactly this value
    pg.t .= t
end

# Helper functions
function full_array(n::Integer, val)
    if isscalar(val)
        return fill(val, n)
    end
    n_here = length(val)
    
    if n_here == 1
        return fill(val[1], n)
    elseif n_here != n
        error("Length mismatch: length(val)=$n_here, but requested n=$n")
    end
    return collect(val)
end

function full_data(data::Dict, exclude=nothing)
    full_data = Dict{String,Any}()
    scalars = Dict{String,Any}()
    
    for (k, v) in data
        if isscalar(v)
            scalars[k] = v
        elseif length(v) == 1
            scalars[k] = v[1]
        else
            # must be array
            full_data[k] = collect(v)
        end
    end

    # Check for single particle
    if isempty(full_data)
        return Dict(k => [v] for (k, v) in scalars)
    end

    # Array data should all have the same length
    nlist = [length(v) for (_, v) in full_data]
    @assert length(unique(nlist)) == 1 "arrays must have the same length. Found lengths: $(Dict(k => length(v) for (k, v) in full_data))"

    for (k, v) in scalars
        full_data[k] = fill(v, nlist[1])
    end

    return full_data
end

# Base methods
Base.length(pg::ParticleGroup) = length(pg.data[pg._settable_array_keys[1]])

function Base.show(io::IO, pg::ParticleGroup)
    print(io, "ParticleGroup with $(pg.n_particle) particles with total charge $(pg.charge) C")
end

function Base.show(io::IO, ::MIME"text/plain", pg::ParticleGroup)
    memloc = string(objectid(pg), base=16)
    print(io, "<ParticleGroup with $(pg.n_particle) particles at 0x$memloc>")
end

# Constructor for single particle
function single_particle(;
    x=0.0, px=0.0, y=0.0, py=0.0, z=0.0, pz=0.0,
    t=0.0, weight=1.0, status=1, species="electron"
)
    data = Dict(
        "x" => x, "px" => px, "y" => y, "py" => py,
        "z" => z, "pz" => pz, "t" => t,
        "weight" => weight, "status" => status,
        "species" => species
    )
    return ParticleGroup(data=data)
end

# Centroid calculation
function centroid(pg::ParticleGroup)
    good = pg.status .== 1
    pg_good = pg[good]
    data = Dict(key => avg(pg_good, key) for key in ["x", "px", "y", "py", "z", "pz", "t"])
    data["species"] = pg.species
    data["weight"] = pg.charge
    data["status"] = 1
    return ParticleGroup(data=data)
end

# Array operations
function Base.getindex(pg::ParticleGroup, x)
    if !isa(x, String)
        return particle_parts(pg, x)
    end

    # Special case for z/c
    if x == "z/c"
        return pg["z"] ./ c_light
    end

    if startswith(x, "cov_")
        subkeys = split(x[5:end], "__")
        @assert length(subkeys) == 2 "Too many properties in covariance request: $x"
        return cov(pg, subkeys...)[1, 2]
    elseif startswith(x, "delta_")
        return delta(pg, x[7:end])
    elseif startswith(x, "sigma_")
        return std(pg, x[7:end])
    elseif startswith(x, "mean_")
        return avg(pg, x[6:end])
    elseif startswith(x, "min_")
        return min(pg, x[5:end])
    elseif startswith(x, "max_")
        return max(pg, x[5:end])
    elseif startswith(x, "ptp_")
        return ptp(pg, x[5:end])
    elseif contains(x, "bunching")
        wavelength = parse_bunching_str(x)
        bunching = bunching(pg, wavelength)  # complex

        # abs or arg (angle):
        if contains(x, "phase_")
            return angle(bunching)
        else
            return abs(bunching)
        end
    else
        return getproperty(pg, Symbol(x))
    end
end

function Base.:+(pg1::ParticleGroup, pg2::ParticleGroup)
    return join_particle_groups(pg1, pg2)
end

function Base.:(==)(pg1::ParticleGroup, pg2::ParticleGroup)
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight", "id"]
        if !all(isapprox.(pg1[key], pg2[key]))
            return false
        end
    end
    return true
end

function Base.in(item, pg::ParticleGroup)
    return item in keys(pg.data)
end

# Helper functions
function split_particles(pg::ParticleGroup, n_chunks=100, key="z")
    # Sorting
    zlist = pg[key]
    iz = sortperm(zlist)

    # Split particles into chunks
    plist = ParticleGroup[]
    for chunk in Iterators.partition(iz, ceil(Int, length(iz)/n_chunks))
        # Prepare data
        data = Dict{String,Any}()
        for k in pg._settable_array_keys
            data[k] = getproperty(pg, Symbol(k))[chunk]
        end
        # These should be scalars
        data["species"] = pg.species

        # New object
        push!(plist, ParticleGroup(data=data))
    end

    return plist
end

function particle_parts(pg::ParticleGroup, x)
    data = Dict{String,Any}()
    for k in pg._settable_array_keys
        data[k] = getproperty(pg, Symbol(k))[x]
    end

    for k in pg._settable_scalar_keys
        data[k] = getproperty(pg, Symbol(k))
    end

    return ParticleGroup(data=data)
end

function join_particle_groups(pgs::ParticleGroup...)
    species = [pg.species for pg in pgs]
    species0 = species[1]
    @assert all(spe == species0 for spe in species) "species must be the same to join"

    data = Dict{String,Any}()
    for key in pgs[1]._settable_array_keys
        data[key] = vcat([pg[key] for pg in pgs]...)
    end

    data["species"] = species0
    data["n_particle"] = sum(pg.n_particle for pg in pgs)

    return ParticleGroup(data=data)
end

# Interface methods
function write_astra(pg::ParticleGroup, filePath::String; verbose=false, probe=false)
    return write_astra(pg, filePath, verbose=verbose, probe=probe)
end

function write_bmad(pg::ParticleGroup, filePath::String; p0c=nothing, t_ref=0, verbose=false)
    return write_bmad(pg, filePath, p0c=p0c, t_ref=t_ref, verbose=verbose)
end

function write_elegant(pg::ParticleGroup, filePath::String; verbose=false)
    return write_elegant(pg, filePath, verbose=verbose)
end

function write_genesis2_beam_file(pg::ParticleGroup, filePath::String; n_slice=nothing, verbose=false)
    # Get beam columns
    beam_columns = genesis2_beam_data(pg, n_slice=n_slice)
    # Actually write the file
    return write_genesis2_beam_file(filePath, beam_columns, verbose=verbose)
end

function write_genesis4_beam(pg::ParticleGroup, filePath::String; n_slice=nothing, return_input_str=false, verbose=false)
    return write_genesis4_beam(pg, filePath, n_slice=n_slice, return_input_str=return_input_str, verbose=verbose)
end

function write_genesis4_distribution(pg::ParticleGroup, filePath::String; verbose=false)
    return write_genesis4_distribution(pg, filePath, verbose=verbose)
end

function write_gpt(pg::ParticleGroup, filePath::String; asci2gdf_bin=nothing, verbose=false)
    return write_gpt(pg, filePath, asci2gdf_bin=asci2gdf_bin, verbose=verbose)
end

function write_impact(pg::ParticleGroup, filePath::String; cathode_kinetic_energy_ref=nothing, include_header=true, verbose=false)
    return write_impact(pg, filePath, cathode_kinetic_energy_ref=cathode_kinetic_energy_ref, include_header=include_header, verbose=verbose)
end

function write_litrack(pg::ParticleGroup, filePath::String; p0c=nothing, verbose=false)
    return write_litrack(pg, filePath, p0c=p0c, verbose=verbose)
end

function write_lucretia(pg::ParticleGroup, filePath::String; ele_name="BEGINNING", t_ref=0, stop_ix=nothing, verbose=false)
    return write_lucretia(pg, filePath, ele_name=ele_name, t_ref=t_ref, stop_ix=stop_ix)
end

function write_simion(pg::ParticleGroup, filePath::String; color=0, flip_z_to_x=true, verbose=false)
    return write_simion(pg, filePath, verbose=verbose, color=color, flip_z_to_x=flip_z_to_x)
end

function write_opal(pg::ParticleGroup, filePath::String; verbose=false, dist_type="emitted")
    return write_opal(pg, filePath, verbose=verbose, dist_type=dist_type)
end

function write(pg::ParticleGroup, h5, name=nothing)
    if isa(h5, Union{String,Path})
        fname = expandvars(h5)
        h5open(fname, "w") do g
            pmd_init(g, basePath="/", particlesPath="particles")
            g = create_group(g, "particles")
            write_pmd_bunch(g, pg, name=name)
        end
    else
        write_pmd_bunch(h5, pg, name=name)
    end
end

# Plotting methods
function plot(pg::ParticleGroup, key1="x", key2=nothing; bins=nothing, xlim=nothing, ylim=nothing, return_figure=false, tex=true, nice=true, ellipse=false, kwargs...)
    if isnothing(key2)
        fig = density_plot(pg, key=key1, bins=bins, xlim=xlim, tex=tex, nice=nice; kwargs...)
    else
        fig = marginal_plot(pg, key1=key1, key2=key2, bins=bins, xlim=xlim, ylim=ylim, tex=tex, nice=nice, ellipse=ellipse; kwargs...)
    end

    if return_figure
        return fig
    end
end

function slice_statistics(pg::ParticleGroup, keys::String...; n_slice=100, slice_key=nothing)
    if isnothing(slice_key)
        if in_t_coordinates(pg)
            slice_key = "z"
        else
            slice_key = "t"
        end
    end

    if slice_key in ("t", "delta_t")
        density_name = "current"
    else
        density_name = "density"
    end

    keys = Set(keys)
    push!(keys, "mean_" * slice_key)
    push!(keys, "ptp_" * slice_key)
    push!(keys, "charge")
    slice_dat = slice_statistics(pg, n_slice=n_slice, slice_key=slice_key, keys=keys)

    slice_dat[density_name] = slice_dat["charge"] / slice_dat["ptp_" * slice_key]

    return slice_dat
end

function slice_plot(pg::ParticleGroup, keys::String...; n_slice=100, slice_key=nothing, tex=true, nice=true, return_figure=false, xlim=nothing, ylim=nothing, kwargs...)
    fig = slice_plot(pg, keys..., n_slice=n_slice, slice_key=slice_key, tex=tex, nice=nice, xlim=xlim, ylim=ylim; kwargs...)

    if return_figure
        return fig
    end
end