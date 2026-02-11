"""
    ParticleGroup{T<:Real}

Represents a collection of particles with their properties and methods for analysis.

# Fields
- `x,y,z::Vector{T}`: positions in [m]
- `px,py,pz::Vector{T}`: momenta in [eV/c]
- `t::Vector{T}`: time in [s]
- `status::Vector{Int}`: particle status (1 for alive, 0 for dead)
- `weight::Vector{T}`: macro-particle charge in [C]
- `species::Species`: particle species
- `id::Vector{Int}`: optional particle IDs
"""

abstract type AbstractParticleGroup{T<:Real, I<:Integer} end

struct ParticleGroup{T<:Real, I<:Integer} <: AbstractParticleGroup{T,I}
    x::Vector{T}
    px::Vector{T}
    y::Vector{T}
    py::Vector{T}
    z::Vector{T}
    pz::Vector{T}
    t::Vector{T}
    status::Vector{I}
    weight::Vector{T}
    species::Species
    id::Vector{Int}
end

# Constructor with optional id
function ParticleGroup(
    x::Vector{T}, px::Vector{T}, y::Vector{T}, py::Vector{T},
    z::Vector{T}, pz::Vector{T}, t::Vector{T}, status::Vector{I},
    weight::Vector{T}, species::Species
) where {T<:Real, I<:Integer}
    n = length(x)
    id = collect(Int,1:n)
    return ParticleGroup(x, px, y, py, z, pz, t, status, weight, species, id)
end

"""
    ParticleGroupView

Non-owning view into a subset of particles from a `ParticleGroup` (or another
`AbstractParticleGroup`). All array fields are `SubArray`s referencing the
parent storage, so mutating a `ParticleGroupView` mutates the original data.
"""
struct ParticleGroupView{T<:Real, I<:Integer,
                         AT<:AbstractVector{T},
                         AI<:AbstractVector{I},
                         AId<:AbstractVector{Int}} <: AbstractParticleGroup{T,I}
    x::AT
    px::AT
    y::AT
    py::AT
    z::AT
    pz::AT
    t::AT
    status::AI
    weight::AT
    species::Species
    id::AId
end

function ParticleGroupView(pg::AbstractParticleGroup{T,I}, Iinds) where {T<:Real, I<:Integer}
    return ParticleGroupView(
        @view(pg.x[Iinds]), @view(pg.px[Iinds]),
        @view(pg.y[Iinds]), @view(pg.py[Iinds]),
        @view(pg.z[Iinds]), @view(pg.pz[Iinds]),
        @view(pg.t[Iinds]), @view(pg.status[Iinds]),
        @view(pg.weight[Iinds]), pg.species,
        @view(pg.id[Iinds]),
    )
end

Base.view(pg::AbstractParticleGroup, Iinds) = ParticleGroupView(pg, Iinds)

# Constructor for single particle
function single_particle(;
    x=0.0, px=0.0, y=0.0, py=0.0, z=0.0, pz=0.0,
    t=0.0, weight=1.0, status=1, species=Species("electron")
)
    return ParticleGroup(
        [x], [px], [y], [py], [z], [pz], [t],
        [status], [weight], species, [1]
    )
end

# Base methods
Base.length(pg::AbstractParticleGroup) = length(pg.x)

# ---- Standalone derived-property functions ----

nalive(pg::AbstractParticleGroup) = count(pg.status .== 1)
ndead(pg::AbstractParticleGroup) = length(pg) - nalive(pg)
mass(pg::AbstractParticleGroup) = massof(pg.species)
species_charge(pg::AbstractParticleGroup) = chargeof(pg.species)
charge(pg::AbstractParticleGroup) = sum(pg.weight)

momentum(pg::AbstractParticleGroup) = sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2)
energy(pg::AbstractParticleGroup) = sqrt.(pg.px.^2 .+ pg.py.^2 .+ pg.pz.^2 .+ mass(pg)^2)
kinetic_energy(pg::AbstractParticleGroup) = energy(pg) .- mass(pg)

xp(pg::AbstractParticleGroup) = pg.px ./ pg.pz
yp(pg::AbstractParticleGroup) = pg.py ./ pg.pz

r(pg::AbstractParticleGroup) = hypot.(pg.x, pg.y)
theta(pg::AbstractParticleGroup) = atan.(pg.y, pg.x)

function pr(pg::AbstractParticleGroup)
    th = theta(pg)
    return pg.px .* cos.(th) .+ pg.py .* sin.(th)
end

function ptheta(pg::AbstractParticleGroup)
    th = theta(pg)
    return -pg.px .* sin.(th) .+ pg.py .* cos.(th)
end

Lz(pg::AbstractParticleGroup) = pg.x .* pg.py .- pg.y .* pg.px

gamma(pg::AbstractParticleGroup) = energy(pg) ./ mass(pg)
beta(pg::AbstractParticleGroup) = momentum(pg) ./ energy(pg)
beta_x(pg::AbstractParticleGroup) = pg.px ./ energy(pg)
beta_y(pg::AbstractParticleGroup) = pg.py ./ energy(pg)
beta_z(pg::AbstractParticleGroup) = pg.pz ./ energy(pg)

# These call normalized_particle_coordinate (defined in statistics.jl, resolved at call time)
x_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, "x")
px_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, "px")
y_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, "y")
py_bar(pg::AbstractParticleGroup) = normalized_particle_coordinate(pg, "py")

Jx(pg::AbstractParticleGroup) = hypot.(x_bar(pg), px_bar(pg))
Jy(pg::AbstractParticleGroup) = hypot.(y_bar(pg), py_bar(pg))

# ---- DERIVED_PROPERTIES dict for getindex resolution ----

const DERIVED_PROPERTIES = Dict{String, Any}(
    "n_particle"     => pg -> length(pg),
    "n_alive"        => pg -> nalive(pg),
    "n_dead"         => pg -> ndead(pg),
    "mass"           => pg -> mass(pg),
    "species_charge" => pg -> species_charge(pg),
    "charge"         => pg -> charge(pg),
    "p"              => pg -> momentum(pg),
    "energy"         => pg -> energy(pg),
    "kinetic_energy" => pg -> kinetic_energy(pg),
    "xp"             => pg -> xp(pg),
    "yp"             => pg -> yp(pg),
    "r"              => pg -> r(pg),
    "theta"          => pg -> theta(pg),
    "pr"             => pg -> pr(pg),
    "ptheta"         => pg -> ptheta(pg),
    "Lz"             => pg -> Lz(pg),
    "gamma"          => pg -> gamma(pg),
    "beta"           => pg -> beta(pg),
    "beta_x"         => pg -> beta_x(pg),
    "beta_y"         => pg -> beta_y(pg),
    "beta_z"         => pg -> beta_z(pg),
    "norm_emit_x"    => pg -> norm_emit_x(pg),
    "norm_emit_y"    => pg -> norm_emit_y(pg),
    "norm_emit_4d"   => pg -> norm_emit_4d(pg),
    "x_bar"          => pg -> x_bar(pg),
    "px_bar"         => pg -> px_bar(pg),
    "y_bar"          => pg -> y_bar(pg),
    "py_bar"         => pg -> py_bar(pg),
    "Jx"             => pg -> Jx(pg),
    "Jy"             => pg -> Jy(pg),
    "average_current"=> pg -> average_current(pg),
)

# ---- set_charge! ----

function set_charge!(pg::AbstractParticleGroup, val)
    val > 0 || throw(ArgumentError("charge must be >0. This is used to weight the particles."))
    pg.weight .*= val / charge(pg)
end

# ---- Display ----

function Base.show(io::IO, pg::ParticleGroup)
    print(io, "ParticleGroup with $(length(pg)) particles (charge = $(charge(pg)) C)")
end

function Base.show(io::IO, ::MIME"text/plain", pg::ParticleGroup)
    print(io, "ParticleGroup with $(length(pg)) particles (charge = $(charge(pg)) C)")
end

# ---- Statistics: proper Base/Statistics extensions ----

function Base.minimum(pg::AbstractParticleGroup, key::AbstractString)
    return minimum(_resolve_key(pg, key))
end

function Base.maximum(pg::AbstractParticleGroup, key::AbstractString)
    return maximum(_resolve_key(pg, key))
end

function Statistics.mean(pg::AbstractParticleGroup, key::AbstractString)
    dat = _resolve_key(pg, key)
    if length(dat) == 1
        return dat[1]
    end
    return mean(dat, weights(pg.weight))
end

function Statistics.std(pg::AbstractParticleGroup, key::AbstractString)
    dat = _resolve_key(pg, key)
    if length(dat) == 1
        return zero(eltype(dat))
    end
    avg_dat = Statistics.mean(pg, key)
    return sqrt(mean((dat .- avg_dat).^2, weights(pg.weight)))
end

function StatsBase.cov(pg::AbstractParticleGroup, keys::AbstractString...)
    dats = hcat([_resolve_key(pg, key) for key in keys]...)
    return StatsBase.cov(dats, weights(pg.weight), 1)
end

function delta(pg::AbstractParticleGroup, key::AbstractString)
    return _resolve_key(pg, key) .- Statistics.mean(pg, key)
end

ptp(x::AbstractArray) = maximum(x) - minimum(x)

function ptp(pg::AbstractParticleGroup, key::AbstractString)
    return ptp(_resolve_key(pg, key))
end

# Internal helper to resolve a key string to data array (for stat methods)
function _resolve_key(pg::AbstractParticleGroup, key::AbstractString)
    # Check struct fields first
    s = Symbol(key)
    if s in fieldnames(ParticleGroup)
        return getfield(pg, s)
    end
    # Check derived properties
    if haskey(DERIVED_PROPERTIES, key)
        return DERIVED_PROPERTIES[key](pg)
    end
    error("Unknown key: $key")
end

# ---- Coordinate checks ----

function in_z_coordinates(pg::AbstractParticleGroup)
    return all(==(pg.z[1]), pg.z)
end

function in_t_coordinates(pg::AbstractParticleGroup)
    return all(==(pg.t[1]), pg.t)
end

# ---- average_current ----

function average_current(pg::AbstractParticleGroup)
    dt = ptp(pg, "t")
    if dt == 0
        dt = ptp(pg, "z") / (Statistics.mean(pg, "beta_z") * C_LIGHT)
    end
    return charge(pg) / dt
end

# ---- Drift methods ----

function drift!(pg::AbstractParticleGroup, delta_t::Number)
    E = energy(pg)
    pg.x .+= (pg.px ./ E) .* C_LIGHT .* delta_t
    pg.y .+= (pg.py ./ E) .* C_LIGHT .* delta_t
    pg.z .+= (pg.pz ./ E) .* C_LIGHT .* delta_t
    pg.t .+= delta_t
end

function drift!(pg::AbstractParticleGroup, delta_t::AbstractVector)
    length(delta_t) == length(pg.x) || throw(ArgumentError("delta_t length must match number of particles"))
    E = energy(pg)
    pg.x .+= (pg.px ./ E) .* C_LIGHT .* delta_t
    pg.y .+= (pg.py ./ E) .* C_LIGHT .* delta_t
    pg.z .+= (pg.pz ./ E) .* C_LIGHT .* delta_t
    pg.t .+= delta_t
end

function drift_to_z!(pg::AbstractParticleGroup, z::Real)
    E = energy(pg)
    dt = (z .- pg.z) ./ ((pg.pz ./ E) .* C_LIGHT)
    drift!(pg, dt)
    pg.z .= z
end

function drift_to_z!(pg::AbstractParticleGroup)
    z = Statistics.mean(pg, "z")
    drift_to_z!(pg, z)
end

function drift_to_t!(pg::AbstractParticleGroup, t::Real)
    dt = t .- pg.t
    drift!(pg, dt)
    pg.t .= t
end

function drift_to_t!(pg::AbstractParticleGroup)
    t = Statistics.mean(pg, "t")
    drift_to_t!(pg, t)
end

# ---- getindex ----

"""
    Base.getindex(pg::ParticleGroup, x)

Returns a property or statistical quantity that can be computed:

- `P["x"]` returns the x array
- `P["sigma_x"]` returns the std(x) scalar
- `P["norm_emit_x"]` returns the norm_emit_x scalar
- `P["cov_x__px"]` returns the covariance of x and px

Parts can also be given. Example: `P[1:10]` returns a new ParticleGroup with the first 10 elements.
"""
function Base.getindex(pg::ParticleGroup, x)
    if !isa(x, AbstractString)
        return particle_parts(pg, x)
    end

    # Special case for z/c
    if x == "z/c"
        return pg.z ./ C_LIGHT
    end

    if startswith(x, "cov_")
        subkeys = split(x[5:end], "__")
        @assert length(subkeys) == 2 "Need 2 properties in covariance request: $x"
        return StatsBase.cov(pg, subkeys...)[1, 2]
    elseif startswith(x, "delta_")
        return delta(pg, x[7:end])
    elseif startswith(x, "sigma_")
        return Statistics.std(pg, x[7:end])
    elseif startswith(x, "mean_")
        return Statistics.mean(pg, x[6:end])
    elseif startswith(x, "min_")
        return Base.minimum(pg, x[5:end])
    elseif startswith(x, "max_")
        return Base.maximum(pg, x[5:end])
    elseif startswith(x, "ptp_")
        return ptp(pg, x[5:end])
    end

    # Check struct fields
    s = Symbol(x)
    if s in fieldnames(ParticleGroup)
        return getfield(pg, s)
    end

    # Check derived properties
    if haskey(DERIVED_PROPERTIES, x)
        return DERIVED_PROPERTIES[x](pg)
    end

    error("Unknown key: $x")
end

# ---- Arithmetic and comparison ----

function Base.:+(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup)
    return join_particle_groups(pg1, pg2)
end

function Base.:(==)(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup)
    if nameof(pg1.species) != nameof(pg2.species)
        return false
    end
    for key in (:x, :px, :y, :py, :z, :pz, :t, :status, :weight, :id)
        if getfield(pg1, key) != getfield(pg2, key)
            return false
        end
    end
    return true
end

function Base.isapprox(pg1::AbstractParticleGroup, pg2::AbstractParticleGroup; kwargs...)
    if nameof(pg1.species) != nameof(pg2.species)
        return false
    end
    for key in (:x, :px, :y, :py, :z, :pz, :t, :weight)
        if !isapprox(getfield(pg1, key), getfield(pg2, key); kwargs...)
            return false
        end
    end
    for key in (:status, :id)
        if getfield(pg1, key) != getfield(pg2, key)
            return false
        end
    end
    return true
end

# ---- Helper functions ----

function split_particles(pg::AbstractParticleGroup; n_chunks=100, key="z")
    zlist = _resolve_key(pg, key)
    iz = sortperm(zlist)
    n = length(iz)
    n < n_chunks && throw(ArgumentError("Cannot split $n particles into $n_chunks chunks"))

    edges = round.(Int, range(1, n + 1, length=n_chunks + 1))
    plist = Vector{ParticleGroup}(undef, n_chunks)
    for i in 1:n_chunks
        chunk = iz[edges[i]:edges[i+1]-1]
        plist[i] = ParticleGroup(
            pg.x[chunk], pg.px[chunk], pg.y[chunk], pg.py[chunk],
            pg.z[chunk], pg.pz[chunk], pg.t[chunk], pg.status[chunk],
            pg.weight[chunk], pg.species, pg.id[chunk]
        )
    end

    return plist
end

function particle_parts(pg::ParticleGroup, x)
    return ParticleGroup(
        pg.x[x], pg.px[x], pg.y[x], pg.py[x],
        pg.z[x], pg.pz[x], pg.t[x], pg.status[x],
        pg.weight[x], pg.species, pg.id[x]
    )
end

function join_particle_groups(pgs::AbstractParticleGroup...)
    species = [pg.species for pg in pgs]
    species0 = species[1]
    @assert all(spe == species0 for spe in species) "species must be the same to join"

    x = reduce(vcat, (pg.x for pg in pgs))
    px = reduce(vcat, (pg.px for pg in pgs))
    y = reduce(vcat, (pg.y for pg in pgs))
    py = reduce(vcat, (pg.py for pg in pgs))
    z = reduce(vcat, (pg.z for pg in pgs))
    pz = reduce(vcat, (pg.pz for pg in pgs))
    t = reduce(vcat, (pg.t for pg in pgs))
    status = reduce(vcat, (pg.status for pg in pgs))
    weight = reduce(vcat, (pg.weight for pg in pgs))
    id = reduce(vcat, (pg.id for pg in pgs))

    return ParticleGroup(
        x, px, y, py, z, pz, t, status, weight, species0, id
    )
end

# ---- HDF5 constructors ----

"""
    ParticleGroup(h5::HDF5.Group)

Load particles into a ParticleGroup from an HDF5 file in openPMD format.
"""
function ParticleGroup(h5::HDF5.Group)

    if !haskey(h5, "position")
        species = keys(h5)
        if length(species) != 1
            error("multiple species in particle paths: $(species)")
        end
        h5 = h5[species[1]]
    end

    attributes = attrs(h5)
    species = Species(attributes["speciesType"])
    n_particle = Int(attributes["numParticles"])
    total_charge = attributes["totalCharge"] * attributes["chargeUnitSI"]

    x = particle_array(h5, "x")
    px = particle_array(h5, "px")
    y = particle_array(h5, "y")
    py = particle_array(h5, "py")
    z = particle_array(h5, "z")
    pz = particle_array(h5, "pz")
    t = particle_array(h5, "t")

    if haskey(h5, "particleStatus")
        status = particle_array(h5, "particleStatus")
    else
        status = fill(1, n_particle)
    end

    if haskey(h5, "weight")
        weight = particle_array(h5, "weight")
    else
        weight = fill(total_charge / n_particle, n_particle)
    end

    return ParticleGroup(x, px, y, py, z, pz, t, status, weight, species)
end

"""
    ParticleGroup(file::String)

Load particles into a ParticleGroup from an HDF5 file in openPMD format.
"""
function ParticleGroup(file::String)
    h5open(file, "r") do h5
        pp = particle_paths(h5)
        @assert length(pp) == 1 "Number of particle paths in $file is $(length(pp))."
        return ParticleGroup(h5[pp[1]])
    end
end
