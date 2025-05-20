"""
    ParticleGroup{T<:Real}

Represents a collection of particles with their properties and methods for analysis.

# Fields
- `x,y,z::Vector{T}`: positions in [m]
- `px,py,pz::Vector{T}`: momenta in [eV/c]
- `t::Vector{T}`: time in [s]
- `status::Vector{Int}`: particle status (1 for alive, 0 for dead)
- `weight::Vector{T}`: macro-particle charge in [C]
- `species::String`: particle species name
- `id::Vector{Int}`: optional particle IDs
"""
struct ParticleGroup{T<:Real}
    x::Vector{T}
    px::Vector{T}
    y::Vector{T}
    py::Vector{T}
    z::Vector{T}
    pz::Vector{T}
    t::Vector{T}
    status::Vector{Int}
    weight::Vector{T}
    species::Species
    id::Vector{Int}
end

# Constructor with optional id
function ParticleGroup(
    x::Vector{T}, px::Vector{T}, y::Vector{T}, py::Vector{T},
    z::Vector{T}, pz::Vector{T}, t::Vector{T}, status::Vector{Int},
    weight::Vector{T}, species::Species
) where T<:Real
    n = length(x)
    id = collect(Int,1:n)
    return ParticleGroup(x, px, y, py, z, pz, t, status, weight, species, id)
end

# Constructor for single particle
function single_particle(;
    x=0.0, px=0.0, y=0.0, py=0.0, z=0.0, pz=0.0,
    t=0.0, weight=1.0, status=1, species=Species("electron")
)
    T = typeof(x)
    return ParticleGroup{T}(
        [x], [px], [y], [py], [z], [pz], [t],
        [status], [weight], species, [1]
    )
end

# Base methods
Base.length(pg::ParticleGroup) = length(pg.x)

function norm_emit_x(pg::ParticleGroup)
    return norm_emit_calc(pg, ["x"])
end

function norm_emit_y(pg::ParticleGroup)
    return norm_emit_calc(pg, ["y"])
end

function norm_emit_4d(pg::ParticleGroup)
    return norm_emit_calc(pg, ["x", "y"])
end

# Derived properties
function Base.getproperty(pg::ParticleGroup, s::Symbol)
    if s == :n
        return length(pg)
    elseif s == :n_alive
        return count(pg.status .== 1)
    elseif s == :n_dead
        return pg.n - pg.n_alive
    elseif s == :mass
        return massof(pg.species)
    elseif s == :species_charge
        return chargeof(pg.species)
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
    elseif s == :norm_emit_x
        return norm_emit_x(pg)
    elseif s == :norm_emit_y
        return norm_emit_y(pg)
    elseif s == :norm_emit_4d
        return norm_emit_4d(pg)
    end
    return getfield(pg, s)
end

function setproperty!(pg::ParticleGroup, s::Symbol, val)
    if key == :charge
        @assert val > 0 "charge must be >0. This is used to weight the particles."
        pg.weight .*= val / pg.charge
    else
        setfield!(pg, s, val)
    end
end

function Base.show(io::IO, pg::ParticleGroup)
    print(io, "ParticleGroup with $(length(pg)) particles with total charge $(pg.charge) C")
end

function Base.show(io::IO, ::MIME"text/plain", pg::ParticleGroup)
    memloc = string(objectid(pg), base=16)
    print(io, "<ParticleGroup with $(length(pg)) particles at 0x$memloc>")
end

function units(pg::ParticleGroup, key::String)
    # to be implemented
end

# Statistical methods
function min(pg::ParticleGroup, key::String)
    return minimum(getproperty(pg, Symbol(key)))
end

function max(pg::ParticleGroup, key::String)
    return maximum(getproperty(pg, Symbol(key)))
end

function ptp(pg::ParticleGroup, key::String)
    return maximum(getproperty(pg, Symbol(key))) - minimum(getproperty(pg, Symbol(key)))
end

function avg(pg::ParticleGroup, key::String)
    dat = getproperty(pg, Symbol(key))
    if length(dat) == 1
        return dat[1]
    end
    return mean(dat, weights(pg.weight))
end

function delta(pg::ParticleGroup, key::String)
    return getproperty(pg, Symbol(key)) .- avg(pg, key)
end

function std(pg::ParticleGroup, key::String)
    dat = getproperty(pg, Symbol(key))
    if length(dat) == 1
        return 0
    end
    avg_dat = avg(pg, key)
    return sqrt(mean((dat .- avg_dat).^2, weights(pg.weight)))
end

function cov(pg::ParticleGroup, keys::AbstractString...)
    dats = hcat([getproperty(pg, Symbol(key)) for key in keys]...)
    return StatsBase.cov(dats, weights(pg.weight), 1)
end

# TODO: multidimensional histograms
# function histogramdd(pg::ParticleGroup, keys::String...; bins=10, range=nothing)

#=
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
=#

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

function drift_to_z!(pg::ParticleGroup, z::Real)
    dt = (z .- pg.z) ./ (pg.beta_z .* c_light)
    drift!(pg, dt)
    # Fix z to be exactly this value
    pg.z .= z
end

function drift_to_z!(pg::ParticleGroup)
    z = avg(pg, "z")
    drift_to_z!(pg, z)
end

function drift_to_t!(pg::ParticleGroup, t::Real)
    dt = t .- pg.t
    drift!(pg, dt)
    # Fix t to be exactly this value
    pg.t .= t
end

function drift_to_t!(pg::ParticleGroup)
    t = avg(pg, "t")
    drift_to_t!(pg, t)
end

"""
    Base.getindex(pg::ParticleGroup, x)

Returns a property or statistical quantity that can be computed:

- `P['x']` returns the x array
- `P['sigmx_x']` returns the std(x) scalar
- `P['norm_emit_x']` returns the norm_emit_x scalar

Parts can also be given. Example: `P[1:10]` returns a new ParticleGroup with the first 10 elements.
"""
function Base.getindex(pg::ParticleGroup, x)
    if !isa(x, String)
        return particle_parts(pg, x)
    end

    # Special case for z/c
    if x == "z/c"
        return pg.z ./ c_light
    end

    if startswith(x, "cov_")
        subkeys = split(x[5:end], "_")
        @assert length(subkeys) == 2 "Need 2 properties in covariance request: $x"
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
    else
        return getproperty(pg, Symbol(x))
    end
end

function Base.:+(pg1::ParticleGroup, pg2::ParticleGroup)
    return join_particle_groups(pg1, pg2)
end

function Base.:(==)(pg1::ParticleGroup, pg2::ParticleGroup)
    for key in ["x", "px", "y", "py", "z", "pz", "t", "status", "weight", "id"]
        if !all(isapprox.(getproperty(pg1, Symbol(key)), getproperty(pg2, Symbol(key))))
            return false
        end
    end
    return true
end

# Helper functions
function split_particles(pg::ParticleGroup; n_chunks=100, key="z")
    # Sorting
    zlist = getproperty(pg, Symbol(key))
    iz = sortperm(zlist)

    # Split particles into chunks
    plist = ParticleGroup[]
    for chunk in Iterators.partition(iz, ceil(Int, length(iz)/n_chunks))
        # Create new particle group with subset of particles
        new_pg = ParticleGroup{eltype(pg.x)}(
            pg.x[chunk], pg.px[chunk], pg.y[chunk], pg.py[chunk],
            pg.z[chunk], pg.pz[chunk], pg.t[chunk], pg.status[chunk],
            pg.weight[chunk], pg.species, pg.id[chunk]
        )
        push!(plist, new_pg)
    end

    return plist
end

function particle_parts(pg::ParticleGroup, x)
    return ParticleGroup{eltype(pg.x)}(
        pg.x[x], pg.px[x], pg.y[x], pg.py[x],
        pg.z[x], pg.pz[x], pg.t[x], pg.status[x],
        pg.weight[x], pg.species, pg.id[x]
    )
end

function join_particle_groups(pgs::ParticleGroup...)
    species = [pg.species for pg in pgs]
    species0 = species[1]
    @assert all(spe == species0 for spe in species) "species must be the same to join"

    # Concatenate arrays
    x = vcat([pg.x for pg in pgs]...)
    px = vcat([pg.px for pg in pgs]...)
    y = vcat([pg.y for pg in pgs]...)
    py = vcat([pg.py for pg in pgs]...)
    z = vcat([pg.z for pg in pgs]...)
    pz = vcat([pg.pz for pg in pgs]...)
    t = vcat([pg.t for pg in pgs]...)
    status = vcat([pg.status for pg in pgs]...)
    weight = vcat([pg.weight for pg in pgs]...)
    id = vcat([pg.id for pg in pgs]...)

    return ParticleGroup{eltype(x)}(
        x, px, y, py, z, pz, t, status, weight, species0, id
    )
end