module BeamPhysicsBeamTrackingExt

using BeamPhysics
using BeamTracking
using AtomicAndPhysicalConstants: C_LIGHT, chargeof, massof

import Statistics: mean

"""
    BeamTracking.Bunch(pg::AbstractParticleGroup; p_over_q_ref=nothing, t_ref=nothing)

Convert a `ParticleGroup` (or `ParticleGroupView`) into a BeamTracking `Bunch`.

Phase space coordinates are mapped as follows:
- `v[:,1]` = `x` [m]
- `v[:,2]` = `px / p_ref` (normalized transverse momentum)
- `v[:,3]` = `y` [m]
- `v[:,4]` = `py / p_ref` (normalized transverse momentum)
- `v[:,5]` = `z` in Bmad coordinates [m]
- `v[:,6]` = `δ = (p - p_ref) / p_ref` (fractional momentum deviation)

where `p_ref` is the reference momentum derived from `p_over_q_ref`.

# Keyword Arguments
- `p_over_q_ref`: Reference momentum / charge / speed of light. Defaults to the
  mean total momentum of alive particles converted via `pc / (c * q)`.
- `t_ref`: Reference time [s]. Defaults to the mean time of alive particles.
"""
function BeamTracking.Bunch(pg::BeamPhysics.AbstractParticleGroup; p_over_q_ref=nothing, t_ref=nothing)
    N = length(pg)
    species = pg.species

    alive = pg.status .== 1

    # Total momentum per particle [eV/c]
    p_total = BeamPhysics.momentum(pg)

    # Reference momentum
    if isnothing(p_over_q_ref)
        p_ref = mean(p_total[alive])
        p_over_q_ref = p_ref / C_LIGHT / chargeof(species)
    else
        p_ref = abs(p_over_q_ref * C_LIGHT * chargeof(species))
    end

    # Reference time
    if isnothing(t_ref)
        t_ref = mean(pg.t[alive])
    end

    # Reference z (mean z of alive particles)
    z_ref = mean(pg.z[alive])

    # Reference beta
    mc2 = massof(species)
    beta_ref = p_ref / sqrt(p_ref^2 + mc2^2)

    # Build N×6 coordinate matrix
    v = similar(pg.x, N, 6)
    v[:, 1] .= pg.x
    v[:, 2] .= pg.px ./ p_ref
    v[:, 3] .= pg.y
    v[:, 4] .= pg.py ./ p_ref
    v[:, 5] .= (pg.z .- z_ref) .- beta_ref .* C_LIGHT .* (pg.t .- t_ref)
    v[:, 6] .= p_total ./ p_ref .- 1

    # Map particle status (reinterpret to UInt8, keeping values as-is)
    state = similar(pg.x, UInt8, N)
    state .= pg.status

    weight = similar(pg.x, N)
    weight .= pg.weight
    coords = BeamTracking.Coords(state, v, nothing, weight)
    return BeamTracking.Bunch(species, p_over_q_ref, t_ref, coords)
end

"""
    BeamTracking.Bunch(file::AbstractString; p_over_q_ref=nothing, t_ref=nothing)

Load particles from an openPMD HDF5 file directly into a BeamTracking `Bunch`.

Equivalent to `Bunch(ParticleGroup(file); p_over_q_ref, t_ref)`.
"""
function BeamTracking.Bunch(file::AbstractString; p_over_q_ref=nothing, t_ref=nothing)
    pg = BeamPhysics.ParticleGroup(file)
    return BeamTracking.Bunch(pg; p_over_q_ref=p_over_q_ref, t_ref=t_ref)
end

end # module
