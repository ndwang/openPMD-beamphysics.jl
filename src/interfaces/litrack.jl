"""
    write_litrack(particle_group, outfile="litrack.zd", p0c=nothing, verbose=false)

LiTrack is a Matlab code for the longitudinal phase space only.

ASCII with columns:
    c*t (mm), relative energy deviation (%).

The head of the bunch has smaller ct.

bunch head is at smaller z. The relative energy deviation unit is really %, which means I time 100 from relative energy spread.

This routine makes ASCII particles, with column labels:
    'x', 'y', 'z', 'GBx', 'GBy', 'GBz', 't', 'q', 'nmacro'
in SI units.

For now, only electrons are supported.
"""
function write_litrack(particle_group, outfile::String="litrack.zd", p0c=nothing, verbose::Bool=false)
    P = particle_group  # convenience

    P.species == "electron" || error("Only electrons supported")  # TODO: add more species
    all(P.weight .> 0) || error("ParticleGroup.weight must be > 0")

    n = particle_group.n_particle
    z = unique(P.z)
    length(z) == 1 || error("All particles must be at the same z. Please call .drift_to_z()")

    if isnothing(p0c)
        p0c = P.mean_p
        verbose && println("Using mean_p as the reference momentum: $p0c eV/c")
    end

    ct = c_light .* P.t
    delta = P.p ./ p0c .- 1.0

    header = """% LiTrack particles
%
% Created using the openPMD-beamphysics Julia package
% https://github.com/ChristopherMayes/openPMD-beamphysics
%
% species: $(P.species)
% n_particle: $n
% total charge: $(P.charge) (C)
% reference momentum p0: $p0c (eV/c)
%
% Columns: ct, delta = p/p0 -1
% Units: mm, percent"""

    outdat = [ct .* 1000 delta .* 100]

    verbose && println("writing $n LiTrack particles to $outfile")

    # Write ASCII
    open(outfile, "w") do io
        println(io, header)
        writedlm(io, outdat, " ", fmt="%20.12e")
    end

    return outfile
end 