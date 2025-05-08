"""
    write_opal(particle_group, outfile, dist_type="emitted", verbose=false)

OPAL's ASCII format is described in:
    https://gitlab.psi.ch/OPAL/Manual-2.2/wikis/distribution

outfile is the name out the ASCII file to be written.

dist_type is one of:
    'emitted' : The longitudinal dimension is
                described by time, and the particles are emitted
                from the cathode plane over a set number of 'emission steps'.

    'injected': the longiudinal dimension is
                described by z, and the particles are born instaneously
                in the simulation at time step 0.
"""
function write_opal(particle_group, outfile::String, dist_type::String="emitted", verbose::Bool=false)
    n = particle_group.n_particle
    x = particle_group.x
    y = particle_group.y

    # Get longitudinal coordinate
    if dist_type == "emitted"
        # Check that z are all the same
        unique_z = unique(particle_group.z)
        length(unique_z) == 1 || error("All particles must be at the same z position")

        z = particle_group.t

    elseif dist_type == "injected"
        # Check that t are all the same
        unique_t = unique(particle_group.t)
        length(unique_t) == 1 || error("All particles must be at the same time")

        z = particle_group.z

    else
        error("unknown dist_type: $dist_type")
    end

    gamma = particle_group.gamma
    GBx = gamma .* particle_group.beta_x
    GBy = gamma .* particle_group.beta_y
    GBz = gamma .* particle_group.beta_z

    header = string(n)
    dat = [x GBx y GBy z GBz]

    verbose && println("writing $dist_type $n particles to $outfile")
    
    open(outfile, "w") do io
        println(io, header)
        writedlm(io, dat, " ", fmt="%20.12e")
    end
end

"""
    opal_to_data(h5)

Converts an OPAL step to the standard data format for openPMD-beamphysics.

In OPAL, the momenta px, py, pz are gamma*beta_x, gamma*beta_y, gamma*beta_z.

TODO: More species.
"""
function opal_to_data(h5)
    D = Dict(h5.attrs)
    mc2 = D["MASS"][1] * 1e9  # GeV -> eV
    charge = D["CHARGE"][1]  # total charge in C
    t = D["TIME"][1]  # s
    ptypes = read(h5["ptype"])  # 0 = electron?

    # Not used: pref = D['RefPartP']*mc2  #
    rref = D["RefPartR"]

    n = length(ptypes)
    all(ptypes .== 0) || error("Only electrons supported")
    species = "electron"
    status = 1
    data = Dict{String,Any}(
        "x" => read(h5["x"]) .+ rref[1],
        "y" => read(h5["y"]) .+ rref[2],
        "z" => read(h5["z"]) .+ rref[3],
        "px" => read(h5["px"]) .* mc2,
        "py" => read(h5["py"]) .* mc2,
        "pz" => read(h5["pz"]) .* mc2,
        "t" => fill(t, n),
        "status" => fill(status, n),
        "species" => species,
        "weight" => fill(abs(charge) / n, n)
    )
    return data
end 