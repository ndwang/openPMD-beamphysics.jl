"""
    read_simion_ION_file(filename)

Read a SIMION *.ion file directly into a dictionary compatible dictionary.

*.ion files are used by SIMION to define input particles for tracking.

The keys of the dictionary are:

- TOB: "time of birth" [microsec]
- MASS: ion mass [atomic mass units]
- CHARGE: ion charge [e], so an electron has CHARGE = -1
- X: x coordinate [mm]  NOTE: the main axis of cylindrical symmetry in SIMION is x, vs. z for particle group
- Y: y coordinate [mm]
- Z: z coordinate [mm]
- AZ: azimuthal angle of ions momentum vector [deg]  NOTE: SIMION's angles are not standard sphereical coordinatre angles
- EL: elevation angle of ions momentum vector [deg]  NOTE: SIMION's angles are not standard sphereical coordinatre angles
- KE: kinetic energy [eV]
- CWF: charge weighting factor (for non-uniform weighting in Space Charge calculations)
- COLOR: int [0, 1, 2...] - defines color of ion track in SIMION plotting
"""
function read_simion_ION_file(filename::String)
    simion_vars = [
        "TOB",
        "MASS",
        "CHARGE",
        "X",
        "Y",
        "Z",
        "AZ",
        "EL",
        "KE",
        "CWF",
        "COLOR"
    ]

    data = readdlm(filename, ',', skipstart=1)
    return Dict(p => data[:, findfirst(isequal(p), simion_vars)] for p in simion_vars)
end

"""
    write_simion_ION_file(particle_group, outfile; verbose=false)

Write a ParticleGroup to a SIMION *.ion file.

Parameters
----------
particle_group : ParticleGroup
    The particle group to write
outfile : String
    The output file name
verbose : Bool, optional
    If true, print status messages
"""
function write_simion_ION_file(particle_group, outfile::String; verbose::Bool=false)
    n = particle_group.n_particle
    
    # Convert to SIMION coordinates
    # SIMION uses x as the main axis, while particle_group uses z
    x = particle_group.z .* 1000  # Convert m to mm
    y = particle_group.y .* 1000  # Convert m to mm
    z = particle_group.x .* 1000  # Convert m to mm
    
    # Calculate angles and energy
    p = sqrt.(particle_group.px.^2 .+ particle_group.py.^2 .+ particle_group.pz.^2)
    KE = (sqrt.(1 .+ (p ./ particle_group.mass).^2) .- 1) .* particle_group.mass  # Kinetic energy in eV
    
    # Calculate angles in degrees
    AZ = atand.(particle_group.py, particle_group.px)  # Azimuthal angle
    EL = atand.(sqrt.(particle_group.px.^2 .+ particle_group.py.^2), particle_group.pz)  # Elevation angle
    
    # Create data matrix
    data = hcat(
        zeros(n),  # TOB
        fill(particle_group.mass, n),  # MASS
        fill(particle_group.species_charge, n),  # CHARGE
        x,  # X
        y,  # Y
        z,  # Z
        AZ,  # AZ
        EL,  # EL
        KE,  # KE
        abs.(particle_group.weight),  # CWF
        zeros(Int, n)  # COLOR
    )
    
    # Write header
    header = "TOB,MASS,CHARGE,X,Y,Z,AZ,EL,KE,CWF,COLOR"
    
    verbose && println("Writing $n particles to $outfile")
    
    # Write file
    open(outfile, "w") do io
        println(io, header)
        writedlm(io, data, ',')
    end
end 