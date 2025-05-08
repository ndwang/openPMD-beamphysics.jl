# Constants
const ASTRA_SPECIES_NAME = Dict(1 => "electron", 2 => "positron", 3 => "proton", 4 => "hydrogen")
const ASTRA_SPECIES_INDEX = Dict(v => k for (k, v) in ASTRA_SPECIES_NAME)

const ASTRA_PARTICLE_STATUS_NAMES = Dict(
    -1 => "standard particle, at the cathode",
    3 => "trajectory probe particle",
    5 => "standard particle"
)

"""
    parse_astra_phase_file(filePath)

Parses ASTRA particle dumps to data dict, that corresponds to the
openPMD-beamphysics ParticleGroup data= input.

Units are in m, s, eV/c

Live particles (status==5) are relabeled as status = 1.
Original status == 2 are relabeled to status = 2 (previously unused by Astra)
"""
function parse_astra_phase_file(filePath::String)
    # Check if file exists
    isfile(filePath) || error("particle file does not exist: $filePath")

    # Read data
    data = readdlm(filePath)
    ref = data[1, :]  # Reference particle

    # Extract data
    x = data[2:end, 1]  # position in m
    y = data[2:end, 2]
    z_rel = data[2:end, 3]
    z_ref = ref[3]
    
    # momenta in eV/c
    px = data[2:end, 4]
    py = data[2:end, 5]
    pz_rel = data[2:end, 6]
    pz_ref = ref[6]

    # Time in seconds
    t_ref = ref[7] * 1e-9
    t_rel = data[2:end, 7] * 1e-9

    # macro charge in Coulomb. The sign doesn't matter, so make positive
    qmacro = abs.(data[2:end, 8] * 1e-9)

    species_index = Int.(data[2:end, 9])
    status = Int.(data[2:end, 10])

    n_particle = length(x)

    # Create data dictionary
    data_dict = Dict{String,Any}(
        "x" => x,
        "y" => y,
        "z" => z_rel .+ z_ref,
        "px" => px,
        "py" => py,
        "pz" => pz_rel .+ pz_ref,
        "t_clock" => t_rel .+ t_ref,
        "t" => t_ref,
        "weight" => qmacro,
        "n_particle" => n_particle
    )

    # Status handling
    # The standard defines 1 as a live particle, but astra uses 1 as a 'passive' particle
    # and 5 as a 'standard' particle. 2 is not used.
    # To preserve this information, make 1->2 and then 5->1
    status[status .== 1] .= 2
    status[status .== 5] .= 1
    data_dict["status"] = status

    # Check species
    unique_species = unique(species_index)
    length(unique_species) == 1 || error("All species must be the same")
    data_dict["species"] = ASTRA_SPECIES_NAME[unique_species[1]]

    return data_dict
end

"""
    write_astra(particle_group, outfile; verbose=false, probe=false)

Writes ASTRA style particles from particle_group type data.

For now, the species must be electrons.

If probe, the six standard probe particles will be written.
"""
function write_astra(particle_group, outfile::String; verbose::Bool=false, probe::Bool=false)
    verbose && println("writing $(particle_group.n_particle) particles to $outfile")

    # number of lines in file
    size = particle_group.n_particle + 1  # Allow one for reference particle
    i_start = 1  # Start for data particles
    if probe
        # Add six probe particles, according to the manual
        size += 6
        i_start += 6
    end

    # Reference particle
    ref_particle = Dict{String,Any}("q" => 0)
    sigma = Dict{String,Float64}()
    for k in ["x", "y", "z", "px", "py", "pz", "t"]
        ref_particle[k] = mean(getproperty(particle_group, Symbol(k)))
        std_val = std(getproperty(particle_group, Symbol(k)))
        sigma[k] = std_val == 0 ? 1e-12 : std_val  # Give some size if zero
    end
    ref_particle["t"] *= 1e9  # s -> nS

    # Create data array
    data = zeros(size, 10)  # 10 columns: x,y,z,px,py,pz,t,q,index,status
    
    # Fill data
    for (i, k) in enumerate(["x", "y", "z", "px", "py", "pz", "t"])
        data[i_start:end, i] = getproperty(particle_group, Symbol(k))
    end
    data[:, 7] .*= 1e9  # s -> nS
    data[i_start:end, 8] = particle_group.weight .* 1e9  # C -> nC

    # Set species index
    data[:, 9] .= ASTRA_SPECIES_INDEX[particle_group.species]

    # Status handling
    status = particle_group.status
    astra_status = copy(status)
    astra_status[status .== 1] .= 5  # Astra normal (alive)
    astra_status[status .== 2] .= 1  # Astra passive
    astra_status[status .== ParticleStatus.CATHODE] .= -1  # Astra cathode
    data[i_start:end, 10] = astra_status

    # Handle reference particle
    ref_particle["status"] = -1 in astra_status ? -1 : 5

    # Subtract off reference z, pz, t
    for (i, k) in enumerate(["z", "pz", "t"])
        data[:, i] .-= ref_particle[k]
    end

    # Put ref particle in first position
    for (i, k) in enumerate(["x", "y", "z", "px", "py", "pz", "t", "q", "status"])
        data[1, i] = ref_particle[k]
    end

    # Optional: probes
    if probe
        data[2, 1] = 0.5 * sigma["x"]
        data[2, 7] = 0.5 * sigma["t"]
        data[3, 2] = 0.5 * sigma["y"]
        data[3, 7] = -0.5 * sigma["t"]
        data[4, 1] = 1.0 * sigma["x"]
        data[4, 7] = sigma["t"]
        data[5, 2] = 1.0 * sigma["y"]
        data[5, 7] = -sigma["t"]
        data[6, 1] = 1.5 * sigma["x"]
        data[6, 7] = 1.5 * sigma["t"]
        data[7, 2] = 1.5 * sigma["y"]
        data[7, 7] = -1.5 * sigma["t"]
        data[2:7, 10] .= -3
        data[2:7, 8] .= 0.5e-5  # ? Seems to be required
        data[2:7, 6] .= 0  # ? This is what the Astra Generator does
    end

    # Save in the 'high_res = T' format
    open(outfile, "w") do io
        for row in eachrow(data)
            println(io, join(@sprintf("%20.12e", x) for x in row[1:8]), " ", join(@sprintf("%4i", x) for x in row[9:10]))
        end
    end
end

"""
    write_astra_1d_fieldmap(fm, filePath)

Writes an ASTRA fieldmap file from a FieldMesh object.

Requires cylindrical geometry for now.

Parameters
----------
filePath: String
    Filename to write to
"""
function write_astra_1d_fieldmap(fm, filePath::String)
    z, fz = astra_1d_fieldmap_data(fm)
    dat = hcat(z, fz)
    writedlm(filePath, dat, "")
end

"""
    astra_1d_fieldmap_data(fm)

ASTRA fieldmap data.

Requires cylindrical geometry for now.

Returns
-------
z: Vector
    z coordinate in meters
fz: Vector
    field amplitude corresponding to z
"""
function astra_1d_fieldmap_data(fm)
    # Implementation depends on FieldMesh structure
    # This is a placeholder that needs to be implemented based on the FieldMesh type
    error("Not implemented yet")
end

"""
    vec_spacing(vec)

Returns the spacing of a vector.
"""
function vec_spacing(vec::Vector)
    return vec[2] - vec[1]
end

"""
    parse_astra_fieldmap_3d(filePath, frequency=0)

Parse ASTRA 3D fieldmap file.
"""
function parse_astra_fieldmap_3d(filePath::String, frequency::Float64=0.0)
    # Implementation depends on fieldmap format
    # This is a placeholder that needs to be implemented
    error("Not implemented yet")
end

"""
    read_astra_3d_fieldmaps(common_filePath, frequency=0)

Read ASTRA 3D fieldmaps.
"""
function read_astra_3d_fieldmaps(common_filePath::String, frequency::Float64=0.0)
    # Implementation depends on fieldmap format
    # This is a placeholder that needs to be implemented
    error("Not implemented yet")
end

"""
    write_astra_3d_fieldmaps(fieldmesh_object, common_filePath)

Write ASTRA 3D fieldmaps.
"""
function write_astra_3d_fieldmaps(fieldmesh_object, common_filePath::String)
    # Implementation depends on fieldmap format
    # This is a placeholder that needs to be implemented
    error("Not implemented yet")
end 