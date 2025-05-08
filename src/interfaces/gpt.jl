"""
    write_gpt(particle_group, outfile; asci2gdf_bin=nothing, verbose=false)

GPT uses a custom binary format GDF for particle data. This can be created with the
asci2gdf utility as:
asci2gdf -o particles.gdf particles.txt

This routine makes ASCII particles, with column labels:
    'x', 'y', 'z', 'GBx', 'GBy', 'GBz', 't', 'q', 'm', 'nmacro'
in SI units.
"""
function write_gpt(particle_group, outfile::String; asci2gdf_bin=nothing, verbose::Bool=false)
    all(particle_group.weight .>= 0) || error("ParticleGroup.weight must be >= 0")

    q = particle_group.species_charge
    mc2 = particle_group.mass  # returns mass_of(particle_group.species) [eV]
    m = mc2 * (e_charge / c_light^2)
    n = particle_group.n_particle
    gamma = particle_group.gamma

    dat = Dict{String,Any}(
        "x" => particle_group.x,
        "y" => particle_group.y,
        "z" => particle_group.z,
        "GBx" => gamma .* particle_group.beta_x,
        "GBy" => gamma .* particle_group.beta_y,
        "GBz" => gamma .* particle_group.beta_z,
        "t" => particle_group.t,
        "q" => fill(q, n),
        "m" => fill(m, n),
        "nmacro" => abs.(particle_group.weight ./ q)
    )

    if hasproperty(particle_group, :id)
        dat["ID"] = particle_group.id
    else
        dat["ID"] = 1:particle_group.n_particle
    end

    header = join(keys(dat), " ")
    outdat = hcat([dat[k] for k in keys(dat)]...)

    verbose && println("writing $n particles to $outfile")

    # Write ASCII
    open(outfile, "w") do io
        println(io, header)
        writedlm(io, outdat, " ")
    end

    if !isnothing(asci2gdf_bin)
        run_asci2gdf(outfile, asci2gdf_bin, verbose=verbose)
    else
        println("ASCII particles written. Convert to GDF using: asci2df -o particles.gdf $outfile")
    end
end

"""
    run_asci2gdf(outfile, asci2gdf_bin; verbose=false)

Helper function to convert ASCII to GDF using asci2gdf
"""
function run_asci2gdf(outfile::String, asci2gdf_bin::String; verbose::Bool=false)
    tempfile = outfile * ".txt"
    mv(outfile, tempfile)

    asci2gdf_bin = expanduser(asci2gdf_bin)
    isfile(asci2gdf_bin) || error("$asci2gdf_bin does not exist")
    cmd = `$asci2gdf_bin -o $outfile $tempfile`
    verbose && println(cmd)
    run(cmd)

    # Cleanup
    rm(tempfile)

    verbose && println("Written GDF file: $outfile")

    return outfile
end

"""
    write_gpt_fieldmap(fm, outfile; asci2gdf_bin=nothing, verbose=false)

Writes a GPT fieldmap file from a FieldMesh object.

Requires cylindrical geometry for now.
"""
function write_gpt_fieldmap(fm, outfile::String; asci2gdf_bin=nothing, verbose::Bool=false)
    if fm.geometry == "cylindrical" && fm.coord_vec("r")[1] == fm.coord_vec("r")[end]
        return write_gpt_1d_fieldmap(fm, outfile, asci2gdf_bin=asci2gdf_bin, verbose=verbose)
    elseif fm.geometry == "cylindrical"
        return write_gpt_2d_fieldmap(fm, outfile, asci2gdf_bin=asci2gdf_bin, verbose=verbose)
    elseif fm.geometry == "rectangular"
        return write_gpt_3d_fieldmap(fm, outfile, asci2gdf_bin=asci2gdf_bin, verbose=verbose)
    else
        error("Unknown geometry $(fm.geometry)")
    end
end

"""
    write_gpt_1d_fieldmap(fm, outfile; asci2gdf_bin=nothing, verbose=false)

Writes a GPT fieldmap file from a FieldMesh object.

Requires cylindrical geometry for now.
"""
function write_gpt_1d_fieldmap(fm, outfile::String; asci2gdf_bin=nothing, verbose::Bool=false)
    fm.geometry == "cylindrical" || error("Geometry: $(fm.geometry) not implemented")
    fm.shape[2] == 1 || error("Cylindrical symmetry required")
    fm.coord_vec("r")[1] == 0 || error("r[1] must equal 0")

    dat = Dict{String,Any}()
    dat["Z"] = fm.coord_vec("z")

    keys = ["Z"]
    if fm.is_static
        if fm.is_pure_magnetic
            keys = ["Z", "Bz"]
            dat["Bz"] = real(fm["Bz"][1, 1, :])
        elseif fm.is_pure_electric
            keys = ["Z", "Ez"]
            dat["Er"] = real(fm["Er"][1, 1, :])
            dat["Ez"] = real(fm["Ez"][1, 1, :])
        else
            error("Mixed static field TODO")
        end
    else
        # Use internal Superfish routine
        keys = ["Z", "Ez"]
        _, dat["Ez"], _, _ = fish_complex_to_real_fields(fm, verbose=verbose)
        dat["Ez"] = dat["Ez"][1, 1, :]
    end

    # Flatten dat
    gptdata = hcat([vec(dat[k]) for k in keys]...)

    # Write file
    open(outfile, "w") do io
        println(io, join(keys, " "))
        writedlm(io, gptdata, " ")
    end

    if !isnothing(asci2gdf_bin)
        run_asci2gdf(outfile, asci2gdf_bin, verbose=verbose)
    elseif verbose
        println("ASCII field data written. Convert to GDF using: asci2df -o field.gdf $outfile")
    end

    return outfile
end

"""
    write_gpt_2d_fieldmap(fm, outfile; asci2gdf_bin=nothing, verbose=false)

Writes a GPT fieldmap file from a FieldMesh object.

Requires cylindrical geometry for now.
"""
function write_gpt_2d_fieldmap(fm, outfile::String; asci2gdf_bin=nothing, verbose::Bool=false)
    fm.geometry == "cylindrical" || error("Geometry: $(fm.geometry) not implemented")
    fm.shape[2] == 1 || error("Cylindrical symmetry required")

    dat = Dict{String,Any}()
    R, Z = meshgrid(fm.coord_vec("r"), fm.coord_vec("z"))
    dat["R"] = R
    dat["Z"] = Z

    keys = ["R", "Z"]
    if fm.is_static
        if fm.is_pure_magnetic
            keys = ["R", "Z", "Br", "Bz"]
            dat["Br"] = real(fm["Br"][:, 1, :])
            dat["Bz"] = real(fm["Bz"][:, 1, :])
        elseif fm.is_pure_electric
            keys = ["R", "Z", "Er", "Ez"]
            dat["Er"] = real(fm["Er"][:, 1, :])
            dat["Ez"] = real(fm["Ez"][:, 1, :])
        else
            error("Mixed static field TODO")
        end
    else
        # Use internal Superfish routine
        keys = ["R", "Z", "Er", "Ez", "Bphi"]
        dat["Er"], dat["Ez"], dat["Bphi"], _ = fish_complex_to_real_fields(fm, verbose=verbose)
    end

    # Flatten dat
    gptdata = hcat([vec(dat[k]) for k in keys]...)

    # Write file
    open(outfile, "w") do io
        println(io, join(keys, " "))
        writedlm(io, gptdata, " ")
    end

    if !isnothing(asci2gdf_bin)
        run_asci2gdf(outfile, asci2gdf_bin, verbose=verbose)
    elseif verbose
        println("ASCII field data written. Convert to GDF using: asci2df -o field.gdf $outfile")
    end

    return outfile
end

"""
    write_gpt_3d_fieldmap(fm, outfile; asci2gdf_bin=nothing, verbose=false)

Writes a 3D GPT fieldmap file from a FieldMesh object.
"""
function write_gpt_3d_fieldmap(fm, outfile::String; asci2gdf_bin=nothing, verbose::Bool=false)
    fm.geometry == "rectangular" || error("Geometry: $(fm.geometry) not implemented")

    dat = Dict{String,Any}()
    X, Y, Z = meshgrid(fm.coord_vec("x"), fm.coord_vec("y"), fm.coord_vec("z"))
    dat["X"] = X
    dat["Y"] = Y
    dat["Z"] = Z

    keys = ["X", "Y", "Z"]
    if fm.is_static
        if fm.is_pure_magnetic
            keys = vcat(keys, ["Bx", "By", "Bz"])
            dat["Bx"] = real(fm["Bx"])
            dat["By"] = real(fm["By"])
            dat["Bz"] = real(fm["Bz"])
        elseif fm.is_pure_electric
            keys = vcat(keys, ["Ex", "Ey", "Ez"])
            dat["Ex"] = real(fm["Ex"])
            dat["Ey"] = real(fm["Ey"])
            dat["Ez"] = real(fm["Ez"])
        else
            error("Mixed static field TODO")
        end
    else
        error("Complex 3D Fields not implement yet!")
        # keys = vcat(['X', 'Y', 'Z'], ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'])
        # dat['Er'], dat['Ez'], dat['Bphi'], _ = fish_complex_to_real_fields(fm, verbose=verbose)
    end

    # Flatten dat
    gptdata = hcat([vec(dat[k]) for k in keys]...)

    # Write file
    open(outfile, "w") do io
        println(io, join(keys, " "))
        writedlm(io, gptdata, " ")
    end

    if !isnothing(asci2gdf_bin)
        run_asci2gdf(outfile, asci2gdf_bin, verbose=verbose)
    elseif verbose
        println("ASCII field data written. Convert to GDF using: asci2df -o field.gdf $outfile")
    end

    return outfile
end

# Helper function for meshgrid
function meshgrid(x, y, z=nothing)
    if isnothing(z)
        X = repeat(x, 1, length(y))
        Y = repeat(y', length(x), 1)
        return X, Y
    else
        X = repeat(repeat(x, 1, length(y)), 1, 1, length(z))
        Y = repeat(repeat(y', length(x), 1), 1, 1, length(z))
        Z = repeat(repeat(z', 1, 1, length(x)), length(y), 1, 1)
        Z = permutedims(Z, (2, 1, 3))
        return X, Y, Z
    end
end 