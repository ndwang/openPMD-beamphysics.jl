"""
    write_elegant(particle_group, outfile; verbose=false)

Elegant uses SDDS files.

Because elegant is an s-based code, particles are drifted to the center.

This routine writes an SDDS1 ASCII file, with a parameter
    Charge
and columns
    't', 'x', 'xp', 'y', 'yp', 'p'
where 'p' is gamma*beta, in units:
    elegant units are:
    s, m, 1, m, 1, 'mass*c'

All weights must be the same.
"""
function write_elegant(particle_group, outfile::String; verbose::Bool=false)
    # Work on a copy, because we will drift
    P = copy(particle_group)

    # Drift to z
    drift_to_z!(P)

    # Form data
    keys = ["t", "x", "xp", "y", "yp", "p"]
    dat = Dict(k => P[k] for k in keys)
    # Correct p, this is really gamma*beta
    dat["p"] ./= P.mass

    verbose && println("writing $(length(P)) particles to $outfile")

    # Correct units for p, depending on the species
    species = particle_group.species
    p_units = if species in ["electron", "positron"]
        "units=\"m\$be\$nc\", "
    elseif species in ["proton"]
        "units=\"m\$bp\$nc\", "
    else
        ""
    end

    # Note that the order of the columns matters below
    header = """SDDS1
!
! Created using the openPMD-beamphysics Julia package
! https://github.com/ChristopherMayes/openPMD-beamphysics
! species: $(P["species"])
!
&parameter name=Charge, type=double, units=C, description="total charge in Coulombs" &end
&column name=t,  type=double, units=s, description="time in seconds" &end
&column name=x,  type=double, units=m, description="x in meters" &end
&column name=xp, type=double, description="px/pz" &end
&column name=y,  type=double, units=m, description="y in meters" &end
&column name=yp, type=double, description="py/pz" &end
&column name=p,  type=double, $(p_units)description="relativistic gamma*beta" &end
&data mode=ascii &end
$(P["charge"])
$(length(P))"""

    # Write ASCII
    outdat = hcat([dat[k] for k in keys]...)
    open(outfile, "w") do io
        println(io, header)
        writedlm(io, outdat, " ")
    end

    return outfile
end

"""
    elegant_h5_to_data(h5, group="page1", species="electron")

Converts elegant data from an h5 handle or file to data for openPMD-beamphysics

Elegant h5 has datasets:

x, xp, y, xp, p=beta*gamma, t
m,  1,   m, 1, 1, s

Momentum are reconstructed as:

pz = p / sqrt(1+xp^2 + yp^2)
px = xp * pz
py = yp * pz

In s-based codes, z=0 by definition.

For now, only electrons are allowed.

Units are checked.

All particles are assumed to be live (status = 1)

TODO: More species.
"""
function elegant_h5_to_data(h5, group::String="page1", species::String="electron")
    # Allow for opening a file
    if isa(h5, String)
        isfile(h5) || error("File does not exist: $h5")
        h5 = h5open(h5, "r")
    end

    g = group == "" ? h5 : h5[group]

    species == "electron" || error("$species not allowed yet. Only electron is implemented.")
    mc2 = mec2

    # These should exist
    col = g["columns"]
    par = g["parameters"]

    p = read(col["p"]) .* mc2
    xp = read(col["xp"])
    yp = read(col["yp"])
    pz = p ./ sqrt.(1 .+ xp.^2 .+ yp.^2)
    px = xp .* pz
    py = yp .* pz

    # Check charge unit
    charge = read(par["Charge"])[1]
    charge_unit = read(attributes(par["Charge"])["units"])
    charge_unit == b"C" || error("Expected C as unit for Charge")

    # Check dataset units
    expected_units = Dict(
        "p" => b"m\$be\$nc",
        "xp" => b"",
        "yp" => b"",
        "x" => b"m",
        "y" => b"m",
        "t" => b"s"
    )
    for (c, v) in expected_units
        u = read(attributes(col[c])["units"])
        u == v || error("Dataset $c units expected to have $v, but have $u")
    end

    # number of particles
    n = length(p)

    status = 1
    data = Dict{String,Any}(
        "x" => read(col["x"]),
        "y" => read(col["y"]),
        "z" => fill(0, n),
        "px" => px,
        "py" => py,
        "pz" => pz,
        "t" => read(col["t"]),
        "status" => fill(status, n),
        "species" => species,
        "weight" => fill(abs(charge) / n, n),
        "id" => read(col["particleID"])
    )
    return data
end

"""
    load_sdds(sddsfile, columns, sdds2plaindata_bin="sdds2plaindata"; verbose=false)

Get tabular data from SDDS file
"""
function load_sdds(sddsfile::String, columns::Vector{String}, sdds2plaindata_bin::String="sdds2plaindata"; verbose::Bool=false)
    outfile = sddsfile * "_table"
    cmd0 = [sdds2plaindata_bin, sddsfile, outfile, "-noRowCount", "-outputMode=ascii"]
    cmd = vcat(cmd0, ["-col=$c" for c in columns], ["-separator= "])

    verbose && println("load_sdds command: ", join(cmd, " "))

    output = read(`$cmd`, String)
    if !isempty(output)
        error("load_sdds error: $output")
    end

    # Read table
    rdat = readdlm(outfile, Float64)

    dat = Dict{String,Vector{Float64}}()
    for (i, key) in enumerate(columns)
        dat[key] = rdat[:, i]
    end

    # Cleanup
    rm(outfile)

    return dat
end

"""
    elegant_to_data(sddsfile, charge=1.0, sdds2plaindata_bin="sdds2plaindata", species="electron"; verbose=false)

Converts elegant SDDS data to data for openPMD-beamphysics.

Similar to elegant_h5_to_data, which is better because it checks units and does not need the conversion to ASCII

Requires the SDDS utility sdds2plaindata

Columns in the file should be:

x, xp, y, xp, p=beta*gamma, t
m,  1,   m, 1, 1, s

Momentum are reconstructed as:

pz = p / sqrt(1+xp^2 + yp^2)
px = xp * pz
py = yp * pz

In s-based codes, z=0 by definition.

For now, only electrons are allowed.

All particles are assumed to be live (status = 1)

TODO: More species.

Also see: elegant_h5_to_data
"""
function elegant_to_data(sddsfile::String, charge::Float64=1.0, sdds2plaindata_bin::String="sdds2plaindata", species::String="electron"; verbose::Bool=false)
    col = load_sdds(
        sddsfile,
        ["x", "xp", "y", "yp", "t", "p"],
        sdds2plaindata_bin=sdds2plaindata_bin,
        verbose=verbose
    )

    species == "electron" || error("$species not allowed yet. Only electron is implemented.")
    mc2 = mec2

    p = col["p"] .* mc2
    xp = col["xp"]
    yp = col["yp"]
    pz = p ./ sqrt.(1 .+ xp.^2 .+ yp.^2)
    px = xp .* pz
    py = yp .* pz

    n = length(p)
    data = Dict{String,Any}(
        "x" => col["x"],
        "y" => col["y"],
        "z" => fill(0, n),
        "px" => px,
        "py" => py,
        "pz" => pz,
        "t" => col["t"],
        "status" => fill(1, n),
        "species" => species,
        "weight" => fill(abs(charge) / n, n)
    )
    return data
end 