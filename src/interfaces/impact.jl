"""
    parse_impact_particles(filePath, names=("x", "GBx", "y", "GBy", "z", "GBz"), skiprows=0)

Parse Impact-T input and output particle data.
Typical filenames: 'partcl.data', 'fort.40', 'fort.50'.

Note that partcl.data has the number of particles in the first line, so skiprows=1 should be used.

Returns a dictionary with the particle data.

Impact-T input/output particles distribions are ASCII files with columns:
x (m)
GBy = gamma*beta_x (dimensionless)
y (m)
GBy = gamma*beta_y (dimensionless)
z (m)
GBz = gamma*beta_z (dimensionless)

Routine from lume-impact:
    https://github.com/ChristopherMayes/lume-impact
"""
function parse_impact_particles(filePath::String, names::Tuple=("x", "GBx", "y", "GBy", "z", "GBz"), skiprows::Int=0)
    data = readdlm(filePath, skipstart=skiprows)
    return Dict(names[i] => data[:,i] for i in 1:length(names))
end

"""
    impact_particles_to_particle_data(tout, mc2=0, species=nothing, time=0, macrocharge=0, cathode_kinetic_energy_ref=nothing, verbose=false)

Convert impact particles to data for ParticleGroup.

particle_charge is the charge in units of |e|

At the cathode, Impact-T translates z to t = z / (beta*c) for emission,
where (beta*c) is the velocity calculated from kinetic energy:
    header['Bkenergy'] in eV.
This is purely a conversion factor.

If cathode_kinetic_energy_ref is given, z will be parsed appropriately to t, and z will be set to 0.

Otherwise, particles will be set to the same time.
"""
function impact_particles_to_particle_data(tout::Dict, mc2::Float64=0.0, species=nothing, time::Float64=0.0, 
                                         macrocharge::Float64=0.0, cathode_kinetic_energy_ref=nothing, verbose::Bool=false)
    mc2 > 0 || error("mc2 must be specified")
    !isnothing(species) || error("species must be specified")

    data = Dict{String,Any}()

    n_particle = length(tout["x"])

    data["x"] = tout["x"]
    data["y"] = tout["y"]
    # data['z'] = tout['z'] will be handled below

    data["px"] = tout["GBx"] .* mc2
    data["py"] = tout["GBy"] .* mc2
    data["pz"] = tout["GBz"] .* mc2

    # Handle z
    if !isnothing(cathode_kinetic_energy_ref)
        # Cathode start
        z = fill(0.0, n_particle)

        # Note that this is purely a conversion factor.
        gamma = 1.0 + cathode_kinetic_energy_ref / mc2
        betac = sqrt(1 - 1 / gamma^2) * c_light

        t = tout["z"] ./ betac
        verbose && println("Converting z to t according to cathode_kinetic_energy_ref = $cathode_kinetic_energy_ref eV")

    else
        # Free space start
        z = tout["z"]
        t = fill(time, n_particle)
    end

    data["z"] = z
    data["t"] = t

    data["status"] = fill(1, n_particle)
    if macrocharge == 0
        weight = 1 / n_particle
    else
        weight = abs(macrocharge)
    end
    data["weight"] = fill(weight, n_particle)

    data["species"] = species
    data["n_particle"] = n_particle
    return data
end

"""
    write_impact(particle_group, outfile, cathode_kinetic_energy_ref=nothing, include_header=true, verbose=false)

Writes Impact-T style particles from particle_group type data.

outfile should ultimately be named 'partcl.data' for Impact-T

For now, the species must be electrons.

If cathode_kinetic_energy_ref is given, t will be used to compute z for cathode emission.

If include_header, the number of particles will be written as the first line. Default is true.

Otherwise, particles must have the same time, and should be started in free space.

A dict is returned with info about emission, for use in Impact-T
"""
function write_impact(particle_group, outfile::String, cathode_kinetic_energy_ref=nothing, 
                     include_header::Bool=true, verbose::Bool=false)
    function vprint(args...; kwargs...)
        verbose && println(args...; kwargs...)
    end

    n_particle = particle_group.n_particle

    vprint("writing $n_particle particles to $outfile")

    mc2 = particle_group.mass

    # Dict output for use in Impact-T
    output = Dict{String,Any}("input_particle_file" => outfile)
    output["Np"] = n_particle

    # Handle z
    if !isnothing(cathode_kinetic_energy_ref)
        # Cathode start
        vprint("Cathode start with cathode_kinetic_energy_ref = $cathode_kinetic_energy_ref eV")

        # Impact-T conversion factor in eV
        output["Bkenergy"] = cathode_kinetic_energy_ref

        # Note that this is purely a conversion factor.
        gamma = 1.0 + cathode_kinetic_energy_ref / mc2
        betac = sqrt(1 - 1 / gamma^2) * c_light

        # z equivalent
        z = -betac .* particle_group.t

        # Get z span
        z_ptp = maximum(z) - minimum(z)
        # Add tiny padding
        z_pad = 1e-20  # Tiny pad

        # Shift all particles, so that z < 0
        z_shift = -(maximum(z) + z_pad)
        z .+= z_shift

        # Starting clock shift
        t_shift = z_shift / betac

        # Suggest an emission time
        output["Temission"] = (z_ptp + 2 * z_pad) / betac

        # Change actual initial time to this shift (just set)
        output["Tini"] = t_shift

        # pz
        pz = particle_group.pz
        # check for zero pz
        all(pz .> 0) || error("pz must be positive")

        # Make sure there as at least some small pz momentum. Simply shift.
        pz_small = 10  # eV/c
        small_pz = pz .< pz_small
        pz[small_pz] .+= pz_small

        gamma_beta_z = pz ./ mc2

    else
        # Free space start
        z = particle_group.z

        t = unique(particle_group.t)
        length(t) == 1 || error("All particles must be at the same time")
        t = t[1]
        output["Tini"] = t
        output["Flagimg"] = 0  # Turn off Cathode start
        gamma_beta_z = particle_group.pz ./ mc2

        vprint("Normal start with at time $t s")
    end

    # Form data table
    dat = [particle_group.x particle_group.px./mc2 particle_group.y particle_group.py./mc2 z gamma_beta_z]

    # Save to ASCII
    open(outfile, "w") do io
        include_header && println(io, n_particle)
        writedlm(io, dat, " ")
    end

    # Return info dict
    return output
end

"""
    riffle(a, b)

Interleave two arrays.
"""
function riffle(a, b)
    return vcat(a', b')[:]
end

"""
    create_fourier_coefficients_via_fft(fz, n_coef=nothing)

Create Fourier coefficients via FFT.
"""
function create_fourier_coefficients_via_fft(fz, n_coef=nothing)
    n = length(fz)
    fcoefs = fft(fz) ./ n

    if !isnothing(n_coef)
        n_coef = min(n_coef, div(n, 2))
        fcoefs = fcoefs[1:n_coef]
    end

    return fcoefs
end

"""
    fourier_field_reconsruction(z, fcoefs, z0=0, zlen=1.0, order=0)

Reconstruct field from Fourier coefficients.
"""
function fourier_field_reconsruction(z, fcoefs, z0::Float64=0.0, zlen::Float64=1.0, order::Int=0)
    n = length(z)
    n_coef = length(fcoefs)
    k = 2π * (0:n_coef-1) ./ zlen

    # Initialize output
    f = zeros(Complex{Float64}, n)

    # Add DC term
    f .+= real(fcoefs[1])

    # Add oscillating terms
    for i in 2:n_coef
        phase = k[i] .* (z .- z0)
        if order == 0
            f .+= 2 .* real(fcoefs[i] .* exp.(im .* phase))
        elseif order == 1
            f .+= 2 .* real(im .* k[i] .* fcoefs[i] .* exp.(im .* phase))
        elseif order == 2
            f .+= 2 .* real(-k[i]^2 .* fcoefs[i] .* exp.(im .* phase))
        else
            error("order must be 0, 1, or 2")
        end
    end

    return f
end

"""
    reconstruction_error(field, fcoefs, n_coef=nothing)

Calculate reconstruction error.
"""
function reconstruction_error(field, fcoefs, n_coef=nothing)
    if !isnothing(n_coef)
        n_coef = min(n_coef, length(fcoefs))
        fcoefs = fcoefs[1:n_coef]
    end

    z = range(0, length=length(field), step=1.0)
    f_recon = fourier_field_reconsruction(z, fcoefs)

    return norm(field - f_recon) / norm(field)
end

"""
    create_fourier_data(field_mesh, component, zmirror=false, n_coef=nothing)

Create Fourier data for a field component.
"""
function create_fourier_data(field_mesh, component::String, zmirror::Bool=false, n_coef=nothing)
    # Get field data
    fz = field_mesh[component][:,1,:]
    z = field_mesh.coord_vec("z")

    # Mirror if requested
    if zmirror
        fz = vcat(reverse(fz, dims=1), fz[2:end,:])
        z = vcat(-reverse(z[2:end]), z)
    end

    # Create Fourier coefficients
    fcoefs = create_fourier_coefficients_via_fft(fz, n_coef)

    return Dict{String,Any}(
        "z" => z,
        "fz" => fz,
        "fcoefs" => fcoefs
    )
end

"""
    create_impact_solrf_fieldmap_fourier(field_mesh, zmirror=false, n_coef=nothing, err_calc=false)

Create Impact-T SOLRF fieldmap using Fourier coefficients.
"""
function create_impact_solrf_fieldmap_fourier(field_mesh, zmirror::Bool=false, n_coef=nothing, err_calc::Bool=false)
    # Get field data
    fdata = create_fourier_data(field_mesh, "Ez", zmirror, n_coef)
    z = fdata["z"]
    fz = fdata["fz"]
    fcoefs = fdata["fcoefs"]

    # Calculate error if requested
    if err_calc
        err = reconstruction_error(fz, fcoefs)
        println("Reconstruction error: $err")
    end

    # Create fieldmap
    fieldmap = Dict{String,Any}(
        "z" => z,
        "fz" => fz,
        "fcoefs" => fcoefs
    )

    return fieldmap
end

"""
    create_impact_solrf_fieldmap_derivatives(field_mesh, method="spline", spline_s=0, spline_k=5)

Create Impact-T SOLRF fieldmap using derivatives.
"""
function create_impact_solrf_fieldmap_derivatives(field_mesh, method::String="spline", spline_s::Float64=0.0, spline_k::Int=5)
    # Get field data
    z = field_mesh.coord_vec("z")
    fz = field_mesh["Ez"][:,1,:]

    # Calculate derivatives
    if method == "spline"
        dfz = spline_derivative_array(z, fz, spline_s, spline_k)
    else
        error("method must be 'spline'")
    end

    # Create fieldmap
    fieldmap = Dict{String,Any}(
        "z" => z,
        "fz" => fz,
        "dfz" => dfz
    )

    return fieldmap
end

"""
    create_impact_solrf_ele(field_mesh; zedge=0, name=nothing, scale=1, phase=0, style="fourier", n_coef=30, zmirror=nothing, spline_s=1e-6, spline_k=5, radius=0.15, x_offset=0, y_offset=0, file_id=666, output_path=nothing)

Create Impact-T SOLRF element.
"""
function create_impact_solrf_ele(field_mesh; zedge::Float64=0.0, name=nothing, scale::Float64=1.0, phase::Float64=0.0, 
                               style::String="fourier", n_coef::Int=30, zmirror=nothing, spline_s::Float64=1e-6, 
                               spline_k::Int=5, radius::Float64=0.15, x_offset::Float64=0.0, y_offset::Float64=0.0, 
                               file_id::Int=666, output_path=nothing)
    # Create fieldmap
    if style == "fourier"
        fieldmap = create_impact_solrf_fieldmap_fourier(field_mesh, zmirror, n_coef)
    elseif style == "derivatives"
        fieldmap = create_impact_solrf_fieldmap_derivatives(field_mesh, "spline", spline_s, spline_k)
    else
        error("style must be 'fourier' or 'derivatives'")
    end

    # Create element
    ele = Dict{String,Any}(
        "zedge" => zedge,
        "name" => name,
        "scale" => scale,
        "phase" => phase,
        "style" => style,
        "radius" => radius,
        "x_offset" => x_offset,
        "y_offset" => y_offset,
        "file_id" => file_id,
        "fieldmap" => fieldmap
    )

    # Write fieldmap if requested
    if !isnothing(output_path)
        if style == "fourier"
            write_impact_solrf_fieldmap_fourier(fieldmap, output_path)
        elseif style == "derivatives"
            write_impact_solrf_fieldmap_derivatives(fieldmap, output_path)
        end
    end

    return ele
end

"""
    parse_impact_emfield_cartesian(filename)

Parse Impact-T Cartesian field data.
"""
function parse_impact_emfield_cartesian(filename::String)
    # Read header
    open(filename, "r") do io
        # Skip first line
        readline(io)
        
        # Read dimensions
        line = readline(io)
        nx, ny, nz = parse.(Int, split(line))
        
        # Read grid
        line = readline(io)
        xmin, xmax = parse.(Float64, split(line))
        line = readline(io)
        ymin, ymax = parse.(Float64, split(line))
        line = readline(io)
        zmin, zmax = parse.(Float64, split(line))
        
        # Read field data
        data = readdlm(io)
        
        # Reshape data
        Ex = reshape(data[:,1], nx, ny, nz)
        Ey = reshape(data[:,2], nx, ny, nz)
        Ez = reshape(data[:,3], nx, ny, nz)
        Bx = reshape(data[:,4], nx, ny, nz)
        By = reshape(data[:,5], nx, ny, nz)
        Bz = reshape(data[:,6], nx, ny, nz)
        
        # Create fieldmap
        fieldmap = Dict{String,Any}(
            "nx" => nx,
            "ny" => ny,
            "nz" => nz,
            "xmin" => xmin,
            "xmax" => xmax,
            "ymin" => ymin,
            "ymax" => ymax,
            "zmin" => zmin,
            "zmax" => zmax,
            "Ex" => Ex,
            "Ey" => Ey,
            "Ez" => Ez,
            "Bx" => Bx,
            "By" => By,
            "Bz" => Bz
        )
        
        return fieldmap
    end
end

"""
    write_impact_emfield_cartesian(fieldmesh, filename)

Write Impact-T Cartesian field data.
"""
function write_impact_emfield_cartesian(fieldmesh, filename::String)
    # Get dimensions
    nx = length(fieldmesh.coord_vec("x"))
    ny = length(fieldmesh.coord_vec("y"))
    nz = length(fieldmesh.coord_vec("z"))
    
    # Get grid
    xmin = minimum(fieldmesh.coord_vec("x"))
    xmax = maximum(fieldmesh.coord_vec("x"))
    ymin = minimum(fieldmesh.coord_vec("y"))
    ymax = maximum(fieldmesh.coord_vec("y"))
    zmin = minimum(fieldmesh.coord_vec("z"))
    zmax = maximum(fieldmesh.coord_vec("z"))
    
    # Get field data
    Ex = fieldmesh["Ex"]
    Ey = fieldmesh["Ey"]
    Ez = fieldmesh["Ez"]
    Bx = fieldmesh["Bx"]
    By = fieldmesh["By"]
    Bz = fieldmesh["Bz"]
    
    # Write file
    open(filename, "w") do io
        println(io, "Impact-T Cartesian field data")
        println(io, "$nx $ny $nz")
        println(io, "$xmin $xmax")
        println(io, "$ymin $ymax")
        println(io, "$zmin $zmax")
        
        # Write field data
        for k in 1:nz, j in 1:ny, i in 1:nx
            println(io, @sprintf("%20.12e %20.12e %20.12e %20.12e %20.12e %20.12e",
                               Ex[i,j,k], Ey[i,j,k], Ez[i,j,k],
                               Bx[i,j,k], By[i,j,k], Bz[i,j,k]))
        end
    end
end

"""
    create_impact_emfield_cartesian_ele(field_mesh; zedge=0, name=nothing, scale=1, phase=0, radius=0.15, x_offset=0, y_offset=0, x_rotation=0, y_rotation=0, z_rotation=0, file_id=666, output_path=nothing)

Create Impact-T Cartesian field element.
"""
function create_impact_emfield_cartesian_ele(field_mesh; zedge::Float64=0.0, name=nothing, scale::Float64=1.0, 
                                           phase::Float64=0.0, radius::Float64=0.15, x_offset::Float64=0.0, 
                                           y_offset::Float64=0.0, x_rotation::Float64=0.0, y_rotation::Float64=0.0, 
                                           z_rotation::Float64=0.0, file_id::Int=666, output_path=nothing)
    # Create element
    ele = Dict{String,Any}(
        "zedge" => zedge,
        "name" => name,
        "scale" => scale,
        "phase" => phase,
        "radius" => radius,
        "x_offset" => x_offset,
        "y_offset" => y_offset,
        "x_rotation" => x_rotation,
        "y_rotation" => y_rotation,
        "z_rotation" => z_rotation,
        "file_id" => file_id
    )

    # Write fieldmap if requested
    if !isnothing(output_path)
        write_impact_emfield_cartesian(field_mesh, output_path)
    end

    return ele
end 