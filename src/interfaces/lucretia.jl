"""
    lucretia_to_data(filename, ele_name="BEGINNING", t_ref=0, exclude_dead_particles=true, verbose=false)

Load one beam in a Lucretia beam file as data for openPMD-beamphysics.

Parameters:
----------
filename : String
           Lucretia '.mat' file name.
ele_name : String
           name of the element at which the beam is located.
           An invalid name results in an error.
           If the beam file has one element, this only one beam is read.
           Default: "BEGINNING"
t_ref : Float64, optional
        reference time of the beam in seconds. Default: 0.
exclude_dead_particles : Bool, optional
                         if true, excludes dead particles.  Default: true.

Returns:
----------
data: Dict
    data for openPMD-beamphysics

Lucretia's format is described in:
    https://www.slac.stanford.edu/accel/ilc/codes/Lucretia/web/beam.html

One Lucretia ".mat" file can include beams at multiple lattice elements.
To find the beam at one element, one has to follow down this order of "fields":
    bstore >> ele_name >> Bunch >> x,

in which x is a 6-to-Np array with:
    Lucretia x  = x in m
    Lucretia px = px/p in radian
    Lucretia y  = y in m
    Lucretia py = py/p in radian
    Lucretia z  = (t - t_ref)*c in m
    Lucretia p  = p in GeV/c

Note that p is the total, not reference, momentum.

To access valid element names in a Lucretia beam file,
use the helper function list_element_names(filename).
"""
function lucretia_to_data(filename::String, ele_name::String="BEGINNING", t_ref::Float64=0.0, 
                         exclude_dead_particles::Bool=true, verbose::Bool=false)
    ele_list = list_element_names(filename)
    verbose && println(length(ele_list), " elements found in the file!")

    # Check if the element exists
    if !(ele_name in ele_list)
        error("The provided element name $ele_name does not exist in the file!")
    elseif length(ele_list) == 1
        ele_name = ele_list[1]
    end

    mdat = matread(filename)

    coords = mdat["bstore"][ele_name][1]["Bunch"][1]["x"][1]
    charges = mdat["bstore"][ele_name][1]["Bunch"][1]["Q"][1]

    Np = size(coords, 2)

    x = coords[1,:]
    px_luc = coords[2,:]  # normalized by total momentum
    y = coords[3,:]
    py_luc = coords[4,:]  # normalized by total momentum
    z_luc = coords[5,:]
    ptot = coords[6,:]  # total momentum in GeV/c

    px = px_luc .* ptot .* 1e9  # in eV/c
    py = py_luc .* ptot .* 1e9
    pz = sqrt.((ptot .* 1e9).^2 .- px.^2 .- py.^2)

    t = z_luc ./ c_light .+ t_ref

    status = fill(1, Np)

    ix = findall(ptot .== 0)
    status[ix] .= 0
    n_dead = length(ix)

    verbose && println("$Np particles detected, $n_dead found dead!")

    data = Dict{String,Any}(
        "x" => x,
        "px" => px,
        "y" => y,
        "py" => py,
        "z" => zeros(Np),
        "pz" => pz,
        "t" => t,
        "status" => status,
        "weight" => charges,
        "species" => "electron"
    )

    if exclude_dead_particles && n_dead > 0
        good = findall(ptot .> 0)

        verbose && println("Excluding $n_dead dead particles")

        for k in keys(data)
            if k != "species"
                data[k] = data[k][good]
            end
        end
    end

    return data
end

"""
    write_lucretia(P, filePath, ele_name="BEGINNING", t_ref=0, stop_ix=nothing, verbose=true)

Write a ParticleGroup beam into a Lucretia beam file.

Parameters:
----------
P: ParticleGroup
           Particles to write.
filePath : String
           Lucretia '.mat' file name to be written to.
ele_name : String
           name of the element at which the beam is located.
           Default: "BEGINNING"
t_ref : Float64, optional
        reference time of the beam in seconds. Default: 0.
stop_ix : Vector{Int}, optional
          If provided, the length must equal to the number of particles.
          See Lucretia website for details.
          Default: nothing
verbose : Bool, optional
          If true, print status messages. Default: true.

Lucretia's format is described in:
    https://www.slac.stanford.edu/accel/ilc/codes/Lucretia/web/beam.html

A general Lucretia beam file can include beams at multiple lattice elements.
This routine only saves one beam at one element (name to be specified).
Contents in the branches of the upper field structures are NOT defined.

Lucretia beam follows:
    Lucretia x  = x in m
    Lucretia px = px/p in radian
    Lucretia y  = y in m
    Lucretia py = py/p in radian
    Lucretia z  = (t - t_ref)*c in m
    Lucretia p  = p in GeV/c

Note that p is the total, not reference, momentum.
"""
function write_lucretia(P, filePath::String, ele_name::String="BEGINNING", t_ref::Float64=0.0, 
                       stop_ix=nothing, verbose::Bool=true)
    Np = P.n_particle

    ptot_nonzero = replace(P.p, 0 => 1)  # this prevents division by zero later
    # The replacement of "1" is dummy
    x_luc = P.x
    px_luc = replace(P.px ./ ptot_nonzero, Inf => 0, NaN => 0)
    y_luc = P.y
    py_luc = replace(P.py ./ ptot_nonzero, Inf => 0, NaN => 0)
    z_luc = (P.t .- t_ref) .* c_light
    ptot = P.p ./ 1e9  # total momentum in GeV/c

    if isnothing(stop_ix)
        stop_ix = zeros(Int, Np)
    elseif all(isnothing, stop_ix)
        stop_ix = zeros(Int, Np)
    else
        if length(stop_ix) != Np
            error("Length of stop_ix must equal to # particles!")
        end
    end

    # Form the lowest level of field structure
    l1 = [x_luc px_luc y_luc py_luc z_luc ptot]
    l2 = [P.weight]
    l3 = [stop_ix]

    # Wrapping upward in level fields
    w1 = [(x=l1, Q=l2, stop=l3)]
    w2 = [(Bunch=w1,)]
    w3 = [(Symbol(ele_name)=w2,)]
    w4 = Dict("bstore" => w3)

    verbose && println("writing $(length(P)) particles in the Lucretia format to $filePath")

    matwrite(filePath, w4)
end

"""
    list_element_names(filePath)

Find the element names in a Lucretia beam file.

Parameters:
----------
filePath : String
           Lucretia '.mat' file name.

Returns:
----------
Vector{String}
    A list with the element names
"""
function list_element_names(filePath::String)
    dat = matread(filePath)
    return string.(keys(dat["bstore"][1]))
end 