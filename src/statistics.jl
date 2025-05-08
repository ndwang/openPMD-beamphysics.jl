"""
    norm_emit_calc(particle_group, planes=["x"])

Calculate 2d, 4d, or 6d normalized emittance.

- `planes = ['x', 'y']` gives the 4d emittance
- `planes = ['x', 'y', 'z']` gives the 6d emittance

Momenta for each plane are taken as p+plane, e.g. 'px' for plane='x'

The normalization factor is (1/mc)^n_planes, so that the units are meters^n_planes
"""
function norm_emit_calc(particle_group, planes=["x"])
    dim = length(planes)
    vars = String[]
    for k in planes
        push!(vars, k)
        push!(vars, "p" * k)
    end

    S = particle_group.cov(vars...)
    mc2 = particle_group.mass
    norm_emit = sqrt(det(S)) / mc2^dim

    return norm_emit
end

"""
    twiss_calc(sigma_mat2)

Calculate Twiss parameters from the 2D sigma matrix (covariance matrix):
```
sigma_mat = <x,x>   <x, p>
            <p, x>  <p, p>
```

This is a simple calculation. Makes no assumptions about units.

- alpha = -<x, p>/emit
- beta  =  <x, x>/emit
- gamma =  <p, p>/emit
- emit = det(sigma_mat)
"""
function twiss_calc(sigma_mat2)
    @assert size(sigma_mat2) == (2, 2) "Bad shape: $(size(sigma_mat2)). This should be (2,2)"
    
    twiss = Dict{String,Float64}()
    emit = sqrt(det(sigma_mat2))
    twiss["alpha"] = -sigma_mat2[1, 2] / emit
    twiss["beta"] = sigma_mat2[1, 1] / emit
    twiss["gamma"] = sigma_mat2[2, 2] / emit
    twiss["emit"] = emit

    return twiss
end

"""
    twiss_ellipse_points(sigma_mat2, n_points=36)

Returns points that will trace the rms ellipse from a 2x2 covariance matrix `sigma_mat2`.

Returns:
- vec: Array with shape (2, n_points) representing x, p for the ellipse points.
"""
function twiss_ellipse_points(sigma_mat2, n_points=36)
    twiss = twiss_calc(sigma_mat2)
    A = A_mat_calc(twiss["beta"], twiss["alpha"])

    theta = range(0, 2π, length=n_points)
    zvec0 = [cos.(theta); sin.(theta)] .* sqrt(2 * twiss["emit"])

    zvec1 = A * zvec0
    return zvec1
end

"""
    twiss_match(x, p, beta0=1, alpha0=0, beta1=1, alpha1=0)

Simple Twiss matching.

Takes positions x and momenta p, and transforms them according to
initial Twiss parameters (beta0, alpha0) into final Twiss parameters (beta1, alpha1).

This is simply the matrix transformation:
```
xnew  = (   sqrt(beta1/beta0)                  0                 ) . ( x )
pnew    (  (alpha0-alpha1)/sqrt(beta0*beta1)   sqrt(beta0/beta1) )   ( p )
```

Returns new x, p
"""
function twiss_match(x, p, beta0=1, alpha0=0, beta1=1, alpha1=0)
    m11 = sqrt(beta1 / beta0)
    m21 = (alpha0 - alpha1) / sqrt(beta0 * beta1)

    xnew = x .* m11
    pnew = x .* m21 .+ p ./ m11

    return xnew, pnew
end

"""
    matched_particles(particle_group, beta=nothing, alpha=nothing, plane="x", p0c=nothing, inplace=false)

Performs simple Twiss 'matching' by applying a linear transformation to
x, px if plane == 'x', or x, py if plane == 'y'

Returns a new ParticleGroup

If inplace, a copy will not be made, and changes will be done in place.
"""
function matched_particles(particle_group, beta=nothing, alpha=nothing, plane="x", p0c=nothing, inplace=false)
    @assert plane in ("x", "y") "Invalid plane: $plane"

    P = inplace ? particle_group : copy(particle_group)

    if isnothing(p0c)
        p0c = P["mean_p"]
    end

    # Use Bmad-style coordinates
    # Get plane
    if plane == "x"
        x = P.x
        p = P.px ./ p0c
    else
        x = P.y
        p = P.py ./ p0c
    end

    # Get current Twiss
    tx = twiss_calc(cov(x, p, weights=P.weight))

    # If not specified, just fill in the current value
    if isnothing(alpha)
        alpha = tx["alpha"]
    end
    if isnothing(beta)
        beta = tx["beta"]
    end

    # New coordinates
    xnew, pnew = twiss_match(x, p, beta0=tx["beta"], alpha0=tx["alpha"], beta1=beta, alpha1=alpha)

    # Set
    if plane == "x"
        P.x = xnew
        P.px = pnew .* p0c
    else
        P.y = xnew
        P.py = pnew .* p0c
    end

    return P
end

"""
    twiss_dispersion_calc(sigma3)

Twiss and Dispersion calculation from a 3x3 sigma (covariance) matrix from particles
x, p, delta

Formulas from:
    https://uspas.fnal.gov/materials/19Knoxville/g-2/creation-and-analysis-of-beam-distributions.html

Returns a dict with:
    alpha
    beta
    gamma
    emit
    eta
    etap
"""
function twiss_dispersion_calc(sigma3)
    # Collect terms
    delta2 = sigma3[3, 3]
    xd = sigma3[1, 3]
    pd = sigma3[2, 3]

    eb = sigma3[1, 1] - xd^2 / delta2
    eg = sigma3[2, 2] - pd^2 / delta2
    ea = -sigma3[1, 2] + xd * pd / delta2

    emit = sqrt(eb * eg - ea^2)

    # Form the output dict
    d = Dict{String,Float64}()
    d["alpha"] = ea / emit
    d["beta"] = eb / emit
    d["gamma"] = eg / emit
    d["emit"] = emit
    d["eta"] = xd / delta2
    d["etap"] = pd / delta2

    return d
end

"""
    particle_twiss_dispersion(particle_group, plane="x", fraction=1, p0c=nothing)

Twiss and Dispersion calc for a ParticleGroup.

Plane must be:
    'x' or 'y'

p0c is the reference momentum. If not given, the mean p will be used.

Returns the same output dict as twiss_dispersion_calc, but with keys suffixed with the plane, i.e.:
    alpha_x
    beta_x
    gamma_x
    emit_x
    eta_x
    etap_x
    norm_emit_x
"""
function particle_twiss_dispersion(particle_group, plane="x", fraction=1, p0c=nothing)
    @assert plane in ["x", "y"]

    P = particle_group  # convenience

    if fraction < 1
        P = P[sortperm(P["J" * plane])][1:Int(fraction * length(P))]
    end

    if isnothing(p0c)
        p0c = P["mean_p"]
    end

    x = P[plane]
    xp = P["p" * plane] ./ p0c
    delta = P["p"] ./ P["mean_p"]  # - 1

    # Form covariance matrix
    sigma3 = cov([x xp delta], weights=P.weight)

    # Actual calc
    twiss = twiss_dispersion_calc(sigma3)

    # Add norm
    twiss["norm_emit"] = twiss["emit"] * P["mean_p"] / P.mass

    # Add suffix
    out = Dict{String,Float64}()
    for (k, v) in twiss
        out[k * "_" * plane] = v
    end

    return out
end

"""
    A_mat_calc(beta, alpha, inverse=false)

Returns the 1D normal form matrix from twiss parameters beta and alpha

    A =   sqrt(beta)         0
         -alpha/sqrt(beta)   1/sqrt(beta)

If inverse, the inverse will be returned:

    A^-1 =  1/sqrt(beta)     0
            alpha/sqrt(beta) sqrt(beta)

This corresponds to the linear normal form decomposition:

    M = A . Rot(theta) . A^-1

with a clockwise rotation matrix:

    Rot(theta) =  cos(theta) sin(theta)
                 -sin(theta) cos(theta)

In the Bmad manual, G_q (Bmad) = A (here) in the Linear Optics chapter.

A^-1 can be used to form normalized coordinates:
    x_bar, px_bar   = A^-1 . (x, px)
"""
function A_mat_calc(beta, alpha, inverse=false)
    a11 = sqrt(beta)
    a22 = 1 / a11
    a21 = -alpha / a11

    if inverse
        return [a22 0; -a21 a11]
    else
        return [a11 0; a21 a22]
    end
end

"""
    amplitude_calc(x, p, beta=1, alpha=0)

Simple amplitude calculation of position and momentum coordinates
relative to twiss beta and alpha.

J = (gamma x^2 + 2 alpha x p + beta p^2)/2

  = (x_bar^2 + px_bar^2)/ 2

where gamma = (1+alpha^2)/beta
"""
function amplitude_calc(x, p, beta=1, alpha=0)
    return (1 + alpha^2) / beta / 2 .* x.^2 .+ alpha .* x .* p .+ beta / 2 .* p.^2
end

"""
    particle_amplitude(particle_group, plane="x", twiss=nothing, mass_normalize=true)

Returns the normalized amplitude array from a ParticleGroup for a given plane.

Plane should be:
    'x' for the x, px plane
    'y' for the y, py plane
Other planes will work, but please check that the units make sense.

If mass_normalize (default=true), the momentum will be divided by the mass, so that the units are sqrt(m).

See: normalized_particle_coordinate
"""
function particle_amplitude(particle_group, plane="x", twiss=nothing, mass_normalize=true)
    x = particle_group[plane]
    key2 = "p" * plane

    if mass_normalize
        # Note: do not do /=, because this will replace the ParticleGroup's internal array!
        p = particle_group[key2] ./ particle_group.mass
    else
        p = particle_group[key2]
    end

    # User could supply twiss
    if isnothing(twiss)
        sigma_mat2 = cov(x, p, weights=particle_group.weight)
        twiss = twiss_calc(sigma_mat2)
    end

    J = amplitude_calc(x, p, beta=twiss["beta"], alpha=twiss["alpha"])

    return J
end

"""
    normalized_particle_coordinate(particle_group, key, twiss=nothing, mass_normalize=true)

Returns a single normalized coordinate array from a ParticleGroup

Position or momentum is determined by the key.
If the key starts with 'p', it is a momentum, else it is a position.

Intended use is for key to be one of:
    x, px, y py

and the corresponding normalized coordinates are named with suffix _bar, i.e.:
    x_bar, px_bar, y_bar, py_bar

If mass_normalize (default=true), the momentum will be divided by the mass, so that the units are sqrt(m).

These are related to action-angle coordinates
    J: amplitude
    phi: phase

    x_bar =  sqrt(2 J) cos(phi)
    px_bar = sqrt(2 J) sin(phi)

So therefore:
    J = (x_bar^2 + px_bar^2)/2
    phi = arctan(px_bar/x_bar)
and:
    <J> = norm_emit_x

Note that the center may need to be subtracted in this case.
"""
function normalized_particle_coordinate(particle_group, key, twiss=nothing, mass_normalize=true)
    # Parse key for position or momentum coordinate
    if startswith(key, "p")
        momentum = true
        key1 = key[2:end]
        key2 = key
    else
        momentum = false
        key1 = key
        key2 = "p" * key
    end

    x = particle_group[key1]

    if mass_normalize
        # Note: do not do /=, because this will replace the ParticleGroup's internal array!
        p = particle_group[key2] ./ particle_group.mass
    else
        p = particle_group[key2]
    end

    # User could supply twiss
    if isnothing(twiss)
        sigma_mat2 = cov(x, p, weights=particle_group.weight)
        twiss = twiss_calc(sigma_mat2)
    end

    A_inv = A_mat_calc(twiss["beta"], twiss["alpha"], inverse=true)

    if momentum
        return A_inv[2, 1] .* x .+ A_inv[2, 2] .* p
    else
        return A_inv[1, 1] .* x
    end
end

"""
    slice_statistics(particle_group, keys=["mean_z"], n_slice=40, slice_key=nothing)

Slices a particle group into n slices and returns statistics from each slice defined in keys.

These statistics should be scalar floats for now.

Any key can be used to slice on.
"""
function slice_statistics(particle_group, keys=["mean_z"], n_slice=40, slice_key=nothing)
    if isnothing(slice_key)
        if particle_group.in_t_coordinates
            slice_key = "z"
        else
            slice_key = "t"
        end
    end

    sdat = Dict{String,Vector{Float64}}()
    twiss_planes = Set{String}()
    twiss = Dict{String,Float64}()

    normal_keys = Set{String}()

    for k in keys
        sdat[k] = zeros(n_slice)
        if startswith(k, "twiss")
            if k == "twiss" || k == "twiss_xy"
                push!(twiss_planes, "x")
                push!(twiss_planes, "y")
            else
                plane = k[end]
                @assert plane in ("x", "y")
                push!(twiss_planes, plane)
            end
        else
            push!(normal_keys, k)
        end
    end

    twiss_plane = join(twiss_planes)  # flatten
    @assert twiss_plane in ("x", "y", "xy", "yx", "")

    for (i, pg) in enumerate(particle_group.split(n_slice, key=slice_key))
        for k in normal_keys
            sdat[k][i] = pg[k]
        end

        # Handle twiss
        if !isempty(twiss_plane)
            twiss = pg.twiss(plane=twiss_plane)
            for (k, v) in twiss
                full_key = "twiss_" * k
                if !haskey(sdat, full_key)
                    sdat[full_key] = zeros(n_slice)
                end
                sdat[full_key][i] = v
            end
        end
    end

    return sdat
end

"""
    resample_particles(particle_group, n=0, equal_weights=false)

Resamples a ParticleGroup randomly.

If n equals particle_group.n_particle or n=0,
particle indices will be scrambled.

Otherwise if weights are equal, a random subset of particles will be selected.

Otherwise if weights are not equal, particles will be sampled according to their weight using a method from Distributions.jl.

Parameters:
- n: int, default = 0
    Number to resample.
    If n = 0, this will use all particles.
- equal_weights: bool, default = false
    If true, will ensure that all particles have equal weights.

Returns:
- data: dict of ParticleGroup data
"""
function resample_particles(particle_group, n=0, equal_weights=false)
    n_old = particle_group.n_particle
    if n == 0
        n = n_old
    end

    if n > n_old
        error("Cannot supersample $n_old to $n")
    end

    weight = particle_group.weight

    # Equal weights
    if length(unique(particle_group.weight)) == 1
        ixlist = randperm(n_old)[1:n]
        weight = fill(particle_group.charge / n, n)
    # variable weights found
    elseif equal_weights || n != n_old
        # Use Distributions.jl for weighted sampling
        pk = weight ./ sum(weight)  # Probabilities
        xk = 1:length(pk)  # index
        dist = Categorical(pk)
        ixlist = rand(dist, n)
        weight = fill(particle_group.charge / n, n)
    else
        @assert n == n_old
        ixlist = randperm(n_old)
        weight = weight[ixlist]  # just scramble
    end

    data = Dict{String,Any}()
    for key in particle_group._settable_array_keys
        data[key] = particle_group[key][ixlist]
    end
    data["species"] = particle_group["species"]
    data["weight"] = weight

    return data
end

"""
    bunching(z, wavelength, weight=nothing)

Calculate the normalized bunching parameter, which is the
complex sum of weighted exponentials.

The formula for bunching is given by:

B(z, λ) = ∑ w_i e^(i k z_i) / ∑ w_i

where:
- z is the position array,
- λ is the wavelength,
- k = 2π/λ is the wave number,
- w_i are the weights.

Parameters:
- z: Array of positions where the bunching parameter is calculated
- wavelength: Wavelength of the wave
- weight: Optional weights for each exponential term. Default is 1 for all terms

Returns:
- Complex: The bunching parameter

Throws:
- ArgumentError: If wavelength is not a positive number
"""
function bunching(z::AbstractArray, wavelength::Real, weight=nothing)
    if wavelength <= 0
        throw(ArgumentError("Wavelength must be a positive number."))
    end

    if isnothing(weight)
        weight = ones(length(z))
    end
    if length(weight) != length(z)
        throw(ArgumentError("Weight array has length $(length(weight)) != length of the z array, $(length(z))"))
    end

    k = 2π / wavelength
    f = exp.(im .* k .* z)
    return sum(weight .* f) / sum(weight)
end
