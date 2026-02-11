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

    S = StatsBase.cov(particle_group, vars...)
    mc2 = mass(particle_group)
    norm_emit = sqrt(det(S)) / mc2^dim

    return norm_emit
end

# Convenience wrappers (moved from particles.jl)
norm_emit_x(pg::ParticleGroup) = norm_emit_calc(pg, ["x"])
norm_emit_y(pg::ParticleGroup) = norm_emit_calc(pg, ["y"])
norm_emit_4d(pg::ParticleGroup) = norm_emit_calc(pg, ["x", "y"])

"""
    twiss_calc(sigma_mat2)

Calculate Twiss parameters from the 2D sigma matrix (covariance matrix):
```
sigma_mat = <x,x>   <x, p>
            <p, x>  <p, p>
```

This is a simple calculation. Makes no assumptions about units.

Returns a NamedTuple with fields: alpha, beta, gamma, emit
"""
function twiss_calc(sigma_mat2::AbstractMatrix{T}) where T
    @assert size(sigma_mat2) == (2, 2) "Bad shape: $(size(sigma_mat2)). This should be (2,2)"

    emit = sqrt(det(sigma_mat2))
    return (;
        alpha = -sigma_mat2[1, 2] / emit,
        beta  =  sigma_mat2[1, 1] / emit,
        gamma =  sigma_mat2[2, 2] / emit,
        emit  =  emit
    )
end

"""
    twiss_dispersion(sigma3)

Twiss and Dispersion calculation from a 3x3 sigma (covariance) matrix from particles
x, p, delta

Returns a NamedTuple with fields: alpha, beta, gamma, emit, eta, etap
"""
function twiss_dispersion(sigma3::Matrix)
    delta2 = sigma3[3, 3]
    xd = sigma3[1, 3]
    pd = sigma3[2, 3]

    eb = sigma3[1, 1] - xd^2 / delta2
    eg = sigma3[2, 2] - pd^2 / delta2
    ea = -sigma3[1, 2] + xd * pd / delta2

    emit = sqrt(eb * eg - ea^2)

    return (;
        alpha = ea / emit,
        beta  = eb / emit,
        gamma = eg / emit,
        emit  = emit,
        eta   = xd / delta2,
        etap  = pd / delta2
    )
end

"""
    twiss_dispersion(particle_group; plane="x", fraction=1, p0c=nothing)

Twiss and Dispersion calc for a ParticleGroup.

Plane must be 'x' or 'y'.

p0c is the reference momentum. If not given, the mean p will be used.

Returns a NamedTuple with fields suffixed by plane:
    alpha, beta, gamma, emit, eta, etap, norm_emit
"""
function twiss_dispersion(particle_group; plane="x", fraction=1, p0c=nothing)
    @assert plane in ["x", "y"]

    P = particle_group

    if fraction < 1
        P = P[sortperm(P["J" * plane])][1:Int(fraction * length(P))]
    end

    if isnothing(p0c)
        p0c = P["mean_p"]
    end

    x = P[plane]
    xp = P["p" * plane] ./ p0c
    d = P["p"] ./ P["mean_p"]  # - 1

    sigma3 = StatsBase.cov([x xp d], weights(P.weight))

    tw = twiss_dispersion(sigma3)

    # Add norm_emit (no plane suffix yet)
    norm_emit = tw.emit * P["mean_p"] / mass(P)

    return (; tw..., norm_emit = norm_emit)
end

"""
    twiss(pg::ParticleGroup; plane="x", fraction=1, p0c=nothing)

Compute Twiss parameters for one or more planes.

Returns a Dict with keys suffixed by plane, e.g. "alpha_x", "beta_x", etc.
"""
function twiss(pg::ParticleGroup; plane="x", fraction=1, p0c=nothing)
    d = Dict{String,Float64}()
    for p in plane
        result = twiss_dispersion(pg, plane=string(p), fraction=fraction, p0c=p0c)
        for (k, v) in pairs(result)
            d[string(k) * "_" * string(p)] = v
        end
    end
    return d
end

"""
    A_mat_calc(beta, alpha)

Returns the 1D normal form matrix from twiss parameters beta and alpha

    A =   sqrt(beta)         0
         -alpha/sqrt(beta)   1/sqrt(beta)
"""
function A_mat_calc(beta, alpha)
    a11 = sqrt(beta)
    a22 = 1 / a11
    a21 = -alpha / a11
    return [a11 0; a21 a22]
end

"""
    A_inverse_mat_calc(beta, alpha)

Returns the inverse of the 1D normal form matrix from twiss parameters beta and alpha

    A^-1 =  1/sqrt(beta)     0
            alpha/sqrt(beta) sqrt(beta)
"""
function A_inverse_mat_calc(beta, alpha)
    a11 = sqrt(beta)
    a22 = 1 / a11
    a21 = -alpha / a11
    return [a22 0; -a21 a11]
end

"""
    amplitude_calc(x, p; beta=1, alpha=0)

Simple amplitude calculation of position and momentum coordinates
relative to twiss beta and alpha.

J = (gamma x^2 + 2 alpha x p + beta p^2)/2
"""
function amplitude_calc(x, p; beta=1, alpha=0)
    return (1 + alpha^2) / beta / 2 .* x.^2 .+ alpha .* x .* p .+ beta / 2 .* p.^2
end

"""
    particle_amplitude(particle_group; plane="x", twiss=nothing, mass_normalize=true)

Returns the normalized amplitude array from a ParticleGroup for a given plane.
"""
function particle_amplitude(particle_group; plane="x", twiss=nothing, mass_normalize=true)
    x = particle_group[plane]
    key2 = "p" * plane

    if mass_normalize
        p = particle_group[key2] ./ mass(particle_group)
    else
        p = particle_group[key2]
    end

    if isnothing(twiss)
        sigma_mat2 = StatsBase.cov(hcat(x, p), weights(particle_group.weight), 1)
        twiss = twiss_calc(sigma_mat2)
    end

    J = amplitude_calc(x, p; beta=twiss.beta, alpha=twiss.alpha)

    return J
end

"""
    normalized_particle_coordinate(particle_group, key; twiss=nothing, mass_normalize=true)

Returns a single normalized coordinate array from a ParticleGroup.
"""
function normalized_particle_coordinate(particle_group, key; twiss=nothing, mass_normalize=true)
    if startswith(key, "p")
        is_momentum = true
        key1 = key[2:end]
        key2 = key
    else
        is_momentum = false
        key1 = key
        key2 = "p" * key
    end

    x = particle_group[key1]

    if mass_normalize
        p = particle_group[key2] ./ mass(particle_group)
    else
        p = particle_group[key2]
    end

    if isnothing(twiss)
        sigma_mat2 = StatsBase.cov(hcat(x, p), weights(particle_group.weight), 1)
        twiss = twiss_calc(sigma_mat2)
    end

    A_inv = A_inverse_mat_calc(twiss.beta, twiss.alpha)

    if is_momentum
        return A_inv[2, 1] .* x .+ A_inv[2, 2] .* p
    else
        return A_inv[1, 1] .* x
    end
end

"""
    slice_statistics(particle_group; keys=["mean_z"], n_slice=40, slice_key=nothing)

Slices a particle group into n slices and returns statistics from each slice.
"""
function slice_statistics(particle_group; keys=["mean_z"], n_slice=40, slice_key=nothing)
    if isnothing(slice_key)
        if in_t_coordinates(particle_group)
            slice_key = "z"
        else
            slice_key = "t"
        end
    end
    slices = split_particles(particle_group; n_chunks=n_slice, key=slice_key)

    statistics = Dict{String,Vector{Float64}}()

    regular_keys = filter(k -> !startswith(k, "twiss"), keys)
    for key in regular_keys
        statistics[key] = zeros(n_slice)
        for (i, slice) in enumerate(slices)
            statistics[key][i] = slice[key]
        end
    end

    twiss_keys = filter(k -> startswith(k, "twiss"), keys)
    if !isempty(twiss_keys)
        planes = Set{String}()
        for key in twiss_keys
            if key == "twiss" || key == "twiss_xy"
                push!(planes, "x", "y")
            else
                p = string(key[end])
                @assert p in ("x", "y") "Invalid Twiss plane: $p"
                push!(planes, p)
            end
        end

        for pl in planes
            for (i, slice) in enumerate(slices)
                result = twiss_dispersion(slice; plane=pl)
                for (k, v) in pairs(result)
                    skey = "twiss_$(k)_$(pl)"
                    if !haskey(statistics, skey)
                        statistics[skey] = zeros(n_slice)
                    end
                    statistics[skey][i] = v
                end
            end
        end
    end

    return statistics
end
