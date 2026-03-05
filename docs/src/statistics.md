# Statistics

## Emittance

```julia
norm_emit_x(pg)                      # 2D x emittance [m]
norm_emit_y(pg)                      # 2D y emittance [m]
norm_emit_4d(pg)                     # 4D (x,y) emittance [m²]
norm_emit_calc(pg, ["x", "y", "z"])  # 6D emittance [m³]
```

See [`norm_emit_calc`](@ref) in the [API Reference](@ref api).

## Twiss Parameters

```julia
# From a ParticleGroup
twiss(pg; plane="x")    # returns Dict with alpha_x, beta_x, emit_x, etc.
twiss(pg; plane="xy")   # both planes

# From a 2×2 covariance matrix
twiss_calc(sigma2)      # returns NamedTuple: alpha, beta, gamma, emit

# With dispersion, from a 3×3 covariance matrix [x, px, delta]
twiss_dispersion(sigma3)               # returns alpha, beta, gamma, emit, eta, etap
twiss_dispersion(pg; plane="x")        # directly from a ParticleGroup
```

See [`twiss`](@ref), [`twiss_calc`](@ref), [`twiss_dispersion`](@ref).

## Action-Angle Variables

```julia
x_bar(pg)   # normalized x per particle  [m^(1/2)]
px_bar(pg)  # normalized px per particle [m^(1/2)]
y_bar(pg)   # normalized y per particle  [m^(1/2)]
py_bar(pg)  # normalized py per particle [m^(1/2)]
Jx(pg)      # x action amplitude per particle [m]
Jy(pg)      # y action amplitude per particle [m]

# Lower-level building blocks
particle_amplitude(pg; plane="x")
normalized_particle_coordinate(pg, "x")
A_mat_calc(beta, alpha)
A_inverse_mat_calc(beta, alpha)
amplitude_calc(x, p; beta=1, alpha=0)
```

See [`particle_amplitude`](@ref), [`normalized_particle_coordinate`](@ref),
[`A_mat_calc`](@ref), [`A_inverse_mat_calc`](@ref), [`amplitude_calc`](@ref).

## Slice Statistics

```julia
stats = slice_statistics(pg; keys=["mean_x", "sigma_x", "norm_emit_x"], n_slice=100)
# Returns a Dict{String, Vector{Float64}} with one value per slice.
```

See [`slice_statistics`](@ref).

## Statistical Helpers

`delta(pg, "x")` returns `x .- mean(x)` per particle.
`ptp(pg, "z")` returns the peak-to-peak range of `z`.
