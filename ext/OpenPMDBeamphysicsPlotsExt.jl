module OpenPMDBeamphysicsPlotsExt

using Plots
using OpenPMDBeamphysics
using OpenPMDBeamphysics: pg_units, plottable_array, nice_array, nice_scale_prefix,
    mathlabel, ptp, limits, in_t_coordinates, slice_statistics

import OpenPMDBeamphysics: slice_plot, density_plot, marginal_plot, density_and_slice_plot
import OpenPMDBeamphysics: fieldmesh_plot, fieldmesh_plot_onaxis
import OpenPMDBeamphysics: FieldMesh, coord_vec, mins, maxs, is_pure_magnetic, is_static, pg_units, nice_array

using StatsBase: weights, fit, Histogram
import StatsBase
import Statistics

# Default colormaps
const CMAP0 = :viridis
const CMAP1 = :plasma

# Sigma ellipse for 2x2 covariance matrix (1-sigma)
function _sigma_ellipse(sigma_mat2; n=100)
    a, b = sigma_mat2[1, 1], sigma_mat2[1, 2]
    d = sigma_mat2[2, 2]
    tr = (a + d) / 2
    disc = sqrt((a - d)^2 / 4 + b^2)
    λ1, λ2 = tr + disc, tr - disc
    angle = atan(2b, a - d) / 2
    t = range(0, 2π, length=n)
    ca, sa = cos(angle), sin(angle)
    r1, r2 = sqrt(max(λ1, 0)), sqrt(max(λ2, 0))
    x = r1 .* ca .* cos.(t) .- r2 .* sa .* sin.(t)
    y = r1 .* sa .* cos.(t) .+ r2 .* ca .* sin.(t)
    return x, y
end

"""
    slice_plot(particle_group, keys...; n_slice=40, slice_key=nothing, xlim=nothing, ylim=nothing, tex=true, kwargs...)

Complete slice plotting routine. Will plot the density of the slice key on the right axis.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `keys`: Keys to calculate the statistics, e.g. `sigma_x`
- `n_slice`: Number of slices (default: 40)
- `slice_key`: Dimension to slice in (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `ylim`: Manual setting of y-axis limits
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function slice_plot(particle_group, keys...; n_slice=40, slice_key=nothing, xlim=nothing, ylim=nothing, tex=true, kwargs...)
    # Determine slice key if not provided
    if isnothing(slice_key)
        if in_t_coordinates(particle_group)
            slice_key = "z"
        else
            slice_key = "t"
        end
    end

    # Convert keys to a vector of strings
    keys_vec = collect(keys)
    push!(keys_vec, "mean_" * slice_key)
    push!(keys_vec, "ptp_" * slice_key)
    push!(keys_vec, "charge")

    # Handle delta prefix
    has_delta_prefix = false
    if startswith(slice_key, "delta_")
        slice_key = slice_key[7:end]
        has_delta_prefix = true
    end

    # Get slice data
    x_key = "mean_" * slice_key
    slice_dat = slice_statistics(particle_group; keys=keys_vec, n_slice=n_slice, slice_key=slice_key)
    slice_dat["density"] = slice_dat["charge"] ./ slice_dat["ptp_" * slice_key]
    y2_key = "density"

    # X-axis
    x = slice_dat["mean_" * slice_key]
    if has_delta_prefix
        x .-= particle_group["mean_" * slice_key]
        slice_key = "delta_" * slice_key
    end

    x, f1, p1, xmin, xmax = plottable_array(x; lim=xlim)
    ux = p1 * string(pg_units(slice_key))

    # Y-axis
    ulist = [pg_units(k) for k in keys]
    uy = ulist[1]
    if !all(u == uy for u in ulist)
        error("Incompatible units: $ulist")
    end
    uy_str = string(uy)

    ymin = maximum(minimum(slice_dat[k]) for k in keys)
    ymax = maximum(maximum(slice_dat[k]) for k in keys)

    _, f2, p2, ymin, ymax = plottable_array([ymin, ymax]; lim=ylim)
    uy_str = p2 * uy_str

    # Calculate plot limits
    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    plot_ylim = isnothing(ylim) ? (ymin/f2, ymax/f2) : ylim

    # Create plot
    p = plot(; xlim=plot_xlim, ylim=plot_ylim, kwargs...)

    # Main curves
    color = length(keys) == 1 ? :black : nothing

    for k in keys
        label = mathlabel(k, units=uy_str, tex=tex)
        plot!(p, x, slice_dat[k] ./ f2, label=label, color=color)
    end

    if length(keys) > 1
        plot!(p, legend=true)
    end

    # Density on right axis
    y2, _, prey2, _, _ = plottable_array(slice_dat[y2_key]; lim=nothing)

    # Convert to Amps if possible
    y2_units = "C/" * string(pg_units(x_key))
    if y2_units == "C/s"
        y2_units = "A"
    end
    y2_units = prey2 * y2_units

    # Labels
    labelx = mathlabel(slice_key, units=ux, tex=tex)
    labely = mathlabel(keys..., units=uy, tex=tex)
    labely2 = mathlabel(y2_key, units=y2_units, tex=tex)

    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Right axis plot
    p2 = twinx(p)
    ylabel!(p2, labely2)
    plot!(p2, x, y2, fillrange=0, fillalpha=0.2, color=:gray, w=0)

    return p
end

"""
    density_plot(particle_group, key="x"; bins=nothing, xlim=nothing, tex=true, kwargs...)

1D density plot.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key`: Key to plot (default: "x")
- `bins`: Number of bins (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function density_plot(particle_group, key="x"; bins=nothing, xlim=nothing, tex=true, kwargs...)
    if isnothing(bins)
        n = length(particle_group)
        bins = Int(sqrt(n/4))
    end

    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(particle_group[key]; lim=xlim)
    w = particle_group["weight"]
    u1 = string(pg_units(key))
    ux = p1 * u1

    # Create label
    labelx = mathlabel(key, units=ux, tex=tex)

    # Create histogram
    hist = fit(Histogram, x, weights(w), nbins=bins)
    hist_x = hist.edges[1][1:end-1] .+ diff(hist.edges[1])./2
    hist_width = diff(hist.edges[1])
    hist_y, hist_f, hist_prefix = nice_array(hist.weights ./ hist_width)

    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    # Create plot
    p = bar(hist_x, hist_y, bar_width=hist_width, lw=0,xlim=plot_xlim, color=:gray; kwargs...)

    # Set labels
    xlabel!(p, labelx)
    if u1 == "s"
        _, hist_prefix = nice_scale_prefix(hist_f/f1)
        ylabel!(p, "$(hist_prefix)A")
    else
        ylabel!(p, "$(hist_prefix)C/$ux")
    end

    return p
end

"""
    marginal_plot(particle_group, key1="t", key2="p"; bins=nothing, xlim=nothing, ylim=nothing, tex=true, ellipse=false, kwargs...)

Density plot and projections.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key1`: Key for x-axis (default: "t")
- `key2`: Key for y-axis (default: "p")
- `bins`: Number of bins (default: automatically determined)
- `xlim`: Manual setting of x-axis limits
- `ylim`: Manual setting of y-axis limits
- `tex`: Use TEX for labels (default: true)
- `ellipse`: Plot sigma matrix ellipse (default: false)

# Returns
- `p`: Plots.Plot object
"""
function marginal_plot(particle_group, key1="t", key2="p"; bins=nothing, xlim=nothing, ylim=nothing, tex=true, ellipse=false, kwargs...)
    if isnothing(bins)
        n = length(particle_group)
        bins = Int(sqrt(n/4))
    end

    # Get data
    x = particle_group[key1]
    y = particle_group[key2]

    if length(x) == 1
        bins = 100

        if isnothing(xlim)
            x0 = x[1]
            if isapprox(x0, 0.0)
                xlim = (-1, 1)
            else
                xlim = sort([0.9*x0, 1.1*x0])
            end
        end
        if isnothing(ylim)
            y0 = y[1]
            if isapprox(y0, 0.0)
                ylim = (-1, 1)
            else
                ylim = sort([0.9*y0, 1.1*y0])
            end
        end
    end

    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(x; lim=xlim)
    y, f2, p2, ymin, ymax = plottable_array(y; lim=ylim)
    w = particle_group["weight"]

    u1 = string(pg_units(key1))
    u2 = string(pg_units(key2))
    ux = p1 * u1
    uy = p2 * u2

    # Create labels
    labelx = mathlabel(key1, units=ux, tex=tex)
    labely = mathlabel(key2, units=uy, tex=tex)

    # Calculate plot limits
    plot_xlim = isnothing(xlim) ? (xmin/f1, xmax/f1) : xlim
    plot_ylim = isnothing(ylim) ? (ymin/f2, ymax/f2) : ylim

    # Create plot
    p = plot(; xlim=plot_xlim, ylim=plot_ylim, kwargs...)

    if all(isnan, x)
        annotate!(p, 0.5, 0.5, text("$key1 is all NaN", :center))
        return p
    end
    if all(isnan, y)
        annotate!(p, 0.5, 0.5, text("$key2 is all NaN", :center))
        return p
    end

    # Main plot
    if length(x) == 1
        scatter!(p, x, y, xlim=plot_xlim, ylim=plot_ylim)
    else
        histogram2d!(p, x, y, weights(w), bins=bins, color=CMAP0, clims=(1e-20, Inf), xlim=plot_xlim, ylim=plot_ylim, colorbar=false)
    end

    # Add ellipse if requested
    if ellipse
        sigma_mat2 = StatsBase.cov(particle_group, key1, key2)
        x_ellipse, y_ellipse = _sigma_ellipse(sigma_mat2)
        x_ellipse .+= Statistics.mean(particle_group, key1)
        y_ellipse .+= Statistics.mean(particle_group, key2)
        plot!(p, x_ellipse ./ f1, y_ellipse ./ f2, color=:red, label="")
    end

    # Add histograms
    hist_x = fit(Histogram, x, weights(w), nbins=bins)
    hist_y = fit(Histogram, y, weights(w), nbins=bins)

    # Top histogram
    hist_x_centers = hist_x.edges[1][1:end-1] .+ diff(hist_x.edges[1])./2
    hist_x_width = diff(hist_x.edges[1])
    hist_x_y, hist_x_f, hist_x_prefix = nice_array(hist_x.weights ./ hist_x_width)

    p_top = bar(hist_x_centers, hist_x_y, bar_width=hist_x_width, lw=0,color=:grey, xlim=plot_xlim,
                xticks=false, xlabel="", legend=false)
    if u1 == "s"
        _, hist_x_prefix = nice_scale_prefix(hist_x_f/f1)
        ylabel!(p_top, "$(hist_x_prefix)A")
    else
        ylabel!(p_top, "$(hist_x_prefix)" * mathlabel("C/$ux"))
    end

    # Side histogram
    hist_y_centers = hist_y.edges[1][1:end-1] .+ diff(hist_y.edges[1])./2
    hist_y_width = diff(hist_y.edges[1])
    hist_y_y, hist_y_f, hist_y_prefix = nice_array(hist_y.weights ./ hist_y_width)

    p_side = bar(hist_y_centers, hist_y_y, bar_width=hist_y_width, lw=0,color=:grey, ylim=plot_ylim, orientation=:h,
                 yticks=false, ylabel="", legend=false)
    xlabel!(p_side, "$(hist_y_prefix)" * mathlabel("C/$uy"))

    # Set labels and limits
    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Combine plots
    layout = @layout [a{0.3h} _; b c{0.3w}]
    return plot(p_top, p, p_side, layout=layout)
end

"""
    density_and_slice_plot(particle_group, key1="t", key2="p"; stat_keys=["norm_emit_x", "norm_emit_y"], bins=100, n_slice=30, tex=true)

Combined density and slice plot.

# Arguments
- `particle_group`: ParticleGroup object to plot
- `key1`: Key for x-axis (default: "t")
- `key2`: Key for y-axis (default: "p")
- `stat_keys`: Statistics keys to plot (default: ["norm_emit_x", "norm_emit_y"])
- `bins`: Number of bins (default: 100)
- `n_slice`: Number of slices (default: 30)
- `tex`: Use TEX for labels (default: true)

# Returns
- `p`: Plots.Plot object
"""
function density_and_slice_plot(particle_group, key1="t", key2="p"; stat_keys=["norm_emit_x", "norm_emit_y"], bins=100, n_slice=30, tex=true)
    # Scale to nice units
    x, f1, p1, xmin, xmax = plottable_array(particle_group[key1])
    y, f2, p2, ymin, ymax = plottable_array(particle_group[key2])
    w = particle_group["weight"]

    u1 = string(pg_units(key1))
    u2 = string(pg_units(key2))
    ux = p1 * u1
    uy = p2 * u2

    # Create labels
    labelx = mathlabel(key1, units=ux, tex=tex)
    labely = mathlabel(key2, units=uy, tex=tex)

    # Create plot
    p = plot()

    xlabel!(p, labelx)
    ylabel!(p, labely)

    # Create 2D histogram
    H = fit(Histogram, (x, y), weights(w), nbins=bins)
    heatmap!(p, H.edges[1], H.edges[2], H.weights', color=CMAP0, clims=(1e-16, Inf))

    # Get slice data
    slice_dat = slice_statistics(
        particle_group,
        n_slice=n_slice,
        slice_key=key1,
        keys=vcat(stat_keys, ["ptp_" * key1, "mean_" * key1, "charge"])
    )

    slice_dat["density"] = slice_dat["charge"] ./ slice_dat["ptp_" * key1]

    # Add slice statistics
    p2 = twinx(p)
    x2 = slice_dat["mean_" * key1] ./ f1
    ulist = [pg_units(k) for k in stat_keys]

    max2 = maximum(ptp(slice_dat[k]) for k in stat_keys)
    f3, p3 = nice_scale_prefix(max2)

    u2 = ulist[1]
    @assert all(u == u2 for u in ulist)
    u2 = p3 * string(u2)
    labely2 = mathlabel(stat_keys..., units=u2, tex=tex)

    for k in stat_keys
        label = mathlabel(k, units=u2, tex=tex)
        plot!(p2, x2, slice_dat[k] ./ f3, label=label)
    end

    plot!(p2, legend=true)
    ylabel!(p2, labely2)

    # Add density
    y2 = slice_dat["density"]
    y2 = y2 .* max2 ./ maximum(y2) ./ f3 ./ 2
    plot!(p2, x2, y2, fillrange=0, fillalpha=0.2, color=:gray, w=0)
    return p
end


# ── FieldMesh plotting ─────────────────────────────────────────────────────

"""
    fieldmesh_plot(fm::FieldMesh; component=nothing, mirror=false, cmap=:plasma, kwargs...)

2D heatmap of a FieldMesh field component.

For cylindrical geometry, plots the (r, z) plane at theta=0. Z is on the x-axis,
r on the y-axis (convention for accelerators).

For rectangular geometry, plots the (x, z) plane at y=0 by default. Pass `y=value`
to select a different slice value.

# Arguments
- `component`: Field component to plot (e.g. `"Bz"`, `"B"`, `"Ez"`).
  Defaults to `"B"` for pure magnetic, `"E"` otherwise.
- `mirror`: For cylindrical only. If `true`, mirrors data to show negative-r side.
  `Br` and `Er` components are sign-flipped on the mirrored side.
- `cmap`: Colormap (default `:plasma`)
- `kwargs...`: Forwarded to `Plots.heatmap`

# Returns
`Plots.Plot`
"""
function fieldmesh_plot(fm::FieldMesh;
    component::Union{String, Nothing} = nothing,
    mirror::Bool = false,
    cmap = :plasma,
    kwargs...
)
    if isnothing(component)
        component = is_pure_magnetic(fm) ? "B" : "E"
    end

    if fm.geometry == "cylindrical"
        return _fieldmesh_plot_cylindrical(fm, component; mirror=mirror, cmap=cmap, kwargs...)
    elseif fm.geometry == "rectangular"
        return _fieldmesh_plot_rectangular(fm, component; cmap=cmap, kwargs...)
    else
        throw(ArgumentError("Unsupported geometry: $(fm.geometry)"))
    end
end

function _fieldmesh_plot_cylindrical(fm::FieldMesh, component::String;
    mirror::Bool = false, cmap = :plasma, kwargs...)

    r = coord_vec(fm, "r")
    z = coord_vec(fm, "z")

    dat = real.(fm[component][:, 1, :])   # shape (nr, nz), real part for display

    if mirror
        if r[1] != 0.0
            throw(ArgumentError("mirror=true only available when r=0 is at grid origin"))
        end
        # Sign convention: radial components flip sign, axial/magnitude do not
        sign = (component in ("Br", "Er")) ? -1 : 1
        dat = vcat(sign .* reverse(dat; dims=1), dat)
        r = vcat(-reverse(r), r)
    end

    u = pg_units(component)
    dat_scaled, fac, prefix = nice_array(dat)
    clabel = mathlabel(component; units=prefix * string(u))
    xlabel = "\$z\$ (m)"
    ylabel = "\$r\$ (m)"

    return heatmap(z, r, dat_scaled;
        xlabel=xlabel, ylabel=ylabel, colorbar_title=clabel,
        color=cmap, kwargs...)
end

function _fieldmesh_plot_rectangular(fm::FieldMesh, component::String;
    cmap = :plasma, kwargs...)

    # Pop slice coordinate from kwargs; default to y=0
    kw = Dict(kwargs)
    slice_axis, slice_val = "y", 0.0
    for ax in fm.axis_labels
        if haskey(kw, Symbol(ax))
            slice_axis = ax
            slice_val = Float64(pop!(kw, Symbol(ax)))
            break
        end
    end

    x = coord_vec(fm, "x")
    y = coord_vec(fm, "y")
    z = coord_vec(fm, "z")

    dat3d = real.(fm[component])   # (nx, ny, nz)

    if slice_axis == "y"
        iy = argmin(abs.(y .- slice_val))
        dat = dat3d[:, iy, :]      # (nx, nz)
        hax, vax = z, x
        xlabel, ylabel = "\$z\$ (m)", "\$x\$ (m)"
    elseif slice_axis == "x"
        ix = argmin(abs.(x .- slice_val))
        dat = dat3d[ix, :, :]      # (ny, nz)
        hax, vax = z, y
        xlabel, ylabel = "\$z\$ (m)", "\$y\$ (m)"
    else  # z
        iz = argmin(abs.(z .- slice_val))
        dat = dat3d[:, :, iz]      # (nx, ny)
        hax, vax = x, y
        xlabel, ylabel = "\$x\$ (m)", "\$y\$ (m)"
    end

    u = pg_units(component)
    dat_scaled, fac, prefix = nice_array(dat)
    clabel = mathlabel(component; units=prefix * string(u))

    return heatmap(hax, vax, dat_scaled;
        xlabel=xlabel, ylabel=ylabel, colorbar_title=clabel,
        color=cmap, kw...)
end

"""
    fieldmesh_plot_onaxis(fm::FieldMesh; kwargs...)

1D on-axis field plot for a FieldMesh.

For cylindrical geometry, plots `Bz` and/or `Ez` at r=0.
For rectangular geometry, plots all present field components along z at x=y=0.

# Arguments
- `kwargs...`: Forwarded to `Plots.plot`

# Returns
`Plots.Plot`
"""
function fieldmesh_plot_onaxis(fm::FieldMesh; kwargs...)
    if fm.geometry == "cylindrical"
        return _fieldmesh_onaxis_cylindrical(fm; kwargs...)
    elseif fm.geometry == "rectangular"
        return _fieldmesh_onaxis_rectangular(fm; kwargs...)
    else
        throw(ArgumentError("Unsupported geometry: $(fm.geometry)"))
    end
end

function _fieldmesh_onaxis_cylindrical(fm::FieldMesh; kwargs...)
    z = coord_vec(fm, "z")
    has_Ez = haskey(fm.components, "electricField/z")
    has_Bz = haskey(fm.components, "magneticField/z")

    p = plot(; xlabel="\$z\$ (m)", kwargs...)

    if has_Ez
        Ez0 = fm["Ez"][1, 1, :]
        Ez_scaled, fac, prefix = nice_array(Ez0)
        plot!(p, z, Ez_scaled;
            label=mathlabel("Ez"; units=prefix * string(pg_units("Ez"))),
            color=:black)
    end
    if has_Bz
        Bz0 = fm["Bz"][1, 1, :]
        Bz_scaled, fac, prefix = nice_array(Bz0)
        if has_Ez
            p2 = twinx(p)
            plot!(p2, z, Bz_scaled;
                label=mathlabel("Bz"; units=prefix * string(pg_units("Bz"))),
                color=:blue, ylabel=mathlabel("Bz"; units=prefix * string(pg_units("Bz"))))
        else
            plot!(p, z, Bz_scaled;
                label=mathlabel("Bz"; units=prefix * string(pg_units("Bz"))),
                color=:blue,
                ylabel=mathlabel("Bz"; units=prefix * string(pg_units("Bz"))))
        end
    end

    return p
end

function _fieldmesh_onaxis_rectangular(fm::FieldMesh; kwargs...)
    z = coord_vec(fm, "z")
    x = coord_vec(fm, "x")
    y = coord_vec(fm, "y")
    ix = argmin(abs.(x))
    iy = argmin(abs.(y))

    p = plot(; xlabel="\$z\$ (m)", kwargs...)

    field_keys = [("electricField/x", "Ex"), ("electricField/y", "Ey"), ("electricField/z", "Ez"),
                  ("magneticField/x", "Bx"), ("magneticField/y", "By"), ("magneticField/z", "Bz")]

    for (comp, alias) in field_keys
        haskey(fm.components, comp) || continue
        fdat = fm[alias][ix, iy, :]
        all(iszero, fdat) && continue
        f_scaled, fac, prefix = nice_array(fdat)
        plot!(p, z, f_scaled;
            label=mathlabel(alias; units=prefix * string(pg_units(alias))))
    end

    return p
end

end # module
