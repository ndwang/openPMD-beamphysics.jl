"""
Label generation utilities for plotting and reporting.
"""

# Dictionary of TeX labels
const TEXLABEL = Dict{String, String}(
    # 'status'
    "t" => "t",
    "energy" => "E",
    "kinetic_energy" => raw"E_\text{kinetic}",
    # 'mass',
    "higher_order_energy_spread" => raw"\sigma_{E_{(2)}}",
    "higher_order_energy" => raw"E_{(2)}",
    "Ex" => "E_x",
    "Ey" => "E_y",
    "Ez" => "E_z",
    "Bx" => "B_x",
    "By" => "B_y",
    "Bz" => "B_z",
    "Etheta" => raw"E_{\theta}",
    "Btheta" => raw"B_{\theta}",
    "px" => "p_x",
    "py" => "p_y",
    "pz" => "p_z",
    "p" => "p",
    "pr" => "p_r",
    "ptheta" => raw"p_{\theta}",
    "x" => "x",
    "y" => "y",
    "z" => "z",
    "r" => "r",
    "Jx" => "J_x",
    "Jy" => "J_y",
    "beta" => raw"\beta",
    "beta_x" => raw"\beta_x",
    "beta_y" => raw"\beta_y",
    "beta_z" => raw"\beta_z",
    "gamma" => raw"\gamma",
    "theta" => raw"\theta",
    "charge" => "Q",
    "twiss_alpha_x" => raw"\text{Twiss}\ \alpha_x",
    "twiss_beta_x" => raw"\text{Twiss}\ \beta_x",
    "twiss_gamma_x" => raw"\text{Twiss}\ \gamma_x",
    "twiss_eta_x" => raw"\text{Twiss}\ \eta_x",
    "twiss_etap_x" => raw"\text{Twiss}\ \eta'_x",
    "twiss_emit_x" => raw"\text{Twiss}\ \epsilon_{x}",
    "twiss_norm_emit_x" => raw"\text{Twiss}\ \epsilon_{n, x}",
    "twiss_alpha_y" => raw"\text{Twiss}\ \alpha_y",
    "twiss_beta_y" => raw"\text{Twiss}\ \beta_y",
    "twiss_gamma_y" => raw"\text{Twiss}\ \gamma_y",
    "twiss_eta_y" => raw"\text{Twiss}\ \eta_y",
    "twiss_etap_y" => raw"\text{Twiss}\ \eta'_y",
    "twiss_emit_y" => raw"\text{Twiss}\ \epsilon_{y}",
    "twiss_norm_emit_y" => raw"\text{Twiss}\ \epsilon_{n, y}",
    # 'species_charge',
    # 'weight',
    "average_current" => raw"I_{av}",
    "norm_emit_x" => raw"\epsilon_{n, x}",
    "norm_emit_y" => raw"\epsilon_{n, y}",
    "norm_emit_4d" => raw"\epsilon_{4D}",
    "Lz" => "L_z",
    "xp" => "x'",
    "yp" => "y'",
    "x_bar" => raw"\overline{x}",
    "px_bar" => raw"\overline{p_x}",
    "y_bar" => raw"\overline{y}",
    "py_bar" => raw"\overline{p_y}",
)

"""
    parse_bunching_str(key::String)

Extract wavelength from a bunching key string, e.g. "bunching_1e-6" → 1e-6.
"""
function parse_bunching_str(key::String)
    m = match(r"bunching_(.+)", key)
    isnothing(m) && error("Cannot parse bunching key: $key")
    return parse(Float64, m[1])
end

"""
    texlabel(key::String)

Returns a tex label from a proper attribute name.

# Arguments
- `key::String`: any pmd_beamphysics attribute

# Returns
- `tex::Union{String, Nothing}`: A TeX string if applicable, otherwise will return nothing

# Examples
```julia
texlabel("cov_x__px")  # returns: "\\left<x, p_x\\right>"
```

# Notes
See matplotlib documentation for more information on mathtext:
    https://matplotlib.org/stable/tutorials/text/mathtext.html
"""
function texlabel(key::AbstractString)
    # Basic cases
    if haskey(TEXLABEL, key)
        return TEXLABEL[key]
    end

    # Operators
    for prefix in ["sigma_", "mean_", "min_", "max_", "ptp_", "delta_"]
        if startswith(key, prefix)
            pre = prefix[1:end-1]
            key0 = key[length(prefix)+1:end]
            tex0 = texlabel(key0)

            if pre == "min"
                return "\\min($tex0)"
            elseif pre == "max"
                return "\\max($tex0)"
            elseif pre == "sigma"
                return "\\sigma_{ $tex0 }"
            elseif pre == "delta"
                return "$tex0 - \\left<$tex0\\right>"
            elseif pre == "mean"
                return "\\left<$tex0\\right>"
            elseif pre == "ptp"
                return "\\Delta $tex0"
            end
        end
    end

    if startswith(key, "cov_")
        subkeys = split(key[5:end], "__")
        tex0 = texlabel(subkeys[1])
        tex1 = texlabel(subkeys[2])
        return "\\left<$tex0, $tex1\\right>"
    end

    if startswith(key, "bunching")
        wavelength = parse_bunching_str(key)
        x, _, prefix = nice_array(wavelength)
        return "\\mathrm{bunching~at}~$x\\mathrm{ $prefix m }"
    end

    return "\\mathrm{ $key }"
end

"""
    mathlabel(keys::String...; units=nothing, tex::Bool=true)

Helper function to return label with optional units from an arbitrary number of keys.

# Arguments
- `keys::String...`: any pmd_beamphysics attributes
- `units`: units to be cast to string (optional)
- `tex::Bool`: if true, a mathtext (TeX) string wrapped in \$ will be returned (default: true)

# Returns
- `label::String`: A TeX string if applicable, otherwise will return a plain string

# Examples
```julia
mathlabel("x_bar", "sigma_x", units="µC")
# returns: "\$\\overline{x}, \\sigma_{ x }~(\\mathrm{ µC } )\$"
```
"""
function mathlabel(keys::String...; units=nothing, tex::Bool=true)
    # Cast to string
    if !isnothing(units)
        units = string(units)
    end

    if tex
        label_list = [texlabel(key) !== nothing ? texlabel(key) : "\\mathrm{ $key }" for key in keys]
        label = join(label_list, ", ")
        if !isnothing(units)
            units = replace(units, "*" => "{\\cdot}")
            label = "$label~(\\mathrm{ $units } )"
        end

        return "\$$label\$"
    else
        label = join(keys, ", ")

        if !isnothing(units)
            label = "$label ($units)"
        end

        return label
    end
end