"""
    expand_1d_static_fieldmap(z0, fz0, spline_s=0)

Expands 1D static fieldmap z, fz into r, z using splines.

Cylindrically symmetric geometry.

This is valid for both electric and magnetic fields.

Returns functions:
    fr(r, z), fz(r, z)
"""
function expand_1d_static_fieldmap(z0, fz0, spline_s=0)
    # Make spline and derivatives
    S = Spline1D(z0, fz0, k=5, s=spline_s)
    Sp = derivative(S)
    Sp2 = derivative(Sp)
    Sp3 = derivative(Sp2)

    function fz(r, z)
        return S(z) - r^2 * Sp2(z) / 4
    end

    function fr(r, z)
        return -r * Sp(z) / 2 + r^3 * Sp3(z) / 16
    end

    return fr, fz
end

"""
    expand_1d_dynamic_fieldmap(z, Ez0, frequency=0, spline_s=0)

Expands 1D dynamic electric fieldmap z, Ez into r, z using splines.

Cylindrically symmetric geometry.

Fields oscillate as:
    Er, Ez ~ sin(wt)  <=> ~  cos(wt)
    Btheta ~ cos(wt)      ~ -sin(wt)
This is the Superfish convention.

Returns real functions:
    Er(r, z), Ez(r, z), Btheta(r, z)
"""
function expand_1d_dynamic_fieldmap(z, Ez0, frequency=0, spline_s=0)
    omega = 2π * frequency

    # Make spline and derivatives
    S = Spline1D(z, Ez0, k=5, s=spline_s)
    Sp = derivative(S)
    Sp2 = derivative(Sp)
    Sp3 = derivative(Sp2)

    function Ez(r, z)
        return S(z) - r^2 / 4 * (Sp2(z) + (omega / c_light)^2 * S(z))  # sin(omega*t)
    end

    function Er(r, z)
        return -r / 2 * Sp(z) + r^3 / 16 * (Sp3(z) + (omega / c_light)^2 * Sp(z))  # sin(omega*t)
    end

    function Btheta(r, z)
        return (r / 2 * S(z) - r^3 / 16 * (Sp2(z) + (omega / c_light)^2 * S(z))) * omega / c_light^2.0  # cos(omega*t)
    end

    return Er, Ez, Btheta
end

"""
    spline_derivative_array(z, fz, s=0, k=5)

Make spline and derivatives
"""
function spline_derivative_array(z, fz, s=0, k=5)
    # Make spline and derivatives
    S = Spline1D(z, fz, k=k, s=s)
    fz1 = derivative(S)(z)
    fz2 = derivative(derivative(S))(z)
    fz3 = derivative(derivative(derivative(S)))(z)

    return hcat(fz, fz1, fz2, fz3)
end

"""
    fft_derivative_array(fz, dz, ncoef=30, max_order=3)

Create derivatives of field data `fz` with regular spacing `dz` using an FFT method.

Parameters:
- fz: array - Field values with regular spacing (assumed to be periodic)
- dz: float - Array spacing
- ncoef: int - Number of Fourier coefficients to keep (default: 30)
- max_order: int - Maximum order of derivatives to compute

Returns:
- array of shape (len(fz), max_order+1) representing the field values and derivatives
"""
function fft_derivative_array(fz, dz, ncoef=30, max_order=3)
    fz = real.(fz)

    n = length(fz)
    L = n * dz

    # Pad odd length
    if isodd(n)
        odd = true
        fz = vcat(fz, 0)
        L += dz
    else
        odd = false
    end

    # Initial FFT
    y = rfft(fz)
    # Cutoff
    y[ncoef+1:end] .= 0
    k = 0:length(y)-1

    derivs = []
    y_current = copy(y)

    for order in 0:max_order
        a = irfft(y_current, n)
        if odd
            # Trim off to return the same length as fz
            a = a[1:end-1]
        end
        push!(derivs, a)
        y_current .*= 2π * im * k / L
    end

    return hcat(derivs...)
end

"""
    expand_radial(r, dfield, frequency=0)

Expand a field at r from its on-axis derivative array.

See, for example, the Astra manual:
https://www.desy.de/~mpyflo/Astra_manual/Astra-Manual_V3.2.pdf
Appendix I: Field expansion formulas

Parameters:
- r: float or array
- dfield: derivative array
- frequency: float - frequency in Hz (default: 0)

Returns:
- fr: array - r field component
- fz: array - z field component
- ftheta: array - theta field component
"""
function expand_radial(r, dfield, frequency=0)
    f0 = dfield[:, 1]  # f
    f1 = dfield[:, 2]  # f'
    f2 = dfield[:, 3]  # f''
    f3 = dfield[:, 4]  # f'''

    omega = 2π * frequency
    ooc2 = (omega / c_light)^2

    if frequency == 0
        fz = f0 - r^2 / 4 * f2
        fr = -r / 2 * f1 + r^3 / 16 * f3
        ftheta = zeros_like(f0)
    else
        fz = f0 - r^2 / 4 * (f2 + ooc2 * f0)  # cos(wt)
        fr = -r / 2 * f1 + r^3 / 16 * (f3 + ooc2 * f0)  # cos(wt)
        ftheta = r / 2 * f0 - r^3 / 16 * (f2 + ooc2 * f0)  # sin(wt) * (-w/c^2 for electric, w for magnetic)
    end

    return fr, fz, ftheta
end

"""
    expand_fieldmesh_from_onaxis(fieldmesh; dr=nothing, nr=10, inplace=false, method="spline", ncoef=nothing, spline_s=0, zmirror="auto")

Create cylindrical FieldMesh data from 1-d on-axis field data.

This uses an FFT method to compute up to third order derivatives,
and then uses field expansion formulas to populate the component data.

Parameters:
- fieldmesh: FieldMesh
- dr: float - radial coordinate spacing (default: same as z spacing)
- nr: int - number of radial coordinates (default: 10)
- frequency: float - frequency in Hz
- method: str - Expansion method to use, one of 'fft' or 'spline' (default: 'spline')
- ncoef: int - Number of Fourier coefficients to use for the expansion (default: nz/4)
- spline_s: float - Spline smoothing factor (default: 0)
- zmirror: str or bool - Mirror the field about the minimum z before the expansion (default: 'auto')
- inplace: bool - If true, modify in-place (default: false)

Returns:
- fieldmesh: FieldMesh
"""
function expand_fieldmesh_from_onaxis(fieldmesh; dr=nothing, nr=10, inplace=false, method="spline", ncoef=nothing, spline_s=0, zmirror="auto")
    if !inplace
        fieldmesh = copy(fieldmesh)
    end

    @assert fieldmesh.geometry == "cylindrical"

    has_Ez = "electricField/z" in fieldmesh.components
    has_Bz = "magneticField/z" in fieldmesh.components

    if has_Bz && has_Ez
        error("Expanding both Ez and Bz not implemented")
    end

    if has_Ez
        fz = fieldmesh["electricField/z"][1, 1, :]
    elseif has_Bz
        fz = fieldmesh["magneticField/z"][1, 1, :]
    else
        error("Neither Ez nor Bz found")
    end

    # Get attrs
    nz = fieldmesh.shape[3]
    dz = fieldmesh.dz
    frequency = fieldmesh.frequency

    # Get real field
    @assert all(isreal.(fz))
    fz = real.(fz)

    zvec = fieldmesh.coord_vec("z")

    # Use the same as z
    if isnothing(dr)
        dr = dz
    end

    # Methods
    if method == "fft"
        # Crude heuristic
        if isnothing(ncoef)
            ncoef = nz ÷ 4
        end

        # Heuristic to mirror field, 1% of field at 0.
        if zmirror == "auto"
            if abs(fz[1] / maximum(abs.(fz))) > 0.01
                zmirror = true
            end
        end
        if zmirror
            fz = vcat(reverse(fz), fz[2:end])
        end

        dfield = fft_derivative_array(fz, dz, ncoef=ncoef)
        # Now strip off the part we want
        if zmirror
            dfield = dfield[end-nz+1:end, :]
        end

        # Collect field derivatives for each r, form large array.
        field = []
        for ir in 1:nr
            r = (ir-1) * dr
            push!(field, expand_radial(r, dfield, frequency=frequency))
        end
        field = cat(field..., dims=3)

        # Extract
        field_r = field[:, :, 1]  #  cos(wt)
        field_z = field[:, :, 2]  #  cos(wt)
        field_theta = field[:, :, 3] * -im * 2π * frequency / c_light^2  # -w/c^2 sin(wt)

    elseif method == "spline"
        rvec = range(0, dr * (nr-1), length=nr)
        RR, ZZ = meshgrid(rvec, zvec)

        if frequency == 0
            Frf, Fzf = expand_1d_static_fieldmap(zvec, fz, spline_s=spline_s)
        else
            Frf, Fzf, Fthetaf = expand_1d_dynamic_fieldmap(zvec, fz, frequency=frequency, spline_s=spline_s)
            field_theta = Fthetaf.(RR, ZZ) * -im
        end

        field_r = Frf.(RR, ZZ)  #  cos(wt)
        field_z = Fzf.(RR, ZZ)  #  cos(wt)

    else
        error("Invalid method: $method, must be one of ['fft', 'spline']")
    end

    # Collect components
    attrs = fieldmesh.attrs
    components = fieldmesh.components
    if has_Ez
        # TM mode
        components["electricField/r"] = reshape(field_r, :, 1, :)  #  cos(wt)
        components["electricField/z"] = reshape(field_z, :, 1, :)  #  cos(wt)

        if frequency != 0
            components["magneticField/theta"] = reshape(field_theta, :, 1, :)  # -w/c^2 sin(wt)
        end

    elseif has_Bz
        # TE mode
        components["magneticField/r"] = reshape(field_r, :, 1, :)
        components["magneticField/z"] = reshape(field_z, :, 1, :)
        if frequency != 0
            components["electricField/theta"] = reshape(field_theta, :, 1, :) * im * 2π * frequency  # w sin(wt)
        end
    end

    # Update attrs
    attrs["gridSize"] = (nr, 1, nz)
    attrs["gridSpacing"] = (dr, 0, dz)

    return fieldmesh 
end