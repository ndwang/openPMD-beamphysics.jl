# ----------------------
# Analysis

"""
    accelerating_voltage_and_phase(z, Ez, frequency)

Computes the accelerating voltage and phase for a v=c positively charged particle in an accelerating cavity field.

    Z = ∫ Ez * e^(-i k z) dz

    where k = omega/c = 2*pi*frequency/c

    voltage = abs(Z)
    phase = arg(Z)

# Arguments
- `z`: z-coordinate array (m)
- `Ez`: On-axis complex Ez field array (V/m), oscillating as exp(-i omega t), with omega = 2*pi*frequency
- `frequency`: RF frequency in Hz

# Returns
- `voltage`: voltage in V
- `phase`: phase in radians
"""
function accelerating_voltage_and_phase(z, Ez, frequency)
    omega = 2π * frequency
    k = omega / c_light
    fz = Ez .* exp.(-1im * k .* z)

    # Integrate using trapezoidal rule
    Z = trapz(z, fz)

    # Max voltage at phase
    voltage = abs(Z)
    phase = angle(Z)

    return voltage, phase
end

"""
    track_field_1d(z, Ez, frequency=0, z0=0, pz0=0, t0=0, mc2=mec2, q0=-1, debug=false, max_step=nothing)

Tracks a particle in a 1d complex electric field Ez, oscillating as Ez * exp(-i omega t)

Uses DifferentialEquations.jl to track the particle.

Equations of motion:

    dz/dt = pc/√((pc)² + m²c⁴) * c
    dp/dt = q E_z
    E_z = Re f(z) exp(-i ω t)

# Arguments
- `z`: positions of the field Ez (m)
- `Ez`: On-axis longitudinal electric field (V/m)
- `frequency`: RF frequency in Hz
- `z0`: initial particle position (m)
- `pz0`: initial particle momentum (eV/c)
- `t0`: initial particle time (s)
- `mc2`: initial particle mass (eV)
- `q0`: initial particle charge (e) (= -1 for electron)
- `max_step`: Maximum timestep for solve_ivp (s)
- `debug`: If true, returns the full solution

# Returns
- `z1`: final z position in (m)
- `pz1`: final particle momentum (eV/c)
- `t1`: final time (s)
"""
function track_field_1d(z, Ez, frequency=0, z0=0, pz0=0, t0=0, mc2=mec2, q0=-1, debug=false, max_step=nothing)
    # Make interpolating function
    field = linear_interpolation(z, Ez * q0 * c_light, extrapolation_bc=Line())
    zmax = maximum(z)
    tmax = 100 / frequency
    omega = 2π * frequency

    # function to integrate
    function fun(t, y, p)
        z = y[1]
        p = y[2]
        zdot = p / hypot(p, mc2) * c_light
        pdot = real(field(z) * exp(-1im * omega * t))
        return [zdot, pdot]
    end

    # Events (stopping conditions)
    function went_backwards(t, y, integrator)
        return y[1]
    end

    function went_max(t, y, integrator)
        return y[1] - zmax
    end

    if isnothing(max_step)
        max_step = 1 / (10 * frequency)
    end

    # Solve
    y0 = [z0, pz0]
    tspan = (t0, tmax)
    prob = ODEProblem(fun, y0, tspan)
    sol = solve(prob, 
                RK4(), 
                dt=1/frequency/1000,
                callback=CallbackSet(
                    ContinuousCallback(went_backwards, terminate!, affect!),
                    ContinuousCallback(went_max, terminate!, affect!)
                ),
                maxstep=max_step)

    if debug
        return sol
    end

    # Final z, p, t
    zf = sol.u[end][1]
    pf = sol.u[end][2]
    tf = sol.t[end]

    return zf, pf, tf
end

"""
    track_field_1df(Ez_f, zstop=0, tmax=0, z0=0, pz0=0, t0=0, mc2=mec2, q0=-1, debug=false, max_step=nothing, method="RK4")

Similar to track_field_1d, except uses a function Ez_f

# Arguments
- `Ez_f`: Ez_f(z, t) callable with two arguments z (m) and t (s)
- `zstop`: z stopping position (m)
- `tmax`: maximum timestep (s)
- `z0`: initial particle position (m)
- `pz0`: initial particle momentum (eV/c)
- `t0`: initial particle time (s)
- `mc2`: initial particle mass (eV)
- `q0`: initial particle charge (e) (= -1 for electron)
- `max_step`: Maximum timestep for solve_ivp (s)
- `debug`: If true, returns the full solution
- `method`: ODE solver method

# Returns
- `z1`: final z position in (m)
- `pz1`: final particle momentum (eV/c)
- `t1`: final time (s)
"""
function track_field_1df(Ez_f, zstop=0, tmax=0, z0=0, pz0=0, t0=0, mc2=mec2, q0=-1, debug=false, max_step=nothing, method="RK4")
    # function to integrate
    function fun(t, y, p)
        z = y[1]
        p = y[2]
        zdot = p / hypot(p, mc2) * c_light
        pdot = Ez_f(z, t) * q0 * c_light
        return [zdot, pdot]
    end

    # Events (stopping conditions)
    function went_backwards(t, y, integrator)
        return y[1]
    end

    function went_max(t, y, integrator)
        return y[1] - zstop
    end

    if isnothing(max_step)
        max_step = tmax / 10
    end

    # Solve
    y0 = [z0, pz0]
    tspan = (t0, tmax)
    prob = ODEProblem(fun, y0, tspan)
    sol = solve(prob, 
                RK4(), 
                dt=tmax/1000,
                callback=CallbackSet(
                    ContinuousCallback(went_backwards, terminate!, affect!),
                    ContinuousCallback(went_max, terminate!, affect!)
                ),
                maxstep=max_step)

    if debug
        return sol
    end

    # Final z, p, t
    zf = sol.u[end][1]
    pf = sol.u[end][2]
    tf = sol.t[end]

    return zf, pf, tf
end

"""
    autophase_field(field_mesh, pz0=0, scale=1, species="electron", tol=1e-9, verbose=false, debug=false)

Finds the maximum accelerating of a FieldMesh by tracking a particle and using Brent's method.

# Arguments
- `field_mesh`: FieldMesh object
- `pz0`: initial particle momentum in the z direction, in eV/c
- `scale`: Additional field scale
- `species`: species to track
- `tol`: Tolerance for brent: Stop if between iteration change is less than tol
- `debug`: If true, returns a function that tracks the field at a given phase in deg
- `verbose`: If true, prints information about the v=c voltage and phase

# Returns
- `phase`: Maximum accelerating phase in deg
- `pz1`: Final particle momentum in the z direction, in eV/c
"""
function autophase_field(field_mesh, pz0=0, scale=1, species="electron", tol=1e-9, verbose=false, debug=false)
    # Get field on-axis
    z = field_mesh.coord_vec("z")
    Ez = field_mesh.Ez[1, 1, :] * scale
    frequency = field_mesh.frequency
    zmin = minimum(z)

    # Get mass and charge state
    mc2 = mass_of(species)
    q0 = charge_state(species)  # -1 for electrons

    # Function for use in brent
    function phase_f(phase_deg)
        zf, pf, _ = track_field_1d(
            z,
            Ez,
            frequency=frequency,
            z0=zmin,
            pz0=pz0,
            t0=phase_deg / 360 / frequency,
            mc2=mc2,
            max_step=1 / frequency / 10,
            q0=q0,
        )
        return pf
    end

    if debug
        return phase_f
    end

    # Get a quick estimate, to use in the bracket
    voltage0, phase0 = accelerating_voltage_and_phase(z, q0 * Ez, frequency)
    phase0_deg = phase0 * 180 / π
    if verbose
        println("v=c voltage: $voltage0 V, phase: $phase0_deg deg")
    end

    alg_sign = -1
    phase_range = [phase0_deg - 90, phase0_deg + 90]
    
    # Use Optim.jl's Brent method
    result = optimize(x -> alg_sign * phase_f(x % 360), phase_range[1], phase_range[2], Brent())
    phase1_deg = result.minimizer
    pz1 = alg_sign * result.minimum

    if verbose
        println("    iterations: $(result.iterations)")
        println("    function calls: $(result.f_calls)")
    end

    return phase1_deg % 360, pz1
end

"""
    autophase_and_scale_field(field_mesh, voltage, pz0=0, species="electron", debug=false, verbose=false)

Finds the maximum accelerating of a FieldMesh.

Uses two iterations of phasing, scaling.

# Arguments
- `field_mesh`: FieldMesh object
- `voltage`: Desired on-crest voltage in V
- `pz0`: initial particle momentum in the z direction, in eV/c
- `species`: species to track
- `debug`: If true, returns a function that tracks the field at a given phase and scale
- `verbose`: If true, prints information about the v=c voltage and phase

# Returns
- `phase`: Maximum accelerating phase in deg
- `scale`: scale factor for the field
"""
function autophase_and_scale_field(field_mesh, voltage, pz0=0, species="electron", debug=false, verbose=false)
    z = field_mesh.coord_vec("z")
    Ez = field_mesh.Ez[1, 1, :]
    frequency = field_mesh.frequency
    zmin = minimum(z)

    # Get mass and charge
    mc2 = mass_of(species)
    q0 = charge_state(species)
    energy0 = hypot(pz0, mc2)

    # Get and initial estimate
    voltage0, phase0 = accelerating_voltage_and_phase(z, q0 * Ez, frequency)
    # convert to deg
    phase0 = phase0 * 180 / π
    scale0 = voltage / voltage0
    if verbose
        println("v=c voltage: $voltage0 V, phase: $phase0 deg")
    end

    function phase_scale_f(phase_deg, scale)
        zf, pf, _ = track_field_1d(
            z,
            Ez * scale,
            frequency=frequency,
            z0=zmin,
            pz0=pz0,
            t0=phase_deg / 360 / frequency,
            mc2=mc2,
            max_step=1 / frequency / 10,
            q0=q0,
        )

        delta_energy = hypot(pf, mc2) - energy0

        return delta_energy
    end

    if debug
        return phase_scale_f
    end

    # Phase 1
    brack = [phase0 - 90, phase0 + 90]
    result1 = optimize(x -> -phase_scale_f(x, scale0), brack[1], brack[2], Brent())
    phase1 = result1.minimizer % 360

    # Scale 1
    s0 = scale0 * 0.9
    s1 = scale0 * 1.1
    result2 = find_zero(x -> phase_scale_f(phase1, x) / voltage - 1.0, (s0, s1), Brent())
    scale1 = result2

    if verbose
        println("    Pass 1 delta energy: $(phase_scale_f(phase1, scale1)) at phase  $phase1 deg")
    end

    # Phase 2
    brack = [phase1 - 10, phase1 + 10]
    result3 = optimize(x -> -phase_scale_f(x, scale1), brack[1], brack[2], Brent())
    phase2 = result3.minimizer % 360

    # Scale 2
    s0 = scale1 * 0.9
    s1 = scale1 * 1.1
    result4 = find_zero(x -> phase_scale_f(phase2, x) / voltage - 1.0, (s0, s1), Brent())
    scale2 = result4

    if verbose
        println("    Pass 2 delta energy: $(phase_scale_f(phase2, scale2)) at phase  $phase2 deg")
    end

    return phase2, scale2
end

# Helper function for trapezoidal integration
function trapz(x, y)
    dx = diff(x)
    return sum(dx .* (y[1:end-1] .+ y[2:end]) ./ 2)
end