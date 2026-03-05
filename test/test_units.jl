using AtomicAndPhysicalConstants: C_LIGHT, E_CHARGE
import OpenPMDBeamphysics:
    KNOWN_UNIT, DIMENSION, DIMENSION_NAME, SI_SYMBOL, SI_NAME,
    PARTICLEGROUP_UNITS, NAMED_UNITS,
    power_unit, sqrt_unit, is_dimensionless, is_identity,
    dimension, dimension_name,
    write_unit_h5, read_unit_h5,
    write_dataset_and_unit_h5, read_dataset_and_unit_h5

using HDF5

# ── PMDUnit struct ──────────────────────────────────────────────────────────

@testset "PMDUnit constructors" begin
    # 3-arg with tuple
    u = PMDUnit("m", 1.0, (1., 0., 0., 0., 0., 0., 0.))
    @test u.symbol == "m"
    @test u.scale  == 1.0
    @test u.dims   == (1., 0., 0., 0., 0., 0., 0.)

    # 3-arg with dimension name
    u2 = PMDUnit("m", 1.0, "length")
    @test u2.dims == DIMENSION["length"]

    # Integer scale is promoted to Float64
    u3 = PMDUnit("test", 2, (0., 0., 0., 0., 0., 0., 0.))
    @test u3.scale === 2.0

    # Integer dims are promoted to Float64
    u4 = PMDUnit("test", 1.0, (1, 0, 0, 0, 0, 0, 0))
    @test u4.dims === (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

@testset "Negative-zero normalization" begin
    # -0.0 in dims should be normalized to +0.0 for consistent hashing
    u1 = PMDUnit("a", 1.0, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    u2 = PMDUnit("a", 1.0, (-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0))
    @test u1 == u2
    @test hash(u1) == hash(u2)
    @test isequal(u1.dims, u2.dims)
end

# ── Dimension dictionaries ─────────────────────────────────────────────────

@testset "Dimension dictionaries" begin
    @test DIMENSION["length"]   == (1., 0., 0., 0., 0., 0., 0.)
    @test DIMENSION["momentum"] == (1., 1., -1., 0., 0., 0., 0.)
    @test DIMENSION["1"]        == (0., 0., 0., 0., 0., 0., 0.)

    # DIMENSION_NAME is the inverse
    for (name, dims) in DIMENSION
        @test DIMENSION_NAME[dims] == name
    end

    @test dimension("energy") == (2., 1., -2., 0., 0., 0., 0.)
    @test_throws ErrorException dimension("nonexistent")

    @test dimension_name((1., 0., 0., 0., 0., 0., 0.)) == "length"
    @test_throws ErrorException dimension_name((9., 9., 9., 9., 9., 9., 9.))
end

@testset "SI symbol dictionaries" begin
    @test SI_SYMBOL["energy"]   == "J"
    @test SI_SYMBOL["momentum"] == "kg⋅m/s"
    @test SI_NAME["m"]          == "length"
    @test SI_NAME["T"]          == "magnetic_field"
end

# ── Predicates & display ───────────────────────────────────────────────────

@testset "Predicates" begin
    @test is_dimensionless(KNOWN_UNIT[""])
    @test is_dimensionless(KNOWN_UNIT["rad"])
    @test is_identity(KNOWN_UNIT[""])
    @test is_identity(KNOWN_UNIT["rad"])
    @test !is_dimensionless(KNOWN_UNIT["m"])
    @test !is_identity(KNOWN_UNIT["eV"])  # scale != 1
end

@testset "show / string" begin
    @test string(KNOWN_UNIT["eV/c"]) == "eV/c"
    @test string(KNOWN_UNIT[""])     == ""
    @test string(KNOWN_UNIT["m"])    == "m"
end

@testset "Equality and hashing" begin
    a = PMDUnit("m", 1.0, "length")
    b = PMDUnit("m", 1.0, (1., 0., 0., 0., 0., 0., 0.))
    @test a == b
    @test hash(a) == hash(b)

    c = PMDUnit("km", 1000.0, "length")
    @test a != c
end

@testset "dimension_name for PMDUnit" begin
    @test dimension_name(KNOWN_UNIT["eV"])   == "energy"
    @test dimension_name(KNOWN_UNIT["eV/c"]) == "momentum"
    @test dimension_name(KNOWN_UNIT["T"])    == "magnetic_field"
    @test dimension_name(KNOWN_UNIT[""])     == "1"
end

# ── Arithmetic ──────────────────────────────────────────────────────────────

@testset "Multiplication" begin
    m   = KNOWN_UNIT["m"]
    evc = KNOWN_UNIT["eV/c"]
    prod = m * evc
    @test prod.symbol == "m*eV/c"
    @test prod.scale  ≈ E_CHARGE / C_LIGHT
    @test prod.dims   == (2., 1., -1., 0., 0., 0., 0.)

    # Same unit squared
    sq = m * m
    @test sq.symbol == "(m)^2"
    @test sq.scale  == 1.0
    @test sq.dims   == (2., 0., 0., 0., 0., 0., 0.)

    # Identity shortcuts
    id = KNOWN_UNIT[""]
    @test (id * m) === m
    @test (m * id) === m
end

@testset "Division" begin
    m = KNOWN_UNIT["m"]
    s = KNOWN_UNIT["s"]
    vel = m / s
    @test vel.symbol == "m/s"
    @test vel.scale  == 1.0
    @test vel.dims   == DIMENSION["velocity"]

    # Same unit → dimensionless identity
    one = m / m
    @test is_identity(one)
    @test one.symbol == ""

    # Division by identity
    @test (m / KNOWN_UNIT[""]) === m
end

@testset "power_unit" begin
    m = KNOWN_UNIT["m"]

    # p=0 → identity
    u0 = power_unit(m, 0)
    @test is_identity(u0)

    # p=1 → same unit
    @test power_unit(m, 1) === m

    # Integer power
    u2 = power_unit(m, 2)
    @test u2.symbol == "m^2"
    @test u2.scale  == 1.0
    @test u2.dims   == (2., 0., 0., 0., 0., 0., 0.)

    # Negative integer power
    um1 = power_unit(m, -1)
    @test um1.symbol == "m^-1"
    @test um1.dims   == (-1., 0., 0., 0., 0., 0., 0.)

    # Fractional power
    uhalf = power_unit(m, 0.5)
    @test uhalf.symbol == "m^(1/2)"
    @test uhalf.dims   == (0.5, 0., 0., 0., 0., 0., 0.)
    @test uhalf.scale  == 1.0

    # Compound symbol gets wrapped in parens
    evc = KNOWN_UNIT["eV/c"]
    evc2 = power_unit(evc, 2)
    @test evc2.symbol == "(eV/c)^2"
    @test evc2.scale  ≈ (E_CHARGE / C_LIGHT)^2

    # Empty symbol stays empty
    id = KNOWN_UNIT[""]
    @test power_unit(id, 3).symbol == ""
end

@testset "sqrt_unit" begin
    m = KNOWN_UNIT["m"]
    u = sqrt_unit(m)
    @test u.symbol == "m^(1/2)"
    @test u.dims == (0.5, 0., 0., 0., 0., 0., 0.)
end

# ── Named-unit registry ────────────────────────────────────────────────────

@testset "KNOWN_UNIT registry" begin
    # All named units are present
    for u in NAMED_UNITS
        @test haskey(KNOWN_UNIT, u.symbol)
        @test KNOWN_UNIT[u.symbol] === u
    end

    # Aliases
    @test KNOWN_UNIT["1"]       === KNOWN_UNIT[""]
    @test KNOWN_UNIT["c_light"] === KNOWN_UNIT["c"]

    # Spot-check scales
    @test KNOWN_UNIT["eV"].scale   == E_CHARGE
    @test KNOWN_UNIT["eV/c"].scale ≈ E_CHARGE / C_LIGHT
    @test KNOWN_UNIT["c"].scale    == C_LIGHT
    @test KNOWN_UNIT["g"].scale    == 1e-3
    @test KNOWN_UNIT["degree"].scale ≈ π / 180

    # Spot-check dims
    @test KNOWN_UNIT["eV"].dims   == DIMENSION["energy"]
    @test KNOWN_UNIT["eV/c"].dims == DIMENSION["momentum"]
    @test KNOWN_UNIT["T"].dims    == DIMENSION["magnetic_field"]
    @test KNOWN_UNIT["V/m"].dims  == DIMENSION["electric_field"]
    @test KNOWN_UNIT["C"].dims    == DIMENSION["charge"]
end

# ── Symbol parser ──────────────────────────────────────────────────────────

@testset "Single-arg PMDUnit (known symbols)" begin
    # Known symbols go through fast path
    @test PMDUnit("eV/c") === KNOWN_UNIT["eV/c"]
    @test PMDUnit("m")    === KNOWN_UNIT["m"]
    @test PMDUnit("")     === KNOWN_UNIT[""]
    @test PMDUnit("1")    === KNOWN_UNIT["1"]
end

@testset "Symbol parsing — compound units" begin
    # kg*m/s = SI momentum
    u = PMDUnit("kg*m/s")
    @test u.symbol == "kg*m/s"
    @test u.scale  == 1.0
    @test u.dims   == DIMENSION["momentum"]

    # W/rad^2 — rad is dimensionless, so same dims as W
    u2 = PMDUnit("W/rad^2")
    @test u2.symbol == "W/rad^2"
    @test u2.scale  == 1.0
    @test u2.dims   == KNOWN_UNIT["W"].dims

    # W/m^2
    u3 = PMDUnit("W/m^2")
    @test u3.symbol == "W/m^2"
    @test u3.scale  == 1.0
    @test u3.dims   == (0., 1., -3., 0., 0., 0., 0.)

    # m*eV/c
    u4 = PMDUnit("m*eV/c")
    @test u4.symbol == "m*eV/c"
    @test u4.scale  ≈ E_CHARGE / C_LIGHT
    @test u4.dims   == (2., 1., -1., 0., 0., 0., 0.)
end

@testset "Symbol parsing — power notation" begin
    u = PMDUnit("m^2")
    @test u.dims == (2., 0., 0., 0., 0., 0., 0.)

    u2 = PMDUnit("m^-1")
    @test u2.dims == (-1., 0., 0., 0., 0., 0., 0.)

    u3 = PMDUnit("m^(1/2)")
    @test u3.dims == (0.5, 0., 0., 0., 0., 0., 0.)

    u4 = PMDUnit("s^(-3/2)")
    @test u4.dims == (0., 0., -1.5, 0., 0., 0., 0.)
    @test u4.scale == 1.0
end

@testset "Symbol parsing — sqrt notation" begin
    u = PMDUnit("sqrt(m)")
    @test u.symbol == "sqrt(m)"
    @test u.scale  == 1.0
    @test u.dims   == (0.5, 0., 0., 0., 0., 0., 0.)
end

@testset "Symbol parsing — errors" begin
    @test_throws ErrorException PMDUnit("nonexistent_unit")
    @test_throws ErrorException PMDUnit("foo*bar")
end

# ── SI prefix / nice_array ─────────────────────────────────────────────────

@testset "nice_scale_prefix" begin
    f, p = nice_scale_prefix(1.5e-6)
    @test f == 1e-6
    @test p == "µ"

    f, p = nice_scale_prefix(0.0)
    @test f == 1.0
    @test p == ""

    f, p = nice_scale_prefix(42.0)
    @test f == 1.0  # within [-1.5, 2] log10 range
    @test p == ""

    f, p = nice_scale_prefix(5000.0)
    @test f == 1e3
    @test p == "k"
end

@testset "nice_array" begin
    a = [1e-3, 2e-3, 3e-3]
    scaled, fac, prefix = nice_array(a)
    @test fac == 1e-3
    @test prefix == "m"
    @test scaled ≈ [1.0, 2.0, 3.0]

    # Scalar version
    val, fac, prefix = nice_array(5e6)
    @test fac == 1e6
    @test prefix == "M"
    @test val ≈ 5.0
end

@testset "limits" begin
    @test limits([3.0, 1.0, 2.0]) == (1.0, 3.0)
end

@testset "plottable_array" begin
    x = [1e-9, 2e-9, 3e-9]
    scaled, factor, prefix, xmin, xmax = plottable_array(x)
    @test xmin == 1e-9
    @test xmax == 3e-9
    @test prefix == "n"
    @test scaled ≈ x ./ factor
end

# ── HDF5 round-trip ────────────────────────────────────────────────────────

@testset "write_unit_h5 / read_unit_h5 round-trip" begin
    tmpfile = tempname() * ".h5"
    try
        # Write several units
        h5open(tmpfile, "w") do h5
            for sym in ["eV/c", "m", "T", ""]
                u = KNOWN_UNIT[sym]
                g = create_group(h5, sym == "" ? "dimensionless" : sym)
                write_unit_h5(g, u)
            end
        end

        # Read back
        h5open(tmpfile, "r") do h5
            for (gname, sym) in [("eV/c", "eV/c"), ("m", "m"), ("T", "T"), ("dimensionless", "")]
                g = h5[gname]
                u = read_unit_h5(g)
                @test u.symbol == sym
                @test u.scale  == KNOWN_UNIT[sym].scale
                @test u.dims   == KNOWN_UNIT[sym].dims
            end
        end
    finally
        rm(tmpfile; force=true)
    end
end

@testset "write/read dataset_and_unit_h5" begin
    tmpfile = tempname() * ".h5"
    try
        data = [1.0, 2.0, 3.0]
        u = KNOWN_UNIT["eV"]

        h5open(tmpfile, "w") do h5
            write_dataset_and_unit_h5(h5, "energy", data; unit=u)
        end

        h5open(tmpfile, "r") do h5
            d, ru = read_dataset_and_unit_h5(h5["energy"])
            @test d ≈ data
            @test ru.symbol == "eV"
            @test ru.scale  == E_CHARGE
        end

        # Read with expected unit (same unit, no conversion needed)
        h5open(tmpfile, "r") do h5
            d, ru = read_dataset_and_unit_h5(h5["energy"], u)
            @test d ≈ data
            @test ru === u
        end

        # Read with expected unit as string
        h5open(tmpfile, "r") do h5
            d, ru = read_dataset_and_unit_h5(h5["energy"], "eV")
            @test d ≈ data
        end

        # Read with conversion to J
        h5open(tmpfile, "w") do h5
            write_dataset_and_unit_h5(h5, "energy_eV", data; unit=u)
        end
        h5open(tmpfile, "r") do h5
            d, ru = read_dataset_and_unit_h5(h5["energy_eV"], KNOWN_UNIT["J"])
            @test ru.symbol == "J"
            # 1 eV = E_CHARGE J, so data in J = data * E_CHARGE
            @test d ≈ data .* E_CHARGE
        end

        # Dimension mismatch error
        h5open(tmpfile, "r") do h5
            @test_throws ErrorException read_dataset_and_unit_h5(h5["energy_eV"], KNOWN_UNIT["m"])
        end
    finally
        rm(tmpfile; force=true)
    end
end

@testset "read_unit_h5 fallback for unknown symbol" begin
    tmpfile = tempname() * ".h5"
    try
        h5open(tmpfile, "w") do h5
            g = create_group(h5, "custom")
            attrs(g)["unitSymbol"]    = "furlong"
            attrs(g)["unitSI"]        = 201.168
            attrs(g)["unitDimension"] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        end

        h5open(tmpfile, "r") do h5
            u = read_unit_h5(h5["custom"])
            @test u.symbol == "furlong"
            @test u.scale  == 201.168
            @test u.dims   == (1., 0., 0., 0., 0., 0., 0.)
        end
    finally
        rm(tmpfile; force=true)
    end
end

# ── ParticleGroup units ────────────────────────────────────────────────────

@testset "pg_units — direct lookup" begin
    @test pg_units("x").symbol  == "m"
    @test pg_units("px").symbol == "eV/c"
    @test pg_units("t").symbol  == "s"
    @test pg_units("energy").symbol == "eV"
    @test pg_units("gamma").symbol  == ""
    @test pg_units("theta").symbol  == "rad"
    @test pg_units("charge").symbol == "C"
    @test pg_units("Bx").symbol     == "T"
    @test pg_units("Ex").symbol     == "V/m"
end

@testset "pg_units — prefix stripping" begin
    @test pg_units("sigma_x").dims  == pg_units("x").dims
    @test pg_units("mean_px").dims  == pg_units("px").dims
    @test pg_units("min_energy").dims == pg_units("energy").dims
    @test pg_units("max_t").dims    == pg_units("t").dims
    @test pg_units("ptp_z").dims    == pg_units("z").dims
    @test pg_units("delta_pz").dims == pg_units("pz").dims
end

@testset "pg_units — covariance" begin
    u = pg_units("cov_x__px")
    @test u.symbol == "m*eV/c"
    @test u.dims   == KNOWN_UNIT["m"].dims .+ KNOWN_UNIT["eV/c"].dims

    u2 = pg_units("cov_x__x")
    @test u2.symbol == "(m)^2"
    @test u2.dims   == (2., 0., 0., 0., 0., 0., 0.)
end

@testset "pg_units — special keys" begin
    # Compound units
    @test pg_units("norm_emit_4d").dims == (2., 0., 0., 0., 0., 0., 0.)
    @test pg_units("Lz").dims == (2., 1., -1., 0., 0., 0., 0.)
    @test pg_units("x_bar").dims == (0.5, 0., 0., 0., 0., 0., 0.)

    # Twiss
    @test pg_units("twiss_beta_x").symbol  == "m"
    @test pg_units("twiss_gamma_x").symbol == "m^-1"
    @test is_identity(pg_units("twiss_alpha_x"))

    # Fields
    @test pg_units("electricField/x").symbol == "V/m"
    @test pg_units("magneticField/z").symbol == "T"

    # Bunching
    @test is_identity(pg_units("bunching"))
    @test pg_units("bunching_phase").symbol == "rad"
end

@testset "pg_units — error on unknown key" begin
    @test_throws ErrorException pg_units("totally_unknown_key")
end

# ── Full ParticleGroup HDF5 round-trip ─────────────────────────────────────

@testset "ParticleGroup HDF5 round-trip" begin
    pg = single_particle(x=1e-3, px=1e6, y=2e-3, py=-5e5,
                         z=0.01, pz=1e9, t=1e-12, weight=1e-15)
    tmpfile = tempname() * ".h5"
    try
        write_particle_group(tmpfile, pg)
        pg2 = ParticleGroup(tmpfile)
        @test length(pg2) == 1
        @test pg2.x  ≈ pg.x   rtol=1e-12
        @test pg2.px ≈ pg.px  rtol=1e-12
        @test pg2.y  ≈ pg.y   rtol=1e-12
        @test pg2.py ≈ pg.py  rtol=1e-12
        @test pg2.z  ≈ pg.z   rtol=1e-12
        @test pg2.pz ≈ pg.pz  rtol=1e-12
        @test pg2.t  ≈ pg.t   rtol=1e-12
        @test pg2.weight ≈ pg.weight rtol=1e-12
    finally
        rm(tmpfile; force=true)
    end
end
