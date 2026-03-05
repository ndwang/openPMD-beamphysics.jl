using HDF5
using AtomicAndPhysicalConstants: Species, C_LIGHT, E_CHARGE
import OpenPMDBeamphysics:
    pmd_init, pmd_field_init,
    write_component_data, write_pmd_bunch, write_pmd_field,
    particle_data_dict, write_particle_group,
    particle_paths, field_paths,
    particle_array, component_data,
    is_constant_component, constant_component_value,
    all_components, component_str, offset_component_name,
    load_field_attrs,
    COMPONENT_FROM_ALIAS, COMPONENT_ALIAS,
    write_unit_h5, read_unit_h5

const _PARTICLE_FILE = joinpath(@__DIR__, "initial_particles.h5")
const _FIELD_FILE    = joinpath(@__DIR__, "lcls2_solenoid_fieldmesh.h5")

@testset "I/O (openPMD)" begin

    # ── pmd_init / pmd_field_init ──────────────────────────────────────────

    @testset "pmd_init / pmd_field_init" begin
        @testset "pmd_init sets required attributes" begin
            tmpfile = tempname() * ".h5"
            try
                h5open(tmpfile, "w") do h5
                    pmd_init(h5)
                end
                h5open(tmpfile, "r") do h5
                    a = attrs(h5)
                    @test a["basePath"]         == "/data/%T/"
                    @test a["dataType"]         == "openPMD"
                    @test a["openPMD"]          == "2.0.0"
                    @test a["openPMDextension"] == "BeamPhysics;SpeciesType"
                    @test a["particlesPath"]    == "./"
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "pmd_init custom basePath" begin
            tmpfile = tempname() * ".h5"
            try
                h5open(tmpfile, "w") do h5
                    pmd_init(h5; basePath="/custom/%T/", particlesPath="particles/")
                end
                h5open(tmpfile, "r") do h5
                    @test attrs(h5)["basePath"]      == "/custom/%T/"
                    @test attrs(h5)["particlesPath"]  == "particles/"
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "pmd_field_init sets required attributes" begin
            tmpfile = tempname() * ".h5"
            try
                h5open(tmpfile, "w") do h5
                    pmd_field_init(h5)
                end
                h5open(tmpfile, "r") do h5
                    a = attrs(h5)
                    @test a["dataType"]           == "openPMD"
                    @test a["openPMD"]            == "2.0.0"
                    @test a["openPMDextension"]   == "BeamPhysics"
                    @test a["externalFieldPath"]  == "/ExternalFieldPath/%T/"
                end
            finally
                rm(tmpfile; force=true)
            end
        end
    end

    # ── write_particle_group + ParticleGroup round-trip ────────────────────

    @testset "ParticleGroup write/read round-trip" begin
        @testset "single particle" begin
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

        @testset "real beam from initial_particles.h5" begin
            pg = ParticleGroup(_PARTICLE_FILE)
            tmpfile = tempname() * ".h5"
            try
                write_particle_group(tmpfile, pg)
                pg2 = ParticleGroup(tmpfile)
                @test length(pg2) == length(pg)
                @test isapprox(pg, pg2; rtol=1e-12)
            finally
                rm(tmpfile; force=true)
            end
        end
    end

    # ── particle_paths / field_paths ───────────────────────────────────────

    @testset "particle_paths" begin
        @testset "written file with %%T pattern" begin
            tmpfile = tempname() * ".h5"
            try
                pg = single_particle(pz=1e9)
                write_particle_group(tmpfile, pg)
                h5open(tmpfile, "r") do h5
                    pp = particle_paths(h5)
                    @test length(pp) == 1
                    @test occursin("data", pp[1])
                    # The resolved path should contain a group with species
                    g = h5[pp[1]]
                    @test haskey(g, "electron")
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "initial_particles.h5" begin
            h5open(_PARTICLE_FILE, "r") do h5
                pp = particle_paths(h5)
                @test length(pp) == 1
                g = h5[pp[1]]
                @test haskey(g, "position") || length(keys(g)) >= 1
            end
        end
    end

    @testset "field_paths" begin
        @testset "file without externalFieldPath returns empty" begin
            tmpfile = tempname() * ".h5"
            try
                h5open(tmpfile, "w") do h5
                    pmd_init(h5)
                end
                h5open(tmpfile, "r") do h5
                    @test field_paths(h5) == String[]
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        if isfile(_FIELD_FILE)
            @testset "solenoid field mesh" begin
                h5open(_FIELD_FILE, "r") do h5
                    fp = field_paths(h5)
                    @test length(fp) >= 1
                    @test occursin("ExternalFieldMesh", fp[1])
                end
            end
        end
    end

    # ── particle_array ─────────────────────────────────────────────────────

    @testset "particle_array" begin
        @testset "reads position, momentum, time, weight, status" begin
            h5open(_PARTICLE_FILE, "r") do h5
                pp = particle_paths(h5)
                pg_h5 = h5[pp[1]]
                # Navigate to species subgroup if needed
                if !haskey(pg_h5, "position")
                    pg_h5 = pg_h5[first(keys(pg_h5))]
                end

                x  = particle_array(pg_h5, "x")
                px = particle_array(pg_h5, "px")
                t  = particle_array(pg_h5, "t")
                @test length(x) > 0
                @test length(px) == length(x)
                @test length(t)  == length(x)
            end
        end

        @testset "aliases resolve correctly" begin
            # Check that alias mapping exists
            @test COMPONENT_FROM_ALIAS["x"]  == "position/x"
            @test COMPONENT_FROM_ALIAS["px"] == "momentum/x"
            @test COMPONENT_FROM_ALIAS["t"]  == "time"
            @test COMPONENT_FROM_ALIAS["status"] == "particleStatus"
            @test COMPONENT_FROM_ALIAS["weight"] == "weight"
        end

        @testset "momentum converted to eV/c" begin
            # Write known momentum, read back, verify unit conversion
            pg = single_particle(px=1e6, pz=1e9)
            tmpfile = tempname() * ".h5"
            try
                write_particle_group(tmpfile, pg)
                h5open(tmpfile, "r") do h5
                    pp = particle_paths(h5)
                    sg = h5[pp[1]]
                    if !haskey(sg, "position")
                        sg = sg[first(keys(sg))]
                    end
                    px_read = particle_array(sg, "px")
                    @test px_read[1] ≈ 1e6 rtol=1e-10
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "include_offset=false" begin
            h5open(_PARTICLE_FILE, "r") do h5
                pp = particle_paths(h5)
                pg_h5 = h5[pp[1]]
                if !haskey(pg_h5, "position")
                    pg_h5 = pg_h5[first(keys(pg_h5))]
                end
                x_with    = particle_array(pg_h5, "x"; include_offset=true)
                x_without = particle_array(pg_h5, "x"; include_offset=false)
                # If no offset group exists, they should be equal
                if !haskey(pg_h5, "positionOffset")
                    @test x_with == x_without
                else
                    @test length(x_with) == length(x_without)
                end
            end
        end
    end

    # ── component_data / is_constant / constant_component_value ────────────

    @testset "component_data and constant components" begin
        @testset "array vs constant components in initial_particles.h5" begin
            h5open(_PARTICLE_FILE, "r") do h5
                pp = particle_paths(h5)
                pg_h5 = h5[pp[1]]
                if !haskey(pg_h5, "position")
                    pg_h5 = pg_h5[first(keys(pg_h5))]
                end

                # position/x should be a dataset (not constant)
                pos_x = pg_h5["position/x"]
                @test !is_constant_component(pos_x)
                dat = component_data(pos_x)
                @test length(dat) > 0

                # position/z should be constant (all zeros)
                pos_z = pg_h5["position/z"]
                @test is_constant_component(pos_z)
                val = constant_component_value(pos_z)
                @test val == 0.0
            end
        end

        @testset "write_component_data constant path" begin
            tmpfile = tempname() * ".h5"
            try
                u = PMDUnit("m", 1.0, "length")
                h5open(tmpfile, "w") do h5
                    # All-same values → should write as constant
                    write_component_data(h5, "const_test", fill(42.0, 10); unit=u)
                end
                h5open(tmpfile, "r") do h5
                    g = h5["const_test"]
                    @test is_constant_component(g)
                    @test constant_component_value(g) == 42.0
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "write_component_data array path" begin
            tmpfile = tempname() * ".h5"
            try
                data = [1.0, 2.0, 3.0]
                h5open(tmpfile, "w") do h5
                    write_component_data(h5, "arr_test", data)
                end
                h5open(tmpfile, "r") do h5
                    ds = h5["arr_test"]
                    @test !is_constant_component(ds)
                    dat = component_data(ds)
                    @test dat ≈ data
                end
            finally
                rm(tmpfile; force=true)
            end
        end

        @testset "component_data with unit_factor" begin
            tmpfile = tempname() * ".h5"
            try
                data = [1.0, 2.0, 3.0]
                u = PMDUnit("m", 1.0, "length")
                h5open(tmpfile, "w") do h5
                    write_component_data(h5, "scaled", data; unit=u)
                end
                h5open(tmpfile, "r") do h5
                    # unit_factor=2 should multiply data by 2
                    dat = component_data(h5["scaled"]; unit_factor=2)
                    @test dat ≈ data .* 2
                end
            finally
                rm(tmpfile; force=true)
            end
        end
    end

    # ── all_components / component_str / offset_component_name ─────────────

    @testset "all_components" begin
        h5open(_PARTICLE_FILE, "r") do h5
            pp = particle_paths(h5)
            pg_h5 = h5[pp[1]]
            if !haskey(pg_h5, "position")
                pg_h5 = pg_h5[first(keys(pg_h5))]
            end
            comps = all_components(pg_h5)
            @test "position/x" in comps
            @test "position/y" in comps
            @test "position/z" in comps
            @test "momentum/x" in comps
            @test "momentum/y" in comps
            @test "momentum/z" in comps
            @test "time"        in comps
            @test "weight"      in comps
        end
    end

    @testset "component_str" begin
        h5open(_PARTICLE_FILE, "r") do h5
            pp = particle_paths(h5)
            pg_h5 = h5[pp[1]]
            if !haskey(pg_h5, "position")
                pg_h5 = pg_h5[first(keys(pg_h5))]
            end
            s = component_str(pg_h5, "position/x")
            @test occursin("position/x", s)
            @test occursin("length", s)

            s_const = component_str(pg_h5, "position/z")
            @test occursin("constant", s_const)
        end
    end

    @testset "offset_component_name" begin
        @test offset_component_name("position/x")  == "positionOffset/x"
        @test offset_component_name("momentum/z")  == "momentumOffset/z"
        @test offset_component_name("time")         == "timeOffset"
        @test offset_component_name("weight")       == "weightOffset"
    end

    # ── write_pmd_bunch / particle_data_dict ───────────────────────────────

    @testset "particle_data_dict" begin
        pg = single_particle(x=1e-3, px=5e5, pz=1e9, weight=1e-15)
        d = particle_data_dict(pg)

        @test d["x"]  == pg.x
        @test d["px"] == pg.px
        @test d["y"]  == pg.y
        @test d["py"] == pg.py
        @test d["z"]  == pg.z
        @test d["pz"] == pg.pz
        @test d["t"]  == pg.t
        @test d["status"]    == pg.status
        @test d["weight"]    == pg.weight
        @test d["id"]        == pg.id
        @test d["species"]   == "electron"
        @test d["n_particle"] == 1
        @test d["charge"]    == charge(pg)
    end

    @testset "write_pmd_bunch round-trip" begin
        pg = single_particle(x=2e-3, px=3e5, pz=2e9, weight=5e-16)
        tmpfile = tempname() * ".h5"
        try
            h5open(tmpfile, "w") do h5
                pmd_init(h5)
                g = create_group(h5, "data/1")
                write_pmd_bunch(g, particle_data_dict(pg))
            end
            pg2 = ParticleGroup(tmpfile)
            @test isapprox(pg, pg2; rtol=1e-12)
        finally
            rm(tmpfile; force=true)
        end
    end

    # ── write_pmd_field ────────────────────────────────────────────────────

    @testset "write_pmd_field" begin
        tmpfile = tempname() * ".h5"
        try
            field_data = Dict{String, Any}(
                "attrs" => Dict{String, Any}(
                    "eleAnchorPt"       => "center",
                    "gridGeometry"      => "rectangular",
                    "axisLabels"        => ["x", "y", "z"],
                    "gridLowerBound"    => [0, 0, 0],
                    "gridOriginOffset"  => [0.0, 0.0, 0.0],
                    "gridSpacing"       => [0.001, 0.001, 0.001],
                    "gridSize"          => [10, 10, 10],
                    "harmonic"          => [0],
                ),
                "components" => Dict{String, Any}(
                    "Ex" => fill(0.0, 10, 10, 10),
                )
            )
            h5open(tmpfile, "w") do h5
                pmd_field_init(h5)
                g = create_group(h5, "ExternalFieldPath/1")
                write_pmd_field(g, field_data)
            end
            h5open(tmpfile, "r") do h5
                fp = field_paths(h5)
                @test length(fp) == 1
                g = h5[fp[1]]
                a, _ = load_field_attrs(attrs(g))
                @test a["gridGeometry"] == "rectangular"
                @test a["eleAnchorPt"]  == "center"
                @test a["gridSize"]     == [10, 10, 10]
                # Ex should be written as constant (all zeros)
                @test haskey(g, "Ex")
            end
        finally
            rm(tmpfile; force=true)
        end
    end

    # ── load_field_attrs with solenoid file ────────────────────────────────

    if isfile(_FIELD_FILE)
        @testset "load_field_attrs on solenoid" begin
            h5open(_FIELD_FILE, "r") do h5
                fp = field_paths(h5)
                g = h5[fp[1]]
                a, other = load_field_attrs(attrs(g))

                @test a["gridGeometry"]    == "cylindrical"
                @test a["eleAnchorPt"]     == "center"
                @test a["axisLabels"]      == ["r", "theta", "z"]
                @test length(a["gridSize"]) == 3
                @test a["gridSize"][1] > 0
                @test haskey(a, "gridSpacing")
                @test haskey(a, "gridOriginOffset")
                @test haskey(a, "gridLowerBound")
                @test haskey(a, "harmonic")

                # Field group should have magneticField
                @test haskey(g, "magneticField")
                mg = g["magneticField"]
                @test "r" in keys(mg)
                @test "z" in keys(mg)
            end
        end
    end

end  # I/O (openPMD)
