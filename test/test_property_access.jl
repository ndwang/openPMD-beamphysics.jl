using AtomicAndPhysicalConstants: Species, massof, chargeof, C_LIGHT
import Statistics: mean, std
import StatsBase

const _PROP_BEAM = ParticleGroup(joinpath(@__DIR__, "initial_particles.h5"))

@testset "Property Access" begin

    @testset "pg[key] direct fields" begin
        pg = _PROP_BEAM

        @test pg["x"]      === pg.x
        @test pg["px"]     === pg.px
        @test pg["y"]      === pg.y
        @test pg["py"]     === pg.py
        @test pg["z"]      === pg.z
        @test pg["pz"]     === pg.pz
        @test pg["t"]      === pg.t
        @test pg["status"] === pg.status
        @test pg["weight"] === pg.weight
        @test pg["id"]     === pg.id
    end

    @testset "pg[\"z/c\"]" begin
        pg = single_particle(z=0.3)
        @test pg["z/c"] ≈ [0.3 / C_LIGHT]
    end

    @testset "Statistical prefix keys" begin
        @testset "sigma_" begin
            pg = _PROP_BEAM
            @test pg["sigma_x"] ≈ std(pg, "x")
            @test pg["sigma_x"] > 0
        end

        @testset "mean_" begin
            pg = _PROP_BEAM
            @test pg["mean_x"] ≈ mean(pg, "x")
        end

        @testset "min_ / max_ / ptp_" begin
            pg = _PROP_BEAM
            @test pg["min_x"]  ≈ minimum(pg.x)
            @test pg["max_x"]  ≈ maximum(pg.x)
            @test pg["ptp_x"]  ≈ maximum(pg.x) - minimum(pg.x)
        end

        @testset "delta_" begin
            pg = _PROP_BEAM
            d = pg["delta_x"]
            @test length(d) == length(pg)
            @test d ≈ pg.x .- mean(pg, "x")
        end

        @testset "cov_" begin
            pg = _PROP_BEAM
            c = pg["cov_x_px"]
            @test isa(c, Number)

            c_xx = pg["cov_x_x"]
            @test c_xx ≈ pg["sigma_x"]^2  rtol=1e-10
        end
    end

    @testset "Standalone functions" begin
        @testset "nalive / ndead" begin
            pg = single_particle() + single_particle() + single_particle()
            @test nalive(pg) == 3
            @test ndead(pg) == 0

            pg.status[1] = 0
            pg.status[2] = 0
            @test nalive(pg) == 1
            @test ndead(pg) == 2
        end

        @testset "mass / species_charge / charge" begin
            pg = _PROP_BEAM
            @test mass(pg) == massof(Species("electron"))
            @test species_charge(pg) == chargeof(Species("electron"))
            @test charge(pg) ≈ sum(pg.weight)
        end

        @testset "energy / momentum / kinetic_energy" begin
            pz = 1e9
            pg = single_particle(pz=pz)
            m = massof(Species("electron"))

            @test momentum(pg)[1] ≈ abs(pz)
            @test energy(pg)[1]   ≈ sqrt(pz^2 + m^2)
            @test kinetic_energy(pg)[1] ≈ energy(pg)[1] - m
        end

        @testset "gamma / beta / beta_x / beta_y / beta_z" begin
            pz = 1e9
            pg = single_particle(px=1e5, py=-2e5, pz=pz)
            m = massof(Species("electron"))
            E = energy(pg)[1]
            p = momentum(pg)[1]

            @test gamma(pg)[1]  ≈ E / m
            @test beta(pg)[1]   ≈ p / E
            @test beta_x(pg)[1] ≈ 1e5 / E
            @test beta_y(pg)[1] ≈ -2e5 / E
            @test beta_z(pg)[1] ≈ pz / E
        end

        @testset "xp / yp" begin
            pg = single_particle(px=1e6, py=-2e6, pz=1e9)
            @test xp(pg)[1] ≈ 1e6 / 1e9
            @test yp(pg)[1] ≈ -2e6 / 1e9
        end

        @testset "r / theta" begin
            pg = single_particle(x=3.0, y=4.0)
            @test r(pg)[1]     ≈ 5.0
            @test theta(pg)[1] ≈ atan(4.0, 3.0)

            pg0 = single_particle(x=1.0, y=0.0)
            @test theta(pg0)[1] ≈ 0.0
        end

        @testset "pr / ptheta" begin
            # Particle on x-axis: pr = px, ptheta = py
            pg = single_particle(x=1.0, y=0.0, px=3e6, py=4e6)
            @test pr(pg)[1]     ≈ 3e6  atol=1e-6
            @test ptheta(pg)[1] ≈ 4e6  atol=1e-6

            # Particle on y-axis: pr = py, ptheta = -px
            pg2 = single_particle(x=0.0, y=1.0, px=3e6, py=4e6)
            @test pr(pg2)[1]     ≈ 4e6  atol=1e-6
            @test ptheta(pg2)[1] ≈ -3e6 atol=1e-6
        end

        @testset "Lz" begin
            pg = single_particle(x=1.0, py=2.0, y=3.0, px=4.0)
            @test Lz(pg)[1] ≈ 1.0 * 2.0 - 3.0 * 4.0
        end

        @testset "set_charge!" begin
            pg = single_particle(weight=1e-15) + single_particle(weight=2e-15)
            new_q = 5e-15
            set_charge!(pg, new_q)
            @test charge(pg) ≈ new_q

            @test_throws ArgumentError set_charge!(pg, -1.0)
        end

        @testset "average_current" begin
            pg = _PROP_BEAM
            I = average_current(pg)
            @test I ≈ charge(pg) / ptp(pg, "t")
        end
    end

    @testset "Emittance and normalized coordinates" begin
        @testset "norm_emit_x / norm_emit_y / norm_emit_4d" begin
            pg = _PROP_BEAM
            ex = norm_emit_x(pg)
            ey = norm_emit_y(pg)
            e4 = norm_emit_4d(pg)
            @test ex > 0
            @test ey > 0
            @test e4 > 0
        end

        @testset "x_bar, px_bar, Jx, Jy" begin
            pg = _PROP_BEAM

            xb  = x_bar(pg)
            pxb = px_bar(pg)
            @test length(xb)  == length(pg)
            @test length(pxb) == length(pg)

            jx = Jx(pg)
            jy = Jy(pg)
            @test all(jx .>= 0)
            @test all(jy .>= 0)
        end
    end

    @testset "pg[key] matches standalone function" begin
        pg = _PROP_BEAM

        @test pg["energy"]         == energy(pg)
        @test pg["kinetic_energy"] == kinetic_energy(pg)
        @test pg["p"]              == momentum(pg)
        @test pg["gamma"]          == gamma(pg)
        @test pg["beta"]           == beta(pg)
        @test pg["beta_x"]         == beta_x(pg)
        @test pg["beta_y"]         == beta_y(pg)
        @test pg["beta_z"]         == beta_z(pg)
        @test pg["xp"]             == xp(pg)
        @test pg["yp"]             == yp(pg)
        @test pg["r"]              == r(pg)
        @test pg["theta"]          == theta(pg)
        @test pg["pr"]             == pr(pg)
        @test pg["ptheta"]         == ptheta(pg)
        @test pg["Lz"]             == Lz(pg)
        @test pg["n_particle"]     == length(pg)
        @test pg["n_alive"]        == nalive(pg)
        @test pg["n_dead"]         == ndead(pg)
        @test pg["mass"]           == mass(pg)
        @test pg["species_charge"] == species_charge(pg)
        @test pg["charge"]         == charge(pg)
    end

    @testset "Coordinate checks" begin
        @test in_z_coordinates(_PROP_BEAM)
        @test !in_t_coordinates(_PROP_BEAM)

        pg_t = single_particle(t=0.0, z=0.1) + single_particle(t=0.0, z=0.2)
        @test in_t_coordinates(pg_t)
        @test !in_z_coordinates(pg_t)
    end

end  # Property Access
