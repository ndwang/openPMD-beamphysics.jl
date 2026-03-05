using AtomicAndPhysicalConstants: Species, massof, chargeof

@testset "ParticleGroup Construction" begin

    @testset "Constructors" begin
        @testset "full constructor (11 args)" begin
            n = 5
            x  = randn(n)
            px = randn(n) .* 1e6
            y  = randn(n)
            py = randn(n) .* 1e6
            z  = randn(n)
            pz = ones(n) .* 1e9
            t  = zeros(n)
            status = ones(Int, n)
            weight = fill(1e-15, n)
            sp = Species("electron")
            id = collect(1:n)

            pg = ParticleGroup(x, px, y, py, z, pz, t, status, weight, sp, id)

            @test length(pg) == n
            @test pg.x === x
            @test pg.px === px
            @test pg.y === y
            @test pg.py === py
            @test pg.z === z
            @test pg.pz === pz
            @test pg.t === t
            @test pg.status === status
            @test pg.weight === weight
            @test pg.species === sp
            @test pg.id === id
        end

        @testset "constructor without id (10 args)" begin
            n = 3
            x  = [0.0, 1.0, 2.0]
            px = [1e6, 2e6, 3e6]
            y  = zeros(n)
            py = zeros(n)
            z  = zeros(n)
            pz = fill(1e9, n)
            t  = zeros(n)
            status = ones(Int, n)
            weight = fill(1e-15, n)
            sp = Species("electron")

            pg = ParticleGroup(x, px, y, py, z, pz, t, status, weight, sp)

            @test length(pg) == n
            @test pg.x == x
            @test pg.id == [1, 2, 3]
        end

        @testset "single_particle with defaults" begin
            pg = single_particle()

            @test length(pg) == 1
            @test pg.x  == [0.0]
            @test pg.px == [0.0]
            @test pg.y  == [0.0]
            @test pg.py == [0.0]
            @test pg.z  == [0.0]
            @test pg.pz == [0.0]
            @test pg.t  == [0.0]
            @test pg.status == [1]
            @test pg.weight == [1.0]
            @test pg.id == [1]
            @test nameof(pg.species) == "electron"
        end

        @testset "single_particle with custom values" begin
            pg = single_particle(
                x=1e-3, px=5e5, y=2e-3, py=-3e5,
                z=0.01, pz=1e9, t=1e-12,
                weight=2e-15, status=0, species=Species("proton")
            )

            @test pg.x  == [1e-3]
            @test pg.px == [5e5]
            @test pg.y  == [2e-3]
            @test pg.py == [-3e5]
            @test pg.z  == [0.01]
            @test pg.pz == [1e9]
            @test pg.t  == [1e-12]
            @test pg.status == [0]
            @test pg.weight == [2e-15]
            @test nameof(pg.species) == "proton"
        end

        @testset "type parameters" begin
            pg64 = single_particle(x=1.0, pz=1e9)
            @test pg64 isa ParticleGroup{Float64,Int64}

            pg32 = ParticleGroup(
                Float32[0, 1], Float32[1e6, 2e6],
                Float32[0, 0], Float32[0, 0],
                Float32[0, 0], Float32[1e9, 1e9],
                Float32[0, 0], Int32[1, 1],
                Float32[1e-15, 1e-15], Species("electron")
            )
            @test pg32 isa ParticleGroup{Float32,Int32}
            @test length(pg32) == 2
        end
    end

    @testset "Species handling" begin
        for name in ("electron", "proton")
            pg = single_particle(species=Species(name))
            @test nameof(pg.species) == name
            @test mass(pg) == massof(Species(name))
            @test species_charge(pg) == chargeof(Species(name))
        end
    end

    @testset "Length and indexing" begin
        @testset "length" begin
            @test length(single_particle()) == 1

            pg3 = single_particle() + single_particle() + single_particle()
            @test length(pg3) == 3
        end

        @testset "integer-range indexing" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)

            sub = pg[1:2]
            @test sub isa ParticleGroup
            @test length(sub) == 2
            @test sub.x == [1.0, 2.0]
            @test sub.pz == pg.pz[1:2]
        end

        @testset "vector indexing" begin
            pg = single_particle(x=10.0) + single_particle(x=20.0) + single_particle(x=30.0)

            sub = pg[[1, 3]]
            @test sub isa ParticleGroup
            @test length(sub) == 2
            @test sub.x == [10.0, 30.0]
        end

        @testset "BitVector indexing" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)

            sub = pg[BitVector([true, false, true])]
            @test length(sub) == 2
            @test sub.x == [1.0, 3.0]
        end
    end

    @testset "Joining (+)" begin
        @testset "two-way addition" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=2.0, pz=2e9)
            pg = pg1 + pg2

            @test length(pg) == 2
            @test pg.x  == [1.0, 2.0]
            @test pg.pz == [1e9, 2e9]
            @test pg.species === pg1.species
        end

        @testset "three-way addition" begin
            pg = single_particle(x=1.0) + single_particle(x=2.0) + single_particle(x=3.0)
            @test length(pg) == 3
            @test pg.x == [1.0, 2.0, 3.0]
        end
    end

    @testset "Equality and comparison" begin
        @testset "==" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=1.0, pz=1e9)
            @test pg1 == pg2

            pg3 = single_particle(x=2.0, pz=1e9)
            @test pg1 != pg3

            pg4 = single_particle(x=1.0, pz=1e9, species=Species("proton"))
            @test pg1 != pg4
        end

        @testset "isapprox" begin
            pg1 = single_particle(x=1.0, pz=1e9)
            pg2 = single_particle(x=1.0 + 1e-14, pz=1e9)
            @test isapprox(pg1, pg2)

            pg3 = single_particle(x=2.0, pz=1e9)
            @test !isapprox(pg1, pg3)
        end
    end

    @testset "show / display" begin
        pg = single_particle(weight=1e-15)
        s = sprint(show, pg)
        @test occursin("ParticleGroup", s)
        @test occursin("1 particles", s)
    end

end  # ParticleGroup Construction
