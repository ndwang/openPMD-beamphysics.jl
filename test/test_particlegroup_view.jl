@testset "ParticleGroupView basic behavior" begin
    # Construct a small ParticleGroup with three particles
    pg1 = single_particle(x = 0.0)
    pg2 = single_particle(x = 1.0)
    pg3 = single_particle(x = 2.0)
    pg = pg1 + pg2 + pg3

    @test length(pg) == 3

    # Create a view onto a subset of particles
    view_inds = 2:3
    pgv = view(pg, view_inds)

    @test length(pgv) == length(view_inds)
    @test all(pgv.x .== pg.x[view_inds])
    @test all(pgv.y .== pg.y[view_inds])
    @test all(pgv.z .== pg.z[view_inds])

    # Mutations through the view should reflect in the parent
    pgv.x .= 42.0
    @test all(pg.x[view_inds] .== 42.0)

    # Derived quantities should be consistent between parent and view
    @test all(momentum(pgv) .== momentum(pg)[view_inds])
    @test all(energy(pgv) .== energy(pg)[view_inds])

    # Drift only the viewed subset and ensure other particles are unchanged
    x_before = copy(pg.x)
    drift!(pgv, 1e-12)
    x_after = pg.x

    other_inds = setdiff(1:length(pg), collect(view_inds))
    @test all(x_after[other_inds] .== x_before[other_inds])
end

