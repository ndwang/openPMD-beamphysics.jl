using Test
using BeamPhysics

@testset "BeamPhysics" begin
    include("test_units.jl")
    include("test_particlegroup_construction.jl")
    include("test_property_access.jl")
    include("test_statistics_beamphysics.jl")
    include("test_particlegroup_view.jl")
    include("test_io_openpmd.jl")
    include("test_fieldmesh.jl")
end

