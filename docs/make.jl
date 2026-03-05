using Documenter
using BeamPhysics
using HDF5
using Plots

const PlotsExt = Base.get_extension(BeamPhysics, :BeamPhysicsPlotsExt)

makedocs(
    sitename = "BeamPhysics.jl",
    authors = "ndwang and contributors",
    modules = [BeamPhysics, PlotsExt],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://ndwang.github.io/BeamPhysics.jl",
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => [
            "Particles" => "particles.md",
            "Field Mesh" => "fieldmesh.md",
            "Statistics" => "statistics.md",
            "I/O" => "io.md",
            "Plotting" => "plotting.md",
        ],
        "API Reference" => "api.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references, :docs_block],
)

deploydocs(
    repo = "github.com/ndwang/BeamPhysics.jl.git",
    devbranch = "main",
)
