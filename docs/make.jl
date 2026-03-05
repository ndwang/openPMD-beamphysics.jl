using Documenter
using OpenPMDBeamphysics
using HDF5
using Plots

const PlotsExt = Base.get_extension(OpenPMDBeamphysics, :OpenPMDBeamphysicsPlotsExt)

makedocs(
    sitename = "OpenPMDBeamphysics.jl",
    authors = "ndwang and contributors",
    modules = [OpenPMDBeamphysics, PlotsExt],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://ndwang.github.io/openPMD-beamphysics",
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
    repo = "github.com/ndwang/openPMD-beamphysics.jl.git",
    devbranch = "main",
)
