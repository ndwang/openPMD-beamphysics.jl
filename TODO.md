# TODO

## P0 — Bugs (wrong results or crashes)

- [x] `units.jl:38` — `electric_potential` dimension tuple has length and mass swapped. Fixed to `(2, 1, -3, -1, 0, 0, 0)`
- [x] `statistics.jl:268` / `particles.jl` — `split_particles` now guarantees exactly `n_chunks` slices; errors if fewer particles than chunks
- [x] `plot.jl:264-269` — Fixed ellipse code: replaced Python method calls with Julia functions, added `_sigma_ellipse` helper
- [x] `plot.jl:260` — Replaced `vmin=1e-20` with `clims=(1e-20, Inf)`
- [x] `plot.jl:348` — Replaced `vmin=1e-16` with `clims=(1e-16, Inf)`
- [x] `units.jl:361-381` — Fixed NoUnits/eV/c round-trip: `write_unit_h5` uses `UNIT_SYMBOL` lookup, `read_unit_h5` uses `KNOWN_UNITS` lookup before falling back to `uparse`
- [x] `units.jl:364` — Fixed by same `KNOWN_UNITS`/`UNIT_SYMBOL` approach above

## P1 — Important improvements

- [ ] No test suite. Add `test/runtests.jl` with at minimum: HDF5 round-trip, unit conversion, statistics correctness, slice statistics
- [x] `Project.toml` — Relaxed compat bounds to accept patch releases (e.g. `HDF5 = "0.17"`). Removed stdlib compat entries
- [x] `Project.toml` — Added `julia = "1.10"` compat bound
- [x] `Project.toml` — Removed unused `Distributions` dependency
- [x] `particles.jl:137-139` — Replaced `@assert` with `throw(ArgumentError(...))`

## P2 — Performance

- [x] `energy()` — compute once in `drift!`, `drift_to_z!` instead of redundantly via `beta_x`/`beta_y`/`beta_z`
- [x] `particles.jl` — `in_z_coordinates`/`in_t_coordinates` now use `all(==(v[1]), v)`: zero allocation, short-circuits
- [x] `particles.jl` — `join_particle_groups` now uses `reduce(vcat, ...)` instead of `vcat([]...)`
- [ ] `particles.jl:358-373` — `split_particles` copies all data per slice. Consider a `ParticleGroupView` wrapper backed by index ranges for read-only use cases like slice statistics

## P3 — Code quality

- [ ] `statistics.jl:114` — Misleading comment `# - 1` on dispersion delta. Remove or clarify that the subtraction is unnecessary because `cov` subtracts the mean
- [ ] `particles.jl:66-67` — `xp`/`yp` divide by `pz` with no zero-check. Will silently produce `Inf`/`NaN`. Document or guard
- [ ] `particles.jl:218-224` — `average_current` can return `Inf` when all particles are at the same position and time. Document or guard
- [ ] `writers.jl:131` — `write_pmd_field` unconditionally converts field data to complex (`complex.(val)`), doubling memory for real-valued fields
- [x] `units.jl` — `write_unit_h5` now uses `UnitfulToOpenPMD(u)` directly (fixed in P0 commit)
- [x] `plot.jl` — Replaced `width=0` with `lw=0` in all `bar` calls
- [ ] `plot.jl:389-630` — 240 lines of commented-out field mesh plotting code with Python-style method calls. Remove from source; track in git history
- [x] `utils.jl` — Moved `ptp` into `particles.jl`, deleted `utils.jl`
- [ ] `particles.jl:99-133` — `DERIVED_PROPERTIES` dict duplicates all standalone functions as lambdas. If a new derived property function is added but not registered in the dict, `pg["prop"]` silently fails. Consider deriving from the function names automatically
- [x] Module renamed from `openPMD_beamphysics` to `OpenPMDBeamphysics`
