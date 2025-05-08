"""
    fieldmesh_rectangular_to_cylindrically_symmetric_data(fieldmesh)

Returns a rectangular fieldmesh from a cartesian one, extracting the y=0 slice.

# Arguments
- `fieldmesh`: FieldMesh object

# Returns
- Dictionary containing `attrs` and `components` for the new fieldmesh
"""
function fieldmesh_rectangular_to_cylindrically_symmetric_data(fieldmesh)
    # Verify geometry
    @assert fieldmesh.geometry == "rectangular" "Input fieldmesh must be rectangular"

    # Find central slice x=0, y=0
    ix = findfirst(x -> x == 0, fieldmesh.coord_vec("x"))
    @assert !isnothing(ix) "Could not find x=0 slice"
    
    iy = findfirst(y -> y == 0, fieldmesh.coord_vec("y"))
    @assert !isnothing(iy) "Could not find y=0 slice"

    # Form components
    components = Dict(
        "electricField/r" => fieldmesh.components["electricField/x"][ix:end, iy:iy, :],
        "magneticField/theta" => fieldmesh.components["magneticField/y"][ix:end, iy:iy, :],
        "electricField/z" => fieldmesh.components["electricField/z"][ix:end, iy:iy, :]
    )

    # Update attrs
    attrs = copy(fieldmesh.attrs)
    attrs["gridGeometry"] = "cylindrical"
    attrs["axisLabels"] = ("r", "theta", "z")
    attrs["gridOriginOffset"] = (0, 0, attrs["gridOriginOffset"][3])
    attrs["gridSize"] = size(components["electricField/r"])

    return Dict("attrs" => attrs, "components" => components)
end 