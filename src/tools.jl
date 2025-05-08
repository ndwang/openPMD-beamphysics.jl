"""
Utility functions for openPMD-beamphysics
"""

"""
    fstr(s)

Makes a fixed string for h5 files by converting to bytes
"""
fstr(s::AbstractString) = Vector{UInt8}(s)

"""
    data_are_equal(d1, d2)

Simple utility to compare data in dicts.
Returns true only if all keys are the same and all data are the same.
"""
function data_are_equal(d1::Dict, d2::Dict)
    if keys(d1) != keys(d2)
        return false
    end
    
    for k in keys(d1)
        if !all(d1[k] .== d2[k])
            return false
        end
    end
    
    return true
end

# -----------------------------------------
# HDF5 utilities

"""
    decode_attr(a)

Decodes:
    ASCII strings and arrays of them to String and arrays of String
    single-length arrays to scalar
"""
function decode_attr(a)
    if a isa Vector{UInt8}
        return String(a)
    end
    
    if a isa Array
        if eltype(a) <: UInt8
            return String.(a)
        end
        if length(a) == 1
            return a[1]
        end
    end
    
    return a
end

"""
    decode_attrs(attrs)

Decodes all attributes in a dictionary
"""
decode_attrs(attrs::Dict) = Dict(k => decode_attr(v) for (k, v) in attrs)

"""
    encode_attr(a)

Encodes attribute for HDF5 storage.
See the inverse function: decode_attr
"""
function encode_attr(a)
    if a isa String
        return fstr(a)
    end
    
    if a isa Array
        if eltype(a) <: String
            return fstr.(a)
        end
    end
    
    return a
end

"""
    encode_attrs(attrs)

Encodes all attributes in a dictionary
"""
encode_attrs(attrs::Dict) = Dict(k => encode_attr(v) for (k, v) in attrs)

"""
    get_version()

Get the installed pmd-beamphysics version.
"""
function get_version()
    return VERSION
end

"""
    current_date_with_tzinfo()

Get current date and time with timezone information
"""
function current_date_with_tzinfo()
    return now(localzone())
end

"""
    pmd_format_date(dt)

Format datetime according to PMD format
"""
function pmd_format_date(dt::DateTime)
    return Dates.format(dt, "yyyy-mm-dd HH:MM:SS zzz")
end 