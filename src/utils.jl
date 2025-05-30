function ptp(x::AbstractArray)
    return maximum(x) - minimum(x)
end
