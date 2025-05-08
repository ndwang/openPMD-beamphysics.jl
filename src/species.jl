# Dictionaries for particle properties
const CHARGE_OF = Dict(
    "electron" => -e_charge,
    "positron" => e_charge,
    "proton" => e_charge,
    "H-" => -e_charge,
    "H2+" => e_charge,
    "muon" => -e_charge,
    "antimuon" => e_charge
)

const CHARGE_STATE = Dict(
    "electron" => -1,
    "positron" => 1,
    "proton" => 1,
    "H-" => -1,
    "H2+" => 1,
    "muon" => -1,
    "antimuon" => 1
)

const MASS_OF = Dict(
    "electron" => mec2,
    "positron" => mec2,
    "proton" => mpc2,
    "H-" => mhmc2,
    "H2+" => mH2pc2,
    "muon" => mmc2,
    "antimuon" => mmc2
)

# Functions
"""
    mass_of(species::String)

Returns the mass energy equivalent (in eV) for the given particle species.
Throws an error if the species is not available.
"""
function mass_of(species::String)
    if haskey(MASS_OF, species)
        return MASS_OF[species]
    end
    throw(ArgumentError("Species not available: $species"))
end

"""
    charge_of(species::String)

Returns the charge (in C) for the given particle species.
Throws an error if the species is not available.
"""
function charge_of(species::String)
    if haskey(CHARGE_OF, species)
        return CHARGE_OF[species]
    end
    throw(ArgumentError("Species not available: $species"))
end

"""
    charge_state(species::String)

Returns the charge state (as an integer) for the given particle species.
Throws an error if the species is not available.
"""
function charge_state(species::String)
    if haskey(CHARGE_STATE, species)
        return CHARGE_STATE[species]
    end
    throw(ArgumentError("Species not available: $species"))
end