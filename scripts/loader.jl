#!/usr/bin/julia


"""
Loader
Module that contains functions for reading and writing Gold system
solutions and parameter files with JSON formatting.
"""
module Loader
using JSON
using LinearAlgebra

"Struct to store relevant system parameters and read/write them in JSON format"
struct GoldParams
    radial_modes::Vector{Float64}
    etas::Vector{Float64}
    participation_factors::Matrix{Float64}
    dac_terms::Dict{String, Float64}

    "Construct GoldParams directly"
    function GoldParams(radial_modes::Vector{Float64},
                        etas::Vector{Float64},
                        participation_factors::Matrix{Float64},
                        dac_terms::Dict{String, Float64})
        # Validate square matrix
        if size(participation_factors, 1) != size(participation_factors, 2)
            throw(ArgumentError("Participation matrix must be square"))
        end
        
        # Validate consistent dimensions
        if !(length(radial_modes) == length(etas) == size(participation_factors, 1))
            throw(ArgumentError("All parameter arrays must have consistent dimensions"))
        end
        
        new(radial_modes, etas, participation_factors, dac_terms)
    end

    "Given a filepath to a Gold Params JSON file, returns a struct containing the information from the file"
    function GoldParams(filepath::String)::GoldParams
        data = JSON.parsefile(filepath)
        # Convert to appropriate types
        radial_modes = Vector{Float64}(data["radial_modes"])
        # vcat takes the transpose, so permutedims avoids the transpose. Enforces type
        part_factors = vcat(permutedims.(Vector{Float64}.(data["participation_factors"]))...) 
        etas = Vector{Float64}(data["lamb_dicke_parameters"])
        dac_terms = Dict(k => v::Float64 for (k, v) in data["dac_terms"])
        
        return GoldParams(radial_modes, etas, part_factors, dac_terms)
    end

end

"Writes the passed GoldParams struct to a JSON file at filepath"
function to_json(gp::GoldParams, filepath::String)::Nothing
    data = Dict(
        "radial_modes" => gp.radial_modes,
        "lamb_dicke_parameters" => gp.etas,
        "participation_factors" => gp.participation_factors,
        "dac_terms" => gp.dac_terms
    )
    
    open(filepath, "w") do f
        JSON3.write(f, data)
    end
    
    return nothing
end

end