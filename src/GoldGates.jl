#

module GoldGates

import ProtoBuf as PB
using JSON
import MSSim: Sequence as Seq

include("gold_gates_pb.jl")

import .gold_gates_pb: ParticipationFactor, XXSolution, SysMetadata, Modes,
    SystemParams, GateSolutionSet

Base.length(f::ParticipationFactor) = length(f.factors)
Base.getindex(f::ParticipationFactor, i) = getindex(f.factors, i)
Base.setindex!(f::ParticipationFactor, v, i) = setindex!(f.factors, v, i)

_load_json(v, ::Type{ParticipationFactor}) = ParticipationFactor(v)
_to_json(v::ParticipationFactor) = v.factors

_load_json(v, ::Type{SysMetadata}) = SysMetadata(get(v, "units", Dict{String,String}()),
                                                 get(v, "structure", ""))
_to_json(v::SysMetadata) = Dict("units"=>v.units, "structure"=>v.structure)

_load_json(v, ::Type{Modes}) = Modes(radial1=v["radial"],
                                     radial2=get(v, "radial2", Float64[]),
                                     axial=get(v, "axial", Float64[]))
_to_json(v::Modes) = Dict("radial"=>v.radial1, "radial2"=>v.radial2, "axial"=>v.axial)

_load_json(v, ::Type{SystemParams}) =
    SystemParams(modes=Modes(radial1=v["radial_modes"]),
                 lamb_dicke_parameters=v["lamb_dicke_parameters"],
                 participation_factors=[_load_json(ve, ParticipationFactor)
                                        for ve in v["participation_factors"]],
                 dac_terms=v["dac_terms"],
                 metadata=_load_json(v["metadata"], SysMetadata))
_to_json(v::SystemParams) = Dict("radial_modes"=>v.modes.radial1,
                                 "lamb_dicke_parameters"=>v.lamb_dicke_parameters,
                                 "participation_factors"=>[_to_json(ve) for ve in
                                                               v.participation_factors],
                                 "dac_terms"=>v.dac_terms,
                                 "metadata"=>_to_json(v.metadata))

_load_json(v, ::Type{XXSolution}) = XXSolution(nsteps=v["nsteps"],
                                               angle_sign=v["angle_sign"],
                                               time=v["time"],
                                               phase=v["phase"],
                                               phase_slope=v["phase_slope"],
                                               amp=v["amp"],
                                               amp_slope=v["amp_slope"])
_to_json(v::XXSolution) = Dict("nsteps"=>v.nsteps,
                               "angle_sign"=>v.angle_sign,
                               "time"=>v.time,
                               "phase"=>v.phase, "phase_slope"=>v.phase_slope,
                               "amp"=>v.amp, "amp_slope"=>v.amp_slope)

_load_json(v, ::Type{GateSolutionSet}) =
    GateSolutionSet(params=SystemParams(modes=_load_json(v["modes"], Modes)),
                    XX=Dict(k=>_load_json(v, XXSolution) for (k, v) in v["XX"]))
_to_json(v::GateSolutionSet) = Dict("modes"=>_to_json(v.params.modes),
                                    "XX"=>Dict(k=>_to_json(v) for (k, v) in v.XX))

function _verify_num_modes(modes::Modes)
    if isempty(modes.radial1)
        throw(ArgumentError("Radial mode missing"))
    end
    nmodes = length(modes.radial1)
    if !isempty(modes.radial2) && length(modes.radial2) == nmodes
        throw(ArgumentError("Mismatch between radial mode number"))
    end
    if !isempty(modes.axial) && length(modes.axial) == nmodes
        throw(ArgumentError("Mismatch between radial and axial mode number"))
    end
    return nmodes
end

function verify(params::SystemParams)
    if params.modes === nothing
        throw(ArgumentError("Mode missing"))
    end
    nmodes = _verify_num_modes(params.modes)
    if !isempty(params.lamb_dicke_parameters) && length(params.lamb_dicke_parameters) != nmodes
        throw(ArgumentError("Lamb Dikie parameters length mismatch"))
    end
    if !isempty(params.participation_factors) && length(params.participation_factors) != nmodes
        throw(ArgumentError("Mode participation factors length mismatch"))
    end
    for factors in params.participation_factors
        if length(factors) != nmodes
            throw(ArgumentError("Mode participation factors ion count mismatch"))
        end
    end
    return params
end

function verify(sol::XXSolution)
    nsteps = sol.nsteps
    if length(sol.time) != nsteps
        throw(ArgumentError("Solution time count mismatch"))
    end
    if length(sol.phase) != nsteps
        throw(ArgumentError("Solution phase count mismatch"))
    end
    if length(sol.phase_slope) != nsteps
        throw(ArgumentError("Solution phase slope count mismatch"))
    end
    if length(sol.amp) != nsteps
        throw(ArgumentError("Solution amp count mismatch"))
    end
    if length(sol.amp_slope) != nsteps
        throw(ArgumentError("Solution amp slope count mismatch"))
    end
    return sol
end

function verify(sol_set::GateSolutionSet)
    if sol_set.params !== nothing
        verify(sol_set.params)
    end
    for (k, sol) in sol_set.XX
        verify(sol)
    end
    return sol_set
end

function Base.read(io::IO, ::Type{SystemParams}; format=:protobuf)
    if format === :json
        return verify(_load_json(JSON.parse(io), SystemParams))
    elseif format === :protobuf
        decoder = PB.ProtoDecoder(io)
        return verify(PB.decode(decoder, SystemParams))
    else
        throw(ArgumentError("Unknown format $(format) for system parameters"))
    end
end

function Base.write(io::IO, v::SystemParams; format=:protobuf)
    if format === :json
        JSON.print(io, _to_json(v), 2)
    elseif format === :protobuf
        encoder = PB.ProtoEncoder(io)
        PB.encode(encoder, v)
    else
        throw(ArgumentError("Unknown format $(format) for system parameters"))
    end
    return
end

function Base.read(io::IO, ::Type{GateSolutionSet}; format=:protobuf)
    if format === :json
        return verify(_load_json(JSON.parse(io), GateSolutionSet))
    elseif format === :protobuf
        decoder = PB.ProtoDecoder(io)
        return verify(PB.decode(decoder, GateSolutionSet))
    else
        throw(ArgumentError("Unknown format $(format) for gate solution set"))
    end
end

function Base.write(io::IO, v::GateSolutionSet; format=:protobuf)
    if format === :json
        JSON.print(io, _to_json(v), 2)
    elseif format === :protobuf
        encoder = PB.ProtoEncoder(io)
        PB.encode(encoder, v)
    else
        throw(ArgumentError("Unknown format $(format) for gate solution set"))
    end
    return
end

function XXSolution(params::Seq.RawParams, angle_sign)
    d = Seq.gate_solution_info(params)
    d["angle_sign"] = sign(angle_sign)
    return verify(_load_json(d, XXSolution))
end

include("thread_utils.jl")
include("solution_candidates.jl")
import .SolutionCandidates: Candidate,
    load_candidates_files, load_candidates_dir, load_candidates_dirs,
    save_candidates

end
