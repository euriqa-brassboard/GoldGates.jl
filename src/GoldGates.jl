#

module GoldGates

import ProtoBuf as PB
using JSON
import MSSim: Sequence as Seq

include("gold_gates_pb.jl")

import .gold_gates_pb: ParticipationFactor, XXSolution, SysMetadata, Modes,
    SystemParams, GateSolutionSet

Base.eltype(f::ParticipationFactor) = eltype(f.factors)
Base.length(f::ParticipationFactor) = length(f.factors)
Base.getindex(f::ParticipationFactor, i) = getindex(f.factors, i)
Base.setindex!(f::ParticipationFactor, v, i) = setindex!(f.factors, v, i)

_load_json(v, ::Type{ParticipationFactor}) = ParticipationFactor(v)
_to_json(v::ParticipationFactor) = v.factors
Base.:(==)(v1::ParticipationFactor, v2::ParticipationFactor) = v1.factors == v2.factors
Base.hash(v::ParticipationFactor, h::UInt) = hash(v.factors, hash(:ParticipationFactor, h))

_load_json(v, ::Type{SysMetadata}) = SysMetadata(get(v, "units", Dict{String,String}()),
                                                 get(v, "structure", ""))
_to_json(v::SysMetadata) = Dict("units"=>v.units, "structure"=>v.structure)
Base.:(==)(v1::SysMetadata, v2::SysMetadata) =
    v1.units == v2.units && v1.structure == v2.structure
Base.hash(v::SysMetadata, h::UInt) =
    hash(v.units, hash(v.structure, hash(:SysMetadata, h)))

_load_json(v, ::Type{Modes}) = Modes(radial1=v["radial"],
                                     radial2=get(v, "radial2", Float64[]),
                                     axial=get(v, "axial", Float64[]))
_to_json(v::Modes) = Dict("radial"=>v.radial1, "radial2"=>v.radial2, "axial"=>v.axial)
Base.:(==)(v1::Modes, v2::Modes) =
    v1.radial1 == v2.radial1 && v1.radial2 == v2.radial2 && v1.axial == v2.axial
Base.hash(v::Modes, h::UInt) =
    hash(v.radial1, hash(v.radial2, hash(v.axial, hash(:Modes, h))))

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
Base.:(==)(v1::SystemParams, v2::SystemParams) =
    (v1.modes == v2.modes && v1.lamb_dicke_parameters == v2.lamb_dicke_parameters &&
    v1.participation_factors == v2.participation_factors &&
    v1.dac_terms == v2.dac_terms && v1.metadata == v2.metadata)
function Base.hash(v::SystemParams, h::UInt)
    h = hash(:SystemParams, h)
    h = hash(v.metadata, h)
    h = hash(v.dac_terms, h)
    h = hash(v.participation_factors, h)
    h = hash(v.lamb_dicke_parameters, h)
    return hash(v.modes, h)
end

_load_json(v, ::Type{XXSolution}) = XXSolution(nsteps=v["nsteps"],
                                               angle_sign=v["angle_sign"],
                                               time=v["time"],
                                               phase=v["phase"],
                                               phase_slope=v["phase_slope"],
                                               amp=v["amp"],
                                               amp_slope=v["amp_slope"],
                                               metadata=get(v, "metadata", ""))
_to_json(v::XXSolution) = Dict("nsteps"=>v.nsteps,
                               "angle_sign"=>v.angle_sign,
                               "time"=>v.time,
                               "phase"=>v.phase, "phase_slope"=>v.phase_slope,
                               "amp"=>v.amp, "amp_slope"=>v.amp_slope,
                               "metadata"=>v.metadata)
Base.:(==)(v1::XXSolution, v2::XXSolution) =
    (v1.nsteps == v2.nsteps && v1.angle_sign == v2.angle_sign && v1.time == v2.time &&
    v1.phase == v2.phase && v1.phase_slope == v2.phase_slope &&
    v1.amp == v2.amp && v1.amp_slope == v2.amp_slope)
function Base.hash(v::XXSolution, h::UInt)
    h = hash(:XXSolution, h)
    h = hash(v.amp_slope, h)
    h = hash(v.amp, h)
    h = hash(v.phase_slope, h)
    h = hash(v.phase, h)
    h = hash(v.time, h)
    h = hash(v.angle_sign, h)
    return hash(v.nsteps, h)
end

_load_json(v, ::Type{GateSolutionSet}) =
    GateSolutionSet(params=SystemParams(modes=_load_json(v["modes"], Modes)),
                    XX=Dict(k=>_load_json(v, XXSolution) for (k, v) in v["XX"]))
_to_json(v::GateSolutionSet) = Dict("modes"=>_to_json(v.params.modes),
                                    "XX"=>Dict(k=>_to_json(v) for (k, v) in v.XX))
Base.:(==)(v1::GateSolutionSet, v2::GateSolutionSet) =
    v1.params == v2.params && v1.XX == v2.XX
Base.hash(v::GateSolutionSet, h::UInt) =
    hash(v.params, hash(v.XX, hash(:GateSolutionSet, h)))

function _verify_num_modes(modes::Modes)
    if isempty(modes.radial1)
        throw(ArgumentError("Radial mode missing"))
    end
    nmodes = length(modes.radial1)
    if !isempty(modes.radial2) && length(modes.radial2) != nmodes
        throw(ArgumentError("Mismatch between radial mode number"))
    end
    if !isempty(modes.axial) && length(modes.axial) != nmodes
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
    if sol.angle_sign != 1 && sol.angle_sign != -1
        throw(ArgumentError("Invalid solution angle sign"))
    end
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

function XXSolution(params::Seq.RawParams, angle_sign; metadata=nothing, kws...)
    d = Seq.gate_solution_info(params; kws...)

    d["angle_sign"] = sign(angle_sign)

    # Optionally attach arbitrary metadata
    if metadata !== nothing
        d["metadata"] = metadata
    end

    return verify(_load_json(d, XXSolution))
end

include("solution_candidates.jl")
import .SolutionCandidates: Candidate,
    load_candidates_files, load_candidates_dir, load_candidates_dirs,
    save_candidates

function vv2m(vv)
    nrow = length(vv)
    if nrow == 0
        return Matrix{eltype(eltype(vv))}(undef, 0, 0)
    end
    ncol = length(vv[1])
    m = Matrix{Float64}(undef, nrow, ncol)
    for i in 1:nrow
        v = vv[i]
        for j in 1:ncol
            e = v[j]
            @inbounds m[i, j] = e
        end
    end
    return m
end
m2vv(m) = [m[i, :] for i in 1:size(m, 1)]

function set_mode_weight!(weights, ηs, bij, ion1, ion2)
    nions = length(ηs)
    for i in 1:nions
        weights[i] = bij[i, ion1] * bij[i, ion2] * ηs[i]^2
    end
    return weights
end
include("optimizers.jl")

include("server.jl")

end
