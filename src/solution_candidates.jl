#

module SolutionCandidates

import MSSim: Sequence as Seq

using JSON
using Printf
import ProtoBuf as PB

using Base.Threads

struct Candidate
    param::Vector{Float64}
    props::Union{Nothing,Seq.SolutionProperties}
end
PB.default_values(::Type{Candidate}) = (;param = Vector{Float64}(), props = nothing)
PB.field_numbers(::Type{Candidate}) = (;param = 1, props = 2)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:Candidate})
    param = PB.BufferedVector{Float64}()
    props = Ref{Union{Nothing,Seq.SolutionProperties}}(nothing)
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            PB.decode!(d, wire_type, param)
        elseif field_number == 2
            PB.decode!(d, props)
        else
            PB.skip(d, wire_type)
        end
    end
    return Candidate(param[], props[])
end

function PB.encode(e::PB.AbstractProtoEncoder, x::Candidate)
    initpos = position(e.io)
    !isempty(x.param) && PB.encode(e, 1, x.param)
    !isnothing(x.props) && PB.encode(e, 2, x.props)
    return position(e.io) - initpos
end
function PB._encoded_size(x::Candidate)
    encoded_size = 0
    !isempty(x.param) && (encoded_size += PB._encoded_size(x.param, 1))
    !isnothing(x.props) && (encoded_size += PB._encoded_size(x.props, 2))
    return encoded_size
end

struct Candidates
    candidates::Vector{Candidate}
    meta::String
end
PB.default_values(::Type{Candidates}) = (;candidates = Vector{Candidate}(), meta = "")
PB.field_numbers(::Type{Candidates}) = (;candidates = 1, meta = 2)

function PB.decode(d::PB.AbstractProtoDecoder, ::Type{<:Candidates})
    candidates = PB.BufferedVector{Candidate}()
    meta = ""
    while !PB.message_done(d)
        field_number, wire_type = PB.decode_tag(d)
        if field_number == 1
            PB.decode!(d, candidates)
        elseif field_number == 2
            meta = PB.decode(d, String)
        else
            PB.skip(d, wire_type)
        end
    end
    return Candidates(candidates[], meta)
end

function PB.encode(e::PB.AbstractProtoEncoder, x::Candidates)
    initpos = position(e.io)
    !isempty(x.candidates) && PB.encode(e, 1, x.candidates)
    !isempty(x.meta) && PB.encode(e, 2, x.meta)
    return position(e.io) - initpos
end
function PB._encoded_size(x::Candidates)
    encoded_size = 0
    !isempty(x.candidates) && (encoded_size += PB._encoded_size(x.candidates, 1))
    !isempty(x.meta) && (encoded_size += PB._encoded_size(x.meta, 2))
    return encoded_size
end

Base.Dict(c::Candidate) = Dict("param"=>c.param, "props"=>Dict(c.props))
Candidate(d::Dict{<:AbstractString}) = Candidate(copy(d["param"]),
                                                 Seq.SolutionProperties(d["props"]))

function load_candidates_json(io::IO)
    data = JSON.parse(io)
    meta = get(data, "meta", nothing)
    return meta, Candidate.(data["candidates"])
end

function load_candidates_pb(io::IO)
    decoder = PB.ProtoDecoder(io)
    candidates = PB.decode(decoder, Candidates)
    if candidates.meta == ""
        meta = nothing
    else
        meta = JSON.parse(candidates.meta)
    end
    return meta, candidates.candidates
end

function load_candidates_files(files; candidates=Candidate[], meta=nothing)
    results = Dict{String,Tuple{Dict,Vector{Candidate}}}()
    lock = ReentrantLock()
    @threads :greedy for f in files
        if endswith(f, ".binpb")
            res = open(load_candidates_pb, f)
        else
            res = open(load_candidates_json, f)
        end
        @lock lock results[f] = res
    end
    for f in files
        file_meta, file_candidates = results[f]
        if meta === nothing
            meta = file_meta
        elseif file_meta != meta
            error("Metadata mismatch")
        end
        append!(candidates, file_candidates)
    end
    return meta, candidates
end

function load_candidates_dir(dir; kwargs...)
    return load_candidates_files(readdir(dir, join=true); kwargs...)
end

function load_candidates_dirs(dirs; kwargs...)
    files = String[]
    for dir in dirs
        append!(files, readdir(dir, join=true))
    end
    return load_candidates_files(files; kwargs...)
end

function save_candidates(prefix, candidates, meta; block_size=2000, format=:protobuf)
    ncandidates = length(candidates)
    meta_str = nothing
    if format === :protobuf && meta !== nothing
        meta_str = JSON.json(meta)
    end
    @threads :greedy for (i, start_idx) in enumerate(1:block_size:ncandidates)
        end_idx = min(start_idx + block_size - 1, ncandidates)
        if format === :protobuf
            open("$(prefix)$(@sprintf("%06d", i)).binpb", "w") do io
                encoder = PB.ProtoEncoder(io)
                PB.encode(encoder, Candidates(candidates[start_idx:end_idx], meta_str))
            end
        elseif format === :json
            open("$(prefix)$(@sprintf("%06d", i)).json", "w") do io
                d = Dict("meta"=>meta,
                         "candidates"=>Dict.(@view candidates[start_idx:end_idx]))
                JSON.print(io, d, 2)
            end
        else
            throw(ArgumentError("Unknown format $(format) for gate solution candidates"))
        end
    end
end

end
