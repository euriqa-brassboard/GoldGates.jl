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
Base.:(==)(v1::Candidate, v2::Candidate) = v1.param == v2.param && v1.props == v2.props
Base.hash(v::Candidate, h::UInt) = hash(v.param, hash(v.props, hash(:Candidate, h)))

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
Base.:(==)(v1::Candidates, v2::Candidates) =
    v1.candidates == v2.candidates && v1.meta == v2.meta
Base.hash(v::Candidates, h::UInt) = hash(v.candidates, hash(v.meta, hash(:Candidates, h)))

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

function Base.Dict(c::Candidate)
    d = Dict{String,Any}("param"=>c.param)
    if c.props !== nothing
        d["props"] = Dict(c.props)
    end
    return d
end
function Candidate(d::Dict{<:AbstractString})
    props = get(d, "props", nothing)
    return Candidate(d["param"],
                     props === nothing ? nothing : Seq.SolutionProperties(props))
end

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
    results = Dict{String,Tuple{Any,Vector{Candidate}}}()
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
        elseif file_meta != meta && file_meta !== nothing
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
    if format === :protobuf
        meta_str = meta === nothing ? "" : JSON.json(meta)
        @threads :greedy for (i, start_idx) in enumerate(1:block_size:ncandidates)
            end_idx = min(start_idx + block_size - 1, ncandidates)
            open("$(prefix)$(@sprintf("%06d", i)).binpb", "w") do io
                encoder = PB.ProtoEncoder(io)
                PB.encode(encoder, Candidates(candidates[start_idx:end_idx], meta_str))
            end
        end
    elseif format === :json
        @threads :greedy for (i, start_idx) in enumerate(1:block_size:ncandidates)
            end_idx = min(start_idx + block_size - 1, ncandidates)
            open("$(prefix)$(@sprintf("%06d", i)).json", "w") do io
                d = Dict("meta"=>meta,
                         "candidates"=>Dict.(@view candidates[start_idx:end_idx]))
                JSON.print(io, d, 2)
            end
        end
    else
        throw(ArgumentError("Unknown format $(format) for gate solution candidates"))
    end
end

end
