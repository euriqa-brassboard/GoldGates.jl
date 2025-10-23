#

module Server

import ..Optimizers: PreOptimizer, opt_all_rounds!
import ..SolutionCandidates: Candidate, load_candidates_files, save_candidates

using AMO.Utils: ThreadObjectPool

using Base.Threads

mutable struct CandidateManager
    const nions::Int
    const datadir::String
    function CandidateManager(nions, datadir)
        mkpath(datadir)
        return new(nions, datadir)
    end
end

function _list_datafile(dir)
    id_map = Dict{Int,String}()
    ids = Int[]
    for f in readdir(dir)
        m = match(r"^data([0-9.]+).binpb$", f)
        if m !== nothing
            id = parse(Int, m[1])
            id_map[id] = joinpath(dir, f)
            push!(ids, id)
        end
    end
    sort!(ids)
    return ids, String[id_map[id] for id in ids]
end

function eachcandidate(@specialize(cb), mgr::CandidateManager)
    for fname in readdir(mgr.datadir)
        m = match(r"^nseg_[0-9]+_ampratio_[0-9.]+$", fname)
        if m !== nothing
            _, files = _list_datafile(joinpath(mgr.datadir, fname))
            load_candidates_files(files) do f, meta, candidates
                cb(meta, candidates)
            end
        end
    end
end

function searchcandidate(mgr::CandidateManager, nseg, nrounds, ωs, ωs2;
                         amp_ratio=0.7, kws...)
    @assert length(ωs) == mgr.nions
    @assert length(ωs2) == mgr.nions * 2
    pre_pool = ThreadObjectPool() do
        return PreOptimizer{nseg}(ωs, ωs2; amp_ratio=amp_ratio, kws...)
    end

    block_size = 2000

    datadir = joinpath(mgr.datadir, "nseg_$(nseg)_ampratio_$(amp_ratio)")
    mkpath(datadir)
    ids, files  = _list_datafile(datadir)

    meta = Dict("nseg"=>nseg, amp_ratio=>amp_ratio)

    candidates = Candidate[]
    if isempty(ids)
        start_block = Ref(1)
    else
        start_block = Ref(ids[end])
        load_candidates_files([files[end]], candidates=candidates, meta=meta)
    end

    lock = ReentrantLock()

    opt_all_rounds!(pre_pool, nrounds) do new_candidates
        @lock lock begin
            append!(candidates, new_candidates)
            ncandidates = length(candidates)
            if ncandidates > 10 * block_size
                nblocks = ncandidates ÷ block_size
                nsaving = nblocks * block_size
                @info "Saving $nblocks blocks of candidates"
                save_candidates(joinpath(datadir, "data"), @view(candidates[1:nsaving]),
                                meta, block_size=block_size, start_block=start_block[],
                                format=:protobuf)
                start_block[] += nblocks
                deleteat!(candidates, 1:nsaving)
            end
        end
        empty!(new_candidates)
    end
end

end
