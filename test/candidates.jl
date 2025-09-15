#

module TestCandidates

include("test_utils.jl")

using Printf

import MSSim: Sequence as Seq

using GoldGates: Candidate, load_candidates_files, load_candidates_dir,
    load_candidates_dirs, save_candidates
using GoldGates.SolutionCandidates: Candidates

@testset "Candidates" begin
    cand = Candidate([0.1, 0.2, 0.3], nothing)
    cand2 = Candidate([0.1, 0.2, 0.3],
                      Seq.SolutionProperties(0.2, [1.2, 3], [0.12, 0.1], [0.0, -0.01],
                                             [1.0, 0.2], [0.1, 2.3], [-0.2, 0.3]))
    @test cand == cand
    @test cand2 == cand2
    @test cand != cand2
    @test hash(cand) != hash(cand2)

    @test PB.default_values(Candidate) == (;param=Float64[], props=nothing)
    @test PB.field_numbers(Candidate) == (;param=1, props=2)

    check_pb(cand)
    check_pb(cand2)

    js = Dict(cand)
    @test js == Dict("param"=>cand.param)
    @test Candidate(js) == cand

    js2 = Dict(cand2)
    @test js2 == Dict("param"=>cand2.param,
                      "props"=>Dict(cand2.props))
    @test Candidate(js2) == cand2


    cs = Candidates([cand], "")
    cs2 = Candidates([cand, cand2], "{\"a\": 2}")
    @test cs == cs
    @test cs2 == cs2
    @test hash(cs) != hash(cs2)

    @test PB.default_values(Candidates) == (;candidates=Candidate[], meta="")
    @test PB.field_numbers(Candidates) == (;candidates=1, meta=2)

    check_pb(cs)
    check_pb(cs2)

    mktempdir() do d
        @test_throws ArgumentError "Unknown format" save_candidates(joinpath(d, "data_"), [cand], nothing, format=:unknown)

        function test_save_dir(name, meta, format)
            suffix = format === :protobuf ? "binpb" : "json"
            path = joinpath(d, name)
            mkdir(path)
            cands = [Candidate(rand(3), Seq.SolutionProperties(rand(), rand(2), rand(2),
                                                               rand(2), rand(2),
                                                               rand(2), rand(2)))
                      for _ in 1:2000]
            save_candidates(joinpath(path, "data_"), cands, meta,
                            block_size=10, format=format)
            @test load_candidates_dir(path) == (meta, cands)
            @test load_candidates_dirs([path]) == (meta, cands)
            @test readdir(path) == ["data_$(@sprintf("%06d", i)).$suffix" for i in 1:200]
            for i in 1:200
                f = "data_$(@sprintf("%06d", i)).$suffix"
                offset = (i - 1) * 10
                @test load_candidates_files([joinpath(path, f)]) == (meta, cands[offset + 1:offset + 10])
            end
            if meta === nothing
                @test load_candidates_dir(path, meta=123) == (123, cands)
            else
                @test_throws ErrorException "Metadata mismatch" load_candidates_dir(path, meta=123)
            end
            cands_out = Candidate[]
            @test load_candidates_dir(path, meta=meta, candidates=cands_out) === (meta, cands_out)
            @test cands_out == cands

            rev_meta, rev_cands =
                load_candidates_files([joinpath(path,
                                                "data_$(@sprintf("%06d", i)).$suffix")
                                       for i in 200:-1:1])
            @test rev_meta == meta
            for i1 in 1:200
                i2 = 201 - i1
                offset1 = (i1 - 1) * 10
                offset2 = (i2 - 1) * 10
                @test rev_cands[offset1 + 1:offset1 + 10] == cands[offset2 + 1:offset2 + 10]
            end
            return path, cands
        end

        path1, cands1 = test_save_dir("save_pb1", nothing, :protobuf)
        path2, cands2 = test_save_dir("save_js1", nothing, :json)
        path3, cands3 = test_save_dir("save_pb2", Dict("a"=>"b", "c"=>[1, 2]), :protobuf)
        path4, cands4 = test_save_dir("save_js2", Dict("a"=>2, "c"=>[2, 3]), :json)
        path5, cands5 = test_save_dir("save_pb3", ["a", "b", "3", [1, 2]], :protobuf)
        path6, cands6 = test_save_dir("save_js3", ["a", 2, "c", [2, 3]], :json)

        @test load_candidates_dirs([path1, path2]) == (nothing, [cands1; cands2])
        @test load_candidates_dirs([path1, path3]) == (Dict("a"=>"b", "c"=>[1, 2]), [cands1; cands3])
        @test load_candidates_dirs([path4, path2]) == (Dict("a"=>2, "c"=>[2, 3]), [cands4; cands2])
        @test_throws ErrorException "Metadata mismatch" load_candidates_dirs([path3, path4])
        @test_throws ErrorException "Metadata mismatch" load_candidates_dirs([path5, path4])
        @test_throws ErrorException "Metadata mismatch" load_candidates_dirs([path3, path6])
    end
end

end
