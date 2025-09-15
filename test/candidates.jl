#

module TestCandidates

include("test_utils.jl")

import MSSim: Sequence as Seq

using GoldGates: Candidate, load_candidates_files, load_candidates_dir,
    load_candidates_dirs, save_candidates

@testset "Candidates" begin
    cand = Candidate([0.1, 0.2, 0.3], nothing)
    cand2 = Candidate([0.1, 0.2, 0.3],
                      Seq.SolutionProperties(0.2, [1.2, 3], [0.12, 0.1], [0.0, -0.01],
                                             [1.0, 0.2], [0.1, 2.3], [-0.2, 0.3]))
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
end

end
