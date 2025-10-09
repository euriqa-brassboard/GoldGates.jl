#

module TestOptimizers

using GoldGates
using GoldGates: Optimizers, Candidate
using AMO.Utils: ThreadObjectPool

import MSSim: Sequence as Seq

using Test

fs = [2.41, 2.50, 2.58]
nseg = 40
amp_ratio = 0.7
pre_pool = ThreadObjectPool() do
    return Optimizers.PreOptimizer{nseg}(
        2π .* fs, 2π .* (fs .+ 0.3); amp_ratio=amp_ratio,
        tmin=100, tmax=300, ntimes=6, ωmin=2π * 2.2, ωmax=2π * 2.7)
end
candidates = Optimizers.opt_all_rounds!(pre_pool, 30)

@testset "PreOptimizer" begin
    @test !isempty(candidates)
end

@testset "PairChecker" begin
    checker = Optimizers.PairChecker([0.1, 0.2], 2, 0.05)
    function dummy_cand(area, areaδ)
        cand = Candidate([0.1, 0.2, 0.3],
                         Seq.SolutionProperties(0.2, [1.2, 3], [1e-3, 1e-8], [1e-2, -2e-7],
                                                [1e-2, 1e-3], area, areaδ))
    end
    @test Optimizers.check(checker, dummy_cand([10.0, 5.1], [0.0, 0.0]))
    @test Optimizers.check(checker, dummy_cand([-10.0, -5.1], [0.0, 0.0]))
    @test !Optimizers.check(checker, dummy_cand([-10.0, 5.1], [0.0, 0.0]))
    @test !Optimizers.check(checker, dummy_cand([10.0, -5.1], [0.0, 0.0]))

    @test !Optimizers.check(checker, dummy_cand([10.0, 5.1], [10.0, 4.9]))
    @test !Optimizers.check(checker, dummy_cand([10.0, 5.1], [-10.0, -4.9]))
    @test Optimizers.check(checker, dummy_cand([10.0, 5.1], [10.0, -4.9]))
    @test Optimizers.check(checker, dummy_cand([10.0, 5.1], [-10.0, 4.9]))
end

end
