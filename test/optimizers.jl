#

module TestOptimizers

using GoldGates
using GoldGates.Optimizers
using AMO.Utils: ThreadObjectPool

using Test

@testset "PreOptimizer" begin
    fs = [2.41, 2.50, 2.58]
    nseg = 40
    amp_ratio = 0.7
    pre_pool = ThreadObjectPool() do
        return Optimizers.PreOptimizer{nseg}(
            2π .* fs, 2π .* (fs .+ 0.3); amp_ratio=amp_ratio,
            tmin=100, tmax=300, ntimes=6, ωmin=2π * 2.2, ωmax=2π * 2.7)
    end
    candidates = Optimizers.opt_all_rounds!(pre_pool, 30, GoldGates.Candidate[])
    @test !isempty(candidates)
end

end
