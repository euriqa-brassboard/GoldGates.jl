#

module GateParams

using GoldGates
using GoldGates: ParticipationFactor, XXSolution, SysMetadata, Modes,
    SystemParams, GateSolutionSet

import ProtoBuf as PB

using Test

@testset "ParticipationFactor" begin
    pf = ParticipationFactor([0.1, 0.2, 0.3])
    pf2 = ParticipationFactor(factors=[0.1, 0.2, 0.3])
    @test pf.factors == pf2.factors
    @test length(pf) == 3
    @test eltype(pf) === Float64
    @test pf[1] == 0.1
    @test pf[2] == 0.2
    @test pf[3] == 0.3
    pf[2] = 0.9
    @test pf[1] == 0.1
    @test pf[2] == 0.9
    @test pf[3] == 0.3

    @test PB.default_values(ParticipationFactor) == (;factors=Float64[])
    @test PB.field_numbers(ParticipationFactor) == (;factors=1)

    io = IOBuffer()
    encoder = PB.ProtoEncoder(io)
    PB.encode(encoder, pf)
    @test PB._encoded_size(pf) == io.size
    seekstart(io)

    decoder = PB.ProtoDecoder(io)
    pf3 = PB.decode(decoder, ParticipationFactor)
    @test pf.factors == pf3.factors

    js = GoldGates._to_json(pf2)
    @test js === pf2.factors
    pf4 = GoldGates._load_json(js, ParticipationFactor)
    @test pf4.factors === pf2.factors
end

end
