#

module GateParams

import GoldGates as GG
using GoldGates: ParticipationFactor, XXSolution, SysMetadata, Modes,
    SystemParams, GateSolutionSet

import ProtoBuf as PB

using Test

function _check_pb(v::T, extra_fld) where T
    io = IOBuffer()
    encoder = PB.ProtoEncoder(io)
    PB.encode(encoder, v)
    if extra_fld
        max_fld = maximum(PB.field_numbers(T))
        PB.encode(encoder, max_fld + 1, 1.2)
    else
        @test PB._encoded_size(v) == io.size
    end
    seekstart(io)

    decoder = PB.ProtoDecoder(io)
    v2 = PB.decode(decoder, T)
    @test v == v2
    @test isbitstype(T) || v !== v2
    @test hash(v) == hash(v2)
end

function check_pb(v)
    _check_pb(v, false)
    _check_pb(v, true)
end

function check_json(v::T) where T
    js = GG._to_json(v)
    v2 = GG._load_json(js, T)
    @test v == v2
    @test hash(v) == hash(v2)
    return js
end

@testset "ParticipationFactor" begin
    pf = ParticipationFactor([0.1, 0.2, 0.3])
    pf2 = ParticipationFactor(factors=[0.1, 0.2, 0.3])
    @test pf == pf2
    @test length(pf) == 3
    @test eltype(pf) === Float64
    @test pf[1] == 0.1
    @test pf[2] == 0.2
    @test pf[3] == 0.3
    pf[2] = 0.9
    @test pf[1] == 0.1
    @test pf[2] == 0.9
    @test pf[3] == 0.3
    @test pf != pf2
    @test hash(pf) != hash(pf2)

    @test PB.default_values(ParticipationFactor) == (;factors=Float64[])
    @test PB.field_numbers(ParticipationFactor) == (;factors=1)

    check_pb(ParticipationFactor())
    check_pb(pf)
    check_pb(pf2)

    @test check_json(pf2) === pf2.factors

    @test GG.vv2m(ParticipationFactor[]) == zeros(0, 0)
    @test GG.vv2m(Vector{Float64}[]) == zeros(0, 0)
    @test GG.vv2m([pf, pf2]) == [0.1 0.9 0.3; 0.1 0.2 0.3]
    @test GG.vv2m([[0.1, 0.9, 0.3], [0.1, 0.2, 0.3]]) == [0.1 0.9 0.3; 0.1 0.2 0.3]
    @test GG.m2vv([0.1 0.9 0.3; 0.1 0.2 0.3]) == [[0.1, 0.9, 0.3], [0.1, 0.2, 0.3]]
end

@testset "SysMetadata" begin
    sysmeta = SysMetadata(Dict("A"=>"B", "C"=>"D"), "abc")
    @test sysmeta.units == Dict("A"=>"B", "C"=>"D")
    @test sysmeta.structure == "abc"
    sysmeta2 = SysMetadata(units=Dict("A"=>"B", "C"=>"D"))
    @test sysmeta2.units == Dict("A"=>"B", "C"=>"D")
    @test sysmeta2.structure == ""
    @test sysmeta != sysmeta2

    @test PB.default_values(SysMetadata) == (;units=Dict{String,String}(), structure="")
    @test PB.field_numbers(SysMetadata) == (;units=1, structure=2)

    check_pb(SysMetadata())
    check_pb(sysmeta)
    check_pb(sysmeta2)

    @test check_json(sysmeta) == Dict("units"=>sysmeta.units,
                                      "structure"=>sysmeta.structure)
end

@testset "Modes" begin
    modes = Modes([1.2, 3.4, 3.4], [0.3, 0.1, 0.4], [1, 2, 3])
    @test modes.radial1 == [1.2, 3.4, 3.4]
    @test modes.radial2 == [0.3, 0.1, 0.4]
    @test modes.axial == [1.0, 2.0, 3.0]
    modes2 = Modes(radial1=[1.2, 3.4, 3.4])
    @test modes2.radial1 == [1.2, 3.4, 3.4]
    @test modes2.radial2 == Float64[]
    @test modes2.axial == Float64[]
    @test modes != modes2

    @test PB.default_values(Modes) == (;radial1=Float64[], radial2=Float64[],
                                       axial=Float64[])
    @test PB.field_numbers(Modes) == (;radial1=1, radial2=2, axial=3)

    check_pb(Modes())
    check_pb(modes)
    check_pb(modes2)

    @test check_json(modes) == Dict("radial"=>modes.radial1,
                                    "radial2"=>modes.radial2,
                                    "axial"=>modes.axial)
    @test GG._load_json(Dict("radial"=>[0.1, 0.2]), Modes) == Modes(radial1=[0.1, 0.2])
end

@testset "SystemParams" begin
    modes = Modes([1.2, 3.4, 3.4], [0.3, 0.1, 0.4], [1, 2, 3])
    pf1 = ParticipationFactor([0.1, 0.2, 0.3])
    pf2 = ParticipationFactor([0.2, 0.1, -0.3])
    pf3 = ParticipationFactor([0.5, -0.1, -0.4])
    sysmeta = SysMetadata(Dict("A"=>"B", "C"=>"D"), "abc")
    params = SystemParams(modes, [0.1, 0.3, 0.2], [pf1, pf2, pf3],
                          Dict("x2"=>0.1, "xy"=>0.2), sysmeta)
    params2 = SystemParams(modes=Modes(radial1=[1.2, 3.4, 3.4]),
                           participation_factors=[pf1, pf2, pf3],
                           lamb_dicke_parameters=[0.1, 0.3, 0.2],
                           dac_terms=Dict("x2"=>0.1, "xy"=>0.2), metadata=sysmeta)
    @test params != params2

    @test PB.default_values(SystemParams) ==
        (;modes = nothing, lamb_dicke_parameters = Float64[],
         participation_factors = ParticipationFactor[],
         dac_terms = Dict{String,Float64}(), metadata = nothing)
    @test PB.field_numbers(SystemParams) == (;modes = 1, lamb_dicke_parameters = 2,
                                             participation_factors = 3, dac_terms = 4,
                                             metadata = 5)

    check_pb(SystemParams())
    check_pb(params)
    check_pb(params2)

    check_json(params2)
    params3 = GG._load_json(GG._to_json(params), SystemParams)
    @test params3 != params
    @test isempty(params3.modes.radial2)
    @test isempty(params3.modes.axial)
end

end
