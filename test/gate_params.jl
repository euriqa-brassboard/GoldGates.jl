#

module GateParams

include("test_utils.jl")

import GoldGates as GG
using GoldGates: ParticipationFactor, XXSolution, SysMetadata, Modes,
    SystemParams, GateSolutionSet

import MSSim: Sequence as Seq

import ProtoBuf as PB

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

function _check_io(v::T, format) where T
    io = IOBuffer()
    write(io, v; format=format)
    seekstart(io)
    v2 = read(io, T; format=format)
    @test v == v2
    @test isbitstype(T) || v !== v2
    @test hash(v) == hash(v2)
end

function check_io(v::T, test_json=true) where T
    _check_io(v, :protobuf)
    if test_json
        _check_io(v, :json)
    end
    @test_throws ArgumentError "Unknown format" write(IOBuffer(), v; format=:unknown)
    @test_throws ArgumentError "Unknown format" read(IOBuffer(), T; format=:unknown)
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

    @test GG._verify_num_modes(modes) == 3
    @test GG._verify_num_modes(modes2) == 3
    @test_throws ArgumentError "Radial mode missing" GG._verify_num_modes(Modes())
    @test_throws(ArgumentError, "Mismatch between radial mode number",
                 GG._verify_num_modes(Modes(radial1=[1, 2, 3], radial2=[1, 2])))
    @test_throws(ArgumentError, "Mismatch between radial and axial mode number",
                 GG._verify_num_modes(Modes(radial1=[1, 2, 3], axial=[1, 0])))

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
    modes2 = Modes(radial1=[1.2, 3.4, 3.4])
    params2 = SystemParams(modes=modes2,
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

    @test GG.verify(params) === params
    @test GG.verify(params2) === params2
    @test_throws ArgumentError "Radial mode missing" GG.verify(SystemParams(modes=Modes()))
    @test_throws ArgumentError "Mode missing" GG.verify(SystemParams())
    @test GG.verify(SystemParams(modes=modes)) == SystemParams(modes=modes)
    @test GG.verify(SystemParams(modes=modes2)) == SystemParams(modes=modes2)
    @test_throws(ArgumentError, "Lamb Dikie parameters length mismatch",
                 GG.verify(SystemParams(modes=modes,
                                        lamb_dicke_parameters=[0.1, 0.3, 0.2, 2])))
    @test_throws(ArgumentError, "Mode participation factors length mismatch",
                 GG.verify(SystemParams(modes=modes, participation_factors=[pf1, pf2])))
    @test_throws(ArgumentError, "Mode participation factors ion count mismatch",
                 GG.verify(SystemParams(modes=modes,
                                        participation_factors=[pf1, pf2,
                                                               ParticipationFactor()])))

    check_json(params2)
    params3 = GG._load_json(GG._to_json(params), SystemParams)
    @test params3 != params
    @test isempty(params3.modes.radial2)
    @test isempty(params3.modes.axial)

    check_io(params, false)
    check_io(params2)
end

@testset "XXSolution" begin
    sol = XXSolution(3, -1, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4], [0.1, -0.1, 0.2],
                     [0.3, 1.0, 0.3], [0.2, -0.2, 0.4])
    sol2 = XXSolution(nsteps=3, angle_sign=-1,
                      time=[1.0, 2.0, 3.0],
                      phase=[0.1, 0.3, 0.4],
                      phase_slope=[0.1, -0.1, 0.2],
                      amp=[0.3, 1.0, 0.3],
                      amp_slope=[0.2, -0.2, 0.4])
    @test sol == sol2
    sol2.time[2] = 0.1
    sol2.phase[3] = 0.1
    sol2.phase_slope[1] = -0.2
    sol2.amp[3] = 0.2
    sol2.amp_slope[1] = -0.1
    @test sol != sol2

    @test PB.default_values(XXSolution) == (;nsteps=0, angle_sign=0, time=Float64[],
                                            phase=Float64[], phase_slope=Float64[],
                                            amp=Float64[], amp_slope=Float64[])
    @test PB.field_numbers(XXSolution) == (;nsteps=1, angle_sign=2, time=3,
                                           phase=4, phase_slope=5, amp=6, amp_slope=7)

    check_pb(XXSolution())
    check_pb(sol)
    check_pb(sol2)

    @test GG.verify(sol) === sol
    @test GG.verify(sol2) === sol2
    @test_throws(ArgumentError, "Invalid solution angle sign",
                 GG.verify(XXSolution(3, 0, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4],
                                      [0.1, -0.1, 0.2], [0.3, 1.0, 0.3],
                                      [0.2, -0.2, 0.4])))
    @test_throws(ArgumentError, "Solution time count mismatch",
                 GG.verify(XXSolution(3, -1, [1.0, 2.0, 3.0, 1.2], [0.1, 0.3, 0.4],
                                      [0.1, -0.1, 0.2], [0.3, 1.0, 0.3],
                                      [0.2, -0.2, 0.4])))
    @test_throws(ArgumentError, "Solution phase count mismatch",
                 GG.verify(XXSolution(3, -1, [1.0, 2.0, 3.0], [0.3, 0.4],
                                      [0.1, -0.1, 0.2], [0.3, 1.0, 0.3],
                                      [0.2, -0.2, 0.4])))
    @test_throws(ArgumentError, "Solution phase slope count mismatch",
                 GG.verify(XXSolution(3, -1, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4],
                                      [0.1, -0.1], [0.3, 1.0, 0.3], [0.2, -0.2, 0.4])))
    @test_throws(ArgumentError, "Solution amp count mismatch",
                 GG.verify(XXSolution(3, -1, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4],
                                      [0.1, -0.1, 0.2], [0.3, 1.0, 0.3, 0.2],
                                      [0.2, -0.2, 0.4])))
    @test_throws(ArgumentError, "Solution amp slope count mismatch",
                 GG.verify(XXSolution(3, -1, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4],
                                      [0.1, -0.1, 0.2], [0.3, 1.0, 0.3], [0.2, -0.2])))

    check_json(sol)
    check_json(sol2)

    params = Seq.RawParams(rand(25))
    time = params.args[1:5:end]
    amp = params.args[2:5:end]
    amp_slope = params.args[3:5:end]
    phase = params.args[4:5:end] ./ 2π
    phase_slope = params.args[5:5:end] ./ 2π
    @test XXSolution(params, 1.2, normalize_amp=false) ==
        XXSolution(nsteps=5, angle_sign=1, time=time, amp=amp, amp_slope=amp_slope,
                   phase=phase, phase_slope=phase_slope)
    @test XXSolution(params, -0.1, normalize_amp=false) ==
        XXSolution(nsteps=5, angle_sign=-1, time=time, amp=amp, amp_slope=amp_slope,
                   phase=phase, phase_slope=phase_slope)
end

@testset "GateSolutionSet" begin
    sol = XXSolution(3, -1, [1.0, 2.0, 3.0], [0.1, 0.3, 0.4], [0.1, -0.1, 0.2],
                     [0.3, 1.0, 0.3], [0.2, -0.2, 0.4])
    sol2 = XXSolution(nsteps=3, angle_sign=1,
                      time=[1.0, 2.0, 1.0],
                      phase=[0.1, 0.3, 0.5],
                      phase_slope=[0.1, -1.1, 0.2],
                      amp=[0.3, 1.0, 0.1],
                      amp_slope=[0.1, -0.2, 0.4])
    modes = Modes([1.2, 3.4, 3.4], [0.3, 0.1, 0.4], [1, 2, 3])
    pf1 = ParticipationFactor([0.1, 0.2, 0.3])
    pf2 = ParticipationFactor([0.2, 0.1, -0.3])
    pf3 = ParticipationFactor([0.5, -0.1, -0.4])
    sysmeta = SysMetadata(Dict("A"=>"B", "C"=>"D"), "abc")
    params = SystemParams(modes, [0.1, 0.3, 0.2], [pf1, pf2, pf3],
                          Dict("x2"=>0.1, "xy"=>0.2), sysmeta)
    params2 = SystemParams(modes=Modes(radial1=[1.2, 3.4, 3.4]))
    solset = GateSolutionSet(params, Dict("1"=>sol, "2"=>sol2))
    solset2 = GateSolutionSet(params2, Dict("2"=>sol, "1"=>sol2))
    @test solset != solset2

    @test PB.default_values(GateSolutionSet) == (;params=nothing,
                                                 XX=Dict{String,XXSolution}())
    @test PB.field_numbers(GateSolutionSet) == (;params=1, XX=2)

    check_pb(GateSolutionSet())
    check_pb(solset)
    check_pb(solset2)

    @test GG.verify(GateSolutionSet()) == GateSolutionSet()
    @test GG.verify(solset) === solset
    @test GG.verify(solset2) === solset2
    @test_throws(ArgumentError, "Solution time count mismatch",
                 GG.verify(GateSolutionSet(XX=Dict("1,2"=>XXSolution(nsteps=3,
                                                                     angle_sign=-1,
                                                                     time=[1.0])))))
    @test_throws ArgumentError "Mode missing" GG.verify(GateSolutionSet(params=SystemParams()))

    check_json(solset2)
    solset3 = GG._load_json(GG._to_json(solset), GateSolutionSet)
    @test solset3 != solset
    @test solset3.XX == solset.XX
    @test solset3.params.modes == solset.params.modes
    @test isempty(solset3.params.lamb_dicke_parameters)
    @test isempty(solset3.params.participation_factors)
    @test isempty(solset3.params.dac_terms)
    @test isnothing(solset3.params.metadata)

    check_io(GateSolutionSet(), false)
    check_io(solset, false)
    check_io(solset2)
end

end
