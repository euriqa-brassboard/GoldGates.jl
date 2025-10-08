#

using Test

import ProtoBuf as PB

macro test_throws(arg1, arg2)
    :(Test.@test_throws $(esc(arg1)) $(esc(arg2)))
end

macro test_throws(arg1, arg2, arg3)
    if VERSION >= v"1.13.0-DEV"
        :(Test.@test_throws $(esc(arg1)) $(esc(arg2)) $(esc(arg3)))
    else
        :(Test.@test_throws $(esc(arg1)) $(esc(arg3)))
    end
end

@test_throws ErrorException error("abc cde")
@test_throws ErrorException "abc" error("abc cde")

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
