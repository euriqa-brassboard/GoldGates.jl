#

using Test

macro test_throws(arg1, arg2)
    :(Test.@test_throws $(esc(arg1)) $(esc(arg2)))
end

macro test_throws(arg1, arg2, arg3)
    if VERSION >= v"1.12"
        :(Test.@test_throws $(esc(arg1)) $(esc(arg2)) $(esc(arg3)))
    else
        :(Test.@test_throws $(esc(arg1)) $(esc(arg3)))
    end
end

@test_throws ErrorException error("abc cde")
@test_throws ErrorException "abc" error("abc cde")
