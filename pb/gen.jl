#

using ProtoBuf

protojl("gold_gates.proto", @__DIR__, joinpath(@__DIR__, "../src/");
        add_kwarg_constructors=true)
