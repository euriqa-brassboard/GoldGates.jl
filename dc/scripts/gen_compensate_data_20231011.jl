#!/usr/bin/julia

using TrapVoltages
using TrapVoltages: Solutions, Potentials, Units
using LinearAlgebra

using HDF5

using JLD2
using JuMP
using HiGHS

using Base.Threads

const potential_file = ARGS[1]
const short_map = Potentials.load_short_map(joinpath(@__DIR__, "../data/electrode_short_202310.csv"))
const potential = Potentials.import_pillbox_64(potential_file,
                                               aliases=short_map, trap="phoenix")
const centers = Solutions.CenterTracker("phoenix")
const fitting = Potentials.Fitting(potential, orders=(8, 2, 2), sizes=(77, 5, 5))

function get_rf_center(xpos_um)
    xidx = Potentials.x_axis_to_index(potential, xpos_um)
    return (xidx, get(centers, xidx)...)
end

const xpos_ums = -3220:1330
const coeff_data = Vector{Dict{String,Any}}(undef, length(xpos_ums))
const coeff_data_nozx = Vector{Dict{String,Any}}(undef, length(xpos_ums))
const coeff_data_x2z = Vector{Dict{String,Any}}(undef, length(xpos_ums))

function compensation_coeff(xpos_um, mask::Val{Terms}) where Terms
    eles, coeff = Solutions.compensate_terms(fitting, get_rf_center(xpos_um),
                                             unit=Units.euriqa, mask=mask,
                                             min_num=20, min_dist=350)
    nterms = count(Terms)
    x0 = coeff \ Matrix(I, nterms, nterms)
    B = qr(coeff').Q[:, (nterms + 1):end]
    return Dict("electrodes"=>eles, "coefficients"=>coeff, "solution"=>x0,
                "free_solution"=>B, "xpos_um"=>xpos_um)

end

@time @threads for i in 1:length(xpos_ums)
    xpos_um = xpos_ums[i]
    coeff_data[i] = compensation_coeff(xpos_um, Solutions.TermMask())
    coeff_data_nozx[i] = compensation_coeff(xpos_um, Solutions.TermMask(zx=false))
    coeff_data_x2z[i] = compensation_coeff(xpos_um, Solutions.TermMask(x2z=true))
end

println("Coefficient generated")

function solve_all(coeff_data, termidx)
    model = Model(HiGHS.Optimizer)
    set_attribute(model, HiGHS.ComputeInfeasibilityCertificate(), false)
    set_attribute(model, "output_flag", false)
    xs = Vector{AffExpr}[]
    maxvs = VariableRef[]
    for data in coeff_data
        solution = data["solution"]
        x0 = @view(solution[:, termidx])
        B = data["free_solution"]
        nx, nt = size(B)
        @assert nx == length(x0)
        t = @variable(model, [1:nt])
        x = @expression(model, B * t .+ x0)
        maxv = @variable(model)
        @constraint(model, maxv .>= x)
        @constraint(model, maxv .>= .-x)
        push!(xs, x)
        push!(maxvs, maxv)
    end
    @objective(model, Min, sum(maxvs))
    JuMP.optimize!(model)
    return [[value(v) for v in x] for x in xs]
end

function pack_data(data, vals)
    solution = data["solution"]
    return Dict("electrodes"=>data["electrodes"], "voltages"=>vals,
                "xpos_um"=>data["xpos_um"])
end

function gen_solution(outputdir, coeff_data, termidx, term_name)
    println("Solving for $(term_name)")
    jldopen(joinpath(outputdir, "$(termidx).jld2"), "w") do fh
        voltages = @time(solve_all(coeff_data, termidx))
        transfer_solutions = [pack_data(data, vals) for (data, vals)
                                  in zip(coeff_data, voltages)]
        write(fh, "electrode_names", potential.electrode_names)
        write(fh, "transfer_solutions", transfer_solutions)
        write(fh, "termidx", termidx)
        write(fh, "termname", term_name)
    end
end

const term_names = ("dx", "dy", "dz", "xy", "yz", "zx", "z2", "x2", "x3", "x4")
const term_names_nozx = ("dx", "dy", "dz", "xy", "yz", "z2", "x2", "x3", "x4")
const term_names_x2z = ("dx", "dy", "dz", "xy", "yz", "zx", "z2", "x2", "x3", "x4", "x2z")

const output_prefix = joinpath(@__DIR__, "../data/compensate_20231011")

const jobs = Tuple{String,Vector{Dict{String,Any}},Int,String}[]
mkpath("$(output_prefix)")
for (termid, term_name) in enumerate(term_names)
    push!(jobs, ("$(output_prefix)", coeff_data, termid, term_name))
end
mkpath("$(output_prefix)_nozx")
for (termid, term_name) in enumerate(term_names_nozx)
    push!(jobs, ("$(output_prefix)_nozx", coeff_data_nozx, termid, term_name))
end
mkpath("$(output_prefix)_x2z")
for (termid, term_name) in enumerate(term_names_x2z)
    push!(jobs, ("$(output_prefix)_x2z", coeff_data_x2z, termid, term_name))
end

@time @threads for (outputdir, coeff_data, termid, term_name) in jobs
    gen_solution(outputdir, coeff_data, termid, term_name)
end
