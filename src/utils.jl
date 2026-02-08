using PyPlot
using MSSim: Sequence as Seq
using NLopt
using JSON
using Dates
using Printf

function amp_base_funcs(n::Integer; atol::Real = 1e-12)
    @assert n ≥ 1 "n must be at least 1"
    m = n + 1
    grid = [((i - 1) / (n / 2) - 1) for i in 1:m]
    fs_vec = [
        let k = k, grid = grid, m = m, atol = atol
            x -> begin
                exclude = k - 1
                lo = 1 + exclude
                hi = m - exclude
                if lo > hi
                    return 0.0
                end
                for j in lo:hi
                    if isapprox(x, grid[j]; atol = atol)
                        return 1.0
                    end
                end
                return 0.0
            end
        end for k in 1:cld(n+1, 2)
    ]
    return tuple(fs_vec...)
end

# Objective function for optimization
function _objfunc(vals)
    dis = vals[1]
    disδ = vals[2]
    area = vals[3]
    areaδ = vals[4]
    τ = vals[5]

    return 5 * dis + disδ / 100 + (abs(area) - π / 2)^2 * 100 + (areaδ / 1e4)^2
end

function setup_modes(sysparams, ion1, ion2, nions)
    modes = Seq.Modes()
    idx1 = ion1 + (nions + 1) ÷ 2
    idx2 = ion2 + (nions + 1) ÷ 2
    for i in 1:nions
        push!(modes,
            2π * sysparams.modes.radial1[i],
            sysparams.participation_factors[i][idx1] * sysparams.participation_factors[i][idx2] * sysparams.lamb_dicke_parameters[i]^2)
    end
    return modes
end

function setup_model(nseg, modes)
    objfunc = Opts.autodiff(_objfunc)
    buf_opt = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_allδ), Val(SS.mask_allδ))
    amp_funcs = amp_base_funcs(nseg)
    nlmodel = Seq.Objective(SL.pmask_full,
        ((:dis2, 0), (:disδ2, 0), (:area, 0),
            (:areaδ, 0), (:τ, 0)),
        objfunc, modes, buf_opt,
        freq=Seq.FreqSpec(false, sym=true),
        amp=Seq.AmpSpec(cb=amp_funcs, sym=true))
    return nlmodel
end

function setup_optimizer(nlmodel, sysparams; pitime=30, τmin=5, τmax=50, maxtime=10, min_mode_index=1, max_mode_index=3)
    Ωmax = π / (2 * pitime)
    ωmin = 2π * sysparams.modes.radial1[min_mode_index]
    ωmax = 2π * sysparams.modes.radial1[max_mode_index]

    nargs = Seq.nparams(nlmodel)
    tracker = Opts.NLVarTracker(nargs)
    for Ω in nlmodel.param.Ωs
        Opts.set_bound!(tracker, Ω, 0, Ωmax)
    end
    Opts.set_bound!(tracker, nlmodel.param.τ, τmin, τmax)
    for ω in nlmodel.param.ωs
        Opts.set_bound!(tracker, ω, ωmin, ωmax)
    end

    opt = NLopt.Opt(:LD_LBFGS, nargs)
    NLopt.min_objective!(opt, nlmodel)
    NLopt.lower_bounds!(opt, Opts.lower_bounds(tracker))
    NLopt.upper_bounds!(opt, Opts.upper_bounds(tracker))
    NLopt.maxtime!(opt, maxtime)

    return opt, tracker
end

function perturb_params(params, tracker, scale)
    result = copy(params)
    lb = Opts.lower_bounds(tracker)
    ub = Opts.upper_bounds(tracker)
    for i in eachindex(result)
        r = (ub[i] - lb[i]) * scale
        result[i] = clamp(params[i] + randn() * r, lb[i], ub[i])
    end
    return result
end

function run_optimization!(opt, tracker, nlmodel; threshold=-Inf,
                           initial_params=nothing, perturbation=0.05)
    iterations = initial_params === nothing ? 300 : 50
    best_obj = Inf
    best_params = nothing

    @time for i in 1:iterations
        if initial_params === nothing
            start_params = Opts.init_vars!(tracker)
        elseif i == 1
            start_params = copy(initial_params)
        else
            start_params = perturb_params(initial_params, tracker, perturbation)
        end
        obj, params, ret = NLopt.optimize(opt, start_params)
        if getfield(NLopt, ret) < 0
            continue
        end
        if obj < best_obj
            best_obj = obj
            area = nlmodel(Val((:area, 0)), params)
            best_status = (
                obj = obj,
                dis = nlmodel(Val((:dis2, 0)), params),
                disδ = nlmodel(Val((:disδ2, 0)), params),
                area = area,
                areaε = abs(area) - π / 2,
                areaδ = nlmodel(Val((:areaδ, 0)), params),
                total_t = nlmodel(Val((:τ, 0)), params),
                Ωmax = params[nlmodel.param.Ωs[1]],
            )
            println(best_status)
            best_params = params
        end
        if best_obj < threshold
            break
        end
    end

    return best_params, best_obj
end

function get_metadata_and_plot(nlmodel, best_params, nseg, sysparams, modes;
                               ion1, ion2, pitime, τmin, τmax, maxtime, min_mode_index, max_mode_index)
    buf_plot = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_full), Val(SS.mask_full));
    kern = SL.Kernel(buf_plot, Val(SL.pmask_full));
    opt_raw_params = Seq.RawParams(nlmodel, best_params)

    fig, axes = subplots(2, 2, figsize=(8, 5.5))

    ts, Ωs = Seq.get_Ωs(opt_raw_params)
    axes[1].plot(ts, Ωs, label="Ω")
    axes[1].set_xlabel(raw"Time ($\mu s$)")
    axes[1].set_ylabel(raw"$\Omega$")
    axes[1].grid(true)
    axes[1].legend()

    plot_δs = range(-1, 1, 10001); # kHz
    axes[2].plot(plot_δs, [Seq.total_dis(kern, Seq.adjust(opt_raw_params, δ=2π * δ / 1000), modes) for δ in plot_δs])
    axes[2].set_xlim([-1, 1])
    axes[2].set_xlabel("Frequency offset (kHz)")
    axes[2].set_ylabel("Total Displacement")
    axes[2].grid(true)

    ts, ωs = Seq.get_ωs(opt_raw_params)
    axes[3].plot(ts, ωs ./ 2π)
    axes[3].set_xlabel(raw"Time ($\mu s$)")
    axes[3].set_ylabel(raw"$\omega$ (MHz)")
    for m in sysparams.modes.radial1
        axes[3].axhline(m, ls="--", color="red", alpha=0.4)
    end

    area0 = Seq.total_area(kern, opt_raw_params, modes)
    axes[4].plot(plot_δs, [Seq.total_area(kern, Seq.adjust(opt_raw_params, δ=2π * δ / 1000), modes) / area0 for δ in plot_δs])
    axes[4].set_xlim([-1, 1])
    axes[4].set_xlabel("Frequency offset (kHz)")
    axes[4].set_ylabel("Total Area")
    axes[4].grid(true)

    total_gate_time = best_params[nlmodel.param.τ] * nseg
    carrier_pi_time = π / maximum(Ωs) / 2
    total_dis = Seq.total_dis(kern, opt_raw_params, modes)
    total_cumdis = Seq.total_cumdis(kern, opt_raw_params, modes)
    total_disδ = Seq.total_disδ(kern, opt_raw_params, modes)
    total_areaδ = Seq.total_areaδ(kern, opt_raw_params, modes)
    dis_at_minus1kHz = Seq.total_dis(kern, Seq.adjust(opt_raw_params, δ=2π * -1 / 1000), modes)
    dis_at_plus1kHz = Seq.total_dis(kern, Seq.adjust(opt_raw_params, δ=2π * 1 / 1000), modes)

    fig.suptitle("Ions ($ion1, $ion2) | Gate time: $(round(total_gate_time, digits=1)) μs | π-time: $(round(carrier_pi_time, digits=1)) μs")
    info_text = join([@sprintf("%-36s %.4g", key, metadata_val) for (key, metadata_val) in [
        ("total_displacement", total_dis),
        ("total_cumulative_displacement", total_cumdis),
        ("gradient_displacement_detuning", total_disδ),
        ("enclosed_area", area0),
        ("gradient_area_detuning", total_areaδ),
    ]], "\n")
    fig.text(0.02, 0.01, info_text, verticalalignment="bottom")
    tight_layout(rect=[0, 0.18, 1, 0.95])

    metadata = Dict(
        "total_gate_time" => total_gate_time,
        "total_displacement" => total_dis,
        "total_cumulative_displacement" => total_cumdis,
        "gradient_displacement_detuning" => total_disδ,
        "displacement_at_-1kHz" => dis_at_minus1kHz,
        "displacement_at_+1kHz" => dis_at_plus1kHz,
        "enclosed_area" => area0,
        "gradient_area_detuning" => total_areaδ,
        "carrier_pi_time_required" => carrier_pi_time,
        "opt_params" => collect(best_params),
        "pitime" => pitime,
        "τmin" => τmin,
        "τmax" => τmax,
        "maxtime" => maxtime,
        "min_mode_index" => min_mode_index,
        "max_mode_index" => max_mode_index,
    )
    return opt_raw_params, metadata
end

function verify_solution(metadata)
    checks = [
        ("total_gate_time",                 v -> v <= 500,   "<= 500"),
        ("displacement_at_+1kHz",           v -> v <= 0.05,  "<= 0.05"),
        ("displacement_at_-1kHz",           v -> v <= 0.05,  "<= 0.05"),
        ("carrier_pi_time_required",        v -> v >= 2.7,   ">= 2.7"),
        ("gradient_displacement_detuning",  v -> v <= 0.002, "<= 0.002"),
        ("total_displacement",              v -> v <= 0.002, "<= 0.002"),
        ("total_cumulative_displacement",   v -> v <= 0.002, "<= 0.002"),
    ]
    failures = String[]
    for (key, test, msg) in checks
        val = metadata[key]
        if !test(val)
            push!(failures, "$key = $(@sprintf("%.4g", val)) ($msg)")
        end
    end
    if isempty(failures)
        println("All verification checks passed.")
    else
        for f in failures
            @warn f
        end
    end
    return failures
end

function save_am_solution(filename, opt_raw_params, metadata, sysparams, ion1, ion2)
    ion_key = join(sort!([ion1, ion2], rev=true), ",")

    new_xx_sol = GoldGates.XXSolution(
        opt_raw_params,
        metadata["enclosed_area"],
        metadata=JSON.json(metadata)
    )

    if isfile(filename)
        existing_file_content = open(filename, "r") do io
            JSON.parse(io)
        end
        xx_dict = get!(existing_file_content, "XX", Dict())
        xx_dict[ion_key] = new_xx_sol
        solution_data = existing_file_content
    else
        solution_data = Dict(
            "XX" => Dict(ion_key => new_xx_sol),
            "modes" => Dict(
                "radial1" => sysparams.modes.radial1,
                "axial" => [],
                "radial2" => []
            )
        )
    end

    open(filename, "w") do io
        JSON.print(io, solution_data, 2)
    end
    println("Saved solution for ions ($ion1, $ion2) to $filename")
end

function load_solution_metadata(filename, ion_pair_key)
    data = open(filename) do io
        JSON.parse(io)
    end
    xx = get(data, "XX", nothing)
    xx === nothing && return nothing
    sol = get(xx, ion_pair_key, nothing)
    sol === nothing && return nothing
    meta_str = get(sol, "metadata", nothing)
    meta_str === nothing && return nothing
    meta = JSON.parse(meta_str)
    if haskey(meta, "opt_params")
        meta["opt_params"] = Float64.(meta["opt_params"])
    end
    return meta
end

function plot_trajectories(opt_raw_params, modes)
    n = length(modes.modes)
    fig, axes = subplots(cld(n, 3), 3, figsize=(8, n))

    for i in 1:n
        _, xs, ys = Seq.get_trajectory(opt_raw_params, 10001; ωm=modes[i][1])
        ax = axes[i]
        ax.plot(xs, ys, lw=1)
        ax.set_title("Mode $i")
        ax.grid(true)
    end

    tight_layout()
    show()
end
