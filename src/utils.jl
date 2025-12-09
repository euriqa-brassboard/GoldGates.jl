using PyPlot
using MSSim: Sequence as Seq 

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

function get_metadata_and_plot(nlmodel, best_params;)
    buf_plot = SL.ComputeBuffer{nseg,Float64}(Val(SS.mask_full), Val(SS.mask_full));
    kern = SL.Kernel(buf_plot, Val(SL.pmask_full));
    opt_raw_params = Seq.RawParams(nlmodel, best_params)

    fig, axes = subplots(2, 2, figsize=(8, 4))

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
    tight_layout()

    total_gate_time = best_params[nlmodel.param.τ] * nseg
    total_dis = Seq.total_dis(kern, opt_raw_params, modes)
    total_cumdis = Seq.total_cumdis(kern, opt_raw_params, modes)
    total_disδ = Seq.total_disδ(kern, opt_raw_params, modes)
    total_areaδ = Seq.total_areaδ(kern, opt_raw_params, modes)
    metadata = Dict(
        "total_gate_time" => total_gate_time,
        "total_displacement" => total_dis,
        "total_cumulative_displacement" => total_cumdis,
        "gradient_displacement_detuning" => total_disδ,
        "enclosed_area" => area0,
        "gradient_area_detuning" => total_areaδ,
        "carrier_pi_time_required" => π/maximum(Ωs)/2,
    )
    return opt_raw_params, metadata
end 