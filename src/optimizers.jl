#

module Optimizers

import MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U

import ..set_mode_weight!, ..vv2m, ..m2vv, ..Candidate

using NLopt

using AMO.Utils: ThreadObjectPool, eachobj

using Base.Threads

export PreOptimizer, update_bounds!, opt_one!, set_mode!, set_time_range!, opt_all_rounds!

struct PreOptimizer{NSeg,PreObj,Sum}
    ωs::Vector{Float64}
    ωs2::Vector{Float64}
    modes::Seq.Modes
    modes2::Seq.Modes

    tmin::Float64
    tmax::Float64
    ntimes::Int

    pre_obj::PreObj
    pre_tracker::Opts.NLVarTracker
    pre_opt::NLopt.Opt

    sum::Sum
    args_buff::Vector{Float64}
    rawparams_buff::Vector{Float64}
    candidates::Vector{Candidate}

    function PreOptimizer{NSeg}(ωs, ωs2; tmin, tmax, ntimes=11, ωmin, ωmax,
                                amp_ratio=0.7, maxiter=2500,
                                disδ_weight=0.01) where NSeg
        nions = length(ωs)
        modes = Seq.Modes()
        for ω in ωs
            push!(modes, ω)
        end
        modes2 = Seq.Modes()
        for ω in ωs
            push!(modes2, ω)
        end
        for ω in ωs2
            push!(modes2, ω)
        end

        freq_spec = Seq.FreqSpec(true, sym=false)
        amp_spec = Seq.AmpSpec(cb=U.BlackmanStartEnd{amp_ratio}())

        pre_obj = Opts.abs_area_obj(NSeg, modes, SL.pmask_tfm,
                                    freq=freq_spec, amp=amp_spec,
                                    dis_weights=fill(1, nions),
                                    disδ_weights=fill(disδ_weight, nions),
                                    area_weights=zeros(nions))
        nargs = Seq.nparams(pre_obj)
        pre_tracker = Opts.NLVarTracker(nargs)
        Opts.set_bound!(pre_tracker, pre_obj.param.Ωs[1], 1, 1)
        for ω in pre_obj.param.ωs
            Opts.set_bound!(pre_tracker, ω, ωmin, ωmax)
        end

        pre_opt = NLopt.Opt(:LD_LBFGS, nargs)
        precompile(pre_obj, (Vector{Float64}, Vector{Float64}))
        NLopt.min_objective!(pre_opt, pre_obj)
        NLopt.maxeval!(pre_opt, maxiter)

        s = Seq.Summarizer{NSeg}()
        return new{NSeg,typeof(pre_obj),typeof(s)}(
            ωs, ωs2, modes, modes2, tmin, tmax, ntimes,
            pre_obj, pre_tracker, pre_opt,
            s, Vector{Float64}(undef, nargs), Vector{Float64}(undef, NSeg * 5),
            Candidate[]
        )
    end
end

function update_bounds!(o::PreOptimizer)
    NLopt.lower_bounds!(o.pre_opt, Opts.lower_bounds(o.pre_tracker))
    NLopt.upper_bounds!(o.pre_opt, Opts.upper_bounds(o.pre_tracker))
    return
end

function opt_one!(o::PreOptimizer)
    objval, args, ret = NLopt.optimize!(o.pre_opt,
                                        Opts.init_vars!(o.pre_tracker, o.args_buff))
    if getfield(NLopt, ret)::NLopt.Result < 0
        return false
    end
    raw_params = Seq.RawParams(o.pre_obj, args, buff=o.rawparams_buff)
    props = get(o.sum, raw_params, o.modes2)
    nions = length(o.ωs)

    dis = 0.0
    disδ = 0.0
    max_area = 0.0
    @inbounds @simd ivdep for i in nions
        dis += abs2(props.dis[i])
        disδ += abs2(props.disδ[i])
        max_area = max(max_area, abs(props.area[i]))
    end
    if dis < 1e-6 * nions && disδ < 1e-4 * nions && max_area >= 100
        push!(o.candidates, Candidate(copy(args), props))
        return true
    end
    return false
end

function set_mode!(o::PreOptimizer, mode_idx)
    nions = length(o.ωs)
    @assert 1 <= mode_idx <= nions
    area_weights = o.pre_obj.obj.area_weights
    for i in 1:nions
        area_weights[i] = i == mode_idx
    end
end

function set_time_range!(o::PreOptimizer, τmin, τmax)
    Opts.set_bound!(o.pre_tracker, o.pre_obj.param.τ, τmin, τmax)
    update_bounds!(o)
end

function opt_all_rounds!(pool::ThreadObjectPool{PreOpt},
                         nrounds, candidates) where {NSeg,PreOpt<:PreOptimizer{NSeg}}
    o0 = get(pool)
    τs = range(o0.tmin / NSeg, o0.tmax / NSeg, o0.ntimes + 1)
    nions = length(o0.ωs)
    ntimes = o0.ntimes
    put!(pool, o0)
    @threads :greedy for (mode_idx, time_idx) in Iterators.product(1:nions, 1:ntimes)
        o = get(pool)
        set_mode!(o, mode_idx)
        set_time_range!(o, τs[time_idx], τs[time_idx + 1])
        for _ in 1:nrounds
            opt_one!(o)
        end
        put!(pool, o)
    end
    for o in eachobj(pool)
        append!(candidates, o.candidates)
        empty!(o.candidates)
    end
    return candidates
end

struct PairOptimizer{NSeg,AreaObj,Sum}
    ωs::Vector{Float64}
    ηs::Vector{Float64}
    bij::Matrix{Float64}
    ωs2::Vector{Float64}
    modes::Seq.Modes

    Ωmax::Float64

    area_obj::AreaObj
    area_tracker::Opts.NLVarTracker
    area_opt::NLopt.Opt

    sum::Sum
    weights_buff::Vector{Float64}
    args_buff::Vector{Float64}
    rawparams_buff::Vector{Float64}

    function PairOptimizer{NSeg}(ωs, ηs, bij, ωs2;
                                 tmin, tmax, ntimes=11,
                                 ωmin, ωmax, δω=2π * 0.0005,
                                 Ωmax=0.4, amp_ratio=0.7,
                                 area_maxiter=10000,
                                 dis2_weight=0.002, dis3_weight=0.0005,
                                 disδ_weight=0.01, disδ2_weight=0,
                                 disδ3_weight=0.0000001) where NSeg
        nions = length(ωs)
        modes = Seq.Modes()
        for ω in ωs
            push!(modes, ω)
        end
        for ω in ωs
            push!(modes, ω - δω)
        end
        for ω in ωs
            push!(modes, ω + δω)
        end
        for ω in ωs2
            push!(modes, ω)
        end

        freq_spec = Seq.FreqSpec(true, sym=false)
        amp_spec = Seq.AmpSpec(cb=U.BlackmanStartEnd{amp_ratio}())

        dis_weights = [fill(1, nions); fill(dis2_weight, nions * 2);
                       fill(dis3_weight, nions)]
        disδ_weights = [fill(disδ_weight, nions); fill(disδ2_weight, nions * 2);
                         fill(disδ3_weight, nions)]
        area_targets = [Opts.AreaTarget(1, area_weights=zeros(nions),
                                        areaδ_weights=zeros(nions)),
                        Opts.AreaTarget(nions + 1, area_weights=zeros(nions)),
                        Opts.AreaTarget(nions * 2 + 1, area_weights=zeros(nions))]

        area_obj = Opts.target_area_obj(NSeg, modes, SL.pmask_full,
                                        freq=freq_spec, amp=amp_spec,
                                        dis_weights=dis_weights,
                                        disδ_weights=disδ_weights,
                                        area_targets=area_targets)
        nargs = Seq.nparams(area_obj)
        area_tracker = Opts.NLVarTracker(nargs)
        Opts.set_bound!(area_tracker, area_obj.param.Ωs[1], 0.01 * Ωmax, Ωmax)
        for ω in area_obj.param.ωs
            Opts.set_bound!(area_tracker, ω, ωmin, ωmax)
        end
        area_opt = NLopt.Opt(:LD_LBFGS, nargs)
        precompile(area_obj, (Vector{Float64}, Vector{Float64}))
        NLopt.min_objective!(area_opt, area_obj)
        NLopt.maxeval!(area_opt, area_maxiter)

        s = Seq.Summarizer{NSeg}()
        return new{NSeg,typeof(area_obj),typeof(s)}(
            ωs, ηs, bij, ωs2, modes, Ωmax,
            area_obj, area_tracker, area_opt,
            s, Vector{Float64}(undef, nions), Vector{Float64}(undef, nargs),
            Vector{Float64}(undef, NSeg * 5)
        )
    end
end

function update_bounds!(o::PairOptimizer)
    NLopt.lower_bounds!(o.area_opt, Opts.lower_bounds(o.area_tracker))
    NLopt.upper_bounds!(o.area_opt, Opts.upper_bounds(o.area_tracker))
    return
end

get_weights!(o::PairOptimizer, ion1, ion2) =
    set_mode_weight!(o.weights_buff, o.ηs, o.bij, ion1, ion2)

function opt_pair!(o::PairOptimizer, args, area, weights;
                   τ_tol=0.05, area_weight=10, area2_weight=0.001, areaδ_weight=0.001)
    args[o.area_obj.param.Ωs[1]] .*= sqrt((π / 2) / area)

    τ0 = args[o.area_obj.param.τ]
    Opts.set_bound!(o.area_tracker, o.area_obj.param.τ,
                    (1 - τ_tol) * τ0, (1 + τ_tol) * τ0)
    update_bounds!(o)

    area_tgt1 = o.area_obj.area_targets[1]
    area_tgt2 = o.area_obj.area_targets[2]
    area_tgt3 = o.area_obj.area_targets[3]

    area_tgt1.target = π / 2 * area_weight
    area_tgt1.area_weights .= weights .* area_weight
    area_tgt1.areaδ_weights .= weights .* (area_weight * areaδ_weight)
    area_tgt3.target = π / 2 * area2_weight
    area_tgt3.area_weights .= weights .* area2_weight
    area_tgt2.target = π / 2 * area2_weight
    area_tgt2.area_weights .= weights .* area2_weight

    objval, args, ret = NLopt.optimize!(o.area_opt, args)
    if getfield(NLopt, ret)::NLopt.Result < 0
        return
    end
    raw_params = Seq.RawParams(o.area_obj, args, buff=o.rawparams_buff)
    props = get(o.sum, raw_params, o.modes)
    return Candidate(args, props)
end

end
