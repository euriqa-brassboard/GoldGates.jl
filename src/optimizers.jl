#

module Optimizers

import MSSim: Optimizers as Opts, SegSeq as SS, SymLinear as SL, Sequence as Seq, Utils as U

import ..set_mode_weight!, ..vv2m, ..m2vv, ..Candidate, ..ThreadObjectPool

using NLopt

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
        @lock pool.lock begin
            append!(candidates, o.candidates)
        end
        empty!(o.candidates)
        put!(pool, o)
    end
    return candidates
end

end
