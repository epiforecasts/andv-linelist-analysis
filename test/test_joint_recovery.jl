## End-to-end simulation-based recovery for the full `joint_model`.
##
## These tests close the loop between data simulation and inference:
## fix the population-level parameters to known truth, simulate
## offspring counts `Zobs` via `rand` on the fixed `joint_model`, then
## re-fit the same model on the simulated data and assert that the 95%
## posterior credible intervals cover the truth.
##
## Two variants:
##
## - Retrospective: `joint_model(d, edges)` with `d.obs_time === nothing`.
##   `truncation_model` is a no-op; the per-case rate is just
##   `R_i = exp(log_R_at(T_onset[i], edges, log_R))`.
##
## - Real-time: `joint_model(d_rt, edges_rt)` with `d_rt.obs_time` set so
##   `truncation_model` fires. This exercises the incubation
##   right-truncation for index cases and the offspring-completeness
##   denominator for sourced cases, plus the per-case `R_eff = R_i · p_i`
##   thinning in `case_model`.
##
## Both tests re-use the bundled Epuyén line list for window structure
## (`onset_lo_day`, `exp_lo_day`, `source_idx`, etc.) and only re-simulate
## the offspring counts `Zobs` from the fixed truth. Latent infection
## and onset times are left as latent parameters in the refit, exactly
## as they are in the production fit.
##
## Scope of coverage assertions: only the parameters whose truth values
## actually propagate into the simulated `Zobs` can be recovered here.
## In `latent_times_model`, `T_onset` and `T_inf` are drawn from
## `Uniform` priors over the per-case windows; the incubation and δ
## log-densities enter via `@addlogprob!`, which does not influence the
## `rand`-time draws. The simulated `Z[i] ~ safe_nb(k, R_i)` therefore
## depends on `log_R_init`, `σ_rw`, `ε`, and `phi_inv_sqrt` only — those
## are the parameters we assert recovery on. Sim-recovery for the delay
## parameters `μ_inc, σ_inc, μ_δ, σ_δ` lives in `test_recovery.jl` /
## `test_submodel_recovery.jl`, which simulate from the marginal delay
## generative models directly. We still `fix` the delay truths during
## simulation so the refit is run on a coherent NamedTuple of truths.

using Random: Random, MersenneTwister
using Statistics: median, quantile
using Dates: Date, Day
using Turing: Turing, DynamicPPL

# Extract simulated offspring counts from the NamedTuple returned by
# `rand(MersenneTwister(seed), fixed_model)`. `case_model` declares
# `Z[i] ~ safe_nb(...)` so each sampled `Z[i]` appears under a VarName
# whose string form is `"Z[i]"`. We match by string to stay decoupled
# from AbstractPPL's internal VarName API.
function _extract_Zobs(sim, N::Integer)
    Z = Vector{Int}(undef, N)
    seen = falses(N)
    for (k, v) in pairs(sim)
        ks = string(k)
        if startswith(ks, "Z[") && endswith(ks, "]")
            idx = parse(Int, ks[3:(end - 1)])
            Z[idx] = Int(v)
            seen[idx] = true
        end
    end
    all(seen) ||
        error("rand() returned only $(count(seen))/$N Z values; " *
              "indices missing: $(findall(!, seen))")
    return Z
end

# Common truth for both retrospective and real-time tests. Modest σ_rw
# so the simulated outbreak doesn't blow up exponentially over the
# weekly knots; phi_inv_sqrt = 1/√4 ⇒ k = 4.
const _JOINT_TRUTH = (;
    μ_inc = 3.0, σ_inc = 0.45,
    μ_δ = 1.5, σ_δ = 1.2,
    log_R_init = log(1.4), σ_rw = 0.15,
    phi_inv_sqrt = 1.0 / sqrt(4.0))

# Per-test sampler budget. Kept modest so each testset fits in ~3 min
# on CI hardware. The 95% CrIs are intentionally wide enough that
# 500 × 4 NUTS draws give reliable coverage despite Monte-Carlo noise.
const _RECOVERY_SAMPLES = 500
const _RECOVERY_CHAINS = 4
const _RECOVERY_SEED = 20260519

# Scalar truth parameters whose values propagate into the simulated
# `Zobs` and so are recoverable from a refit on the simulated data.
# See the module-level note above for why the delay parameters are not
# included.
const _RECOVERY_PARAMS = (:log_R_init, :σ_rw, :phi_inv_sqrt)

# Build a data tuple with `Zobs` swapped for the simulated counts.
# `merge` preserves all the window structure (`onset_lo_day`, etc.)
# from the original line list, so only the observed offspring counts
# change relative to the production fit.
_with_Zobs(d, Z_sim) = merge(d, (; Zobs = Z_sim))

# Replace Zobs with a `missing` vector so `case_model` treats `Z[i]`
# as latent and `rand` simulates from the NB likelihood.
function _with_missing_Zobs(d)
    merge(d, (; Zobs = Vector{Union{Missing, Int}}(missing, d.N)))
end

@testset "joint_model: retrospective sim-then-recover" begin
    ll = TransmissionLinelist.load_linelist()
    t0 = minimum(ll.onset_date) - Day(60)
    d_obs = TransmissionLinelist.build_data(ll; t0 = t0)
    edges = TransmissionLinelist.prepare_rt_edges(t0)

    # Step 1: simulate Zobs from the fixed truth.
    sim_model = TransmissionLinelist.joint_model(
        _with_missing_Zobs(d_obs), edges)
    fixed = DynamicPPL.fix(sim_model, _JOINT_TRUTH)
    sim = rand(MersenneTwister(_RECOVERY_SEED), fixed)
    Z_sim = _extract_Zobs(sim, d_obs.N)

    # Step 2: refit on the simulated counts and check coverage.
    d_sim = _with_Zobs(d_obs, Z_sim)
    fit_model = TransmissionLinelist.joint_model(d_sim, edges)
    chn = TransmissionLinelist.sample_fit(fit_model;
        samples = _RECOVERY_SAMPLES,
        chains = _RECOVERY_CHAINS,
        seed = _RECOVERY_SEED,
        progress = false)

    for p in _RECOVERY_PARAMS
        draws = _chain_vec(chn, p)
        lo, hi = quantile(draws, [0.025, 0.975])
        truth_val = getfield(_JOINT_TRUTH, p)
        @test lo <= truth_val <= hi
    end
end

@testset "joint_model: real-time sim-then-recover" begin
    ll = TransmissionLinelist.load_linelist()
    obs_date = Date("2019-01-07")
    t0 = minimum(ll.onset_date) - Day(60)
    ll_rt = TransmissionLinelist.filter_realtime(ll, obs_date)
    d_obs = TransmissionLinelist.build_data(ll_rt;
        obs_time = obs_date, t0 = t0)
    edges = TransmissionLinelist.prepare_rt_edges(t0;
        obs_time = obs_date)

    # Step 1: simulate Zobs from the fixed truth, with `obs_time` set
    # so `truncation_model` is engaged. Note that `latent_times_model`
    # and `truncation_model` are exercised during simulation (latent
    # times are drawn from their priors; truncation adds log-prob
    # terms but does not affect the sampled values of `Z`).
    sim_model = TransmissionLinelist.joint_model(
        _with_missing_Zobs(d_obs), edges)
    fixed = DynamicPPL.fix(sim_model, _JOINT_TRUTH)
    sim = rand(MersenneTwister(_RECOVERY_SEED), fixed)
    Z_sim = _extract_Zobs(sim, d_obs.N)

    # Step 2: refit on the simulated counts. `obs_time` is preserved
    # through `merge`, so the real-time truncation fires during the
    # refit as well.
    d_sim = _with_Zobs(d_obs, Z_sim)
    fit_model = TransmissionLinelist.joint_model(d_sim, edges)
    chn = TransmissionLinelist.sample_fit(fit_model;
        samples = _RECOVERY_SAMPLES,
        chains = _RECOVERY_CHAINS,
        seed = _RECOVERY_SEED,
        progress = false)

    for p in _RECOVERY_PARAMS
        draws = _chain_vec(chn, p)
        lo, hi = quantile(draws, [0.025, 0.975])
        truth_val = getfield(_JOINT_TRUTH, p)
        @test lo <= truth_val <= hi
    end
end
