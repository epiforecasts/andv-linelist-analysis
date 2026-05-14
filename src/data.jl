## Line-list loading, exposure/onset encoding, and R(t) bin definitions.

const LINELIST_PATH = joinpath(pkgdir(@__MODULE__), "data", "linelist.csv")
const OUTPUT_DIR = joinpath(pkgdir(@__MODULE__), "output")
const FIGURES_DIR = joinpath(pkgdir(@__MODULE__), "figures")

# Weekly R(t) bin edges spanning the outbreak.
const BIN_EDGES = collect(Date("2018-11-12"):Day(7):Date("2019-02-04"))

"""
$(TYPEDSIGNATURES)

Load and clean the Epuyén line list from a CSV file.

Reads the CSV at `path`, drops any duplicated rows with patient IDs ending
in `_alt`, parses `exposure_lower`, `exposure_upper`, `onset_date`,
`onset_lower`, and `onset_upper` as `Date`s (defaulting `onset_lower` and
`onset_upper` to `onset_date` when absent), and sorts the rows by
integer-valued `patient_id`. Used by [`build_data`](@ref) to produce the
model input tuple.

# Arguments
- `path`: path to a line-list CSV. Defaults to the bundled
  `data/linelist.csv` shipped with the package.

# Returns
A `DataFrame` with one row per case and parsed date columns ready for
[`build_data`](@ref).

# Examples
```julia
using TransmissionLinelist
ll = load_linelist()
first(ll, 3)
```
"""
function load_linelist(path = LINELIST_PATH)
    ll = CSV.read(path, DataFrame; missingstring = ["NA"],
        types = Dict(:patient_id => String))
    ll = filter(r -> !endswith(r.patient_id, "_alt"), ll)
    ll.exposure_lower = passmissing(Date).(ll.exposure_lower)
    ll.exposure_upper = passmissing(Date).(ll.exposure_upper)
    ll.onset_date = Date.(ll.onset_date)
    # Default onset_lower / onset_upper to onset_date if not present, allowing
    # the model to support multi-day onset uncertainty when the data has it.
    if !hasproperty(ll, :onset_lower)
        ;
        ll.onset_lower = copy(ll.onset_date);
    end
    if !hasproperty(ll, :onset_upper)
        ;
        ll.onset_upper = copy(ll.onset_date);
    end
    ll.onset_lower = Date.(ll.onset_lower)
    ll.onset_upper = Date.(ll.onset_upper)
    sort!(ll, :patient_id, by = x -> parse(Int, x))
    return ll
end

# Source attributions: NA / "index" → no source; bare ID → that patient;
# "i/j/..." → multi-source, take the first per the paper's tree.
function _parse_source(s)
    (ismissing(s) || s in ("NA", "index")) && return missing
    return parse(Int, occursin("/", s) ? split(s, "/")[1] : s)
end

"""
$(TYPEDSIGNATURES)

Build the model input from a line-list `DataFrame`.

Anchors all times in days relative to `t0` (defaults to
`minimum(onset_date) - 60 d`) and encodes interval-censored onset and
exposure windows as `[lo, hi)` pairs. Resolves the `source_case` column
to integer indices into the line list (`0` denotes a zoonotic index case
with no human source) and reads observed offspring counts from the `Z`
column. When `obs_time` is supplied, encodes the per-case real-time
cut-off in days relative to `t0`.

# Arguments
- `ll`: a `DataFrame` returned by [`load_linelist`](@ref).

# Keyword Arguments
- `obs_time`: optional real-time observation cut-off. Either a single
  `Date` (broadcast to all cases) or an `AbstractVector` of `Date`s with
  one entry per case. `nothing` (the default) collapses the model to its
  retrospective form.
- `t0`: optional explicit time origin. Defaults to
  `minimum(onset_date) - 60 d`.

# Returns
A named tuple `(t0, onset_lo_day, onset_hi_day, exp_lo_day, exp_hi_day,
source_idx, Zobs, N, obs_time)` ready to pass to [`joint_model_def`](@ref).

# Examples
```julia
ll = load_linelist()
d  = build_data(ll)
d.N
```
"""
function build_data(
        ll; obs_time::Union{Nothing, Date, AbstractVector} = nothing,
        t0::Union{Nothing, Date} = nothing)
    t0 = t0 === nothing ? minimum(ll.onset_date) - Day(60) : t0

    onset_lo_day = Float64.(Dates.value.(ll.onset_lower .- t0))
    onset_hi_day = Float64.(Dates.value.(ll.onset_upper .- t0)) .+ 1.0
    exp_lo_day = [ismissing(d) ? missing : Float64(Dates.value(d - t0))
                  for d in ll.exposure_lower]
    exp_hi_day = [ismissing(d) ? missing : Float64(Dates.value(d - t0)) + 1
                  for d in ll.exposure_upper]

    source_id = passmissing(_parse_source).(ll.source_case)
    id_to_idx = Dict(r.patient_id => i for (i, r) in enumerate(eachrow(ll)))
    source_idx = [ismissing(s) ? 0 : id_to_idx[string(s)] for s in source_id]

    obs_time_day = _encode_obs_time(obs_time, t0, nrow(ll), exp_hi_day)

    # Zobs[i] is the observed offspring count of case i — the number of
    # secondaries in the line list attributed to i as their source. Read
    # directly from the Z column of the line list.
    return (; t0, onset_lo_day, onset_hi_day, exp_lo_day, exp_hi_day,
        source_idx, Zobs = Int.(ll.Z), N = nrow(ll),
        obs_time = obs_time_day)
end

# Encode the optional per-case observation cut-off into days since `t0`.
# Returns nothing for a missing obs_time so the model collapses to the
# retrospective form. A scalar Date broadcasts to all cases.
function _encode_obs_time(::Nothing, t0, N, exp_hi_day)
    return nothing
end
function _encode_obs_time(d::Date, t0, N, exp_hi_day)
    return _encode_obs_time(fill(d, N), t0, N, exp_hi_day)
end
function _encode_obs_time(dates::AbstractVector, t0, N, exp_hi_day)
    length(dates) == N || error("obs_time length $(length(dates)) ≠ N=$N")
    out = Vector{Float64}(undef, N)
    for i in 1:N
        d = dates[i]
        ismissing(d) &&
            error("obs_time[$i] is missing; supply a Date for every case")
        out[i] = Float64(Dates.value(d - t0))
        if !ismissing(exp_hi_day[i]) && out[i] < exp_hi_day[i]
            error("obs_time[$i] ($(d)) precedes the upper exposure bound for that case")
        end
    end
    return out
end

"""
$(TYPEDSIGNATURES)

Return the weekly R(t) knot dates expressed as days relative to `t0`.

The knots span the outbreak in weekly steps; combined with [`log_R_at`](@ref)
this defines the piecewise-linear log R(t) trajectory used by
[`joint_model`](@ref).

# Arguments
- `t0`: the model's time origin (the `t0` field of the tuple returned by
  [`build_data`](@ref)).

# Returns
A `Vector{Float64}` of length `length(BIN_EDGES)` giving the knot positions
in days.
"""
bin_edges_day(t0) = Float64[Dates.value(d - t0) for d in BIN_EDGES]

"""
$(TYPEDSIGNATURES)

Build the joint model from a line-list `ll`.

Wraps [`build_data`](@ref), [`bin_edges_day`](@ref), and
[`joint_model_def`](@ref) into a single call so the analysis walkthrough
and CLI share the same model construction code path.

# Arguments
- `ll`: a line-list `DataFrame` as returned by [`load_linelist`](@ref).

# Returns
A NamedTuple `(; model, d, edges)`: the Turing model, the augmented data
named tuple from [`build_data`](@ref), and the weekly knot edges from
[`bin_edges_day`](@ref).
"""
function joint_model(ll)
    d = build_data(ll)
    edges = bin_edges_day(d.t0)
    return (; model = joint_model_def(d, edges), d, edges)
end

# Piecewise-linear interpolation: log_R[b] is the value at knot b, with
# linear interpolation inside the knot range and clamping outside.
"""
$(TYPEDSIGNATURES)

Piecewise-linear interpolation of `log R(t)` between weekly knots.

Linearly interpolates `log_R` against `knots` at the time `t`, clamping to
the endpoint values outside the knot range.

# Arguments
- `t`: time (in days from `t0`) at which to evaluate log R.
- `knots`: knot positions in days, as returned by [`bin_edges_day`](@ref).
- `log_R`: vector of log R values at each knot.

# Returns
The interpolated log R value at `t`.
"""
function log_R_at(t::Real, knots::AbstractVector{<:Real}, log_R)
    t <= knots[1] && return log_R[1]
    t >= knots[end] && return log_R[end]
    b = searchsortedlast(knots, t)
    w = (t - knots[b]) / (knots[b + 1] - knots[b])
    return (1 - w) * log_R[b] + w * log_R[b + 1]
end

"""
$(TYPEDSIGNATURES)

Return string labels for the weekly R(t) knots.

One entry per `log_R` element. Used to label plots and posterior summaries
produced by [`summarise`](@ref) and [`plot_rt`](@ref).

# Returns
A `Vector{String}` of ISO-format knot dates.
"""
bin_labels() = string.(BIN_EDGES)
