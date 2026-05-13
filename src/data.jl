## Line-list loading, exposure/onset encoding, and R(t) bin definitions.

const LINELIST_PATH = joinpath(pkgdir(@__MODULE__), "data", "linelist.csv")
const OUTPUT_DIR    = joinpath(pkgdir(@__MODULE__), "output")
const FIGURES_DIR   = joinpath(pkgdir(@__MODULE__), "figures")

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
using Hantavirus
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
    ll.onset_date     = Date.(ll.onset_date)
    # Default onset_lower / onset_upper to onset_date if not present, allowing
    # the model to support multi-day onset uncertainty when the data has it.
    if !hasproperty(ll, :onset_lower); ll.onset_lower = copy(ll.onset_date); end
    if !hasproperty(ll, :onset_upper); ll.onset_upper = copy(ll.onset_date); end
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

Anchors all times in days relative to `t0 = minimum(onset_date) - 60 d` and
encodes interval-censored onset and exposure windows as `[lo, hi)` pairs.
Resolves the `source_case` column to integer indices into the line list
(`0` denotes a zoonotic index case with no human source) and reads
observed offspring counts from the `Z` column.

# Arguments
- `ll`: a `DataFrame` returned by [`load_linelist`](@ref).

# Returns
A named tuple `(t0, onset_lo_day, onset_hi_day, exp_lo_day, exp_hi_day,
source_idx, Zobs, N)` ready to pass to [`joint_model`](@ref).

# Examples
```julia
ll = load_linelist()
d  = build_data(ll)
d.N
```
"""
function build_data(ll)
    t0 = minimum(ll.onset_date) - Day(60)

    onset_lo_day = Float64.(Dates.value.(ll.onset_lower .- t0))
    onset_hi_day = Float64.(Dates.value.(ll.onset_upper .- t0)) .+ 1.0
    exp_lo_day = [ismissing(d) ? missing : Float64(Dates.value(d - t0))     for d in ll.exposure_lower]
    exp_hi_day = [ismissing(d) ? missing : Float64(Dates.value(d - t0)) + 1 for d in ll.exposure_upper]

    source_id  = passmissing(_parse_source).(ll.source_case)
    id_to_idx  = Dict(r.patient_id => i for (i, r) in enumerate(eachrow(ll)))
    source_idx = [ismissing(s) ? 0 : id_to_idx[string(s)] for s in source_id]

    # Zobs[i] is the observed offspring count of case i — the number of
    # secondaries in the line list attributed to i as their source. Read
    # directly from the Z column of the line list.
    return (; t0, onset_lo_day, onset_hi_day, exp_lo_day, exp_hi_day,
            source_idx, Zobs = Int.(ll.Z), N = nrow(ll))
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
[`joint_model`](@ref) into a single call so the analysis walkthrough
and CLI share the same model construction code path.

# Arguments
- `ll`: a line-list `DataFrame` as returned by [`load_linelist`](@ref).

# Returns
A 3-tuple `(model, d, edges)`: the Turing model, the augmented data
named tuple from [`build_data`](@ref), and the weekly knot edges from
[`bin_edges_day`](@ref).
"""
function prepare_model(ll)
    d     = build_data(ll)
    edges = bin_edges_day(d.t0)
    model = joint_model(d, edges)
    return (model, d, edges)
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
    t <= knots[1]   && return log_R[1]
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
