## Line-list loading, exposure/onset encoding, and R(t) bin definitions.

using CSV, DataFrames, Dates

const PROJECT_ROOT  = normpath(joinpath(@__DIR__, ".."))
const LINELIST_PATH = joinpath(PROJECT_ROOT, "data", "linelist.csv")
const OUTPUT_DIR    = joinpath(PROJECT_ROOT, "output")

# Weekly R(t) bin edges spanning the outbreak.
const BIN_EDGES = collect(Date("2018-11-12"):Day(7):Date("2019-02-04"))

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

bin_edges_day(t0) = Float64[Dates.value(d - t0) for d in BIN_EDGES]

function which_bin(t::Real, edges::Vector{Float64})
    for (b, e) in enumerate(edges)
        t < e && return b
    end
    return length(edges) + 1
end

# Pretty labels for the R(t) bins, used in summaries and posterior output.
function bin_labels()
    labels = String[string("≤ ", BIN_EDGES[1])]
    for i in 2:length(BIN_EDGES)
        push!(labels, string(BIN_EDGES[i-1], " – ", BIN_EDGES[i]))
    end
    push!(labels, string("> ", BIN_EDGES[end]))
    return labels
end
