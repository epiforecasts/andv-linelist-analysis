Baseline regression-summary CSVs lived alongside the PR that introduced the feature.

# regression-baseline

Per-walkthrough CSVs produced by `TransmissionLinelist.save_regression_summary`,
treated as the **expected** posterior summary for each docs walkthrough.

The Documenter CI step `scripts/regression_diff.py` compares the freshly
built `output/regression/*.csv` against the matching baseline here and
posts a sticky PR comment with any drift outside the tolerances.

## Re-baselining

When a PR intentionally changes the model (prior, likelihood, refactor
with substantive effects) the baseline needs to move with it. Workflow:

1. Land the model change on a branch.
2. Locally run the docs build (`julia --project=docs docs/make.jl`),
   which writes the new CSVs to `output/regression/`.
3. Copy them into `regression-baseline/` and commit alongside the
   model change. Note in the PR description **why** the baseline moved
   (e.g. "tightened σ_R prior; lifts late-knot R(t) ≈ +0.6σ").
4. Reviewers can sanity-check the diff against the previous baseline
   in the same PR (the comment will compare PR build vs the
   pre-update baseline file as long as you don't update it until last).

A baseline file is **just a CSV** — no schema beyond the four columns
`section, quantity, median, lower_95, upper_95`. Walkthroughs that
don't yet emit a regression summary simply have no baseline file and
the CI step skips them.
