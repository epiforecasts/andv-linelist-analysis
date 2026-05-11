## Run with:  julia --project=. -t auto scripts/run.jl

using Hantavirus

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
