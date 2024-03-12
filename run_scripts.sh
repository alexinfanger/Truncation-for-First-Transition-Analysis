#!/bin/bash
SCRIPTS=("scripts/MM1_discounted_rewards.jl" \
         "scripts/TandemQueue_hitting_time.jl" 
         )

for script in "${SCRIPTS[@]}"
do
  julia --project=. -e "using Pkg; Pkg.instantiate(); include(\"scripts/startup_fta.jl\"); include(\"$script\")"
done

echo "All scripts have been run."