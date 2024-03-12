### Code for "A Posteriori Error Bounds for Truncated Markov Chain Linear Systems Arising from First Transition Analysis"

This repository contains the code to reproduce the plots in the paper "A Posteriori Error Bounds for Truncated Markov Chain Linear Systems Arising from First Transition Analysis" by A. Infanger and P. W. Glynn.

We ran this code on Julia v1.10.2. The versions of the packages we used can be found
in the project.toml and manifest.toml files, and the correctly versioned packages can be instantiated using Julia's package manager (see code below).

To reproduce our results, clone the repository and then run the following code to instantiate the correct packages and load the necessary startup files:

```
> julia --project="." 

julia> using Pkg
julia> Pkg.instantiate()
julia> include("scripts/startup_fta.jl")
```

Then you can run the script for the $M/M/1$ queue discounted rewards example with:

```
include("scripts/MM1_discounted_rewards.jl")
```

And you can run the script for the Tandem Queue hitting time example with:

```
include("scripts/TandemQueue_hitting_time.jl")
```
Alternatively, you can run the following bash script to reproduce the plots for both examples (they will be saved in the "assets" folder):

```
> ./run_scripts.sh
```
