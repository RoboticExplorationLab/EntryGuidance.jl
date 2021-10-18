# EntryGuidance.jl

[![Build Status](https://travis-ci.com/RoboticExplorationLab/EntryGuidance.jl.svg?token=SqgAVz1CAUik7HkKWcXq&branch=master)](https://travis-ci.com/RoboticExplorationLab/EntryGuidance.jl)

This repository accompanies the following paper, [CPEG: A Predictor-corrector Entry Guidance Algorithm](https://github.com/RoboticExplorationLab/EntryGuidance.jl/blob/master/cpeg_paper.pdf), submitted to the [2022 IEEE Aerospace Conference](https://aeroconf.org/).

<!-- The examples from the paper can be run in the following manner: -->
### Running the examples from the paper

1. Download and start Julia https://julialang.org/downloads/ (1.6.3)
2. Navigate to the `EntryGuidance.jl` directory. 
3. run `using Pkg; Pkg.activate("."); Pkg.instantiate();`

From here, the following three examples can be run:
- `include("cpeg_examples/bank_angle/run_L1_example.jl")`
- `include("cpeg_examples/bank_angle/run_quad_example.jl")`
- `include("cpeg_examples/full_lift_control/run_example.jl")`

### Dependencies
All dependencies will be automatically configured with `Pkg.instantiate()`, except for the following:
- MATLAB (for plotting)
- Mosek (conic solver)
