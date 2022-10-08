# Real Time Simulations of Quantum Spin Chains: Density-of-States and Reweighting approaches
### by Pavel Buividovich and Johann Ostmeyer

This repository contains the scripts and data required to reproduce the results presented in [arXiv:2209.13970 [cond-mat.stat-mech]](https://arxiv.org/abs/2209.13970).

For questions concerning the code contact [J.Ostmeyer@liverpool.ac.uk](mailto:J.Ostmeyer@liverpool.ac.uk).

## Simulations
The `simulation` directory contains all the important scripts and data. The Monte Carlo simulations are implemented with an `R` front-end allowing a simple usage as well as quick testing and a `C` back-end yielding high performance. The code is parallelised for single-node multi-core usage with R's native [mclapply](https://www.rdocumentation.org/packages/parallel/versions/3.4.0/topics/mclapply) and [OpenMP](https://www.openmp.org/) in the C-part.

### Compile
You might have to update the path to your `R/include` library in the `Makefile`. Find the path by executing `Sys.getenv("R_INCLUDE_DIR")` in your R environment.

After adjusting the Makefile according to your needs (possibly switching to another compiler), simply type `make` to compile all the `C` files.

### Run Code
Individual functions can be executed by running `source("ising_dos.R")` in an `R` environment and then simply using the pre-defined functions.

The data used in the paper can be reproduced with the help of the `llr_vs_reweighting.R` script (e.g. run with `Rscript llr_vs_reweighting.R` in the console). Parameters like length `L` might have to be adjusted.

Depending on the computer at hand it might be advantageous to submit the program as `slurm` jobs. See the scripts `submit.sh` and `job_script.sh` for this.

## Data

Most of the data is stored in the `simulation/R.Data` directory using the `.RData` file format that allows to `load()` the complete R objects/variables into the current environment.

The randomly chosen field configurations are stored (as plain text) in the respective subdirectory of `production`. In the case of length `L=16` we additionally provide the eigenvalues. We did not add the `RData` files for `L=6` because they are rather large and can be reproduced quickly.

In case of interest we can provide the `L=6` data as well as more results from exact diagonalisation.

## Plots
Plots are performed using `R` as well. All the plots provided in the paper have been generated with the `ising_plots.R` script and different settings for the length `L`.

## Mathematica notebooks
`partition-sums.nb` contains the error estimations leading to the scaling results of Appendix B.

`transverse-ising.nb` implements the Hamiltonian and some simple routines for its exact diagonalisation as well as a simple test example.

## License
Code and data follow the `mit` license.
