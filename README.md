# Fast and forward stable randomized algorithms for linear least-squares problems

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=eepperly/Randomized-Least-Squares-Solvers)

This repository contains code for the paper _Fast and forward stable randomized algorithms for linear least-squares problems_ by [Ethan N. Epperly](https://www.ethanepperly.com).

## Experiments 

Code to reproduce the numerical experiments from the paper are found in `experiments/`. 

- Figure 1 (`experiments/test_iterative_sketching.m`): Computes the forward, residual, and backward errors for iterative sketching for different condition numbers and residual norms.
- Figure 2 (`experiments/test_variants.m`): Compares the performance iterative sketching (including versions with damping and momentum) to sketch-and-precondition (with both the zero and sketch-and-solve initializations).
- Figure 3 (`experiments/test_bad_iterative_sketching.m`): Compare the stable implementation of iterative sketching (in `code/iterative_sketching.m`) to three "bad" implementations.
- Figure 4 (`experiments/test_timing.m`): Compares the runtime of iterative sketching (including versions with damping and momentum) to MATLAB's `mldivide` on a simplified kernel regression task.
