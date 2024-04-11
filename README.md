# Fast and forward stable randomized algorithms for linear least-squares problems

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=eepperly/Randomized-Least-Squares-Solvers)

This repository contains code for the paper [_Fast and forward stable randomized algorithms for linear least-squares problems_](https://arxiv.org/abs/2311.04362) by [Ethan N. Epperly](https://www.ethanepperly.com).

## Experiments 

Code to reproduce the numerical experiments from the paper are found in `experiments/`. 

- Figure 1 (`experiments/test_iterative_sketching.m`): Computes the forward, residual, and backward errors for iterative sketching for different condition numbers and residual norms.
- Figure 2 (`experiments/test_variants.m`): Compares the performance iterative sketching (including versions with damping and momentum) to sketch-and-precondition (with both the zero and sketch-and-solve initializations).
- Figure 3 (`experiments/test_bad_iterative_sketching.m`): Compare the stable implementation of iterative sketching (in `code/iterative_sketching.m`) to three "bad" implementations.
- Figure 4 (`experiments/test_timing.m`): Compares the runtime of iterative sketching (including versions with damping and momentum) to MATLAB's `mldivide` on a simplified kernel regression task.
- Figure 5 (`experiments/test_sparse.m`): Compares the runtime of iterative sketching to MATLAB's `mldivide` for solving random sparse least-squares problems with three nonzeros per row.

## Randomized least-squares solvers

This repository contains code for [iterative sketching](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS4) and [sketch-and-precondition](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS5).
Code for these methods is found in the `code/` directory.

### Sparse sign embeddings

The core ingredient for both iterative sketching and sketch-and-precondition is a [_fast random embedding_](https://ar5iv.labs.arxiv.org/html/2002.01387#S9).
Based on the empirical comparison in [this paper (see Fig. 2)](https://arxiv.org/abs/2104.05877) and our own testing, our preferred embedding is the [_sparse sign embedding_](https://ar5iv.labs.arxiv.org/html/2002.01387#S9.SS2).

To generate a sparse sign embedding in our code, first build the [mex file](https://www.mathworks.com/help/matlab/ref/mex.html) using the following command:

```
mex sparsesign.c
```

Then, to generate a `d` by `m` sparse sign embedding with `zeta` nonzeros per column, run

```
S = sparsesign(d, m, zeta);
```

### Iterative sketching

_Iterative sketching_ is a randomized iterative method for solving $Ax = b$ in the least-squares sense.
The first step of the algorithm is to collect a sketch $SA$ of the matrix $A$ and compute its [QR factorization](https://en.m.wikipedia.org/wiki/QR_decomposition) $SA = QR$.
We use a sparse sign embedding for the embedding matrix $S$.
After which, iterative sketching produces a sequence of better-and-better least-squares solutions using the iteration

$$
x_{i+1} = x_i + R^{-1} ( R^{-\top} (A^\top(b - Ax_i))).
$$

The main result of [my paper](https://arxiv.org/abs/2311.04362) is that iterative sketching is [(forward) stable](https://nhigham.com/2020/08/04/what-is-numerical-stability/): If you run it for sufficiently many iterations, the (forward) error $\lVert x - x_i \rVert$ and residual error $\lVert (b-Ax) - (b-Ax_i) \rVert$ are roughly as small as for a standard direct method like (Householder) QR factorization.

Iterative sketching can optionally be accelerated by incorporating _damping_ and _momentum_, resulting in a modified iteration

$$
x_{i+1} = x_i + \alpha R^{-1} ( R^{-\top} (A^\top(b - Ax_i))) + \beta (x_i - x_{i-1}).
$$

We call $\alpha$ and $\beta$ the _damping parameter_ and _momentum parameter_ respectively.
The optimal damping and momentum parameters were computed in [these](https://web.stanford.edu/~pilanci/papers/IHSMomentum18.pdf) [papers](https://arxiv.org/abs/1911.02675).

To run iterative sketching using our code, the command is

```
[x, stats] = iterative_sketching(A, b, [d, q, summary, verbose, damping, momentum])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: see paper)
- `q`: number of iterations (if `q>=1`) or tolerance (if `q<1`). If `q>=1`, iterative sketching will be run for `q` iterations. Otherwise, iterative sketching is run for an adaptive number of steps until the norm change in residual is less than `q*(Anorm * norm(x) + 0.01*Acond*norm(r))`. Here, `Anorm` and `Acond` are estimates of the norm and condition number of `A`. (_Default value_: `eps`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `damping`: damping coefficient. To use the optimal value, set `damping` to `'optimal'`. (_Default value_: 1, i.e., no damping).
- `momentum`: momentum coefficient. To use the optimal value, set both `damping` and `momentum` to `'optimal'`. (_Default value_: 0, i.e., no momentum).
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.

### Sketch-and-precondition

While it is not the main focus of our paper, we also provide an implementation of the sketch-and-precondition method.

Sketch-and-precondition is also based on a QR factorization $SA = QR$ of a sketch of the matrix $A$.
It then uses $R$ as a preconditioner for solving $Ax = b$ using the [LSQR](https://stanford.edu/group/SOL/software/lsqr/lsqr-toms82a.pdf) or [CGNE](https://en.wikipedia.org/wiki/Conjugate_gradient_method#Conjugate_gradient_on_the_normal_equations) method.
To call sketch-and-precondition, use the following command

```
[x, stats] = sketch_and_precondition(A, b, [d, q, summary, verbose, opts])
```

All but the first two inputs are optional.
The optional inputs are as follows:

- `d`: sketching dimension. (_Default value_: `2*size(A,2)`)
- `q`: number of iterations. (_Default value_: `100`)
- `summary`: a function of `x` to be recorded at every iteration. The results will be outputted in the optional output argument `stats`. (_Default value_: None)
- `verbose`: if true, output at every iteration. (_Default value_: `false`)
- `opts`: specifies the initial iterate $x_0$ and iterative method (LSQR or CGNE). If `'cold'` is a substring of `opts`, then the initial iterate is chosen to be $x_0 = 0$. Otherwise, we use a warm start and choose $x_0$ to be the [sketch-and-solve solution](https://ar5iv.labs.arxiv.org/html/2002.01387#S10.SS3). If `cgne` is a substring of `opts`, then we solve $Ax = b$ using CGNE; otherwise, we use LSQR. (_Default value_: `''`)
- `reproducible`: if true, use a slow, but reproducible implementation of sparse sign embeddings. (_Default value_: `false`)

Inputting a value of `[]` for an optional argument results in the default value.
