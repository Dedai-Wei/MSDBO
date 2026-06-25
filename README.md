# MSDBO-CEC2017

Minimal MATLAB code for the MSDBO optimizer and CEC2017 benchmark testing.

This repository keeps only the optimizer, the CEC2017 evaluation interface, and two simple running examples. The earlier KELM, bankruptcy prediction data, comparison optimizers, generated tables, generated figures, and long reproducibility notes are intentionally removed.

## Repository structure

```text
MSDBO-CEC2017/
├── src/
│   └── MSDBO.m
├── benchmarks/
│   ├── cec2017_eval.m
│   └── cec2017_basic.m
├── examples/
│   ├── run_single_function.m
│   └── run_cec2017_demo.m
├── .gitignore
└── README.md
```

## Requirements

- MATLAB R2020a or newer is recommended.
- No Python dependency is required.
- No MATLAB toolbox is required for the basic examples.

## Quick start

Clone the repository, open MATLAB, and run:

```matlab
cd examples
run_single_function
```

To run MSDBO on functions F1-F30:

```matlab
cd examples
run_cec2017_demo
```

The demo writes the following files to the repository root:

```text
msdbo_cec2017_results.csv
msdbo_cec2017_curves.mat
```

These generated files are ignored by Git.

## MSDBO usage

```matlab
addpath('src')

popSize = 30;
maxIter = 500;
dim = 30;
lb = -100;
ub = 100;
fobj = @(x) sum(x.^2);

[bestScore, bestPosition, curve] = MSDBO(popSize, maxIter, lb, ub, dim, fobj);
```

Input arguments:

| Name | Meaning |
|---|---|
| `popSize` | Population size. Use `popSize >= 20`. |
| `maxIter` | Maximum number of iterations. |
| `lb` | Lower bound. Scalar or `1 x dim` vector. |
| `ub` | Upper bound. Scalar or `1 x dim` vector. |
| `dim` | Dimension of the search space. |
| `fobj` | Objective function handle. The input is a row vector. |

Output arguments:

| Name | Meaning |
|---|---|
| `bestScore` | Best objective value found. |
| `bestPosition` | Best solution found. |
| `curve` | Best objective value at initialization and each iteration. |

## CEC2017 interface

Use `cec2017_eval` as the objective function:

```matlab
addpath('src')
addpath('benchmarks')

funcNum = 1;
fobj = @(x) cec2017_eval(x, funcNum);
[bestScore, bestPosition, curve] = MSDBO(30, 500, -100, 100, 30, fobj);
```

`cec2017_eval.m` first checks whether an official CEC2017 evaluator named `cec17_func` is available on the MATLAB path. If it is available, the wrapper calls it. If it is not available, the wrapper uses `cec2017_basic.m`, a lightweight built-in benchmark file.

For formal CEC2017 reporting, use the official shifted and rotated CEC2017 benchmark files. The built-in `cec2017_basic.m` file is included only to make the repository runnable without external files.

## Notes

- The code is intentionally small.
- The repository does not include KELM, bankruptcy data, other optimizers, or generated experiment results.
- The examples are for checking that MSDBO runs correctly on benchmark functions.
