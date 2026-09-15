# GEMSTAT

This code implements the method from the paper

- Thermodynamics-based models of transcriptional regulation by enhancers: the roles of synergistic activation, cooperative binding and short-range repression
- Author: Xin He <xinhe2@illinois.edu>

## Authors

See the file [AUTHORS](./AUTHORS) for information about additional authors and maintainers of this software.

## Documentation

See the file named [README.txt](./README.txt) for documentation about building and running the program.

## building

After checking out the repo, do:

```bash
git submodule init
git submodule update --recursive
mkdir build
cd build
cmake ..
make
```

## Testing

The build has a golden-output regression suite driven by ctest. Every
directory under `tests/cases/` is one case: the arguments to run
`seq2expr` with and the expected outputs. Cases cover every model,
objective and interaction option, the error paths, short training runs
and a check of the automatic-differentiation gradient against central
differences.

```bash
cd build
ctest -j4                  # everything (about a minute)
ctest -j4 -L fast          # skip the training runs
ctest -R gradcheck         # only the gradient checks
```

After an intentional change of behaviour, regenerate the expected files
with `make golden-update` and review the diff under `tests/cases/*/expected/`
before committing. See `tests/run_case.py` for the case format.

## Options added since the paper

- `--threads N`: predict sequences in parallel (OpenMP; default all cores).
- `--gradient ad|fd`: exact reverse-mode automatic-differentiation gradient
  (default) or the old forward-difference gradient.
- `--check_gradient`: compare the two at the initial parameters and exit.
- `--random_starts N`: N random restarts in one process.
- `--seed S`: reproducible runs.

See `README.txt` for all options.

### Maintainers

See the file named [MAINTAINER_README](./MAINTAINER_README) for more hints on how to maintain this software.

### Builds

Continuous integration runs on GitHub Actions (`.github/workflows/ci.yml`):
Release and Debug builds, then the full test suite.
