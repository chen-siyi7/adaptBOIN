# Contributing to adaptBOIN

Bug reports and focused feature requests are welcome through the GitHub issue
tracker. Please include a minimal reproducible example and the output of
`sessionInfo()`.

For code contributions:

1. Create a branch from `main`.
2. Add or update tests under `tests/testthat/`.
3. Run `R CMD check --no-manual` locally.
4. Open a pull request describing the change and its tests.

Please do not commit compiled files from `src/` or generated simulation output.
