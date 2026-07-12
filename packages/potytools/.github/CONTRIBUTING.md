# Contributing to potytools

This outlines how to propose a change to potytools.

## Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.
Edit the roxygen comment in the `.R` file below `R/`, not the generated `.Rd`
file below `man/`.

## Bigger changes

If you want to make a bigger change, it's a good idea to first
[open an issue](https://github.com/MintyHybrid/project-hail-mary/issues) and
file it before you make the change, so we can discuss whether it's needed.

### Pull request process

1. Fork the package and clone onto your computer.
2. Create a Git branch for your pull request (PR).
3. Make your changes, and update `NEWS.md` with a bullet describing the change.
   Add your GitHub username, and the PR number.
4. Documentation is generated with [roxygen2](https://roxygen2.r-lib.org):
   edit the roxygen comments in the `.R` files and run
   `roxygen2::roxygenise()` (or `devtools::document()`) to regenerate
   `NAMESPACE` and `man/`.
5. Add tests under `tests/testthat/`. Network-dependent tests should be recorded
   with [vcr](https://docs.ropensci.org/vcr/) so they run offline.
6. Run `devtools::check()` (or `R CMD check`) and ensure it passes cleanly.
7. Open the PR. The PR description should describe what the change does and
   reference the issue it addresses.

### Code style

- New code should follow the tidyverse [style guide](https://style.tidyverse.org).
  You can use [styler](https://styler.r-lib.org) to apply styling and
  [lintr](https://lintr.r-lib.org) to check it.
- Use `snake_case` for function and argument names.

## Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
