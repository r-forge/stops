## Test environments
* local Linux Mint 19.2 install, R 4.2.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check())
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was the following NOTE:

Possibly misspelled words in DESCRIPTION:
  Isomap (8:668)
  Jaakola (8:933)
  Luus (8:928)
  MDS (8:167, 8:369, 8:540, 8:609, 8:649, 8:660)
  SMACOF (8:437, 8:488, 8:506)
  Sammon (8:445, 8:574)
  Torgerson (8:406)
  hyperparameters (8:736, 8:774)
  structuredness (8:1076)

I believe these are false positives (mainly own names, abbreviations and technical terms).

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Prior submission

> Please do not start the description with "This package", package name, title or similar.

Description now starts with "A collection of methods that fit nonlinear distance transformations"

> Please always explain all acronyms in the description text. -> 'SMACOF'

Done.

> If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> [...]

Done.

> Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation.> Please write about the structure of the output (class) and also what the output means. (If a function does not re> turn a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
> Missing Rd-tags in up to 36 .Rd files, e.g.:
>     bcStressMin.Rd: \value
>     c_association.Rd: \value
>     c_clumpiness.Rd: \value
>     c_clusteredness.Rd: \value
>     c_complexity.Rd: \value
>     c_convexity.Rd: \value
>     ...

Added @return statement to all the roxygen documentations of exported functions.

> "Using foo:::f instead of foo::f allows access to unexported objects. This is generally not recommended, as the semantics of unexported objects may be changed by the package author in routine maintenance."
> Used ::: in documentation:
>      man/c_regularity.Rd:
>        hpts <- sp:::genHexGrid(dx = 0.9, ll = c(-2, -2), ur = c(2, 2))
> Please omit one colon.

Done.

> We see:
> Warning in tryError(code[which(nok)], Rdnames = R[which(nok)]) :
> 
> Unexecutable code in tests/testthat/test-cops.R:
>   unexpected '{': test_that("cops ndim argument"{
>
> Warning in tryError(code[which(nok)], Rdnames = R[which(nok)]) :
>
> Unexecutable code in tests/testthat/test-cordillera.R:
>   unexpected symbol: test_that
>
> Please make sure that your tests run smoothly.

The test were included in error; removed from the package via .Rbuildignore.

> Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.
> e.g.:
>
> ...
> old <- options() # code line i
> on.exit(options(old)) # code line i+1
> ...
> options(locatorBell = FALSE) # somewhere after
> ...
> e.g.: R/extras.R
> If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.

Removed options(locatorBell=FALSE).

