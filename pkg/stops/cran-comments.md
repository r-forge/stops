## Test environments
* local Linux Mint 21.3 install, R 4.4.0 and R 4.4.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check() in devtools_2.4.5)
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was the following NOTE:

* checking CRAN incoming feasibility ... [15s] NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

Possibly misspelled words in DESCRIPTION:
  ALSCAL (6:913)
  CLCA (6:1304, 6:1491)
  CLDA (6:1342)
  curvilinear (6:1272, 6:1311)
  smacofx (6:629, 6:1469)
  sparsified (6:1352)

I believe these are false positives as they are abbreviations with full names given. smacofx is a package. curvilinear and sparsified are technical terms.

## Downstream dependencies
There are currently no downstream dependencies for this package.

