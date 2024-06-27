## Test environments
* local Linux Mint 21.3 install, R 4.4.0 and R 4.4.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check() in devtools_2.4.5)
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was one NOTE for r-release, r-devel and r-oldrel.

      * checking CRAN incoming feasibility ... NOTE
      	Maintainer: ‘Thomas Rusch <thomas.rusch@wu.ac.at>’

	Possibly misspelled words in DESCRIPTION:
  	hyperparameters (6:738)
  	smacofx (6:973)

The first flagged words seems to be a false positive. smacofx is a package name. 

## Downstream dependencies
There are currently no downstream dependencies for this package.