## Test environments
* local Linux Mint 17.1 64-bit install, R 3.4.3
* win-builder (R-devel, R-release, R-oldrelease)

## R CMD check --as-cran --run-donttest results
There were no ERRORs or WARNINGs or NOTEs. 

## Win-Builder

There was 1 NOTE:

      * Possibly mis-spelled words in DESCRIPTION:
      	Hornik (8:273)
  	Mair (8:281)
  	Rusch (8:266)

	These are own names of the paper authors and spelled correctly.

## System requirements 
For e_optics() and e_cordillera() to work there needs to be a recent version of ELKI installed (https://elki-project.github.io/). I tested the package functionality with ELKI 0.6.1 release.

## Downstream dependencies
There are currently no downstream dependencies for this package.
       