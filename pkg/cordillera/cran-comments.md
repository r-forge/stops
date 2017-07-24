## Test environments
* local Linux Mint 17.1 64-bit install, R 3.4.0
* local Ubuntu 14.04, R 3.4.0
* win-builder (devel and release)

## R CMD check --as-cran --run-donttest results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

      * checking CRAN incoming feasibility ... NOTE
      Maintainer: ‘Thomas Rusch <thomas.rusch@wu.ac.at>’

      New submission

This is my first package on CRAN.      

## System requirements 
For e_optics() and e_cordillera() to work there needs to be a recent version of ELKI installed (https://elki-project.github.io/). I tested the package functionality with ELKI 0.6 release.

## Downstream dependencies
There are currently no downstream dependencies for this package
       