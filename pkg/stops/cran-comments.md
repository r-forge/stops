## Test environments
* local Linux Mint 19.2 install, R 4.2.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check())
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There were the following NOTEs:

Possibly misspelled words in DESCRIPTION:
  LJ (8:886)
  MDS (8:167, 8:369, 8:540, 8:609)
  SMACOF (8:437, 8:488, 8:506)
  STructure (2:8)
  Sammon (8:445, 8:574)
  Torgerson (8:406)
  boxcox (8:615)
  hyperparameter (8:750)
  hyperparameters (8:712)
  structuredness (8:1001)

I believe these are false positives (mainly own names and technical terms).


## Downstream dependencies
There are currently no downstream dependencies for this package.