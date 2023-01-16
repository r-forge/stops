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