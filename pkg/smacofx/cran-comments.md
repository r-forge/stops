## Test environments
* local Linux Mint 22.1 install, R 4.5.0
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check())
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was 1 NOTE.

Version jumps in minor (submitted: 1.20.1, existing: 1.6.1)

Possibly misspelled words in DESCRIPTION:
  nonmetric (10:889, 10:2111)

The version jump is correct. The misspelled word is a false positive.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of smacofx. All packages passed.
