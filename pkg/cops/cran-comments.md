## Test environments
* local Linux Mint 19.2. install, R 3.6.1
* Ubuntu 18.04, R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

New submission
 
Possibly mis-spelled words in DESCRIPTION:
  LJ (7:1162)
  MDS (7:93, 7:471, 7:680, 7:747, 7:900, 7:1040)
  SMACOF (7:797, 7:848, 7:866)
  Sammon (7:805, 7:934)
  Torgerson (7:778)
  clusteredness (7:295)
  proximities (7:211)

  proximities and clusteredness are flagged but I think they are false positives. All the other words are own names or abbreviations and not misspelled.

## Downstream dependencies
There are currently no downstream dependencies for this package.