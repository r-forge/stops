## Fixes after first CRAN submission 

# Please always explain all acronyms in the description text.
The acronyms are now explained in the DESCRIPTION.

# Please add references describing the methods in your package to the description field of your DESCRIPTION file in the form
# authors (year) <doi:...>
# authors (year) <arXiv:...>
# authors (year, ISBN:...)
# or only if none those are available:  <https:...>
# with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
I have now done this. 

# Please add small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.
All exported functions (apart from standard S3 methods like print or summary) now have a small example.  

# I think there is a "," missing in test-cops.R in
# test_that("cops ndim argument"{
# ?
Thanks, fixed. 

# Please make sure that you do not change the user's options, par or working directory. If you really have to do so, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited. e.g.:
# ...
# oldoptions <- options(SCIPEN = 100)   # code line i
# on.exit(options(oldoptions))          # code line i+1
# ...
# e.g.:plot3dstatic.cmdscale(),...
The options call has been removed. 

# Please add \value to .Rd files regarding methods and explain the functions results in the documentation.
# e.g. pdist.Rd,... 
I added these to all functions (except the S3 methods) 

## Test environments
* local Linux Mint 19.2. install, R 3.6.1
* Ubuntu 18.04, R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

## Win-builder results
## r-relase and r-devel
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

New submission
 
Possibly mis-spelled words in DESCRIPTION:
  Buja (7:1357)
  Groenen (7:549, 7:1284)
  Hornik (7:346, 7:1503, 7:1756)
  Jakola (7:1723)
  Leeuw (7:1002, 7:1225, 7:1277)
  Luus (7:1718)
  MDS (7:93, 7:537, 7:793, 7:860, 7:1599)
  Mair (7:355, 7:1294, 7:1496, 7:1749)
  Rusch (7:339, 7:1489, 7:1742)
  SMACOF (7:991)
  Sammon (7:1061, 7:1077, 7:1443)
  Swayne (7:1364)
  Takane (7:1206)
  Torgerson (7:891, 7:910)
  clusteredness (7:295)
  de (7:999, 7:1222, 7:1274)
  majorizing (7:960)
  proximities (7:211)

majorizing, proximities and clusteredness are flagged but I think they are false positives. The acronyms MDS and SMACOF are explained in the text. All the other words are own names and not misspelled.


## Downstream dependencies
There are currently no downstream dependencies for this package.