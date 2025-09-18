## Test environments
* local Linux Mint 22.1 install, R 4.5.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check())
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs or NOTEs.

There was 1 NOTE for r-oldrel

* checking DESCRIPTION meta-information ... NOTE
Author field differs from that derived from Authors@R
  Author:    'Thomas Rusch [aut, cre] (ORCID: <https://orcid.org/0000-0002-7773-2096>), Jan de Leeuw [aut], Lisha Chen [aut], Patrick Mair [aut] (ORCID: <https://orcid.org/0000-0003-0100-6511>)'
  Authors@R: 'Thomas Rusch [aut, cre] (<https://orcid.org/0000-0002-7773-2096>), Jan de Leeuw [aut], Lisha Chen [aut], Patrick Mair [aut] (<https://orcid.org/0000-0003-0100-6511>)'

I'm not sure about this NOTE, there is only an Authors@R field in the DESCRIPTION file.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of smacofx. All packages passed.
