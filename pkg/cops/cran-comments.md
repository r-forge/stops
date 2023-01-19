## Test environments
* local Linux Mint 19.2 install, R 4.2.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check())
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was one NOTE for r-release and r-devel.

      * checking CRAN incoming feasibility ... NOTE
      	Maintainer: ‘Thomas Rusch <thomas.rusch@wu.ac.at>’

       	Possibly misspelled words in DESCRIPTION:
        Jaakola (7:1808)

The flagged word is an own name.

There was one NOTE for r-release.

      * checking CRAN incoming feasibility ... NOTE
      	Maintainer: ‘Thomas Rusch <thomas.rusch@wu.ac.at>’

       	Possibly misspelled words in DESCRIPTION:
        Jaakola (7:1808)

	Found the following (possibly) invalid DOIs:
  	DOI: 10.1111/j.2044-8317.1966.tb00367.x
    	From: DESCRIPTION
    	Status: Service Unavailable
    	Message: 503

The flagged word is an own name. The invalid DOI NOTE seems to be a false positive on CRAN's side. The DOI presents as valid and leads to the article in question, see https://bpspsychub.onlinelibrary.wiley.com/doi/10.1111/j.2044-8317.1966.tb00367.x 

There was one NOTE for r-oldrel.

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

Possibly mis-spelled words in DESCRIPTION:
  Buja (7:1413)
  Groenen (7:602, 7:1340)
  Hornik (7:191, 7:413, 7:1585, 7:1843)
  Jaakola (7:1808)
  Leeuw (7:1055, 7:1280, 7:1333)
  Luus (7:1803)
  MDS (7:93, 7:590, 7:846, 7:913, 7:1683)
  Mair (7:184, 7:422, 7:1350, 7:1578, 7:1836)
  Rusch (7:177, 7:406, 7:1571, 7:1829)
  SMACOF (7:1044)
  Sammon (7:1115, 7:1131, 7:1547)
  Swayne (7:1420)
  Takane (7:1261)
  Torgerson (7:944, 7:963)
  clusteredness (7:362)
  de (7:1052, 7:1277, 7:1330)
  majorizing (7:1013)
  proximities (7:278)


The flagged words are own names, abbreviations (that are expanded) or technical terms. Proximities seems to be a false positive.

## Downstream dependencies
There are currently no downstream dependencies for this package.