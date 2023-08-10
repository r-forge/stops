## Test environments
* local Linux Mint 19.2 install, R 4.3.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check() in devtools_2.4.5)
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs.

There was one NOTE for r-devel and r-oldrel.

      * checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

New submission

Possibly misspelled words in DESCRIPTION:
  ALSCAL (16:6)
  Buja (20:16, 24:25, 25:23)
  De (16:30, 18:36)
  Demartines (26:6)
  Groenen (18:46)
  Herault (26:19)
  Hornik (21:44, 23:20)
  Leeuw (16:33, 18:39)
  Lendasse (27:11)
  MDS (10:49, 12:99, 15:16, 15:111, 18:31, 19:59, 20:11, 24:13, 25:11, 28:30, 28:54)
  Mair (18:56, 21:37, 23:13)
  Majorization (31:5)
  Multiscale (15:5)
  Rusch (21:30, 23:6)
  Sammon (14:5, 14:21)
  Swayne (20:23)
  Takane (16:14)
  Torgerson (13:42, 13:61)
  Verleysen (27:22)
  curvilinear (25:66, 26:63)
  majorization (10:91)
  proximities (29:93, 30:73)
  sparsified (28:19, 28:38)


The flagged words are all own names or acronyms that are explained in the text. Curvilinear, sparsified and proximities seem to be false positives. 

There was one NOTE for r-release.

      * checking CRAN incoming feasibility ... NOTE
Maintainer: 'Thomas Rusch <thomas.rusch@wu.ac.at>'

New submission

Possibly misspelled words in DESCRIPTION:
  ALSCAL (16:6)
  Buja (20:16, 24:25, 25:23)
  De (16:30, 18:36)
  Demartines (26:6)
  Groenen (18:46)
  Herault (26:19)
  Hornik (21:44, 23:20)
  Leeuw (16:33, 18:39)
  Lendasse (27:11)
  MDS (10:49, 12:99, 15:16, 15:111, 18:31, 19:59, 20:11, 24:13, 25:11, 28:30, 28:54)
  Mair (18:56, 21:37, 23:13)
  Majorization (31:5)
  Multiscale (15:5)
  Rusch (21:30, 23:6)
  Sammon (14:5, 14:21)
  Swayne (20:23)
  Takane (16:14)
  Torgerson (13:42, 13:61)
  Verleysen (27:22)
  curvilinear (25:66, 26:63)
  majorization (10:91)
  proximities (29:93, 30:73)
  sparsified (28:19, 28:38)

Found the following (possibly) invalid DOIs:
  DOI: 10.1080/10618600.2020.1869027
    From: DESCRIPTION
          DESCRIPTION
    Status: Forbidden
    Message: 403
  DOI: 10.1111/j.2044-8317.1966.tb00367.x
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
  DOI: 10.1198/jasa.2009.0111
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

The flagged words are all own names or acronyms that are explained in the text. Curvilinear, sparsified and proximities seem to be false positives. The invalid DOI NOTE seem to be false positive on CRAN's side. The DOI presents as valid and leads to the article in question, see the URLs https://www.tandfonline.com/doi/full/10.1080/10618600.2020.1869027 https://bpspsychub.onlinelibrary.wiley.com/doi/abs/10.1111/j.2044-8317.1966.tb00367.x https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.0111.

## Downstream dependencies
There are currently no downstream dependencies for this package.