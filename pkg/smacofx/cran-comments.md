## Test environments
* local Linux Mint 21.3 install, R 4.4.0 and R 4.4.1
* win-builder (devel, release, oldrel)

## R CMD check results (via devtools::check() in devtools_2.4.5)
There were no ERRORs or WARNINGs or NOTEs. 

## Win-builder results
## r-release, r-devel and r-oldrel
There were no ERRORs or WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Fixes upon Resubmission
Flavor: r-devel-linux-x86_64-debian-gcc
Check: Rd cross-references, Result: NOTE
  Found the following Rd file(s) with Rd \link{} targets missing package
  anchors:
    alscal.Rd: smacofSym
    apStressMin.Rd: smacofSym
    elscal.Rd: smacofSym
    multiscale.Rd: smacofSym
    powerStressFast.Rd: smacofSym
    powerStressMin.Rd: smacofSym
    rStressMin.Rd: smacofSym
    rpowerStressMin.Rd: smacofSym
    sammonmap.Rd: smacofSym
    spmdda.Rd: smacofSym
    spmds.Rd: smacofSym
  Please provide package anchors for all Rd \link{} targets not in the
  package itself and the base packages.


I changed the offending expression to \code{\link[smacof][smacofSym}} in all thee .Rd