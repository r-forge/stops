#' onload function creates the registry object of c-structuredness indices called struc_reg
#' @importFrom registry registry
#' @param libname libraryname dummy
#' @param pkgname pkgname dummy
.onLoad <- function(libname,pkgname)
    {
     match_partial_ignorecase_nopunct <- function(lookup,entry,...)
     {
      lookup <- gsub("[[:punct:]]", "", lookup)    
      !is.na(pmatch(tolower(lookup), tolower(entry), ...))
     }
     ## check for c-structuredness; must only return a numeric scalar
     #checkReturn  <- function(f) {
     #                       tmpmat <- matrix(rnorm(10),ncol=2)
     #                       res <- f(tmpmat)  
     #                       stopifnot(length(res)==1 && is.numeric(res))
     #}
     R <- registry::registry(registry_class="structureindices")#, validity_FUN=checkReturn)
     # matches partially, ignores case and  punctuation
     #### setting up the registry   
     R$set_field("name",type="character",is_key=TRUE,index_FUN=match_partial_ignorecase_nopunct)
     R$set_field("index",type="function")#,validity_FUN=checkReturn)
     
     
     #### create the registry entries
     ## c inequality
     R$set_entry(name="cinequality",index=c_inequality)
     ## c linearity
     R$set_entry(name="clinearity",index=c_linearity)
     ## c dependence
     R$set_entry(name="cdependence",index=c_dependence)
     ## c manifoldness 
     R$set_entry(name="cmanifoldness",index=c_manifoldness)
     ## c association
     R$set_entry(name="cassociation",index=c_association)
     ## c nonmonotonicity 
     R$set_entry(name="cnonmonotonicity",index=c_nonmonotonicity)
     ## c functionality
     R$set_entry(name="cfunctionality",index=c_functionality)
     ## c complexity
     R$set_entry(name="ccomplexity",index=c_complexity)
     ## c faithfulness
     R$set_entry(name="cfaithfulness",index=c_faithfulness)
     ## c clusteredness
     R$set_entry(name="cclusteredness",index=c_clusteredness)
     ## c regularity
     R$set_entry(name="cregularity",index=c_regularity)
     ## c hierarchy
     R$set_entry(name="chierarchy",index=c_hierarchy)
     ## c outlying
     R$set_entry(name="coutlying",index=c_outlying)
     ## c convexity
     R$set_entry(name="cconvexity",index=c_convexity)
     ## c skinniness
     R$set_entry(name="cskinniness",index=c_skinniness)
     ## c stringiness
     R$set_entry(name="cstringiness",index=c_stringiness)
     ## c sparsity
     R$set_entry(name="csparsity",index=c_sparsity)
     ## c clumpiness
     R$set_entry(name="cclumpiness",index=c_clumpiness)
     ## c striatedness
     R$set_entry(name="cstriatedness",index=c_striatedness)
     ## c shepardness
     R$set_entry(name="cshepardness",index=c_shepardness)
        
  #### Assign the registry to variable struc__reg in the topenv()       
  assign("struc_reg", R, envir = parent.env(environment()))
  }
#struc_reg <- .onLoad()
