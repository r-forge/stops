        
#### Assign the registry to variable struc__reg in the topenv()       
#  assign("struc_reg", R, envir = parent.env(environment()))
#  }
#struc_reg <- .onLoad()
#' function for lookup that partially matched and ignores cases and punctuation
#' @param lookup the lookup string
#' @param entry the registry entry
#' @param ... additional arguments to pmatch 
match_partial_ignorecase_nopunct <- function(lookup,entry,...)
     {
      lookup <- gsub("[[:punct:]]", "", lookup)    
      !is.na(pmatch(tolower(lookup), tolower(entry), ...))
     }

#' @importFrom registry registry
struc_reg <- registry::registry(registry_class="structureindices")#, validity_FUN=checkstruc_regeturn)

#### setting up the registry   
struc_reg$set_field("name",type="character",is_key=TRUE,index_FUN=match_partial_ignorecase_nopunct)
struc_reg$set_field("index",type="function")#,validity_FUN=checkReturn)
     
#### create the registry entries
## c inequality
     struc_reg$set_entry(name="cinequality",index=c_inequality)
     ## c linearity
     struc_reg$set_entry(name="clinearity",index=c_linearity)
     ## c dependence
     struc_reg$set_entry(name="cdependence",index=c_dependence)
     ## c manifoldness 
     struc_reg$set_entry(name="cmanifoldness",index=c_manifoldness)
     ## c association
     struc_reg$set_entry(name="cassociation",index=c_association)
     ## c nonmonotonicity 
     struc_reg$set_entry(name="cnonmonotonicity",index=c_nonmonotonicity)
     ## c functionality
     struc_reg$set_entry(name="cfunctionality",index=c_functionality)
     ## c complexity
     struc_reg$set_entry(name="ccomplexity",index=c_complexity)
     ## c faithfulness
     struc_reg$set_entry(name="cfaithfulness",index=c_faithfulness)
     ## c clusteredness
     struc_reg$set_entry(name="cclusteredness",index=c_clusteredness)
     ## c regularity
     struc_reg$set_entry(name="cregularity",index=c_regularity)
     ## c hierarchy
     struc_reg$set_entry(name="chierarchy",index=c_hierarchy)
     ## c outlying
     struc_reg$set_entry(name="coutlying",index=c_outlying)
     ## c convexity
     struc_reg$set_entry(name="cconvexity",index=c_convexity)
     ## c skinniness
     struc_reg$set_entry(name="cskinniness",index=c_skinniness)
     ## c stringiness
     struc_reg$set_entry(name="cstringiness",index=c_stringiness)
     ## c sparsity
     struc_reg$set_entry(name="csparsity",index=c_sparsity)
     ## c clumpiness
     struc_reg$set_entry(name="cclumpiness",index=c_clumpiness)
     ## c striatedness
     struc_reg$set_entry(name="cstriatedness",index=c_striatedness)
     ## c shepardness
     struc_reg$set_entry(name="cshepardness",index=c_shepardness)
        
  #### Assign the registry to variable struc__reg in the topenv()       
#  assign("struc_reg", R, envir = parent.env(environment()))
#  }
#struc_reg <- .onLoad()

 
