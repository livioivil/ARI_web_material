#' @description  Valid Circular Inference for Brain Imaging. As a main feature, 
#' it estimate lower bounds for the proportion of active voxels in a set of clusters as, for example, given by a cluster-wise analysis.
#' @author all of us
#' @docType package
#' @name ARI-package
#' @title Valid Circular Inference
#' @import hommel
#' @examples 
#' pvalue_name <- system.file("extdata", "pvalue.nii.gz", package="ARI")
#' cluster_name <- system.file("extdata", "cluster_th_3.2.nii.gz", package="ARI")
#' zstat_name <- system.file("extdata", "zstat.nii.gz", package="ARI")
#' mask_name <- system.file("extdata", "mask.nii.gz", package="ARI")
#' ARI(Pmap = pvalue_name, clusters= cluster_name, 
#'     mask=mask_name, Statmap = ztat_name)
#' 
NULL