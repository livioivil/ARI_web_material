compute_stats <- function(clstr,disco,id_signal,
                          clusters_FSL,m1,id_H1_clusters,id_H0_clusters){
  # clstr=c(0,0,0,1,1,0,0,2,2,2)
  # cls1 not signif
  
  clstr_signif=force_zero_these_clusters(clstr,
                                         clusters_FSL$`Cluster.Index`[disco==0])
  # clstr_signif=c(0,0,0,0,0,0,0,2,2,2)
  
  results=list()
  if(length(id_signal)>0){
    results$ave_power=sum(clstr_signif[id_signal]>0)/m1
  } else {
    results$ave_power=NA
    }
  results$any_power=(results$ave_power>0)*1
  results$all_power=(results$ave_power==1)*1
  
  #### proportion of non active voxels selected
  if(length(id_signal)>0){
    results$FFP_voxels=sum(clstr_signif[-id_signal]>0)/
    max(1,sum(clstr_signif>0))
  } else {
    results$FFP_voxels=any(clstr_signif>0)
    }
  
  #######select clusters without any active voxels and check if rejected
  clstr_signif_H0=force_zero_these_clusters(clstr_signif,id_H1_clusters)
  
  if(length(id_signal)>0){
    results$any_false_disc=any(clstr_signif_H0[-id_signal]>0)
  } else {
    results$any_false_disc=any(clstr_signif_H0>0)
  }
  
  #- N_cluster TRUE returned by FLS and signif (rejected H0)
  results$N_cluster_H0_rejected = sum(clusters_FSL$P[id_H0_clusters]<=alpha)
  #- N_cluster FALSE returned by FLS and signif (rejected H1)
  if(length(id_signal)>0){
    results$N_cluster_H1_rejected = sum(clusters_FSL$P[id_H1_clusters]<=alpha)
  } else {
    results$N_cluster_H1_rejected = NA
    }
  
  results$ANY_cluster_H0_rejected=(results$N_cluster_H0_rejected>0)*1
  results$cluster_power=results$N_cluster_H1_rejected/
    max(results$N_cluster_H1,1)
  results
}
