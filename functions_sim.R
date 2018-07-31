######UTILS
source('myggplot.R')
source('functions_output.R')
### take the map of cluster, set to 0 the non signif ones.
force_zero_these_clusters<- function(map, cluster_id){
  for(ii in cluster_id){
    map[map==ii]=0
  }
  map
}

###########
# take output and transform to data.frame
convert_FSL_table_to_data_frame <- function(out){
  temp=matrix(unlist(strsplit(out,"\\t")),nrow = length(out),byrow = TRUE)
  colnames(temp)=temp[1,]
  temp=temp[-1,,drop=FALSE]
  rownames(temp)=temp[,'Cluster Index']
  if(nrow(temp)>0)
    temp=apply(temp,2,as.numeric) else 
      temp=t(aaply(temp,2,as.numeric)) 
  if(is.null(dim(temp))){
    temp=as.data.frame(t(temp))
    } else{
      temp=data.frame(temp)
    }
  names(temp)=gsub("Cluster.Index","Cluster.Index",names(temp))
  
  if(nrow(temp)>0) temp=temp[nrow(temp):1,,drop=FALSE]
  # names(temp) = gsub(' ','',names(temp))
  temp
}
############################
generate_signal <- function(effsize=1,radius=10,
                            coordList = list(c(32,32,1))){
  signalstrength <- list(effsize) #list of signal strengths (.5, 1, 2) (SNR = signalstrength / SD(noise)). SD(noise) = 1.
   #coordList: list of coordinates of spheres
  sphereradius <- list(radius) #list of radii of spheres
  
  #make signal spatially
  signalVolume <- array(0,dim = dims)
  for(i in 1:length(coordList)) {
    signalVolume <- signalVolume + 
      specifyregion(dim = dims,
                    coord = coordList[[i]],
                    radius = sphereradius[[1]],
                    form = 'sphere',
                    fading = 0)*signalstrength[[1]]
  }
  signalVolume
}

######################
sim_function <- function(j,signal4D,t_threshold) {
  set.seed(j)
  filenamekey=gsub("\\W","",paste0(gsub(" ", "", date()),runif(1)))
  #empty results list
  results=list()
  
  signal=as.vector(signal4D>0)
  id_signal=which(signal)
  id_NOTsignal=which(!signal)
  m1=length(id_signal)
  m0=length(id_NOTsignal)
  
  
  
  #set-up spatial properties for noise and make noise volume (generate z-scores directly)
  noise4D <- array(spatialnoise(dim = dims,
                                sigma = 1,
                                nscan = 1,
                                method = 'gaussRF',
                                FWHM = FWHM),dim = dims)
  
  
  
  # add spatial noise
  z4D <- abs(signal4D + noise4D)
  
  #make into nifti and save
  dat <- new('fmri.data')
  dat@dims = c(4,dims,1,rep(0,3))
  dat@datavec = as.vector(z4D)
  dat@filename = paste0('testdata',filenamekey)
  dat@fullpath = data_folder
  dat@descrip = 'ARI - Simulation data with Neurosim'
  writeData(dat)
  
  # convert Zs to p-values
  p4D=pnorm(-abs(z4D))*2
  
  #### FSL analysis
  
out=system(paste0('. /etc/fsl/5.0/fsl.sh
                    cluster -i ',data_folder,'/testdata',filenamekey,' -t ',t_threshold,' -p 1 -d ',DLH,' --volume=',VOLUME,' --mm  -o ',data_folder,'/cluster_sim',filenamekey,'_th',t_threshold,'.nii'),intern=TRUE)
  
  #capture output and transform to data.frame
  clusters_FSL=try(convert_FSL_table_to_data_frame(out),silent = TRUE)
  if(is(clusters_FSL,"try-error")) browser()
 
    clusters_FSL$P[clusters_FSL$P<.0001]=.0001
    n_clusters=nrow(clusters_FSL)
    # read FLS clusters
    clstr = readNifti(paste0(data_folder,'/cluster_sim',filenamekey,'_th',t_threshold,'.nii'))
  
    file.remove(paste0(data_folder,'/testdata',filenamekey,'.nii.gz'))
    file.remove(paste0(data_folder,'/cluster_sim',filenamekey,'_th',t_threshold,'.nii.gz'))
    
  
    tab= table(clstr*signal)[-1] #is the signal in the clusters
    results$N_cluster_H1= length(tab) #returned by FLS TRUE
    results$N_cluster_H0= n_clusters-results$N_cluster_H1 # returned by FLS FALSE 
    results$N_suprathreshold=sum(clstr!=0)  
    results$N_ACTIVES_in_suprathreshold_set=sum(clstr[signal]>0)
  
    results$TRUE_PTD_in_suprathreshold_set=
      results$N_ACTIVES_in_suprathreshold_set/
      max(1,results$N_suprathreshold,na.rm=TRUE)
    
    if(results$N_cluster_H0<0){
      cat("\nBUG HERE\n")
      print(results$N_cluster_H0)
    }
    
    id_H1_clusters=as.numeric(names(tab))
    id_H0_clusters=setdiff(1:n_clusters,id_H1_clusters)
    
    
  

  #############
  ## for FSL:
  disco_fsl=(clusters_FSL$P<=alpha)*1
    # browser()
    # if(nrow(clusters_FSL)==0) {
      # if(sum(id_signal)>0) browser()  } 
    
  results_fsl=compute_stats(clstr = clstr,disco = disco_fsl,id_signal,
                            clusters_FSL,m1,id_H1_clusters,id_H0_clusters)
  names(results_fsl)=paste0("FSL_",names(results_fsl))
  
  
  #######cherry
  hF<-hommel(as.vector(p4D))
  # str(hF)
  # sum(hF@p<alpha)
  # sum(hF@adjusted<alpha)
  # cherry_discoveries <- discoveries(hF, alpha=alpha)
  
  # table(as.vector(clstr))
  
  disco_ari=laply(clusters_FSL$`Cluster.Index`,function(h)discoveries(hF, ix = clstr==h,alpha=alpha))
  # disco_ari
  
  results_ari=compute_stats(clstr = clstr,disco = disco_ari,id_signal,
                            clusters_FSL,m1,id_H1_clusters,id_H0_clusters)
  # results_ari$any_power_any <- (discoveries(hF)>0)*1
  results_ari$EST_discovs_in_supraThreshold=sum(disco_ari)
  results_ari$voxel_power_in_clusters=
    results_ari$EST_discovs_in_supraThreshold/
    max(1,results$N_ACTIVES_in_suprathreshold_set)
  results_ari$PTD_in_suprathreshold_set=
    results_ari$EST_discovs_in_supraThreshold/
    max(1,results$N_suprathreshold)
  
    names(results_ari)=paste0("ARI_",names(results_ari))
  

  # PTD active set: true one
  if(m1>0){
    results_ari$ARI_PTD_active_set <- tdp(hF, alpha=alpha,ix = signal)
    results_ari$ARI_PTD_active_set_median <- tdp(hF, alpha=.5,ix = signal)
    
  } else {
    results_ari$ARI_PTD_active_set <- 0
    results_ari$ARI_PTD_active_set_median <- 0
    
  } 
  
  # larger than 0:
  results_ari$ARI_anyPTD_active_set=results_ari$ARI_PTD_active_set>0
  
  # PTD in entire brain: true one	
  results_ari$ARI_PTD_entire_brain <- tdp(hF, alpha=alpha)
  results_ari$ARI_PTD_entire_brain_median <- tdp(hF, alpha=.5)

  results_ari$ARI_PTD_ave_power <- results_ari$ARI_PTD_entire_brain/(m1/(m1+m0))
  results_ari$ARI_PTD_median_ave_power <- results_ari$ARI_PTD_entire_brain_median/(m1/(m1+m0))
  
  # larger than 0:
  results_ari$ARI_anyPTD_entire_brain=results_ari$ARI_PTD_entire_brain>0
  
  # % of time signal found in non-active set: true one (PTD inactive set)	
  results_ari$ARI_anyFP_null_set <- (discoveries(hF, alpha=alpha,ix = !signal)>0)*1
  
  # % of time signal found somewhere in the brain	: see:
  # results$ARI_any_power

  results=c(results,results_ari,results_fsl)
  ######### voxel-wise
  
  if(m1>0){
    # voxels found
    p.adj.holm=p.adjust(hF@p,"holm")
    results$HOLM_ave_power=sum(p.adj.holm[id_signal]<alpha)/m1
    # p.adj.hommel=p.adjust(hF@p,"hommel")
    results$HOMMEL_ave_power_base=NA#sum(p.adj.hommel[id_signal]<alpha)/m1
    results$HOMMEL_ave_power=sum(hF@adjusted[id_signal]<alpha)/m1
  
    p.adj.bh=p.adjust(hF@p,"BH")
    results$BH_ave_power=sum(p.adj.bh[id_signal]<alpha)/m1
  } else {
    results$BH_ave_power <- results$HOLM_ave_power <-
      results$HOMMEL_ave_power <-results$HOMMEL_ave_power_base <- NA
  }
  
  unlist(results)
}


###################
run_sims_all_effsizes <- function(effsizes,thrshld=2.3,radius=5){
  res_all=laply(effsizes, function(effsize,thrshld2=thrshld){
    signal4D <- generate_signal(effsize,radius = radius,coordList)
    out=laply(1:B,sim_function,signal4D=signal4D,t_threshold=thrshld2)
    out
  })
  # str(res)
  dimnames(res_all)[[1]]=effsizes
  res_all
}
