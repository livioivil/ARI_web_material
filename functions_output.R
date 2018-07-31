merge_tables_print<- function(res_all_sim1,res_all_sim4,res_all_sim16,options,...){
  
  temp=array(c(res_all_sim1,res_all_sim4,res_all_sim16),c(dim(res_all_sim1),3))
  temp_names=dimnames(res_all_sim1)
  temp_names$nClsuters=c("1","4","16")
  dimnames(temp)=temp_names
  temp=aperm(temp,c(4,1,2,3))
  # rm(res_all_sim)
  names(dimnames(temp))[1]="#clusters"
  names(dimnames(temp))[2]="Radius/sqrt(#clusters)"
  dimnames(temp)[[2]]=paste(dimnames(temp)[[2]]," $(\\pi_1=",round(as.numeric(dimnames(temp)[[2]])^2*pi/(64^2),3),sep="",")$")
  
  
  x=ftable(as.table(temp))
  colnames(x)=dimnames(temp)[[4]]
  
  knitr::asis_output(
    format_html(x,
                digits=if(length(options$ftable.digits))
                  options$ftable.digits
                else 0,
                ...))
}

merge_tables_print_fixed_effsize<- function(effsize=1,res_all_sim1,res_all_sim4,res_all_sim16,options,...){
  
  temp=array(c(res_all_sim1,res_all_sim4,res_all_sim16),c(dim(res_all_sim1),3))

    temp_names=dimnames(res_all_sim1)
  temp_names$nClsuters=c("1","4","16")
  dimnames(temp)=temp_names
  temp=aperm(temp,c(4,1,2,3))
  temp=temp[,,effsize,]
  
  # rm(res_all_sim)
  names(dimnames(temp))[1]="#clusters"
  names(dimnames(temp))[2]="Radius/sqrt(#clusters)"
  # dimnames(temp)[[2]]=paste(dimnames(temp)[[2]]," $(\\pi_1=",round(as.numeric(dimnames(temp)[[2]])^2*pi/(64^2),3),sep="",")$")

  dnames=dimnames(temp)
  temp=array(c(round(as.vector(t(matrix(pis,4,3))),4),temp),dim(temp)+c(0,0,1))
  dnames$method=c("$\\pi_1$",dnames$method)
  
  dimnames(temp)=dnames
  
  x=ftable(as.table(temp))
  colnames(x)=dimnames(temp)[[3]]
  
  knitr::asis_output(
    format_html(x,
                digits=if(length(options$ftable.digits))
                  options$ftable.digits
                else 0,
                ...))
}

reshape_table<- function(res_all_sim){
  res=aaply(res_all_sim,c(1,2,3,5),mean)
  res=res[,,,c("FSL_any_power",
               "ARI_anyPTD_entire_brain")]#,"ARI_PTD_median_ave_power")]
  res=aperm(res,c(2,1,3,4))
  # res=res[,,-1,]
  temp1=aperm(res[,,,1],c(2,3,1))
  temp2=array(res[1,,,2],c(dim(res[1,,,2]),1))
  res2=array(c(temp1,temp2),dim(temp1)+c(0,0,1))
  dimnames(res2) =list(radius=dimnames(res)[[2]],effectsize=dimnames(res)[[3]],method=c("FSL(2.3th)","FSL(3.2th)","ARI"))
  
  res2
}

# make table knitr
library(memisc)
library(knitr)
knit_print.ftable_matrix <-function(res_all_sim,options,...){
  res=aaply(res_all_sim,c(1,2,3,5),mean)
  res=res[,,,c("FSL_any_power",
               "ARI_anyPTD_entire_brain")]#,"ARI_PTD_median_ave_power")]
  res=aperm(res,c(2,1,3,4))
  res=res[,,-1,]
  
  
  x=ftable(as.table(res))
  colnames(x)=dimnames(res)[[4]]
  
  knitr::asis_output(
    format_html(x,
                digits=if(length(options$ftable.digits))
                  options$ftable.digits
                else 0,
                ...))
}