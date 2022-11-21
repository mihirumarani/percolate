library(tidyverse)
library(RANN)
library(furrr)
library(magrittr)
library(parallel)

perc.up=function(dat,d_cutoff=seq(1,100,by=5)){#dat is a full grid of 0s and 1s.
  
  dat=dat%>%
      as_tibble()%>%
      mutate(y=1:nrow(dat))%>%
      pivot_longer(cols="1":"1000",names_to = "x",values_to = "pres")%>%
      mutate(x=as.numeric(x),
             y=as.numeric(y))
  
  result=NULL    
  
  for(i in d_cutoff){
  
  cutx=seq(0, Lx, i)
  cuty=seq(0,Ly, i)
  
   dat1=dat%>%
    mutate(x = cut(x, cutx, labels = FALSE,include.lowest=TRUE),
          y = cut(y, cuty, labels = FALSE,include.lowest=TRUE))%>%
     group_by(x,y)%>%
     summarize(pres=sum(pres))%>%
     ungroup()%>%
     pivot_wider(names_from = x,values_from = pres)%>%
     select(-y)%>%
     as.matrix()
  
    dat1[dat1!=0]=1
   
  clust=matrix(0,nrow=nrow(dat1),ncol=ncol(dat1))
  
  largest=0
  
  lbl=0:length(clust)
  
  for(y in seq(ncol(dat1))){
    for(x in seq(nrow(dat1))){
      
      if((dat1[x,y]==1)){
        
        left=ifelse(length(clust[x,y-1])==0,0,clust[x,y-1])
        above=ifelse(length(clust[x-1,y])==0,0,clust[x-1,y])
        
        if((left==0) & (above==0)){
          
          largest=largest+1
          clust[x,y]=largest
          lbl[10*(y-1)+x]=largest
          
        }else{
          
          if((left!=0) & (above==0)){
            
            clust[x,y]=clust[x,y-1]
            
            lbl[10*(y-1)+x]=lbl[10*(y-2)+x]
            
          }else{
            
            if(left==0 & above!=0){
              
              clust[x,y]=clust[x-1,y]
              
              lbl[10*(y-1)+x]=lbl[10*(y-1)+x-1]
              
            }else{
              
              newlab=min(left,above)
              
              lbl[which(lbl==left | lbl==above)]=newlab
              
              clust[x,y]=newlab
              
              clust[which(clust==max(left,above))]=newlab
              
            }
          }
        }
      }
    }
  }
  
  #Get the #clusters and the size of biggest cluster
  clustab=as.vector(clust)%>%as_tibble%>%count(value)%>%filter(value!=0)
  
  result=bind_rows(
          result,
          tibble(
            d_cutoff=i,
            clust=clustab$value,
            freq=clustab$n
          )
        )

  }
  
  return(result)
}

Lx=1000
Ly=500

data=readRDS('bci_raw.rds')%>%
  select(census,treeID,sp,gx,gy,dbh,status)

dat=data%>%
  filter(census==7,status=='A')%>%
  count(sp)%>%
  filter(n>200)%>%
  slice(sample(length(n),5))%>%
  left_join(data)%>%
  filter(status=="A")

dat=dat%>%
  inner_join(
    dat%>%
      group_by(sp)%>%
      mutate(baldeck=quantile(dbh,0.56,na.rm=TRUE))
  )%>%
  filter(dbh>baldeck)

splist=dat%>%count(sp)




params=expand.grid(
  d_cutoff=seq(10,150,by=10),
  census=1:8,
  sps=splist$sp
)

cores=detectCores()
plan(multisession,workers=cores)

   result=params%>%
     future_pmap_dfr(
       .f=function(
            d_cutoff,
            census,
            sps
            
       ){
         datsp1=dat%>%
           filter(sp==sps,census==census)%>%
           drop_na()%>%
           mutate(
             x = cut(gx, seq(0, Lx, d_cutoff), labels = FALSE,include.lowest=TRUE),
             y = cut(gy, seq(0,Ly, d_cutoff), labels = FALSE,include.lowest=TRUE)
           )%>%
           count(x,y)%>%
           complete(x,y)%>%
           pivot_wider(names_from = y,values_from = n)%>%
           select(-x)%>%
           as.matrix()
         
         datsp1[is.na(datsp1)]=0
         datsp1[datsp1!=0]=1
         
         res=perc.up(datsp1)
         
         res=tibble(
           census=census,
           sp=sps,
           d_cutoff=d_cutoff,
           nclust=res$nclust,
           clust.size=res$clust.size
         )%>%return()
         
         
       },
       .options = furrr_options(seed = TRUE)
     )
    
    result%>%
      ggplot(aes(d_cutoff,clust.size,col=census))+geom_line()+facet_wrap(~sp)
    
    
  
  




