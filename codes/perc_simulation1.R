
library(tidyverse)
library(magrittr)
library(RandomFields)
library(furrr)
library(parallel)
library(starsExtra)

ind = commandArgs(trailingOnly = TRUE)

initset=read.csv("initset.csv")

reps=initset[ind,1]

rangepars=initset[ind,2]

generate_landscape = function(
    Lx,
    Ly,
    quadrat_length,
    rangepar,
    sillpar,
    nuggetpar,
    seed
) {
  stopifnot(nuggetpar >= 0 & nuggetpar <= 1)
  RFoptions(seed = seed)
  
  stress = 
    RFsimulate(
      RMgauss(
        scale = rangepar + 1e-16,
        var = 2 * sillpar * (1 - nuggetpar)
      ) + 
        RMtrend(mean = 0) + 
        RMnugget(var = 2 * sillpar * nuggetpar),
      x = seq(Lx / quadrat_length),
      y = seq(Ly / quadrat_length)
    )@data$variable1
  
  
  stress=stress-min(stress)
  stress=stress/max(stress)
  
  landscape = 
    expand_grid(
      y = seq(Ly / quadrat_length),
      x = seq(Lx / quadrat_length)
    ) %>%
    mutate(
      soiltype=stress
    )
  
  return(as_tibble(landscape))
}



#generate a simulation
#Inputs- Lx,Ly,quadrat length, dispersal kernel, rangepar of habitat autocorrelation, population size, mortality=fecundity,seed

Lx=1000
Ly=500
abun=1000
quadrat_length=1 
timepoints=10000
nbirth=ndeath=floor(abun*0.1)



landscape =  
  generate_landscape(
    Lx = Lx, 
    Ly = Ly, 
    quadrat_length = 1, 
    rangepar = rangepars, 
    sillpar = 1, 
    nuggetpar = 0.01, 
    seed = reps
  )


initial=
  landscape%>%
  mutate(pres=0)

initial$pres[sample(seq(nrow(landscape)),abun)]=1



disps=c(6,10,20,50,100)
methods=c("neutral","positive","negative")


params=expand.grid(disps=disps,
                   methods=methods
)

saveRDS(tibble(reps=reps,
               rangepars=rangepars,
               initial%>%filter(pres==1)%>%select(x,y,pres)),
        paste0("initial_",ind,".rds")
)

cores=detectCores()
plan(multisession,workers=cores)


result=params%>%
  future_pmap_dfr(
    .f=function(
    disps,
    methods){
      

      popmat=initial%>% 
        select(x,y,pres)%>%
        pivot_wider(names_from = x, values_from = pres)%>%
        select(-y)%>%
        as.matrix()
      
      maxnb=((disps-1)^2-1)
      
      neigh=log(seq(maxnb),base=maxnb)
      
      habmat=initial%>%
        select(x,y,soiltype)%>%
        pivot_wider(names_from = x, values_from = soiltype)%>%
        select(-y)%>%  
        as.matrix()
      
      
      nbirth=ndeath=floor(abun*0.1)
      
      samples=floor(seq(1,10000,length.out=15))
      res=NULL
      
      for(i in 1:timepoints){
      
        
        samp=which(popmat==1)
        samp=sample(samp,ndeath)
        
        popmat[samp]=0
        
        j=0
        
        
        while(j<nbirth){
          
          #Randomly select a cell
          
          samp=sample(1:(Lx*Ly),1)
          
          if(popmat[samp]==0){ #If it is empty
            
            x=samp%%Ly
            
            x=ifelse(x==0,Ly,x)
            
            y=((samp-x)/Ly)+1
            
            nb=sum(matrix_get_neighbors(popmat,c(x,y),disps-1))
            
              if(methods=="neutral"){
                disp.prob=ifelse(nb>0,1,0)
              }else{
                if(methods=="positive"){
                  disp.prob=ifelse(nb>0,0.5+neigh[nb],0)
                }else{
                  if(methods=="negative"){
                    disp.prob=ifelse(nb>0,1-(neigh[nb]),0)
                  }
                }
              }
    

            
            if(habmat[x,y]>0 & runif(1)<=disp.prob) { #If it is a suitable cell
              
                popmat[x,y]=1
                
                j=j+1
              
            }
            
          }
        }
        
        
      }
      
      

      dat=popmat%>%
        as_tibble()%>%
        mutate(y=1:nrow(popmat))%>%
        pivot_longer(cols="1":"1000",names_to = "x",values_to = "pres")%>%
        mutate(x=as.numeric(x),
               y=as.numeric(y))%>%
        filter(pres==1)
      
      
      return(tibble(
        disp=disps,
        rangepar=rangepars,
        method=methods,
        rep=reps,
        dat
      ))
      
    },
    .options = furrr_options(seed = TRUE)
    
  )


saveRDS(result,paste0("percsim_",ind,".rds"))



result%>%ggplot(aes(x,y))+geom_point()+facet_wrap(disp~method)

perc.up=function(dat,d_cutoff=seq(1,100,by=5)){
  
  a=tibble(expand.grid(x=1:Lx,y=1:Ly))
  
  dat=dat%>%
    as_tibble()%>%
    right_join(a)%>%
    mutate(pres=replace_na(pres,0))%>%
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


pars=result%>%select(disp,method)%>%unique()

plan(multisession,workers=cores)



pars=result%>%select(disp,method)%>%unique()%>%rename(disp1=disp,method1=method)

plan(multisession,workers=cores)
plot_dat=pars%>%
        future_pmap_dfr(
          .f=function(disp1,method1){
            
            dat=result%>%
              filter(disp==disp1,method==method1)
            
            res1=perc.up(dat)
            
            res1=res1%>%
              group_by(d_cutoff)%>%
              summarize(freq=max(freq))%>%
              ungroup()
            
            return(tibble(disp=disp1,
                   method=method1,
                   res1))
          }
        )

plot_dat%>%ggplot(aes(d_cutoff,freq,col=method))+geom_line()+facet_wrap(vars(disp))


init=readRDS('initial_1.rds')
init.res=perc.up(init%>%select(x,y,pres))%>%ggplot(aes(d_cutoff,freq))+geom_line()

init.res%>%group_by(d_cutoff)%>%summarize(freq=max(freq))%>%
  ggplot(aes(d_cutoff,freq))+geom_line()
