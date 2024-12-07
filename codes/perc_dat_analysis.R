library(tidyverse)
library(RANN)
library(furrr)
library(magrittr)
library(parallel)


theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)



#Generate landscapes
############################################################################

generate_landscape = function(
    Lx,
    Ly,
    quadrat_length,
    rangepar, 
    sillpar = 1, 
    nuggetpar = 0.01,
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

#Params to vary

reps=1:100
rangepars=c(1,8,24)
pars=expand.grid(reps=reps,rangepars=rangepars)

#Fixed params

quadrat_length=1  
sillpar = 1       #Set as default in a function
nuggetpar = 0.01  #Set as default in a function

wd1=getwd()

setwd("C:/Users/Mihir/Documents/landscapes128")

pars%>%
  future_pmap_dfr(
    .f=function(
    reps,
    rangepars){
      
      landscape =  
        generate_landscape(
          Lx = Lx, 
          Ly = Ly, 
          quadrat_length = 1, 
          rangepar = rangepars, 
          sillpar = 1, 
          nuggetpar = 0.01, 
          seed = reps
        )%>%mutate(rep=reps,
                   rangepar=rangepars)
      
      write.csv(landscape,paste0("land_",rangepars,"_",reps,".csv"))
      
    }
  )

setwd(wd1)

###############################################################################################

#Sample simulation



divisors=function(x){
  y=seq_len(x)
  y[ x%%y ==0 ]
}

fisher.fun=function(a,b){
  a=log(a,base=2)
  b=log(b+1,base=2)
  
  res=lm(b~a)
  
  return(res$coefficients[2])
  
}

std.error <- function(x) sd(x)/sqrt(length(x))


d_cutoffs=divisors(min(Lx,Ly))

d_cutoffs=sort(union(seq(2,50,by=2),d_cutoffs))

perc.up=function(dat,d_cutoff=d_cutoffs,unit){#dat is a full grid of 0s and 1s.
  
  options(dplyr.summarise.inform = FALSE)

   coords=expand.grid(x=1:Lx,y=1:Ly)
  
  dat=dat%>%
    as_tibble()%>%
    select(x,y,pres)%>%
    right_join(coords,by=c('x','y'))%>%
    replace_na(list(pres=0))
  
  
  result=NULL    
  
  for(i in d_cutoff){
    
    cutx=unique(c(seq(0, Lx, i),Lx))
    cuty=unique(c(seq(0, Ly, i),Ly))
    
    if(unit=='tree'){
      
      dat1=dat%>%
        mutate(x = cut(x, cutx, labels = FALSE,include.lowest=TRUE),
               y = cut(y, cuty, labels = FALSE,include.lowest=TRUE))%>%
        group_by(x,y)%>%
        summarize(pres=(sum(pres)))%>%
        ungroup()%>%
        pivot_wider(names_from = x,values_from = pres)%>%
        select(-y)%>%
        as.matrix()
    }
    
    if(unit=='cell'){
      
      dat1=dat%>%
        mutate(x = cut(x, cutx, labels = FALSE,include.lowest=TRUE),
               y = cut(y, cuty, labels = FALSE,include.lowest=TRUE))%>%
        group_by(x,y)%>%
        summarize(pres=sign(sum(pres)))%>%
        ungroup()%>%
        pivot_wider(names_from = x,values_from = pres)%>%
        select(-y)%>%
        as.matrix()
    }
    
      
      clust=matrix(0,nrow=nrow(dat1),ncol=ncol(dat1))
      
      largest=0
      
      lbl=0:length(clust)
      
      for(y in seq(ncol(dat1))){
        for(x in seq(nrow(dat1))){
          
          if((dat1[x,y]>0)){
            
            left=ifelse(length(clust[x,y-1])==0,0,clust[x,y-1])
            above=ifelse(length(clust[x-1,y])==0,0,clust[x-1,y])
            
            if((left==0) & (above==0)){
              
              largest=largest+1
              clust[x,y]=largest
              lbl[(nrow(clust)*(y-1))+x]=largest
              
            }else{
              
              if((left!=0) & (above==0)){
                
                clust[x,y]=clust[x,y-1]
                
                lbl[(nrow(clust)*(y-1))+x]=lbl[(nrow(clust)*(y-2))+x]
                
              }else{
                
                if(left==0 & above!=0){
                  
                  clust[x,y]=clust[x-1,y]
                  
                  lbl[(nrow(clust)*(y-1))+x]=lbl[(nrow(clust)*(y-1))+x-1]
                  
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
     
      dat2=dat1%>%
        as_tibble()%>%
        mutate(y=1:nrow(dat1))%>%
        pivot_longer("1":as.character(ncol(dat1)),names_to="x",values_to = "pres")
      
      colnames(clust)=1:ncol(clust)
      
      clust2=clust%>%
        as_tibble()%>%
        mutate(y=1:nrow(clust))%>%
        pivot_longer("1":as.character(ncol(clust)),names_to="x",values_to = "cluster")
      
      
      clustab=dat2%>%
        inner_join(clust2,by=c("x","y"))%>%
        filter(pres!=0)%>%
        select(pres,cluster)%>%
        group_by(cluster)%>%
        summarize(freq=sum(pres))%>%
        ungroup()%>%
        count(freq)
      
      result=bind_rows(
        result,
        tibble(
          d_cutoff=i,
          clust=clustab$freq,
          freq=clustab$n
        )
      )
       
    }
  
  return(result)
  }


setwd("C:/Users/mihir/Documents/landscapes128")

files=list.files("C:/Users/mihir/Documents/landscapes128")

landscape.files=files[grep("landscape",files)]
final.files=files[grep("final",files)]


sets=1:length(final.files)

cores=detectCores()
plan(multisession,workers=cores)

results=NULL

for(j in 1:length(sets)){
      
      finset=read.csv(paste0("C:/Users/mihir/Documents/landscapes128/",final.files[j]))%>%
            as_tibble()%>%
            mutate(comp_dist=rep(rep(c(3,6),each=400),3))
      
      rep=finset$rep[1]
      rangepar=finset$rangepar[1]
      
      params=finset%>%count(disp,method,comp_dist)%>%select(-n)%>%as.data.frame()
      
      fin.res=NULL
      
      for(i in 1:nrow(params)){
        
        final=finset%>%
          filter(disp==params[i,1],
                 method==params[i,2],
                 comp_dist==params[i,3])%>%
          select(x,y,pres)
        
        final.res=perc.up(final,unit='tree')%>%
          mutate(time=750,
                 rep=rep,
                 rangepar=rangepar,
                 disp=params[i,1],
                 method=params[i,2],
                 comp_dist=params[i,3])
        
        results=bind_rows(
          results,final.res
        )
      }
}
      

write.csv(results,"perc.res128.csv")
  result=read.csv("perc_res128.csv")



sets=1:length(final.files)

cores=detectCores()
plan(multisession,workers=cores)

results=sets%>%
  future_map_dfr(
    .f=function(sets){
      
      initial=read.csv(paste0("C:/Users/mihir/Documents/landscapes/","initial_",114,".csv"))
      
      rep=initial[1,1]
      rangepar=initial[1,2]
      
      init.res=perc.up(initial%>%
                         select(x,y,pres))%>%
        mutate(time=0,
               rep=rep,
               rangepar=rangepar)%>%
        select(rep,rangepar,time,d_cutoff,clust,freq)
      
      
      return(init.res)
    })

saveRDS(results,'perc.res.rds')


results=read.csv('perc_res128.csv')%>%as_tibble()

#Summarize the results

pars=results%>%
      select(rep,rangepar,disp,method,comp_dist)%>%
      unique()%>%
      as.data.frame()




clusdat=NULL
for(i in 1:nrow(pars)){
  
  dat=results%>%
      filter(rep==pars[i,1],
             rangepar==pars[i,2],
             disp==pars[i,3],
             method==pars[i,4],
             comp_dist==pars[i,5])%>%
      select(d_cutoff,clust,freq)
  
  
  dat.fisher1=as_tibble(expand.grid(d_cutoff=d_cutoffs,clust=1:100))%>%
             full_join(dat)%>%
              replace_na(list(freq=0))%>%
              group_by(d_cutoff)%>%
              mutate(cumfreq=cumsum(freq))%>%
              ungroup()%>%
              arrange(d_cutoff)%>%
              filter(clust%in%c(2,4,8,16,32,64,100))%>%
              group_by(d_cutoff)%>%
              mutate(freq1=c(cumfreq[1],diff(cumfreq)))%>%
              ungroup()%>%
              group_by(d_cutoff)%>%
              summarize(f1=fisher.fun(clust,freq1))%>%
              ungroup()
  
  dat.fisher=as_tibble(expand.grid(d_cutoff=d_cutoffs,clust=1:100))%>%
            full_join(dat)%>%
            replace_na(list(freq=0))%>%
            arrange(d_cutoff)%>%
            group_by(d_cutoff)%>%
            summarize(f=fisher.fun(clust,freq))%>%
            ungroup()%>%
            inner_join(dat.fisher1)
  
  clusmax=dat%>%
          #mutate(clust=clust*(d_cutoff^2)/(Lx*Ly))%>%
          group_by(d_cutoff)%>%
          slice_max(clust)%>%
          ungroup()%>%
          rename(clusmax=clust)%>%
          select(-freq)
  
  #x=clusmax$d_cutoff
  #y=clusmax$clust
  
  #fit=nls(y~1/(1+exp((xmid-x)/scal)),data.frame(x,y),start=list(xmid=30,scal=5))
  
fun1=function(i){
  
  res=NULL
  
  for(j in 1:length(i)){
  
  if(i[j]==2){res=c(res,"Facilitation")}
  if(i[j]==3){res=c(res,"Competition")}
  }
    
    return(res)
}
  
  
  
  clusdat=bind_rows(
                clusdat,
                dat.fisher%>%
                  inner_join(clusmax)%>%
                  mutate(rep=pars[i,1],
                         rangepar=pars[i,2],
                         disp=pars[i,3],
                         method=pars[i,4],
                         comp_dist=pars[i,5]))

}

write.csv(clusdat,"clusdat128.csv")
clusdat=read.csv("clusdat128.csv")


clusdat.sum=clusdat%>%
  group_by(rangepar,disp,method,comp_dist,d_cutoff)%>%
  summarize(f.m=mean(f),
            f.sd=std.error(f),
            f1.m=mean(f1),
            f1.sd=std.error(f1),
            clus.m=mean(clusmax),
            clus.sd=std.error(clusmax))%>%
  ungroup()

clusdat.sum%>%
  mutate(method=as.factor(method))%>%
  filter(rangepar==2)%>%
  ggplot(aes(d_cutoff,f.m,col=method))+
  geom_line()+
  facet_grid(disp~comp_dist)


#Calculate the inflection points
#Use piecewise regresssion model to get two inflection points
pardat=clusdat%>%select(rep,rangepar,disp,comp_dist,method)%>%unique()

infl.dat=NULL

library(segmented)

for(i in 1:nrow(pardat)){
  
  samps=clusdat%>%
        filter(rep==pardat[i,1],
               rangepar==pardat[i,2],
               disp==pardat[i,3],
               comp_dist==pardat[i,4],
               method==pardat[i,5])
  
  x=samps$d_cutoff
  y=samps$f1
  z=samps$clusmax
  
  #fill data
  x=c(x,x+c(diff(x)/2,0))
  y=c(y,y+c(diff(y)/2,0))
  z=c(z,z+c(diff(z)/2,0))
  
  lm.y=lm(y~x)
  seg.y=summary(segmented(lm.y,seg.Z=~x,npsi=2))$psi[,2]
  
  lm.z=lm(z~x)
  seg.z=summary(segmented(lm.z,seg.Z=~x,npsi=2))$psi[,2]
  
  
  infl.dat=bind_rows(infl.dat,
                     tibble(
                       rep=pardat[i,1],
                       rangepar=pardat[i,2],
                       disp=pardat[i,3],
                       comp_dist=pardat[i,4],
                       method=pardat[i,5],
                       #f1=seg.y[1],
                       #cm1=seg.z[1],
                       f2.1=seg.y[1],
                       f2.2=seg.y[2],
                       cm2.1=seg.z[1],
                       cm2.2=seg.z[2]
                     ))
}


write.csv(infl.dat,"infl1.csv")
infl1.dat=read.csv("infl1.csv")

  
infl1.sum=infl1.dat%>%
          group_by(rangepar,disp,comp_dist,method)%>%
          summarize(f=mean(f1),
                    f.sd=std.error(f1),
                    cm=mean(cm1),
                    cm.sd=std.error(cm1))

infl1.sum%>%
  mutate(method=as.factor(method),
         comp_dist=as.factor(comp_dist),
         rangepar=as.factor(rangepar))%>%
  ggplot(aes(disp,f,col=rangepar,fill=rangepar))+
  geom_line()+
  #geom_ribbon(aes(ymax=f+f.sd,ymin=f-f.sd),alpha=0.2)+
  facet_grid(method~comp_dist)+
  ylab("Inflection point for Fisher's tau")+
  xlab("Maximum Dispersal range")+
  ggtitle("Column panels: Range of interaction; Row panels:Mode of interaction")


plotdat=infl1.sum%>%
  filter(comp_dist==3)%>%
  mutate(method=as.factor(method),
         rangepar=as.factor(rangepar))

plotdat$method=recode(plotdat$method,"1"="Neutral","2"="Positive","3"="Negative")

  plotdat%>%ggplot(aes(disp,cm,col=rangepar,fill=rangepar))+
  geom_line(size=2)+
  #geom_ribbon(aes(ymax=f+f.sd,ymin=f-f.sd),alpha=0.2)+
  facet_wrap(vars(method))+
  ylab("Inflection point for the biggest cluster size")+
  xlab("Maximum Dispersal range")

infl1.sum%>%
  mutate(method=as.factor(method),
         comp_dist=as.factor(comp_dist))%>%
  ggplot(aes(disp,cm,col=method))+
  geom_line()+
  facet_grid(rangepar~comp_dist)+
  ylab("Inflection point for the biggest cluster size")+
  xlab("Maximum Dispersal range")+
  ggtitle("Column panels: Range of interaction; Row panels:Scale of habitat spatial autocorrelation")

infl1.dat%>%
  mutate(method=as.factor(method),
         rangepar=as.factor(rangepar),
         disp=as.factor(disp),
         comp_dist=as.factor(comp_dist))%>%
  ggplot(aes(f1,cm1,col=disp))+
  geom_point()

infl2.sum=infl1.dat%>%
  mutate(dist=disp/comp_dist)%>%
  group_by(rangepar,dist,method)%>%
  summarize(f=mean(f1),
            f.sd=std.error(f1),
            cm=mean(cm1),
            cm.sd=std.error(cm1))

infl2.sum%>%
  mutate(method=as.factor(method),
         rangepar=as.factor(rangepar))%>%
  ggplot(aes(dist,f,col=rangepar))+
  geom_line()+
  geom_errorbar(aes(ymin=f-f.sd,ymax=f+f.sd))+
  facet_wrap(vars(method))+
  xlab("Dispersal range/Interaction range")

infl2.sum%>%
  mutate(method=as.factor(method),
         rangepar=as.factor(rangepar))%>%
  ggplot(aes(dist,cm,col=rangepar))+
  geom_line()+
  #geom_errorbar(aes(ymin=f-f.sd,ymax=f+f.sd))+
  facet_wrap(vars(method))+
  xlab("Dispersal range/Interaction range")

#Can the summary statistics inform the param values?

#cursory C5.0 analysis 
#Test if the algorithm can identify different values of individual parameters (one at a time)



######################################################
#Apply the percolation analysis to BCI data

perc.up1=function(dat,d_cutoff=d_cutoffs,unit){#dat is a full grid of 0s and 1s.
  
  options(dplyr.summarise.inform = FALSE)
  
  coords=expand.grid(x=1:Lx,y=1:Ly)
  
  dat=dat%>%
    as_tibble()%>%
    select(x,y,pres)%>%
    #mutate(x=ceiling(x),y=ceiling(y))%>%
    right_join(coords,by=c('x','y'))%>%
    replace_na(list(pres=0))
  
  
  result=NULL    
  
  for(i in d_cutoff){
    
    cutx=unique(c(seq(0, Lx, i),Lx))
    cuty=unique(c(seq(0, Ly, i),Ly))
    
    if(unit=='tree'){
      
      dat1=dat%>%
        mutate(x = cut(x, cutx, labels = FALSE,include.lowest=TRUE),
               y = cut(y, cuty, labels = FALSE,include.lowest=TRUE))%>%
        group_by(x,y)%>%
        summarize(pres=(sum(pres)))%>%
        ungroup()%>%
        pivot_wider(names_from = x,values_from = pres)%>%
        select(-y)%>%
        as.matrix()
    }
    
    if(unit=='cell'){
      
      dat1=dat%>%
        mutate(x = cut(x, cutx, labels = FALSE,include.lowest=TRUE),
               y = cut(y, cuty, labels = FALSE,include.lowest=TRUE))%>%
        group_by(x,y)%>%
        summarize(pres=sign(sum(pres)))%>%
        ungroup()%>%
        pivot_wider(names_from = x,values_from = pres)%>%
        select(-y)%>%
        as.matrix()
    }
    
    
    clust=matrix(0,nrow=nrow(dat1),ncol=ncol(dat1))
    
    largest=0
    
    lbl=0:length(clust)
    
    for(y in seq(ncol(dat1))){
      for(x in seq(nrow(dat1))){
        
        if((dat1[x,y]>0)){
          
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
    
    dat2=dat1%>%
      as_tibble()%>%
      mutate(y=1:nrow(dat1))%>%
      pivot_longer("1":as.character(ncol(dat1)),names_to="x",values_to = "pres")
    
    colnames(clust)=1:ncol(clust)
    
    clust2=clust%>%
      as_tibble()%>%
      mutate(y=1:nrow(clust))%>%
      pivot_longer("1":as.character(ncol(clust)),names_to="x",values_to = "cluster")
    
    
    clustab=dat2%>%
      inner_join(clust2,by=c("x","y"))%>%
      filter(pres!=0)%>%
      select(pres,cluster)%>%
      group_by(cluster)%>%
      summarize(freq=sum(pres))%>%
      ungroup()%>%
      count(freq)
    
    result=bind_rows(
      result,
      tibble(
        d_cutoff=i,
        clust=clustab$freq,
        freq=clustab$n
      )
    )
    
  }
  
  return(result)
}

bci_spnames=read_tsv('C:/Users/Mihir/Downloads/Condit_ValidBCITaxa.tsv')%>%
            select(Latin,Mnemonic)%>%
            rename(sp=Mnemonic)


mylist=c('Alchornea costaricensis',"Alseis blackiana",
         "Casearia arborea","Cecropia insignis" ,"Cordia alliodora" ,
         "Croton billbergianus","Jacaranda copaia" ,"Luehea seemannii",
         "Miconia argentea","Palicourea guianensis" ,"Terminalia amazonia",
         "Zanthoxylum acuminatum")


splist=bci_spnames%>%
  filter(Latin%in%mylist)%>%
  select(sp)%>%pull()


bci=readRDS("bci_raw.rds")%>%
  select(census,treeID,sp,status,dbh,gx,gy)%>%
  filter(sp%in%splist,status %in% c("A","D"))

demodat=NULL

for(i in 1:length(splist)){
  
  for(j in 2:8){
    
    prevdat=bci%>%
            filter(sp==splist[i],
                   census==j-1,
                   status=='A')
    
    currdat=bci%>%
      filter(sp==splist[i],
             census==j)
    
    deaths=currdat%>%filter(status=='D')%>%nrow()
    
    births=nrow(currdat)-length(intersect(currdat$treeID,prevdat$treeID))
      
    demodat=bind_rows(
            demodat,
            tibble(
              sp=splist[i],
              census=j-1,
              births=births,
              deaths=deaths
            )
    )
    
  }
}

demodat=demodat%>%
  mutate(births=births/10,
         deaths=deaths/10)%>%
  group_by(sp)%>%
  summarize(
    births.m=ceiling(mean(births)),
    deaths.m=ceiling(mean(deaths)),
    births.sd=ceiling(sd(births)),
    deaths.sd=ceiling(sd(deaths)))%>%
  ungroup()
  

write.csv(demodat,'perc.demodata.bci.csv')


bci_spnames%<>%
  filter(Latin%in%mylist) 


dat=bci%>%
    filter(sp%in%splist,
           census==8)

dat=dat%>%
    inner_join(
      dat%>%
      group_by(sp)%>%
      summarize(baldeck=quantile(dbh,0.56,na.rm=TRUE))%>%
      ungroup()
    )%>%
    filter(dbh>baldeck)%>%
    select(sp,status,gx,gy)
    


#Plot distributions of each species
dat%>%
  ggplot(aes(gx,gy))+
  geom_point()+
  facet_wrap(vars(sp))+
  theme(aspect.ratio = 0.5)



Lx=1000
Ly=500

d_cutoffs=divisors(min(Lx,Ly))
d_cutoffs=sort(union(seq(5,75,by=5),d_cutoffs))

#Percolation curves for the observed empirical data
bci_res=NULL

for(i in splist){
  
  dat1=dat%>%filter(sp==i)%>%
        mutate(pres=1)%>%
        rename(x=gx,y=gy)%>%
        mutate(x=ceiling(x),y=ceiling(y))%>%
        select(x,y,pres)
  
  res=perc.up1(dat1,d_cutoffs,unit="tree")
  
  bci_res=bind_rows(bci_res,
                    res%>%
                      mutate(sp=i))
}

write.csv(bci_res,"bci_perc.csv")

bci_res=read.csv("bci_perc.csv")%>%
        as_tibble()


#Percolation analysis for the simulated species data under different mechanisms
files=list.files("C:/Users/Mihir/Documents/landscapes")
files=files[-grep('landscape',files)]
files=files[grep('.csv',files)]


popdat=tibble(
  sp=c("alchco", "alsebl", "casear", "cecrin", "cordal", "crotbi", "jac1co", 
       "luehse", "micoar", "paligu", "termam", "zantpr"),
  pop=c(271,9913,108,1381,223,634,309,199,652,1326,58,72))


bci_perc_null=NULL

for(i in 1:length(files)){

  dat=read.csv(paste0('C:/Users/Mihir/Documents/landscapes/',files[i]))%>%
      as_tibble()
  
  pop=popdat%>%filter(sp==dat$sp[1])%>%pull(pop)
  
  dat=dat%>%
      mutate(comp_dist=rep(rep(c(3,6),each=pop),30))
  
  pars=dat%>%
      select(rep,rangepar,method,comp_dist)%>%
      unique()
  
  for(j in 1:nrow(pars)){
    
    dat1=dat%>%filter(rep==pars$rep[j],
                      rangepar==pars$rangepar[j],
                      method==pars$method[j],
                      comp_dist==pars$comp_dist[j])%>%
        select(x,y,pres)
    
    res=perc.up1(dat1,d_cutoffs,unit="tree")
    
    bci_perc_null=bind_rows(
      bci_perc_null,
      tibble(
        res%>%
          mutate(
            sp=dat$sp[i],
            rep=pars$rep[j],
            rangepar=pars$rangepar[j],
            method=pars$method[j],
            comp_dist=pars$comp_dist[j]
            )
      )
    )
  }

}

write.csv(bci_perc_null,'bci_perc_null.csv')

bci_null_res=read.csv('bci_perc_null.csv')%>%as_tibble()

#Calculate the inflection points
#Use piecewise regresssion model to get two inflection points

splist=intersect(bci_res$sp,popdat$sp)

infl_bci =NULL

library(segmented)

for(i in splist){
  
  samps=bci_res%>%
        filter(sp==i)
  
  pop=samps%>%
      slice_min(d_cutoff)%>%
      summarize(sum(clust*freq))%>%
      pull()
  
  dat.fisher=as_tibble(expand.grid(d_cutoff=d_cutoffs,clust=1:pop))%>%
    full_join(samp1)%>%
    replace_na(list(freq=0))%>%
    arrange(d_cutoff)%>%
    group_by(d_cutoff)%>%
    summarize(f=fisher.fun(clust,freq))%>%
    ungroup()%>%
    filter(d_cutoff%in%unique(samp1$d_cutoff))
  
  clusmax=samps%>%
    group_by(d_cutoff)%>%
    slice_max(clust)%>%
    ungroup()%>%
    rename(clusmax=clust)
  
  x=dat.fisher$d_cutoff
  y=dat.fisher$f
  z=clusmax$clusmax
  
  
  lm.y=lm(y~x)
  seg.y=summary(segmented(lm.y,seg.Z=~x,npsi=1))$psi[,2]
  
  lm.z=lm(z~x)
  seg.z=summary(segmented(lm.z,seg.Z=~x,npsi=1))$psi[,2]
  
  
  infl_bci=bind_rows(infl_bci,
                     tibble(sp=i,
                            pop=pop,
                       f1=seg.y,
                       cm1=seg.z,
                       #f2.1=seg.y[1],
                       #f2.2=seg.y[2],
                       #cm2.1=seg.z[1],
                       #cm2.2=seg.z[2]
                     ))
}


write.csv(infl_bci,"infl_bci.csv")
infl_bci=read.csv("infl_bci.csv")%>%
        filter(sp%in%splist)

#Get the inflection points for the null data
infl_bci_null =NULL

pars=bci_null_res%>%
  count(rangepar,method,comp_dist)%>%
  select(rangepar,method,comp_dist)

for(i in unique(bci_null_res$sp)){
  
  samps=bci_null_res%>%
    filter(sp==i)
  
  pop=popdat%>%filter(sp==i)%>%pull(pop)
  
  for(j in 1:nrow(pars)){
    
    samp1=samps%>%
          filter(rangepar==pars$rangepar[j],
                 method==pars$method[j],
                 comp_dist==pars$comp_dist[j])
    
    
    dat.fisher=as_tibble(expand.grid(d_cutoff=d_cutoffs,clust=1:pop))%>%
      full_join(samp1)%>%
      replace_na(list(freq=0))%>%
      arrange(d_cutoff)%>%
      group_by(d_cutoff)%>%
      summarize(f=fisher.fun(clust,freq))%>%
      ungroup()%>%
      filter(d_cutoff%in%unique(samp1$d_cutoff))
    
    clusmax=samp1%>%
      group_by(d_cutoff)%>%
      slice_max(clust)%>%
      ungroup()%>%
      rename(clusmax=clust)
    
    x=dat.fisher$d_cutoff
    y=dat.fisher$f
    z=clusmax$clusmax
    
    
    lm.y=lm(y~x)
    seg.y=summary(segmented(lm.y,seg.Z=~x,npsi=1))$psi[,2]
    
    lm.z=lm(z~x)
    seg.z=summary(segmented(lm.z,seg.Z=~x,npsi=1))$psi[,2]
    
    infl_bci_null=bind_rows(infl_bci_null,
                       tibble(sp=i,
                              rangepar=pars$rangepar[j],
                              method=pars$method[j],
                              comp_dist=pars$comp_dist[j],
                              pop=pop,
                              f1=seg.y,
                              cm1=seg.z
                              #f2.1=seg.y[1],
                              #f2.2=seg.y[2],
                              #cm2.1=seg.z[1],
                              #cm2.2=seg.z[2]
                       ))
    
  }
  
}


write.csv(infl_bci_null,"infl_bci_null.csv")
infl_bci_null=read.csv("infl_bci_null.csv")%>%
  mutate(treat='real')
  
infl_bci=
  infl_bci%>%
  rename("f1.st"=f1,
         "cm1.st"=cm1)

infl_bci_null=infl_bci_null%>%
  inner_join(infl_bci%>%select(sp,f1.st,cm1.st))%>%
  mutate(f1=f1.st-f1,
         cm1=cm1.st-cm1)



infl_bci_null%>%
  filter(comp_dist==6)%>%
  mutate(method=as.factor(method))%>%
  ggplot(aes(rangepar,f1,col=method))+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  facet_wrap(vars(sp))

infl_bci_null%>%
  filter(comp_dist==6)%>%
  mutate(method=as.factor(method))%>%
  ggplot(aes(rangepar,cm1,col=method))+
  geom_point()+
  geom_hline(aes(yintercept=0))+
  facet_wrap(vars(sp))


infl1.sum=infl1.dat%>%
  group_by(rangepar,disp,comp_dist,method)%>%
  summarize(f=mean(f1),
            f.sd=std.error(f1),
            cm=mean(cm1),
            cm.sd=std.error(cm1))

infl1.sum%>%
  mutate(method=as.factor(method),
         comp_dist=as.factor(comp_dist),
         rangepar=as.factor(rangepar))%>%
  ggplot(aes(disp,f,col=rangepar))+
  geom_line()+
  facet_grid(method~comp_dist)+
  ylab("Inflection point for Fisher's tau")+
  xlab("Maximum Dispersal range")+
  ggtitle("Column panels: Range of interaction; Row panels:Mode of interaction")

infl1.sum%>%
  mutate(method=as.factor(method),
         comp_dist=as.factor(comp_dist))%>%
  ggplot(aes(disp,cm,col=method))+
  geom_line()+
  facet_grid(rangepar~comp_dist)+
  ylab("Inflection point for the biggest cluster size")+
  xlab("Maximum Dispersal range")+
  ggtitle("Column panels: Range of interaction; Row panels:Scale of habitat spatial autocorrelation")


#Test if the summary statistics strongly associate with param combinations
infl2=infl1.dat%>%
      pivot_longer(cols=c('rangepar','disp','method','comp_dist'),
                   names_to='param',
                   values_to='vals')%>% 
      pivot_longer(cols=c('f1','cm1'),
                   names_to='res',
                   values_to='values')%>%
      group_by(param)%>%
      summarize()

infl2%>%
  ggplot(aes(x=vals))+
  geom_point(aes(vals,values,col=res))+
  facet_grid(param~res)
  


testdat=infl1.sum%>%
  as_tibble()%>%
  pivot_longer(
    cols=c('rangepar','disp','comp_dist','method'),
           names_to='param',
           values_to='value')
  

testdat%>%
  filter(param=="rangepar")%>%
  mutate(f=as.factor(f),
         cm=as.factor(cm),
         value=as.numeric(value))%>%
  ggplot(aes(f,cm,col='value',fill='value'))+
  geom_tile()

library(C50)
library(caret)

c5dat=infl1.dat%>%
      filter(rangepar%in%c(5,20),
             disp%in%c(3,10))%>%
      mutate(pars=paste0(rangepar,"_",disp,"-",comp_dist,"_",method))%>%
      mutate(pars=as.factor(pars))

rp_res=C5.0(c5dat%>%select(f1,cm1),
     c5dat$pars,
     rule=TRUE,
     control=C5.0Control(winnow=TRUE),
     trials=1)


c5dat.rp=infl1.dat%>%
      select(rangepar,f1,cm1)

c5.rangepar=train(
  pars~.,
  data=c5dat%>%select(pars,f1,cm1),
  method="C5.0",
  trcontrol=trainControl(method='repeatedcv',repeats=10),
  metric='Kappa'
)


C5_model = 
  train(
    cata ~ ., 
    data = c5dat, 
    method = 'C5.0',
    trControl = trainControl(method = 'repeatedcv', repeats = 10),
    metric = 'Kappa'
  )

  
#######################################################
results=as_tibble(read.csv("clustdat.csv"))
results=results%>%
        mutate(total=cluster*counts)%>%
        group_by(rep,rangepar,disp,method,d_cutoff)%>%
        mutate(clust=cluster/sum(total))%>%
        ungroup()

pars=results%>%
    select(rep,rangepar,disp,d_cutoff)%>%
    unique()

pars%>%
  future_pmap_dfr(
    .f=function(reps,rangepars,disps,d_cutoffs){
      
      dat=results%>%
          filter(rep==reps,
                 rangepar==rangepars,
                 disp==disps,
                 d_cutoff==d_cutoffs)
      
      samp=sample(dat$clust,sum(dat$total),prob=dat$counts,replace=TRUE)
      
      fit=MASS::fitdistr(samp,"beta",list(shape1=1,shape2=3))
      
      return(tibble(
        rep=reps,
        rangepar=rangepars,
        disp=disps,
        d_cutoff=d_cutoffs,
        shape1=fit$estimate[1],
        shape2=fit$estimate[2]
      ))
      
    }
  )

####################################################################

files=list.files("C:/Users/mihir/Documents/landscapes2")
finals=files[grep("final",files)]
Lx=128
Ly=128
abun=100
#Get null expectations
cores=detectCores()
plan(multisession,workers=cores)

fileno=1:300

nullres=fileno%>%
  future_map_dfr(
    .f=function(fileno){
      
      land=read.csv(paste0("C:/Users/mihir/Documents/landscapes/landscape_",fileno,".csv"))%>%
        as_tibble()%>%
        filter(x<129,y<129)
      
      mid=quantile(land$soiltype,0.5)
      rangepar=land$rangepars[1]
      
      land=land%>%
            filter(soiltype>mid)%>%
            mutate(pres=1)%>%
            select(x,y,pres)
      
      res1=NULL
      
      for(i in 1:5){
        land1=land%>%
              slice(sample((1:nrow(land)),abun))
        
        res.tree=perc.up(land1,unit="tree")
        
        res1=bind_rows(res1,
                       tibble(rep=i,
                              rangepar=rangepar,
                              d_cutoff=res.tree$d_cutoff,
                              clust=res.tree$clust,
                              freq=res.tree$freq))
        
      }
      
      return(res1)
      
    }
  )

write.csv(nullres,"clustdat_null.csv")

nullres=read.csv("clustdat_null.csv")%>%as_tibble()

nullres_sum=nullres%>%
  group_by(rep,rangepar,d_cutoff)%>%
  summarize(means=sum((clust*freq))/sum(freq),
            max=max(clust),
            mode=clust[which.max(freq)],
            peak=max(freq))%>%
  ungroup()


min_res=function(data,par){
  with(data,sum((y-(1/(1+exp(-(x-par[1])/par[2]))))^2))
}



meandat.null=NULL

parms=nullres_sum%>%select(rep,rangepar)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=nullres_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i])%>%
    select(d_cutoff,means)%>%
    rename(x=d_cutoff,y=means)%>%
    mutate(y=y/abun)%>%
    as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  meandat.null=bind_rows(meandat.null,
                    tibble(res=parms$rep[i],
                           rangepar=parms$rangepar[i],
                           mid=foo.res$par[1],
                           scale=foo.res$par[2]))
}



maxdat.null=NULL

parms=results_sum%>%select(rep,rangepar)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i])%>%
    select(d_cutoff,max)%>%
    rename(x=d_cutoff,y=max)%>%
    mutate(y=y/abun)%>%
    as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  maxdat.null=bind_rows(maxdat.null,
                   tibble(res=parms$rep[i],
                          rangepar=parms$rangepar[i],
                          mid=foo.res$par[1],
                          scale=foo.res$par[2]))
}



  #reps=2
#rangepars=5
#disps=10
#methods='positive' 

cores=detectCores()
plan(multisession,workers=cores)

fileno=1:length(finals)
results_raw=fileno%>%
  future_map_dfr(
    .f=function(fileno){
      
      fin1=read.csv(paste0("C:/Users/mihir/Documents/landscapes2/final_",fileno,".csv"))%>%as_tibble()
      
      pars=fin1%>%select(rep,rangepar,disp,method)%>%unique()
      
      res1=NULL
      
      for(i in 1:nrow(pars)){
        
        dat=fin1%>%
          filter(rep==pars$rep[i],
                 rangepar==pars$rangepar[i],
                 disp==pars$disp[i],
                 method==pars$method[i])
        
        res=perc.up(dat,unit='tree')
        
        res1=bind_rows(res1,
                       tibble(rep=pars$rep[i],
                              rangepar=pars$rangepar[i],
                              disp=pars$disp[i],
                              method=pars$method[i],
                              d_cutoff=res$d_cutoff,
                              clust=res$clust,
                              freq=res$freq))
        
      }
      return(res1)
    }
  )
 write.csv(results_raw,"clustdat_raw_tree.csv")
 
results_raw=read.csv("clustdat_raw_tree.csv")%>%as_tibble()

results_sum=results_raw%>%
        group_by(rep,rangepar,disp,method,d_cutoff)%>%
        summarize(means=sum((clust*freq))/sum(freq),
                  var=sd(rep(clust,freq)),
                  max=max(clust),
                  mode=clust[which.max(freq)],
                  peak=max(freq))%>%
        ungroup()

results_sum1=results_sum%>%
            group_by(rangepar,disp,method,d_cutoff)%>%
            summarize(means.m=mean(means),
                      means.sd=sd(means),
                      var.m=mean(var),
                      var.sd=sd(var),
                      max.m=mean(max),
                      max.sd=sd(max),
                      mode.m=mean(mode),
                      mode.sd=sd(mode),
                      peak.m=mean(peak),
                      peak.sd=sd(peak))%>%
            ungroup()

results_sum1%>%
  ggplot(aes(d_cutoff,means.m,col=method,fill=method))+
  geom_point()+
  #geom_ribbon(aes(ymin=means.m-means.sd,ymax=means.m+means.sd),alpha=0.2)+
  facet_grid(rangepar~disp)


results_sum1%>%
  ggplot(aes(d_cutoff,mode.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=mode.m-mode.sd,ymax=mode.m+mode.sd),alpha=0.2)+
  facet_grid(rangepar~disp)




meandat=NULL

parms=results_sum%>%select(rep,rangepar,disp,method)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i],
           disp==parms$disp[i],
           method==parms$method[i])%>%
    select(d_cutoff,means)%>%
    rename(x=d_cutoff,y=means)%>%
    mutate(y=y/abun)%>%
    as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  meandat=bind_rows(meandat,
                    tibble(res=parms$rep[i],
                           rangepar=parms$rangepar[i],
                           disp=parms$disp[i],
                           method=parms$method[i],
                           mid=foo.res$par[1],
                           scale=foo.res$par[2]))
}
 
meandat=meandat%>%
        rename(mean.mid=mid,mean.scale=scale)

meandat_sum=meandat%>%
            group_by(rangepar,disp,method)%>%
            summarize(mid.m=mean(mid),
                      mid.sd=sd(mid),
                      scale.m=mean(scale),
                      scale.sd=sd(scale))%>%
            ungroup()


meandat_sum%>%
  ggplot(aes(disp,mid.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  #geom_ribbon(aes(ymin=mid.m-mid.sd,ymax=mid.m+mid.sd))+
  facet_wrap(vars(rangepar))
  
meandat_sum%>%
  ggplot(aes(disp,scale.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  #geom_ribbon(aes(ymin=scale.m-scale.sd,ymax=scale.m+scale.sd))+
  facet_wrap(vars(rangepar))

maxdat=NULL

parms=results_sum%>%select(rep,rangepar,disp,method)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i],
           disp==parms$disp[i],
           method==parms$method[i])%>%
    select(d_cutoff,max)%>%
    rename(x=d_cutoff,y=max)%>%
    mutate(y=y/abun)%>%
    as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  maxdat=bind_rows(maxdat,
                    tibble(res=parms$rep[i],
                           rangepar=parms$rangepar[i],
                           disp=parms$disp[i],
                           method=parms$method[i],
                           mid=foo.res$par[1],
                           scale=foo.res$par[2]))
}


maxdat=maxdat%>%
  rename(max.mid=mid,max.scale=scale)
  
modedat=NULL

parms=results_sum%>%select(rep,rangepar,disp,method)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i],
           disp==parms$disp[i],
           method==parms$method[i])%>%
    select(d_cutoff,mode)%>%
    rename(x=d_cutoff,y=mode)%>%
    mutate(y=y/abun)%>%
    as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  modedat=bind_rows(modedat,
                   tibble(res=parms$rep[i],
                          rangepar=parms$rangepar[i],
                          disp=parms$disp[i],
                          method=parms$method[i],
                          mid=foo.res$par[1],
                          scale=foo.res$par[2]))
}

modedat=modedat%>%
  rename(mode.mid=mid,mode.scale=scale)


 results_sum=results_sum%>%
          group_by(rangepar,disp,method,d_cutoff)%>%
          summarize(means.m=mean(means),
                    means.sd=sd(means),
                    max.m=mean(max),
                    max.sd=sd(max),
                    mode.m=mean(mode),
                    mode.sd=sd(mode),
                    peak.m=mean(peak),
                    peak.sd=sd(peak))%>%
          ungroup()
 
 
 
 result_summary=meandat%>%
                inner_join(maxdat)%>%
                inner_join(modedat)
 
 
 c5dat=result_summary%>%
        filter(disp!=5,rangepar!=30,method!="neutral")%>%
        mutate(category=paste0(rangepar,"_",disp,"_",method))
        
refs=c5dat%>%
      select(category)%>%
      unique()%>%
      mutate(cata=seq(nrow(refs)))
 
c5dat=c5dat%>%
      inner_join(refs)%>%
      select(cata,mean.mid:mode.scale)%>%
      mutate(cata=as.factor(cata))%>%
      as.data.frame()

C5_model = 
  train(
    cata ~ ., 
    data = c5dat, 
    method = 'C5.0',
    trControl = trainControl(method = 'repeatedcv', repeats = 10),
    metric = 'Kappa'
  )


results_sum%>%
  ggplot(aes((d_cutoff),means.m,col=method,fill=method))+
  geom_line()+
  #geom_ribbon(aes(ymin=means.m-means.sd,ymax=means.m+means.sd),alpha=0.2)+
  facet_grid(disp~rangepar)


results_sum%>%
  ggplot(aes((d_cutoff),max.m,col=method,fill=method))+
  geom_line()+
  #geom_ribbon(aes(ymin=max.m-max.sd,ymax=max.m+max.sd),alpha=0.2)+
  facet_grid(disp~rangepar)

results_sum%>%
  ggplot(aes((d_cutoff),mode.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=mode.m-mode.sd,ymax=mode.m+mode.sd),alpha=0.2)+
  facet_grid(disp~rangepar)



min_res=function(data,par){
  with(data,sum((y-(1/(1+exp(-(x-par[1])/par[2]))))^2))
}


optim(par=c(20, 0.1), fn=min_res, data=data)



#Plots
clus_null=read.csv("clustdat_null.csv")%>%as_tibble()
clus_tree=read.csv("clustdat_raw.csv")%>%as_tibble()
clus_cell=read.csv("clustdat_raw_cell.csv")%>%as_tibble()


#Null data
clus_null%>%
  mutate(rangepar=as.factor(rangepar))%>%
  filter(unit=="tree")%>%
  ggplot(aes(x=d_cutoff,y=clust))+
  geom_point()+
  facet_wrap(vars(rangepar))



clus_cell_sum=clus_cell%>%
              group_by(rep,rangepar,disp,method,d_cutoff)%>%
              summarize(tot=sum(clust*freq))%>%
              ungroup()

clus_cell_sum=clus_cell%>%
  group_by(rep,rangepar,disp,method,d_cutoff)%>%
  #mutate(clust=clust/max(clust))%>%
  summarize(max=max(clust),
            peak=max(freq),
            mean=sum(clust*freq)/length(clust),
            unique=length(clust))%>%
  ungroup()%>%
  group_by(rangepar,disp,method,d_cutoff)%>%
  summarize(max.m=mean(max),
            max.sd=sd(max),
            peak.m=mean(peak),
            peak.sd=sd(peak),
            mean.m=mean(mean),
            mean.sd=sd(mean),
            unique.m=mean(unique),
            unique.sd=sd(unique))%>%
  ungroup()

clus_cell_sum%>%
  mutate(rangepar=as.factor(rangepar))%>%
  ggplot(aes(d_cutoff,max.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=max.m-max.sd,ymax=max.m+max.sd),alpha=0.2)+
  facet_grid(rangepar~disp)


clus_cell_sum%>%
  ggplot(aes(d_cutoff,peak.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=peak.m-peak.sd,ymax=peak.m+peak.sd),alpha=0.2)+
  facet_grid(rangepar~disp)

clus_cell_sum%>%
  ggplot(aes(d_cutoff,mean.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd),alpha=0.2)+
  facet_grid(rangepar~disp)

clus_cell_sum%>%
  ggplot(aes(d_cutoff,unique.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=unique.m-unique.sd,ymax=unique.m+unique.sd),alpha=0.2)+
  facet_grid(rangepar~disp)



clus_stat=clus_tree%>%
  group_by(rep,rangepar,disp,method,d_cutoff)%>%
  summarize(max=max(clust),
            peak=max(freq),
            mean=sum(clust*freq)/length(clust),
            unique=length(clust))%>%
  ungroup() 


clus_tree_thr=clus_stat%>%
              group_by(rep,rangepar,disp,method)%>%
              slice_max(unique)%>%
              ungroup()%>%
              group_by(rangepar,disp,method)%>%
              summarize(thr=mean(d_cutoff),
                        thr.sd=sd(d_cutoff))%>%
              ungroup()

clus_tree_thr%>%
  ggplot(aes(disp,thr,col=method))+
  geom_point()+
  geom_errorbar(aes(ymin=thr-thr.sd,ymax=thr+thr.sd))+
  facet_wrap(vars(rangepar))

clus_tree_thr%>%
  ggplot(aes(x=disp,y=method,col=thr,fill=thr))+
  geom_tile()+
  facet_wrap(vars(rangepar))




clus_tree_sum=clus_tree%>%
group_by(rep,rangepar,disp,method,d_cutoff)%>%
  summarize(max=max(clust),
            peak=max(freq),
            mean=sum(clust*freq)/length(clust),
            unique=length(clust))%>%
  ungroup()%>%
  group_by(rangepar,disp,method,d_cutoff)%>%
  summarize(max.m=mean(max),
            max.sd=sd(max),
            peak.m=mean(peak),
            peak.sd=sd(peak),
            mean.m=mean(mean),
            mean.sd=sd(mean),
            unique.m=mean(unique),
            unique.sd=sd(unique))%>%
  ungroup()
  
  
clus_tree_sum%>%
  mutate(rangepar=as.factor(rangepar))%>%
  ggplot(aes(d_cutoff,max.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=max.m-max.sd,ymax=max.m+max.sd),alpha=0.2)+
  facet_grid(rangepar~disp)

 
clus_tree_sum%>%
  ggplot(aes(d_cutoff,peak.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=peak.m-peak.sd,ymax=peak.m+peak.sd),alpha=0.2)+
  facet_grid(rangepar~disp)

clus_tree_sum%>%
  ggplot(aes(d_cutoff,mean.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd),alpha=0.2)+
  facet_grid(rangepar~disp)

clus_tree_sum%>%
  ggplot(aes(d_cutoff,unique.m,col=method,fill=method))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=unique.m-unique.sd,ymax=unique.m+unique.sd),alpha=0.2)+
  facet_grid(rangepar~disp)



####################################################################3
#Fit weibull distribution to cluster size data

library(tidyverse)
library(parallel)
library(furrr)



cores = detectCores() - 2
plan(multisession, workers = cores)

#Load data
clusdat=read.csv("C:/Users/mihir/Documents/percolate/codes/clusdat.csv")%>%as_tibble()%>%
        inner_join(clusdat%>%
                     group_by(sp,rep,comp_dist,disp,k,d_cutoff)%>%
                     summarize(pop=sum(cluster*freq))%>%
                     ungroup())%>%
        mutate(cluster=cluster/pop)%>%
        mutate(cluster=seq(0,1,0.05)[cut(cluster,seq(0,1,0.05),label=FALSE,include.lowest=TRUE)])




pars=clusdat%>%
      select(sp,rep,comp_dist,disp,k,k,d_cutoff)%>%
      unique()%>%
      as.data.frame()

names(pars)=c("sps","reps","comp_dists","disps","ks","d_cutoffs")


dat1=clusdat%>%
    filter(sp=="alchco",
           rep==1,
           comp_dist==15,
           disp==2,
           k==-20,
           d_cutoff==1)


res=NULL


for(i in 1:nrow(pars)){
  
      dat=clusdat%>%
          filter(sp==pars$sps[i],
                 rep==pars$reps[i],
                 comp_dist==pars$comp_dists[i],
                 disp==pars$disps[i],
                 k==pars$ks[i],
                 d_cutoff==pars$d_cutoffs[i])
      
      dat$cluster=dat$cluster/maxpop
      
      tryCatch({
        shape1=coef(MASS::fitdistr(sample(dat$cluster,500,replace=TRUE,prob=dat$freq),"Weibull"))[1]
      }, error=function(e){})
        
        
      res=bind_rows(res,
                    tibble(
                      sp=pars$sps[i],
                      rep=pars$reps[i],
                      comp_dist=pars$comp_dists[i],
                      disp=pars$disps[i],
                      k=pars$ks[i],
                      d_cutoff=pars$d_cutoffs[i],
                      shape=shape1))
}


p1=res%>%
  group_by(pop,comp_dist,d_cutoff)%>%
  summarize(shape.m=mean(shape),
            shape.sd=sd(shape))%>%
  ungroup()%>%
  mutate(ld=log(d_cutoff))%>%
  mutate(comp_dist=as.factor(comp_dist))%>%
  ggplot(aes(log(d_cutoff),(shape.m),col=comp_dist,fill=comp_dist))+
  geom_line()+
  #geom_ribbon(aes(ymax=shape.m+shape.sd,ymin=shape.m-shape.sd,alpha=0.2))+
  facet_wrap(vars(pop),scales="free")

p2=res%>%
  group_by(pop,comp_dist,d_cutoff)%>%
  summarize(shape.m=mean(shape),
            shape.sd=sd(shape))%>%
  ungroup()%>%
  mutate(ld=log(d_cutoff))%>%
  mutate(comp_dist=as.factor(comp_dist))%>%
  ggplot(aes((d_cutoff),(shape.m),col=comp_dist,fill=comp_dist))+
  geom_line()+
  #geom_ribbon(aes(ymax=shape.m+shape.sd,ymin=shape.m-shape.sd,alpha=0.2))+
  facet_wrap(vars(pop),scales="free")


library(patchwork)
           
p1+p2

res%>%
  group_by(pop,k,d_cutoff)%>%
  summarize(shape.m=mean(shape),
            shape.sd=sd(shape))%>%
  ungroup()%>%
  mutate(k=as.factor(k))%>%
  ggplot(aes((d_cutoff),log(shape.m),col=k))+
  geom_line()+
  #geom_ribbon(aes(ymax=shape.m+shape.sd,ymin=shape.m-shape.sd,alpha=0.2))+
  facet_wrap(vars(pop))

res%>%
  group_by(pop,disp,d_cutoff)%>%
  summarize(shape.m=mean(shape),
            shape.sd=sd(shape))%>%
  ungroup()%>%
  mutate(disp=as.factor(disp))%>%
  ggplot(aes((d_cutoff),log(shape.m),col=disp))+
  geom_line()+
  #geom_ribbon(aes(ymax=shape.m+shape.sd,ymin=shape.m-shape.sd,alpha=0.2))+
  facet_wrap(vars(pop))
 

res%>%filter(
  sp==pars$sp[1],
  #rep==pars$rep[1],
  comp_dist==pars$comp_dist[1],
  disp==pars$disp[1],
  k==pars$k[1])%>%
  mutate(rep=as.factor(rep))%>%
  group_by(d_cutoff)%>%
  summarize(shape=mean(shape))%>%
  ungroup()%>%
  ggplot(aes(d_cutoff,shape))+
  geom_line()
  
spdata=tibble(sp=c( "alchco", "alsebl", "casear", "cecrin", "cordal", "crotbi", "jac1co",
          "luehse", "micoar", "paligu", "termam", "zantpr"),pop=c(271,9913,108,1381,223,634,309,
                                                             199,652,1326,58,72))
res%<>%inner_join(spdata)

