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

Lx=128
Ly=128

divisors=function(x){
  y=seq_len(x)
  y[ x%%y ==0 ]
}

d_cutoffs=divisors(min(Lx,Ly))

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



setwd("C:/Users/mihir/Documents")

files=list.files("C:/Users/mihir/Documents/landscapes")

init.files=files[grep("initial",files)]
landscape.files=files[grep("landscape",files)]
final.files=files[grep("final",files)]


sets=1:length(final.files)

cores=detectCores()
plan(multisession,workers=cores)

results=sets%>%
  future_map_dfr(
    .f=function(sets){
      
      finset=read.csv(paste0("C:/Users/mihir/Documents/landscapes/","final_",sets,".csv"))
      rep=finset[1,4]
      rangepar=finset[1,5]
      
      params=finset%>%count(disp,method)%>%select(-n)
      
      fin.res=NULL
      
      for(i in 1:nrow(params)){
        
        final=finset%>%
          filter(disp==params[i,1],
                 method==params[i,2])%>%
          select(x,y,pres)
        
        final.res=perc.up(final)%>%
          mutate(time=5000,
                 rep=rep,
                 rangepar=rangepar,
                 disp=params[i,1],
                 method=params[i,2])
        
        fin.res=bind_rows(
          fin.res,final.res
        )
        
        
      }
      
      return(fin.res)
    })
      
      
  result=readRDS("perc.res.rds")



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


result=readRDS('perc.res.rds')
#Summarize the results
pars=results%>%
      select(rep,rangepar,disp,method)%>%
      unique()%>%
      as.data.frame()

clusdat=NULL
for(i in 1:nrow(pars)){
  
  dat=results%>%
      filter(rep==pars[i,1],
             rangepar==pars[i,2],
             disp==pars[i,3],
             method==pars[i,4])%>%
      select(d_cutoff,clust,freq)
  
  clusmax=dat%>%
          mutate(clust=clust*(d_cutoff^2)/(Lx*Ly))%>%
          group_by(d_cutoff)%>%
          slice_max(clust)%>%
          ungroup()
  
  x=clusmax$d_cutoff
  y=clusmax$clust
  
  fit=nls(y~1/(1+exp((xmid-x)/scal)),data.frame(x,y),start=list(xmid=30,scal=5))
  

  clusdat=bind_rows(
                clusdat,
                tibble(
                  rep=pars[i,1],
                  rangepar=pars[i,2],
                  disp=pars[i,3],
                  method=pars[i,4],
                  start=coef(fit)[2],
                  infl=coef(fit)[1],
                ))

 }


rangepar.labs <- c("Spatial Autocorrelation: low", 
                   "Spatial Autocorrelation: medium",
                   "Spatial Autocorrelation: high")
names(rangepar.labs) <- c("5", "15","30")


clusdat%>%
  select(rep,rangepar,disp,method,infl)%>%
  mutate(rangepar=as.factor(rangepar),
         disp=as.factor(disp),
         method=as.factor(method))%>%
  group_by(rangepar,disp,method)%>%
  summarize(infl.mean=mean(infl),
            infl.se=(sd(infl)/sqrt(length(infl))))%>%
  ungroup()%>%
  mutate(infl.upper=infl.mean+infl.se,
         infl.lower=infl.mean-infl.se)%>%
  ggplot(aes(x=disp,y=infl.mean,col=method))+
  geom_line()+
  geom_pointrange(aes(ymin=infl.lower, ymax=infl.upper))+
  facet_grid(~rangepar,labeller = labeller(rangepar=rangepar.labs))+
  xlab("Dispersal range")+
  ylab("Inflection point")



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

fileno=c(1:50,151:200)

nullres=fileno%>%
  future_map_dfr(
    .f=function(fileno){
      
      land=read.csv(paste0("C:/Users/mihir/Documents/landscapes/landscape_",fileno,".csv"))%>%
        as_tibble()
      
      mid=quantile(land$soiltype,0.5)
      rangepar=land$rangepars[1]
      
      land=land%>%
            filter(soiltype>mid)%>%
            mutate(pres=1)%>%
            select(x,y,pres)
      
      res1=NULL
      
      for(i in 1:5){
        land1=land%>%
              slice(sample((1:nrow(land)),1000))
        
        res.cell=perc.up(land1,unit="cell")
        res.tree=perc.up(land1,unit="tree")
        
        res1=bind_rows(res1,
                       tibble(rangepar=rangepar,
                              unit="cell",
                              d_cutoff=res.cell$d_cutoff,
                              clust=res.cell$clust,
                              freq=res.cell$freq),
                       tibble(rangepar=rangepar,
                              unit="tree",
                              d_cutoff=res.tree$d_cutoff,
                              clust=res.tree$clust,
                              freq=res.tree$freq))
        
      }
      
      return(res1)
      
    }
  )

write.csv(nullres,"clustdat_null.csv")

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
                  max=max(clust),
                  mode=clust[which.max(freq)],
                  peak=max(freq))%>%
        ungroup()


min_res=function(data,par){
  with(data,sum((y-(1/(1+exp(-(x-par[1])/par[2]))))^2))
}




meandat=NULL

parms=results_sum%>%select(rep,rangepar,disp,method)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i],
           disp==parms$disp[i],
           method==parms$method[i])%>%
    select(d_cutoff,means)%>%as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  meandat=bind_rows(meandat,
                    tibble(res=parms$rep[i],
                           rangepar=parms$rangepar[i],
                           disp=parms$disp[i],
                           method=parms$method[i],
                           mid=foo.res$par[1],
                           scale=foo.res$par[2]))
}
  
  

maxdat=NULL

parms=results_sum%>%select(rep,rangepar,disp,method)%>%unique()

for(i in 1:nrow(parms)){
  
  dat=results_sum%>%
    filter(rep==parms$rep[i],
           rangepar==parms$rangepar[i],
           disp==parms$disp[i],
           method==parms$method[i])%>%
    select(d_cutoff,max)%>%as.data.frame()
  
  foo.res=optim(par=c(20,0.1),fn=min_res,data=dat)
  
  maxdat=bind_rows(maxdat,
                    tibble(res=parms$rep[i],
                           rangepar=parms$rangepar[i],
                           disp=parms$disp[i],
                           method=parms$method[i],
                           mid=foo.res$par[1],
                           scale=foo.res$par[2]))
}

  



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


clus_tree_max=


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



############################################################################################
#Spectral density profiles

library(RPANDA)
library(igraph)
library(mclust)

theme_set(theme_bw())
theme_update(
  panel.grid = element_blank(),
  strip.background = element_rect(fill = 'orange')
)

sigma = 0.1
gKernel <- function(x) 1/(sigma * sqrt(2 * pi)) * exp(-(x^2)/2 * 
                                                        sigma^2)
kernelG <- function(x, mean = 0, sd = 1) dnorm(x, mean = mean, 
                                               sd = sd)

dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096, 
                 from = min(x) - 3 * sd, to = max(x) + 3 * sd, adjust = 1, 
                 ...) {
  if (has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if (length(x) == 0) 
      stop("no finite or non-missing data!")
  }
  sd <- (if (is.numeric(bw)) 
    bw[1]
    else bw(x)) * adjust
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernel, sd = sd, ...)
  structure(list(x = X, y = rowMeans(M), bw = sd, call = match.call(), 
                 n = length(x), data.name = deparse(substitute(x)), 
                 has.na = has.na), class = "density")
}

integr <- function(x, f) {
  if (!is.numeric(x)) {
    stop("The variable of integration \"x\" is not numeric.")
  }
  if (!is.numeric(f)) {
    stop("The integrand \"f\" is not numeric.")
  }
  if (length(x) != length(f)) {
    stop("The lengths of the variable of integration and the integrand do not match.")
  }
  n = length(x)
  integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + 
                                                    f[1:(n - 1)]))
  return(integral)
}




files=list.files("C:/Users/mihir/Documents/landscapes")
finals=files[grep("final",files)]

Lx=1000
Ly=500

results=NULL

cores=detectCores()
plan(multisession,workers=cores)

#reps=1
#rangepars=5
#disps=10
#methods='positive'

for(j in sample((1:length(finals)),50)){
  
  fin1=read.csv(paste0("C:/Users/mihir/Documents/landscapes/",finals[j]))%>%as_tibble()
  
  pars=fin1%>%select(rep,rangepar,disp,method)%>%unique()
  
  res2=pars%>%
    future_pmap_dfr(
      .f=function(reps,rangepars,disps,methods){
        
        dat=fin1%>%
          filter(rep==reps,
                 rangepar==rangepars,
                 disp==disps,
                 method==methods
          )
        
        adj=as.matrix(dist(dat[,c(1,2)],method='euclidean'))
        e = eigen(graph.laplacian(graph.adjacency(adj,weighted = T), normalized = T), symmetric = T, only.values = T)
        x = subset(e$values, e$values >= 0)
        d = dens(log(x))
        
        dsc = d$y/(integr(d$x, d$y))
        principal_eigenvalue <- max(x)
        skew <- moments::skewness(x)
        peak_height <- max(dsc)
        gaps <- abs(diff(x))
        gapMat <- as.matrix(gaps)
        modalities <- c(1:length(gapMat))
        gapMatCol <- cbind(modalities, gapMat)
        eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[,2]))
        
        return(tibble(rep=reps,
               rangepar=rangepars,
               disp=disps,
               method=methods,
               max_eigen=principal_eigenvalue,
               skewness=skew,
               peak=peak_height,
               eigengap=eigenGap[1]))
        
      },
      .options=furrr_options(seed=TRUE)
    )
  
  results=bind_rows(results,res2)
}

result1=results%>%
  group_by(rangepar,disp,method)%>%
  summarize(eigenmax.m=mean(max_eigen),
            eigenmax.sd=sd(max_eigen),
            skewness.m=mean(skewness),
            skewness.sd=sd(skewness),
            peak.m=mean(peak),
            peak.sd=sd(peak),
            gap.m=mean(eigengap),
            gap.sd=sd(eigengap))%>%
  ungroup()
            

result1%>%
  ggplot(aes(disp,eigenmax.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=eigenmax.m-eigenmax.sd,ymax=eigenmax.m+eigenmax.sd),alpha=0.2)+
  #geom_pointrange(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd))+
  facet_grid(~rangepar)

result1%>%
  ggplot(aes(disp,skewness.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=skewness.m-skewness.sd,ymax=skewness.m+skewness.sd),alpha=0.2)+
  #geom_pointrange(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd))+
  facet_grid(~rangepar)

result1%>%
  ggplot(aes(disp,peak.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=peak.m-peak.sd,ymax=peak.m+peak.sd),alpha=0.2)+
  #geom_pointrange(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd))+
  facet_grid(~rangepar)

result1%>%
  ggplot(aes(disp,gap.m,col=method,fill=method))+
  geom_line()+
  geom_ribbon(aes(ymin=gap.m-gap.sd,ymax=gap.m+gap.sd),alpha=0.2)+
  #geom_pointrange(aes(ymin=mean.m-mean.sd,ymax=mean.m+mean.sd))+
  facet_grid(~rangepar)


    j=1

fin1=read.csv(paste0("C:/Users/mihir/Documents/landscapes/",finals[j]))%>%as_tibble()

fin1%>%
  ggplot(aes(x,y))+
  geom_point()+
  facet_grid(disp~method)

pars=fin1%>%select(rep,rangepar,disp,method)%>%unique()

reps=1
rangepars=5
disps=10
method='positive'

dat=fin1%>%
    filter(rep==reps,
           rangepar==rangepars,
           disp==disps,
           method==methods)

adj=as.matrix(dist(dat[,c(1,2)],method='euclidean'))
e = eigen(graph.laplacian(graph.adjacency(adj,weighted = T), normalized = T), symmetric = T, only.values = T)
x = subset(e$values, e$values >= 0)
d = dens(log(x))

dsc = d$y/(integr(d$x, d$y))
principal_eigenvalue <- max(x)
skewness <- skewness(x)
peak_height <- max(dsc)
gaps <- abs(diff(x))
gapMat <- as.matrix(gaps)
modalities <- c(1:length(gapMat))
gapMatCol <- cbind(modalities, gapMat)
eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
                                                              2]))
res <- list(eigenvalues = x, principal_eigenvalue = principal_eigenvalue, 
            asymmetry = skewness, peakedness = peak_height, eigengap = eigenGap[, 
                                                                                1])

plot(d)




res=read.csv("result1.csv")

res%>%
  ggplot(aes(x,y))+
  geom_point()+
  facet_grid(disp~method)


clustdat=NULL

for(i in 1:nrow(pars)){
  
  dat=res%>%filter(disp==pars$disp[i], method==pars$method[i])
  
  dat1=perc.up(dat,unit="tree")
  
  clustdat=bind_rows(
    clustdat,
    tibble(disp=pars$disp[i],
           method=pars$method[i],
           d_cutoff=dat1$d_cutoff,
           clust=dat1$clust,
           freq=dat1$freq)
  )
  
  
  
}
