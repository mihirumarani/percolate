using CSV, DataFrames, Random, LinearAlgebra, Distributions, Distances, 
CategoricalArrays, GLM, Plots, Optim, EcologicalNetworks, Graphs, NearestNeighbors

#Declare the functions


#Functions to calculate clusters and further process the simulated data
function divisors(Lx::Int64,Ly::Int64)
    
    d=collect(1:min(Lx,Ly))
    d=d[mod.(min(Lx,Ly),d) .==0]
    
    return sort(union(d, collect(2:2:ceil(Int64,min(Lx,Ly)/4))))
end

function perc_up(dat::DataFrame,Lx::Int64,Ly::Int64)
    
    d=divisors(Lx,Ly)
    
    a=DataFrame(vec(collect(Iterators.product(1:Lx,1:Ly))))
    
    select!(a,"1"=>"x","2"=>"y")
    
    dat=outerjoin(dat,a, on= [:x=>:x,:y=>:y],order=:right)
    
    replace!(dat.pres,missing=>0)
    
    result=DataFrame()
    
    for i in d
        
        cutx=ceil.(Int64,collect((1:Lx) ./ i))
        cuty=ceil.(Int64,collect((1:Ly) ./ i))
        
        dat1=deepcopy(dat)
        
        dat1[!,:x1] =cutx[dat1.x]
        dat1[!,:y1] = cuty[dat1.y]
        
        dat1=combine(groupby(dat1,[:x1,:y1]),:pres=>sum)
        
        mat1=unstack(dat1,:y1,:x1,:pres_sum)
        
        mat1=Matrix{Int64}(select(mat1,Not(:y1)))
        
        clust=zeros(Int64,size(mat1)[1],size(mat1)[2])
        
        largest=0
        
        lbl=collect(0:length(clust))
        
        rw=size(mat1)[1]
        cl=size(mat1)[2]
        
        for i2 in 1:cl 
            for i1 in 1:rw
                
                if mat1[i1,i2] > 0
                    
                    left= (i2-1)==0 ? 0 : clust[i1,i2-1]
                    above= (i1-1)==0 ? 0 : clust[i1-1,i2]
                    
                    if (left==0) && (above==0)
                        
                        largest+=1
                        clust[i1,i2]=largest
                        lbl[(rw*(i2-1))+i1]=largest
                        
                    else
                        
                        if (left!=0) && (above==0)
                            
                            clust[i1,i2]=clust[i1,i2-1]
                            lbl[(rw*(i2-1))+i1]=lbl[(rw*(i2-2))+i1]
                            
                        else
                            
                            if (left==0) && (above!=0)
                                
                                clust[i1,i2]=clust[i1-1,i2]
                                lbl[(rw*(i2-1))+i1]=lbl[(rw*(i2-1))+i1-1]
                                
                            else
                                
                                newlab=min(left,above)
                                lbl[findall(x-> x==left || x==above, lbl)] .=newlab
                                clust[i1,i2]=newlab
                                clust[findall(x-> x==max(left,above), clust)].=newlab
                            end
                        end
                    end
                end
            end
        end
        
        res=DataFrame(pres=reshape(mat1,length(mat1)),cluster=reshape(clust,length(clust)))
        res=combine(groupby(res,:cluster),:pres=>sum)
        res=res[res.cluster .!=0,:]

        freq=combine(groupby(res,:pres_sum),:pres_sum=>length)

        append!(result,DataFrame(d_cutoff=i,cluster=freq.pres_sum,freq=freq.pres_sum_length)) 
        
    end
    
    return result
end


#############################################################################################
#Function to get a matrix of pairwise distances (shortest and between centroids) 
Lx=24
Ly=24
abun=trunc(Int,Lx*Ly*0.1)

dat=DataFrame(x=sample(1:Ly,abun,replace=true),y=sample(1:Lx,abun,replace=true),pres=1)

d=divisors(Lx,Ly)
    
    a=DataFrame(vec(collect(Iterators.product(1:Lx,1:Ly))))
    
    select!(a,"1"=>"x","2"=>"y")
    
    dat=outerjoin(dat,a, on= [:x=>:x,:y=>:y],order=:right)
    
    replace!(dat.pres,missing=>0)
    
    result=DataFrame()
    
    i=2
        
        cutx=ceil.(Int64,collect((1:Lx) ./ i))
        cuty=ceil.(Int64,collect((1:Ly) ./ i))
        
        dat1=deepcopy(dat)
        
        dat1[!,:x1] =cutx[dat1.x]
        dat1[!,:y1] = cuty[dat1.y]
        
        dat1=combine(groupby(dat1,[:x1,:y1]),:pres=>sum)
        
        mat1=unstack(dat1,:y1,:x1,:pres_sum)
        
        mat1=Matrix{Int64}(select(mat1,Not(:y1)))
        
        clust=zeros(Int64,size(mat1)[1],size(mat1)[2])
        
        largest=0
        
        lbl=collect(0:length(clust))
        
        rw=size(mat1)[1]
        cl=size(mat1)[2]
        
        for i2 in 1:cl 
            for i1 in 1:rw
                
                if mat1[i1,i2] > 0
                    
                    left= (i2-1)==0 ? 0 : clust[i1,i2-1]
                    above= (i1-1)==0 ? 0 : clust[i1-1,i2]
                    
                    if (left==0) && (above==0)
                        
                        largest+=1
                        clust[i1,i2]=largest
                        lbl[(rw*(i2-1))+i1]=largest
                        
                    else
                        
                        if (left!=0) && (above==0)
                            
                            clust[i1,i2]=clust[i1,i2-1]
                            lbl[(rw*(i2-1))+i1]=lbl[(rw*(i2-2))+i1]
                            
                        else
                            
                            if (left==0) && (above!=0)
                                
                                clust[i1,i2]=clust[i1-1,i2]
                                lbl[(rw*(i2-1))+i1]=lbl[(rw*(i2-1))+i1-1]
                                
                            else
                                
                                newlab=min(left,above)
                                lbl[findall(x-> x==left || x==above, lbl)] .=newlab
                                clust[i1,i2]=newlab
                                clust[findall(x-> x==max(left,above), clust)].=newlab
                            end
                        end
                    end
                end
            end
        end

heatmap(clust)


clusid=sort(unique(clust))
deleteat!(clusid,1)
xc=zeros(Float64,length(clusid))
yc=deepcopy(xc)
clusind=DataFrame()

for j in 1:(length(clusid))
    inds=findall(clust .==clusid[j])
    xs=[inds[x][1] for x in 1:length(inds)]
    ys=[inds[x][2] for x in 1:length(inds)]
    append!(clusind,DataFrame(id=clusid[j],x=xs,y=ys))
    xc[j]=mean(xs)
    yc[j]=mean(ys)
end

distmat1=pairwise(Euclidean(),[xc yc]')


#Calculate shortest distances between pairs of clusters
#Brute force method
distmat2=zeros(Float64, size(distmat1))
numclus=DataFrame(clus=unique(clusind.id),id=1:length(unique(clusind.id)))

for i in 1:(length(numclus.id)-1), j in (i+1):length(numclus.id)

    cl1=clusind[clusind.id .== numclus.clus[i],:]
    cl1=Matrix([cl1.x cl1.y])

    cl2=clusind[clusind.id .== numclus.clus[j],:]
    cl2=Matrix([cl2.x cl2.y])

    dist=Inf
  
    for i1 in 1:size(cl1)[1], j1 in 1:size(cl2)[1]

        dist=min(dist,sqrt((cl1[i1,1]-cl2[j1,1])^2 + (cl1[i1,2]-cl2[j1,2])^2))
    end

    distmat2[i,j]=dist
end

h1=heatmap(distmat1)
h2=heatmap(distmat2)

plot(h1,h2)

p1=histogram(unique(distmat1))
p2=histogram(unique(distmat2))
plot(p1,p2)

StatsPlots.density(unique(distmat1))
StatsPlots.density!(unique(distmat2))
#############################################################################################
 
#Kill #ndeath individuals with the prob equal to habitat suitability 
#A is a map with a presence/absence data and B is a map with habitat suitability values
#Ndeath=no. of inds to kill
function kill!(A::Matrix{Int8}, B::Matrix{Float64},ndeath::Int64)
    
    samp=findall(A.==1)
    A[wsample(samp,B[samp],ndeath,replace=false)].=0
    return nothing

end


#Calculate recruitment probability for a cell as a function of dispersal kernel and neighborhood 
#crowding effect 
#a: Coordinates of the cell, A:presence/absence matrix,disps: 2*st. dev.of gaussian dispersal
#kernel, d:Max distance over which inds interact, comp_list: prob.vs. distance function for 
#interaction effects.
function recprob(a::CartesianIndex{2},A::Matrix{Int8},disps::Int64,d::Int64,method::Int64,
    comp_list::Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}})

xs=((a[1]-disps):(a[1]+disps))
xs=xs[0 .<xs .<=Ly]

ys=(a[2]-disps):(a[2]+disps)
ys=ys[0 .<ys .<=Lx]

prob=0

if sum(A[xs,ys])>0 

parents=findall(popmat[xs,ys].==1)

for i in 1:length(parents)
prob+=pdf(Normal(0,disps/2),euclidean([a[1],a[2]],[xs[parents[i][1]],ys[parents[i][2]]]))
end

xc=((a[1]-d):(a[1]+d))
xc=xc[0 .<xc .<=Ly]

yc=(a[2]-d):(a[2]+d)
yc=yc[0 .<yc .<=Lx]

prob *= comp_list[method][sum(A[xc,yc])+1]
end

return prob

end

#Function to introduce next generation on a map using the rules of habitat filtering,
#dispersal kernels and local interactions
    
function spawns!(A::Matrix{Int8},disps::Int64,comp_dist::Int64,method::Int64,
    comp_list::Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}})  

samp2=shuffle(findall(A.==0))

probs=map(x->recprob(x,A,disps,comp_dist,method,comp_list),samp2)

A[wsample(samp2,probs,nbirth,replace=false)].=1

return nothing

end 


#Run a population simulation of birth/death for 100 time steps
function sampsim!(A::Matrix{Int8},B::Matrix{Float64},disps::Int64,comp_dist::Int64,method::Int64,
    comp_list::Tuple{Vector{Int64}, Vector{Float64}, Vector{Float64}}) 

for i in 1:100
    
    kill!(A,B)
    spawns!(A,disps,comp_dist,method,comp_list)
    
end

return nothing
end


#Simulation starts here------

Lx=128 
Ly=128
abun=100
nbirth=ndeath=convert(Int64,floor(abun*0.15))

#randomly introduce #abun individuals as an initial population



