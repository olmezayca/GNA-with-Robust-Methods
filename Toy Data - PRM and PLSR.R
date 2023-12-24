

source("My PRM.R")
my_prm
source("NA and PCA.R")

index_95_percent <- which(cumulative_variance >= 95)[1]
a<-index_95_percent

data<- datExpr

##PRM

part1and2<-matrix(nrow=ncol(data)-1,ncol = ncol(data))

for( gene in 1:(ncol(data))) {
  
  prm_model<-my_prm(data[,-(gene)],data[,gene],a=a)
  
  part1ve2[,gene]<-prm_model$yscores %*% t(prm_model$xweights)
  
}
part1and2<-apply(part1and2,2,normalized)

s_ik<-matrix(nrow=ncol(data),ncol = ncol(data))
for (gene in 1:(ncol(data)-1) ) {
  for (i in gene:(ncol(data)-1)) {
    
    s_ik[gene,i+1]<-part1and2[i,gene]+part1and2[gene,(i+1)]
    
  } 
} 

s_ik<- s_ik/2
s_ik[is.na(s_ik)]<-0
s_ik<-t(s_ik)+s_ik

for(i in 1:ncol(s_ik)) {
  s_ik[i,i]<-1
}

write.csv2(s_ik,"s_ik_prm.csv")  

prm_s_ik <- s_ik


##PLSR

data<- datExpr

part1and2<-matrix(nrow=ncol(data)-1,ncol = ncol(data))

for( gene in 1:(ncol(data))) {
  
  model_plsr<-pls1_nipals(data[,-(gene)],data[,gene],a=a)
  
  part1and2[,gene]<-model_plsr$C %*% t(model_plsr$W)
  
}


part1and2<-apply(part1and2,2,normalized)

s_ik<-matrix(nrow=ncol(data),ncol = ncol(data))
for (gene in 1:(ncol(data)-1) ) {
  for (i in gene:(ncol(data)-1)) {
    
    s_ik[gene,i+1]<-part1and2[i,gene]+part1and2[gene,(i+1)]
    
  } 
} 

s_ik<- s_ik/2
s_ik[is.na(s_ik)]<-0
s_ik<-t(s_ik)+s_ik

for(i in 1:ncol(s_ik)) {
  s_ik[i,i]<-1
}

write.csv2(s_ik,"s_ik_plsr.csv")  

plsr_s_ik <- s_ik






