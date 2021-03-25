#####################################################
#This R code tests whether cell type composition significantly mediates the effects of breeding status on gene expression
#The code is modified from Snyder-Mackler et al 2016 Science: https://github.com/nsmackler/status_genome_2016/blob/master/cell%20specific%20expression%20analysis.R#
#####################################################
require(limma);require(edgeR); require(EMMREML);require(RColorBrewer);require(ggplot2)
require(gridExtra); require(grid); require(doParallel)

#Create residuals dataframe of 329 breeder status-associated genes 
res<-read.table('residuals_propinfeat_litter.txt',header=T,sep='\t') #See script Analysis_code_gene_expression_linear_modeling.R for how this file had been generated
modelgenes=read.table('model_bone_bonequ_bonetreat_dayscult_EMMREML_regrlitterpropinfeatvoom.txt',sep='\t',header=T)
head(modelgenes)
modelgenes<-subset(modelgenes,fdr_Queen_at_lb<=0.2|fdr_Queen_at_vert<=0.2)
res2=res[modelgenes$Row.names,]
write.table(res2,file='residuals_propinfeat_litter_fdr0.2.txt',sep='\t',quote=F)

#Read in data files. We'll only analyze the 329 genes that were significantly differentially expressed with breeder status
residuals <- lapply("~/Desktop/hiseq_allbmscs/residuals_propinfeat_litter_fdr0.2.txt", read.table) 
names(residuals)=c("regrproplitter")
info_cell <- lapply("~/Desktop/github_dmr_bone/gene_expression_linear_modeling/traits_47bmscs.txt", read.table,header=T) 
names(info_cell)=c("regrproplitter")
rownames(info_cell[[1]])<-info_cell[[1]]$sample
job13<-read.table("CIBERSORT.Output_Job13.txt",sep="\t",header=T)
j13<-job13;rownames(j13)<-j13[,1];j13=j13[,-c(1,29:34)]

#Obtain the first three principle components of the cell type composition, and add to the traits
pca2<-prcomp(j13,scale=T,center=T)
load2<-pca2$x[as.character(traits[,"sample"]),]
load2<-as.data.frame(load2)
info_cell[[1]]$ciberpc1<-load2$PC1
info_cell[[1]]$ciberpc2<-load2$PC2
info_cell[[1]]$ciberpc3<-load2$PC3

#Format traits
info_cell[[1]]$queen<-as.factor(info_cell[[1]]$queen)
info_cell[[1]]$bone<-as.factor(info_cell[[1]]$bone)
info_cell[[1]]$treatment<-as.factor(info_cell[[1]]$treatment)
info_cell[[1]]$treatment <- relevel(info_cell[[1]]$treatment, ref = "null")
info_cell[[1]]$ciberpc1<-sample(info_cell[[1]]$ciberpc1,replace=F)
traits<-info_cell[[1]]
traits$animal<-as.factor(traits$animal)

#Read in relatedness matrix, K
K<-read.table("~/Desktop/github_dmr_bone/gene_expression_linear_modeling/relatedness.txt",sep="\t",header=T)
K<-as.matrix(K[levels(traits$animal),levels(traits$animal)])
kin<-K

#Create matrix Z, which is a matrix mapping samples to animals
Z=matrix(rep(0,nrow(traits)*ncol(K)),nrow=nrow(traits),ncol=ncol(K)) 
rownames(Z)=rownames(traits)
colnames(Z)=colnames(K)
for(i in 1:ncol(Z))
{
  set=which(traits$animal == colnames(Z)[i])
  Z[set,i]=1
}
#################################################################################################
#mediation effects of ciberpc1
#################################################################################################
## First, run the model with ciberpc1 (or permutated ciberpc1 if running permutations)
## ciberpc1 model
design = model.matrix(~bone+bone:queen+bone:treatment+daysinculture+ciberpc1,data=info_cell[[1]])
EMMA_ciberpc1=lapply(c("regrproplitter"),function(x){
  tmp=t(apply(residuals[[x]],1,function(y){
    design = model.matrix(~bone+bone:queen+bone:treatment+daysinculture+ciberpc1,data=info_cell[[x]])
    emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    p=emma$pvalbeta[,"none"]
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p))
  }));
  colnames(tmp)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(tmp)[(ncol(design)+1):(2*ncol(design))]=paste0("se_",colnames(design))
  colnames(tmp)[((2*ncol(design))+1):(3*ncol(design))]=paste0("pval_",colnames(design));
  return(tmp)}); names(EMMA_ciberpc1)=c("regrproplitter")
#################################################################################################
## Next, bootstrap the data to calculate confidence intervals of the indirect effect
clus <- makeCluster(10)
registerDoParallel(cores=10)  

niter=1000  ## NOTE: This will take a long time. We found that it was best to run these in groups of 10 and to parallelize across mutliple machines with multiple nodes each. Make sure that you set a different seed, using set.seed(), for each new job. Otherwise you will end up with the same data from each run...
set.seed(12) #seed 12
system.time(for (t in 1:1000){
  print(t)
  if (t==1) {bootstrap_results1=lapply(names(info_cell), function(cell){
    C=info_cell[[cell]] 
    n=sample(rownames(C),replace=T) ## sample with replacement
    C=C[n,]  ## covariate matrix
    C$queen<-as.factor(C$queen); C$bone<-as.factor(C$bone); C$treatment<-as.factor(C$treatment); C$treatment <- relevel(C$treatment, ref = "null")
    E=residuals[[cell]][,n] ## expression matrix
    Zmat=matrix(nrow=nrow(C),ncol=length(unique(C$animal))) 
    colnames(Zmat)=unique(C$animal)
    rownames(Zmat)=C$animal
    for (r in 1:nrow(Zmat)) {
      for (c in 1:ncol(Zmat)) {
        if (rownames(Zmat)[r]==colnames(Zmat)[c]) {Zmat[r,c]=1}else {Zmat[r,c]=0}} }
    K=kin[unique(C$animal),unique(C$animal)]
    ## run models
    clusterExport(clus,varlist=c("C","E","Zmat","K"),envir=environment())
    ## model 1: model
    m1=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(terms(~C$bone+C$daysinculture+C$bone:C$queen+C$bone:C$treatment,keep.order=T)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T); return(c(emma$betahat[4:5],emma$varbetahat[4:5]))}))
    ## model 2: model + ciberpc1
    m2=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(terms(~C$bone+C$daysinculture+C$bone:C$queen+C$bone:C$treatment+C$ciberpc1,keep.order=T)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[c(4,5,8)],emma$varbetahat[c(4,5,8)]))}))
    ## model a1: #ciberpc1 as a function of bone and queen in each bone
    a1=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[3]
    var_a1=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[3]
    a2=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[4]
    var_a2=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[4]
    ## combine the data into a dataframe that will be returned from the lapply function
    data=cbind(m1,m2,rep(a1,nrow(m1)),rep(a2,nrow(m1)),rep(var_a1,nrow(m1)),rep(var_a2,nrow(m1)))
    colnames(data)=c("c_qlb","c_qvert","var_cqlb","var_cqvert","qlb_prime_ciberpc1_SM","qvert_prime_ciberpc1_SM","b1_ciberpc1_SM","var_cqlb_prime_ciberpc1_SM","var_cqvert_prime_ciberpc1_SM","var_b1_ciberpc1_SM","a1_qlb","a2_qvert","var_a1_qlb","var_a2qvert")
    return(data)
  })
  names(bootstrap_results1)=names(info_cell)
  bootstrap_results=bootstrap_results1}
  if (t!=1) {
    bootstrap_results1=lapply(names(info_cell), function(cell){
      C=info_cell[[cell]]
      n=sample(rownames(C),replace=T) ## sample with replacement
      C=C[n,]  ## covariate matrix
      E=residuals[[cell]][,n] ## expression matrix
      Zmat=matrix(nrow=nrow(C),ncol=length(unique(C$animal)))
      colnames(Zmat)=unique(C$animal)
      rownames(Zmat)=C$animal
      for (r in 1:nrow(Zmat)) {
        for (c in 1:ncol(Zmat)) {
          if (rownames(Zmat)[r]==colnames(Zmat)[c]) {Zmat[r,c]=1}else {Zmat[r,c]=0}} }
      K=kin[unique(C$animal),unique(C$animal)]
      ## run models
      clusterExport(clus,varlist=c("C","E","Zmat","K"),envir=environment())
      ## model 1: model
      m1=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(terms(~C$bone+C$daysinculture+C$bone:C$queen+C$bone:C$treatment,keep.order=T)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T); return(c(emma$betahat[4:5],emma$varbetahat[4:5]))}))
      ## model 2: model + ciberpc1
      m2=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(terms(~C$bone+C$daysinculture+C$bone:C$queen+C$bone:C$treatment+C$ciberpc1,keep.order=T)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[c(4,5,8)],emma$varbetahat[c(4,5,8)]))}))
      ## model a1: #ciberpc1 as a function of bone and queen
      a1=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[3]
      var_a1=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[3]
      a2=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[4]
      var_a2=emmreml(y=C$ciberpc1,X=model.matrix(~C$bone+C$bone:C$queen),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[4]
      ## combine the data into a dataframe that will be returned from the lapply function
      data=cbind(m1,m2,rep(a1,nrow(m1)),rep(a2,nrow(m1)),rep(var_a1,nrow(m1)),rep(var_a2,nrow(m1)))
      colnames(data)=c("c_qlb","c_qvert","var_cqlb","var_cqvert","qlb_prime_ciberpc1_SM","qvert_prime_ciberpc1_SM","b1_ciberpc1_SM","var_cqlb_prime_ciberpc1_SM","var_cqvert_prime_ciberpc1_SM","var_b1_ciberpc1_SM","a1_qlb","a2_qvert","var_a1_qlb","var_a2qvert")
      return(data)
    })
    names(bootstrap_results1)=names(info_cell)
    bootstrap_results2=lapply(names(bootstrap_results),function(cell){cbind(bootstrap_results[[cell]],bootstrap_results1[[cell]])}); names(bootstrap_results2)=names(bootstrap_results);bootstrap_results=bootstrap_results2}
})


### Write the data to a text file for posterity (and to combine if you parallelized across multiple machines)
lapply(names(bootstrap_results), function(x){write.table(as.matrix(bootstrap_results[[x]]),col.names=T,quote=F,row.names=T,file=paste0(x,"ciberpc1_mediation.txt"))})  

## calculate the median and 95% CIs of the indirect effects (animalE) of ciberpc1
## Queen in lb
#niter=9044 #specify if niter stopped prematurely
ciberpc1_qlb_animalE=lapply(bootstrap_results,function(x){t(apply(x,1,function(z){
  y=as.numeric(z)
  a=y[seq(1,niter*14,by=14)]/sqrt(y[seq(3,niter*14,by=14)]) ## standardized effect in model 1 
  b=y[seq(5,niter*14,by=14)]/sqrt(y[seq(8,niter*14,by=14)]) ## standardized effect of in model 2 
  ciberpc1_qlb_animalE=a-b
  q1=quantile(ciberpc1_qlb_animalE,c(0.025,0.5,0.975))
  return(q1)
}))})
#Queen in vert
ciberpc1_qvert_animalE=lapply(bootstrap_results,function(x){t(apply(x,1,function(z){
  y=as.numeric(z)
  a=y[seq(2,niter*14,by=14)]/sqrt(y[seq(4,niter*14,by=14)]) ## standardized effect in model 1 
  b=y[seq(6,niter*14,by=14)]/sqrt(y[seq(9,niter*14,by=14)]) ## standardized effect in model 2 
  ciberpc1_qvert_animalE=a-b
  q1=quantile(ciberpc1_qvert_animalE,c(0.025,0.5,0.975))## 
  return(q1)
}))})

## Count the number of genes that with a bootstrapped 95% CI that does not include 0 (so are significant at p<0.05)
lapply(ciberpc1_qlb_animalE,function(x){sum(x[,3]*x[,1]>0)})
lapply(ciberpc1_qvert_animalE,function(x){sum(x[,3]*x[,1]>0)})

## NOTE: be sure to filter these to only include those genes that are signficantly affected breeding status (as determined from the permutation-based FDR above)
qlb<-ciberpc1_qlb_animalE[[1]]
real=read.delim("model_bone_bonequ_bonetreat_dayscult_EMMREML_regrlitterpropinfeatvoom.txt",header=T)

par(mfrow=c(1,2))
for (i in c(0.1,0.2)) {
  real2<-subset(real,fdr_Queen_at_lb<=i)
  dim(real2)
  qlb2<-qlb[as.character(real2$Row.names),]
  dim(qlb2)
  sum(qlb2[,3]*qlb2[,1]>0)
  aa<-qlb2[which(qlb2[,3]*qlb2[,1]>0),]
  r<-subset(real2,Row.names==as.character(rownames(aa))
            q<-bootstrap_results[[1]]
            q[,1:14]
            
            
            realpc1<-read.delim("model_bone_bonequ_bonetreat_dayscult_ciberpc1_EMMREML_regrlitterpropfeat_voom.txt",header=T)
            dim(real)
            dim(realpc1)
            head(realpc1)
            realpc2<-realpc1[as.character(real2$Row.names),]
            
            plot(density(real2$beta_bone0.queen1),col='steelblue', xlab='Queen in lb betas',lwd=3,main='',cex.lab=1.3,cex.axis=1.3,ylim=c(0,0.42))
            lines(density(realpc2$beta_bone0.queen1),col='darkred',lwd=3)
            legend('topleft',c('With ciberpc1','Without ciberpc1'),col=c('darkred','steelblue'),lwd=c(3,3),cex=1.3,bty='n')
            abline(v=0,lty="dashed",col="darkgray")
}

####Save PDF as 12x6
