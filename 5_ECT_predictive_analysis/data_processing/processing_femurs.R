library(fBasics)
library(kernlab)
directory="femora/"
R<-list.files(path=directory)
t=length(R)
Names<-c()

integrate= function(v){
  mean=mean(v)
  smooth_v=rep(0,length(v))
  v_cent=v-mean(v)
  sum=0
  for (i in 1:length(v)) {
    sum=sum+v_cent[i]
    smooth_v[i]=sum
  }
  return(smooth_v)
}


# In computing the ECT for femurs, we deliberately chose oversampled the heights so that every part of all of the shapes are well-represented.
# However, it is not computationally feasible to fit a model with all of these features, so we want to subsample these features by thinning.
# In this oversampling procedure, there are heights that contain no information: those that do not even intersect any of the shapes, and those that contain all of every shape.
# The purpose of this loop is to discard those features before thinning, so that we have the best possible representation that we can achieve with a given number of covariates.
# We need to loop through all of the shapes to make sure that these features are indeed redundant.

start=2500
end=2500
indices=data.frame(matrix(rep(2500,2*162),ncol=2,nrow=162))
Features=data.frame(matrix(NA, ncol=810000,nrow=t))
for (i in 1:t) {
  A=R[i] #strsplit(R[i],".off") #If we want file names without the .off extension
  name<-toString(A[[1]][1])
  Names<-c(Names,name)
  V=c()
  lines<-readLines(paste(directory,R[i], sep=""))
  runner=1
  for (line in lines) {
    start=2500
    end=2500
    #To get of the rid of the extra comma at the end of each row
    vline=substr(line,1,nchar(line)-1)
    v=as.numeric(unlist((strsplit(vline, ","))))
    tempstart=which(v!=0)[1]
    start=min(start,tempstart)
    indices[runner,1]=min(start,indices[runner,1])
    ENDID=v[-5000]-v[-1]
    tempend=max(which(ENDID!=0))
    end=max(tempend,end)
    indices[runner,2]=max(end,indices[runner,2])
    V=c(V,v)
    runner=runner+1
  }
  Features[i,]<-V
}


# Now that we know what the redundancies are, we can read in the features

sum=0
for (i in 1:dim(indices)[1]) {
  sum=sum+indices[i,2]-indices[i,1]+1
}

Features=data.frame(matrix(NA,ncol=sum,nrow=t))
for (i in 1:t) {
  A=R[i]
  name<-toString(A[[1]][1])
  V=c()
  runner=1
  lines<-readLines(paste(directory,R[i], sep=""))
  for (line in lines) {
    #To get of the rid of the extra comma at the end of each row
    start=indices[runner,1]
    end=indices[runner,2]
    vline=substr(line,1,nchar(line)-1)
    v=as.numeric(unlist((strsplit(vline, ","))))
    v=v[start:end]
    V=c(V,v)
    runner=runner+1
  }
  Features[i,]<-V
}


Names=Names[1:37]
# Thin the features:
Features<-Features[,seq(1, length(Features), 50)]
Data<-cbind(Names,Features)
names(Data)[1]<-"filename"
Labels<-read.csv("molerat-traits-femur-all.csv")
Merged<-merge(Labels,Data, by="filename")
# Drop observation Z3F022 From the analysis
Merged<-Merged[-37,]
New<-Merged[,-c(1,3,4,5)]

saveRDS("Reproduced_femurdata",New)
