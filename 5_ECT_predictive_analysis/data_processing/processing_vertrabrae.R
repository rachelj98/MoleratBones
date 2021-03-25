library(fBasics)
library(kernlab)
directory="vertabrae/"
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

Features=data.frame(matrix(NA, ncol=16200,nrow=t))
for (i in 1:t) {
  A=R[i] #strsplit(R[i],".off")
  name<-toString(A[[1]][1])
  Names<-c(Names,name)
  V=c()
  lines<-readLines(paste(directory,R[i], sep=""))
  for (line in lines) {
    #To get of the rid of the extra comma at the end of each row
    vline=substr(line,1,nchar(line)-1)
    v=as.numeric(unlist((strsplit(vline, ","))))
    integrate(v)
    V=c(V,v)
  }
  Features[i,]<-V
}


Names=Names[1:37]
Features<-Features[,seq(1, length(Features), 50)]
Data<-cbind(Names,Features)
names(Data)[1]<-"filename"
Labels<-read.csv("molerat-traits.csv")
Merged<-merge(Labels,Data, by="filename")
# Drop observation Z3F022 From the analysis
Merged<-Merged[-37,]
New<-Merged[,-c(1,3)]
saveRDS("Reproduced_LumbarData",New)
