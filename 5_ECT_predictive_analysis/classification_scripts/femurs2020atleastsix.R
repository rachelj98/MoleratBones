setwd(getSrcDirectory()[1])
New<-readRDS("femurdata")
#Omit the molerats that have at most 5 offsprings
New<-New[-c(5,8,10,18,24,34,37),]
library(e1071)

#Balance
costs <- table(New[,1]) 
costs[1]<-1
costs[2]<-100000

# Classification rate:

iter=dim(New)[1]
CR=rep(0,dim(New)[1])
predictions=rep(0,iter)
for (k in 1:dim(New)[1]) {
  trainingdata<-New[-k,]
  testdata<-New[k,]
  model<-svm(queenstatus~., data=trainingdata,type="C-classification",kernel="linear",cost=10,scale=FALSE, class.weights=costs)
  preds<-predict(model,testdata[,-1])
  predictions[k]<-preds
  CR[k]<-sum((preds==testdata[,1]))/length(testdata[,1])  
}


# permutation tests for the classification rate
# Note: This computes the misclassification rates for a random permutation, and then stores it. The p-value is then obtained by looking at the quantile at which the original classification rate ranks among these randomly permutated ones, as mentioned in the paper


# helper function to fit the model
onego=function(){
  iter=dim(New)[1]
  CR=rep(0,dim(New)[1])
  predictions=rep(0,iter)
  for (k in 1:dim(New)[1]) {
    trainingdata<-New[-k,]
    testdata<-New[k,]
    model<-svm(queenstatus~., data=trainingdata,type="C-classification",kernel="linear",cost=10,scale=FALSE, class.weights=costs)
    preds<-predict(model,testdata[,-1])
    predictions[k]<-preds
    CR[k]<-sum((preds==testdata[,1]))/length(testdata[,1])  
  }
  return(1-sum(CR)/length(CR))
}


iter=100
saveRDS(CR,"atleastsixCR")
saveRDS(predictions,"atleastsixpreds")
size=length(New[,1])

for (i in 1:iter) {
  permutation=sample(1:size)
  New[,1]=New[permutation,1]
  cr=onego()
  saveRDS(cr, file=paste("femuratleastsix/cr",i,sep=""))
}
