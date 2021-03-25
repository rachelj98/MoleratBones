#Select one of the following cases by commenting/uncommenting the pertinent lines.

#Case 1: Analyzing all the femurs
#trueCR=readRDS("femurfullCR")
#directory='femurs/'

#Case 2: Analyzing femurs of queens with at least 6 offsprings
#trueCR=readRDS("atleastsixCR")
#directory='femursatleastsix/'

#Case 3: Analyzing vertabrae:
trueCR=readRDS("lumbarCR")
directory='vertabrae/'

trueMCR=1-sum(trueCR)/length(trueCR)
directory='femurs/'
stat=rep(0,100)
#stat[1]=trueMCR

R<-list.files(path=directory)
for (i in 1:length(R)) {
  stat[i]=readRDS(paste(directory,"cr",i,sep=''))
}

sum(trueMCR>=stat)/length(stat)
