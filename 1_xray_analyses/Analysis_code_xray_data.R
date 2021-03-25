##############################
#Test whether animals in queen/nonbreeder conditions differed in age or weight at experimental *START*.
##############################
#Read in and format info file
info<-read.table('~/Desktop/github_dmr_bone/Johnston_et_al_Supplementary_Table_S1.txt',sep='\t')
info[2,] <- gsub(' ', '_', info[2,])
colnames(info)<-info[2,]
info<-info[-c(1,2),]
info$Imin_at_femoral_midshaft<-as.numeric(info$Imin_at_femoral_midshaft)
info$Weight_in_grams_at_experimental_start<-as.numeric(info$Weight_in_grams_at_experimental_start)
info$Age_in_months_at_experimental_start<-as.numeric(info$Age_in_months_at_experimental_start)
info$LV5_length_at_experimental_start<-as.numeric(info$LV5_length_at_experimental_start)
info<-info[which(info$Experimental_status1=='experimental animal'),]

q<-info[which(info$Breeding_status_at_sample_collection=='queen'),]
h<-info[which(info$Breeding_status_at_sample_collection=='helper'),]
s<-info[which(info$Breeding_status_at_sample_collection=='solitaire'),]

#Body mass
t.test(q$Weight_in_grams_at_experimental_start,h$Weight_in_grams_at_experimental_start,paired=F)
t.test(q$Weight_in_grams_at_experimental_start,s$Weight_in_grams_at_experimental_start,paired=F)
t.test(s$Weight_in_grams_at_experimental_start,h$Weight_in_grams_at_experimental_start,paired=F)

#Age at experimental start
t.test(q$Age_in_months_at_experimental_start,h$Age_in_months_at_experimental_start,paired=F)
t.test(q$Age_in_months_at_experimental_start,s$Age_in_months_at_experimental_start,paired=F)
t.test(s$Age_in_months_at_experimental_start,h$Age_in_months_at_experimental_start,paired=F)
#LV5 at experimental start
t.test(q$LV5_length_at_experimental_start,h$LV5_length_at_experimental_start,paired=F)
t.test(q$LV5_length_at_experimental_start,s$LV5_length_at_experimental_start,paired=F)
t.test(s$LV5_length_at_experimental_start,h$LV5_length_at_experimental_start,paired=F)

#Plot the results as boxplots, Supplementary Figure 1.
library(beeswarm)
par(mfrow=c(1,3))
#Body mass
boxplot(info$Weight_in_grams_at_experimental_start~info$Breeding_status_at_sample_collection,cex.lab=1.4,cex.axis=1.4,outline = FALSE,xlab="",ylab='Mass (g) at experimental start');beeswarm(info[,26]~info[,6], pch = 16, col=4, add = TRUE)
#Age at experimental start
boxplot(info$Age_in_months_at_experimental_start~info$Breeding_status_at_sample_collection,cex.lab=1.4,cex.axis=1.4,outline = FALSE,xlab="",ylab='Age (months) at experimental start');beeswarm(info[,9]~info[,6], pch = 16, col=4, add = TRUE)
#LV5 at experimental start
boxplot(info$LV5_length_at_experimental_start~info$Breeding_status_at_sample_collection,cex.lab=1.4,cex.axis=1.4,outline = FALSE,xlab="",ylab='LV5 length (mm) at experimental start');beeswarm(info[,24]~info[,6], pch = 16, col=4, add = TRUE)
#saved as 4x9 pdf 

#################################
#Test for interaction effects between breeding status and post-pairing time point on LV5 length, 4 months relative to 0 months, and all subsequent time intervals
#################################
library(lme4)
library(lmerTest)

#Read in data
setwd("~/Desktop/github_dmr_bone")
df<-read.table('~/Desktop/github_dmr_bone/Molerat_xray_measures.txt',header=T,sep="\t")
df1<-subset(df,LHAM_mo_experiment=='4'|LHAM_mo_experiment=='0')
df2<-subset(df,LHAM_mo_experiment=='8'|LHAM_mo_experiment=='4')
df3<-subset(df,LHAM_mo_experiment=='12'|LHAM_mo_experiment=='8')
df4<-subset(df,LHAM_mo_experiment=='16'|LHAM_mo_experiment=='12')
df5<-subset(df,LHAM_mo_experiment=='22'|LHAM_mo_experiment=='16')

#Run mixed models, with treatment (Treatment) and months of treatment (LHAM_mo_experiment) as fixed effects, and animal ID as random effect. The 1 indicates that an intercept is to be fitted for each animal ID.
model_4_0 = lmer(LV5 ~ LHAM_mo_experiment + Treatment + LHAM_mo_experiment*Treatment + (1|AnimalID),data=df1,REML=TRUE) 
model_8_4 = lmer(LV5 ~ LHAM_mo_experiment + Treatment + LHAM_mo_experiment*Treatment + (1|AnimalID),data=df2,REML=TRUE) 
model_12_8 = lmer(LV5 ~ LHAM_mo_experiment + Treatment + LHAM_mo_experiment*Treatment + (1|AnimalID),data=df3,REML=TRUE) 
model_16_12 = lmer(LV5 ~ LHAM_mo_experiment + Treatment + LHAM_mo_experiment*Treatment + (1|AnimalID),data=df4,REML=TRUE) 
model_22_16 = lmer(LV5 ~ LHAM_mo_experiment + Treatment + LHAM_mo_experiment*Treatment + (1|AnimalID),data=df5,REML=TRUE) 
summary(model_4_0)
summary(model_8_4)
summary(model_12_8)
summary(model_16_12)
summary(model_22_16)

#####################################################
#Test whether LV5 length change from month 0 to month 4 differs between 1) nonpregnant queens vs nonbreeders, and 2) pregnant queens vs nonbreeders
#####################################################
library("beeswarm")

#Read in data
data<-read.table('~/Desktop/github_dmr_bone/xray_change_4month_vs_0month.txt',sep='\t',header=T)
nonbreeders<-subset(data,Treatment2=='0BREEDER')
breeders<-subset(data,Treatment2=='BREEDER')
notpreg<-subset(breeders,pregnantby4mo==0)
preg<-subset(breeders,pregnantby4mo==1)

# 1) Unpaired t-test, queens not yet pregnant vs nonbreeders:
t.test(notpreg$LV5_change04,nonbreeders$LV5_change04, paired=F) 

# 2) Unpaired t-test, queens that experienced pregnancy vs nonbreeders:
t.test(preg$LV5_change04,nonbreeders$LV5_change04, paired=F)

#####################################################
#Test whether LV5 or total lumbar vertebrae lengths differ between breeders and nonbreeders at 12 months of experimental treatment
#####################################################
m<-read.table('Molerat_xray_measures.txt',sep='\t',header=T)
m<-subset(m,experimental_animal==1)
m$LHAM_mo_experiment<-as.numeric(m$LHAM_mo_experiment)
m<-subset(m,LHAM_mo_experiment==12)
nb<-subset(m,Treatment=="0BREEDER")
b<-subset(m,Treatment=="BREEDER")

#Test whether LV5 length differs with breeder status at 12 months of experimental treatment:
aa<-t.test(b$LV5,nb$LV5,paired=F)
aa
#Calculate % change in mean length in breeders relative to nonbreeders:
(aa$estimate[1]-aa$estimate[2])/aa$estimate[2]*100

#Test whether total length of lumbar vertebral column differs with breeder status at 12 months of experimental treatment:
bb<-t.test(b$L_vert_1.7length,nb$L_vert_1.7length)
bb
#Calculate % change in mean length in breeders relative to nonbreeders:
(bb$estimate[1]-bb$estimate[2])/bb$estimate[2]*100

#Test whether femur length differs with breeder status at 12 months of experimental treatment:
cc<-t.test(b$Femur,nb$Femur,paired=F)
cc

#Test whether tibia length differs with breeder status at 12 months of experimental treatment:
dd<-t.test(b$R_tibia,nb$R_tibia,paired=F)
dd

#Controlling for zygomatic arch, test whether LV5 length differs with breeder status at 12 months of experimental treatment:
aa<-t.test(b$LV5/b$ZygoArchW,nb$LV5/nb$ZygoArchW,paired=F)
aa
#Calculate % change in mean length in breeders relative to nonbreeders:
(aa$estimate[1]-aa$estimate[2])/aa$estimate[2]*100

#Controlling for zygomatic arch, test whether total length of lumbar vertebral column differs with breeder status at 12 months of experimental treatment:
bb<-t.test(b$L_vert_1.7length/b$ZygoArchW,nb$L_vert_1.7length/nb$ZygoArchW)
bb
#Calculate % change in mean length in breeders relative to nonbreeders:
(bb$estimate[1]-bb$estimate[2])/bb$estimate[2]*100

#####################################################
#Create Figure 1C showing vertebral length differences between breeders and nonbreeders at 12 months of experimental treatment
#####################################################
library(ggplot2)
library(cowplot)

#Read in data for 0 months experimental treatment
m<-read.table('Molerat_xray_measures.txt',sep='\t',header=T)
m<-subset(m,experimental_animal==1)
m$LHAM_mo_experiment<-as.numeric(m$LHAM_mo_experiment)
m<-subset(m,LHAM_mo_experiment==0)
nb<-subset(m,Treatment=="0BREEDER")
nb<-nb[,c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')]
b<-subset(m,Treatment=="BREEDER")
b<-b[,c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')]

#Create dataframe "result" that stores values for plotting length differences between breeder and nonbreeder lumbar vertebrae:
result<-as.data.frame(matrix(nrow=7,ncol=7))
colnames(result)<-c('beta','lower','upper','t','df','p','vert')
result$vert<-c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')

for (p in 1:7) {
  aa<-t.test(b[,p],nb[,p])
  result[p,1]<-aa$estimate[1]-aa$estimate[2]
  result[p,2]<-aa$conf.int[1]
  result[p,3]<-aa$conf.int[2]
  result[p,4]<-aa$statistic
  result[p,5]<-aa$parameter
  result[p,6]<-aa$p.value
  
}
a<-result
a$decile<-1:7

#Create 0 months plot with ggplot2
p <- ggplot(a, aes(x = decile, y = beta)) 
plot1<- p + 
  geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed",col="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), data=a, width = .2, color = "black",size=1) +
  geom_point(data=a,size = 3, bg = "black",pch=21,col="black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("B - NB length (mm)") +
  theme(axis.text=element_text(size=30,color='black')) +
  theme(axis.title=element_text(size=30,color='black')) +
  xlab('')+ylim(-0.42,0.55)+
  scale_x_continuous(breaks=c(1:7),labels=a$vert)

#Read in data for 12 months experimental treatment
m<-read.table('Molerat_xray_measures.txt',sep='\t',header=T)
m<-subset(m,experimental_animal==1)
m$LHAM_mo_experiment<-as.numeric(m$LHAM_mo_experiment)
m<-subset(m,LHAM_mo_experiment==12)
nb<-subset(m,Treatment=="0BREEDER")
nb<-nb[,c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')]
b<-subset(m,Treatment=="BREEDER")
b<-b[,c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')]

#Create dataframe "result" that stores values for plotting length differences between breeder and nonbreeder lumbar vertebrae:
result<-as.data.frame(matrix(nrow=7,ncol=7))
colnames(result)<-c('beta','lower','upper','t','df','p','vert')
result$vert<-c('LV7','LV6','LV5','LV4','LV3','LV2','LV1')

for (p in 1:7) {
  aa<-t.test(b[,p],nb[,p])
  result[p,1]<-aa$estimate[1]-aa$estimate[2]
  result[p,2]<-aa$conf.int[1]
  result[p,3]<-aa$conf.int[2]
  result[p,4]<-aa$statistic
  result[p,5]<-aa$parameter
  result[p,6]<-aa$p.value
  
}
a<-result
a$decile<-1:7

#Create 12 months plot with ggplot2
p <- ggplot(a, aes(x = decile, y = beta)) 
plot2<- p + 
  geom_hline(aes(yintercept = 0), size = 0.5, linetype = "dashed",col="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), data=a, width = .2, color = "black",size=1) +
  geom_point(data=a,size = 3, bg = "black",pch=21,col="black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("B - NB length (mm)") +ylim(-0.42,0.55)+
  theme(axis.text=element_text(size=30,color='black')) +
  theme(axis.title=element_text(size=30,color='black')) +
  xlab('')+
  scale_x_continuous(breaks=c(1:7),labels=a$vert)


plot_grid(plot1+theme(legend.position="none"),align='v',
          plot2+theme(legend.position="none"),
          ncol=2,nrow=1) +
  theme(plot.margin = unit(c(0,0,0,2), "cm")) 
#save panel as PDF

##############################
#Model litter size and pup mass as function of mother's body length. Create plots for Figure 1B and 1D.
###############################
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(cowplot)

df<-read.delim('pupweightinfo30nov2019_outliersremoved.txt')
df<-df[,c('Weight','BodyLength','flag','Number_total','MotherID','LitterRef','DateOfBirth')]

#Remove outliers for birth weight:
df<-subset(df,Weight<=16 & Weight>=6) 
#Remove rows with missing data
df<-na.omit(df)
#Calculate pup weight controlling for litter size 
df$a<-df$Weight/df$Number_total

#Reduce data to litter for litter model
out2<- df %>% distinct(LitterRef, .keep_all = TRUE) #each row corresponds to a litter

#Model litter size controlling for whether litter was mother's first litter (flag)
model2 = lmer(Number_total ~ BodyLength + flag + (1|MotherID),data=out2,REML=TRUE) 
summary(model2)
# #Regress covariates to plot litter size residuals 
# model2b = lmer(Number_total ~ flag + (1|MotherID),data=out2,REML=TRUE) #each row corresponds to a litter

#Model pup mass, with mother's body length as fixed effect, controlling for litter size (Number_total) and if pup was born in the mother's first litter (flag). Mother id and litter are included as random effects.
model1 = lmer(Weight ~ BodyLength + flag + Number_total + (1|MotherID) + (1|LitterRef),data=df,REML=TRUE) #the 1 indicates that an intercept is to be fitted for each level of the random variable.
summary(model1)

#Estimate the % weight gained per 1 cm increase in mother's mass 
a<-summary(model1)
b<-a$coefficients[2,1] #Increase in weight (grams) per 1 cm increase in mother's mass
c<-mean(df$Weight)
b/c

#Plot 
library(beeswarm)
library(ggplot2)
library(cowplot)

#Create color for plotting
r<-adjustcolor( '#41006A', alpha.f = 0.8)

#Create Main Figure panel 1D
p <- ggplot(out2, aes(x = BodyLength, y = Number_total)) 
plot1D<- p + 
  geom_point(data=out2,size = 4, pch=21,col='black',bg=r,) +
  geom_smooth(method='lm', size = 0.5, linetype = "solid",col="black",level=0)+
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab('Litter size') +
  theme(axis.text=element_text(size=30,color='black')) +
  theme(axis.title=element_text(size=30,color='black')) +
  xlab('Body length (cm)')+
  scale_y_continuous(breaks=seq(2,12,by=2),labels=seq(2,12,by=2))+
  scale_x_continuous(breaks=c(16:21),labels=c(16:21))

#Read in data for Figure panel 1B
this<-read.table('~/Desktop/github_dmr_bone/Molerat_xray_measures.txt',header=T,sep="\t")

#Create dataframe for queen data
m<-'LV5'
result<-as.data.frame(matrix(nrow=6,ncol=2));rownames(result)<-c(0,4,8,12,16,22);colnames(result)<-c('mean','se')
mo<-c(0,4,8,12,16,22)
for (p in 1:nrow(result)) {
  mon<-mo[p]
  t<-subset(this,LHAM_mo_experiment==mon & Treatment=='BREEDER')
  result[p,1]<-mean(t[,m])
  result[p,2]<-sd(t[,m], na.rm=TRUE)/sqrt(length(t[,m][!is.na(t[,m])]))
}
a<-result
a$month<-as.numeric(rownames(a))
a$upper<-a$mean+a$se
a$lower<-a$mean-a$se

#Create dataframe for nonbreeders data
result2<-as.data.frame(matrix(nrow=6,ncol=2));rownames(result2)<-c(0,4,8,12,16,22);colnames(result2)<-c('mean','se')
mo<-c(0,4,8,12,16,22)
for (p in 1:nrow(result2)) {
  mon<-mo[p]
  t<-subset(this,LHAM_mo_experiment==mon & Treatment=='0BREEDER')
  result2[p,1]<-mean(t[,m])
  result2[p,2]<-sd(t[,m], na.rm=TRUE)/sqrt(length(t[,m][!is.na(t[,m])]))
}
b<-result2
b$month<-as.numeric(rownames(b))
b$upper<-b$mean+b$se
b$lower<-b$mean-b$se

#Create Main Figure panel 1B
p <- ggplot(a, aes(x = month, y = mean)) 
plot1B<- p + 
  geom_errorbar(aes(ymax=upper, ymin=lower), data=b, width = .7, color = "#009193",size=1) +
  geom_point(data=b,size = 4, bg = "#009193",pch=23,col="black") +
  geom_errorbar(aes(ymax=upper, ymin=lower), data=a, width = .7, color = "#41006A",size=1) +
  geom_point(data=a,size = 4, bg = "#41006A",pch=21,col="black") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab('LV5') +
  theme(axis.text=element_text(size=30,color='black')) +
  theme(axis.title=element_text(size=30,color='black')) +
  xlab('Months of treatment')+
  scale_x_continuous(breaks=seq(0,20,by=4),labels=seq(0,20,by=4))


plot_grid(plot1B+theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
          plot1D+theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
          label_size = 32,
          nrow=2,ncol=1,align='v')
#saved as 12x6 PDF 
