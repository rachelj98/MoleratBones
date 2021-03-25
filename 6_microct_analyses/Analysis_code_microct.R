#######################################################################
#Assess effects of breeder status and number offspring on growth plate fusion and chondrocyte column density in tibia and LV7
#######################################################################

setwd('~/Desktop/github_dmr_bone/microct_analyses/')

#Read in tibia growth plate data
this<-read.table('Table_S9_growthplate.txt',header=T,sep='\t')
th<-subset(this,Proportion_growth_plate_fusion_tibia!='NA');th<-droplevels(th)
q<-subset(th,Breeding_status=='breeder')
n<-subset(th,Breeding_status=='non-breeder')
aa<-length(which(q$Proportion_growth_plate_fusion_tibia==0));aa #queen successes
bb<-length(which(q$Proportion_growth_plate_fusion_tibia!=0));bb #queen failures
cc<-length(which(n$Proportion_growth_plate_fusion_tibia==0));cc #expected successes (non-breeders)
dd<-length(which(n$Proportion_growth_plate_fusion_tibia!=0));dd #expected failures (non-breeders)
pp<-cc/(cc+dd);pp

# binom.test(x, p = 0.5, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
# x: vector of 2 specifying number successes, number failures
# p: hypothesized probability of successes

#Here, "success" is unfused tibia growth plate. 1 queen had success, 6 queens had failure. In nonbreeders, probability of success was 3/5=0.6:
binom.test(c(aa,bb), p = pp, alternative = c("two.sided"),conf.level = 0.95)

#Are ages different between queens and helpers for tibia growth plate dataset?
t.test(q$Age_in_years,n$Age_in_years,alternative='two.sided')

#Test whether proportion tibia growth plate fusion is explained by number offspring, controlling for age
l<-lm(th$Proportion_growth_plate_fusion_tibia~th$Number_offspring_born+th$Age_in_years)
summary(l)

#Test whether number chondrocyte columns per mm tibia growth plate is explained by number offspring, controlling for age
th<-subset(this,Chondrocyte_columns_per_mm_tibia!='NA');th<-droplevels(th)
l<-lm(th$Chondrocyte_columns_per_mm_tibia~th$Number_offspring_born+th$Age_in_years)
summary(l)

#Now test whether growth plate differs with breeding status in LV7:
this<-read.table('Table_S9_growthplate.txt',header=T,sep='\t')
th<-subset(this,Proportion_growth_plate_fusion_LV7_average!='NA');th<-droplevels(th)
head(th)
q<-subset(th,Breeding_status=='breeder')
n<-subset(th,Breeding_status=='non-breeder')
aa<-length(which(q$Proportion_growth_plate_fusion_LV7_average==0));aa #queen successes
bb<-length(which(q$Proportion_growth_plate_fusion_LV7_average!=0));bb #queen failures
cc<-length(which(n$Proportion_growth_plate_fusion_LV7_average==0));cc #expected successes (non-breeders)
dd<-length(which(n$Proportion_growth_plate_fusion_LV7_average!=0));dd #expected failures (non-breeders)
pp<-cc/(cc+dd);pp

binom.test(c(aa,bb), p = pp, alternative = c("two.sided"),conf.level = 0.95)

#Are ages different between queens and helpers for LV7 growth plate dataset?
t.test(q$Age_in_years,n$Age_in_years,alternative='two.sided')

#########################################
#Figure 3c and 3d. Plot growth plate fusion and chondrocyte activity as function of number offspring
#########################################
library(beeswarm)
this<-read.table('Table_S9_growthplate.txt',header=T,sep='\t')
head(this)
this$queen2<-ifelse(this$Breeding_status=='non-breeder','Nonbreeder','Queen')
this$queen2 <- factor(this$queen2, levels = c("Nonbreeder", "Queen"))
r<-adjustcolor( '#009193', alpha.f = 0.8)
k<-adjustcolor( '#41006A', alpha.f = 0.8)

#Create grayscale according to number offspring
this$gray<-100-this$Number_offspring_born*5
this$gray2<-paste('gray',this$gray,sep='')

th<-subset(this,Proportion_growth_plate_fusion_tibia!='NA');th<-droplevels(th)

par(mfrow=c(1,2))
boxplot(th$Proportion_growth_plate_fusion_tibia*100~th$queen2,cex.axis=1.4,cex.lab=1.4,ylab="% growth plate fusion",border=c('#009193','#41016A'))
beeswarm(th$Proportion_growth_plate_fusion_tibia*100~th$queen2, pch = 21, pwbg=th$gray2,cex=1.5,add = TRUE)

#Chondrocyte columns/mm gp tibia
th<-subset(this,Chondrocyte_columns_per_mm_tibia!='NA');th<-droplevels(th)
boxplot(th$Chondrocyte_columns_per_mm_tibia~th$queen2,cex.axis=1.4,cex.lab=1.4,ylab="Chondrocyte columns/mm gp",border=c('#009193','#41016A'))
beeswarm(th$Chondrocyte_columns_per_mm_tibia~th$queen2,pch = 21,pwbg=th$gray2,cex=1.5,add = TRUE)

b<-subset(th,queen2=='Queen')
nb<-subset(th,queen2=='Nonbreeder')
#saved as 5x8 PDF

#########################################
# Supplementary Figure. Plot LV7 growth plate fusion and chondrocyte activity as function of number offspring
#########################################
library(beeswarm)
this<-read.table('Table_S9_growthplate.txt',header=T,sep='\t')
head(this)
this$queen2<-ifelse(this$Breeding_status=='non-breeder','Nonbreeder','Queen')
this$queen2 <- factor(this$queen2, levels = c("Nonbreeder", "Queen"))
r<-adjustcolor( '#009193', alpha.f = 0.8)
k<-adjustcolor( '#41006A', alpha.f = 0.8)

#Create grayscale according to number offspring
this$gray<-100-this$Number_offspring_born*5
this$gray2<-paste('gray',this$gray,sep='')

th<-subset(this,Proportion_growth_plate_fusion_LV7_average!='NA');th<-droplevels(th)

par(mfrow=c(1,2))
boxplot(th$Proportion_growth_plate_fusion_LV7_average*100~th$queen2,cex.axis=1.4,cex.lab=1.4,ylab="% growth plate fusion",border=c('#009193','#41016A'))
beeswarm(th$Proportion_growth_plate_fusion_LV7_average*100~th$queen2, pch = 21, pwbg=th$gray2,cex=1.5,add = TRUE)

#Chondrocyte columns/mm gp tibia
this$columns_lv7<-(this$Chondrocyte_columns_per_mm_LV7_caudal+this$Chondrocyte_columns_per_mm_LV7_cranial)/2
th<-subset(this,columns_lv7!='NA');th<-droplevels(th)
th$columns_lv7<-as.numeric(as.character(th$columns_lv7))
boxplot(th$columns_lv7~th$queen2,cex.axis=1.4,cex.lab=1.4,ylab="Chondrocyte columns/mm gp",border=c('#009193','#41016A'))
beeswarm(th$columns_lv7~th$queen2,pch = 21,pwbg=th$gray2,cex=1.5,add = TRUE)
#saved as 5x8 PDF

###########################################
#Paired t-tests to test effect of breeder status on trabecular bone volume/total volume in femur, tibia, lv6, lv7
###########################################
#Run the following python script on the PDFs that were output from the microCT scanner. This extracts the values from the PDFs and saves them to a text file.
#We provide an example input PDF output from the vivaCT scanner ('C0000351_G10F023_FEM_MORPHO.PDF'), and the output files for femur, tibia, lv6, and lv7 (eg, 'bvtv_femur_30oct2019.txt').
python bvtv_extraction_2.py #Make sure script is in same directory as target files

#Perform paired t-tests in R
setwd('~/Desktop/github_dmr_bone/microct_analyses')
info<-read.table('paired_animals_breeder_status.txt',sep='\t',header=T)
rownames(info)<-info$Name
out<-as.data.frame(matrix(ncol=4,nrow=1))
rownames(out)<-colnames(queen[c(10)])
colnames(out)<-c('femur','tibia','lv6','lv7')

for (p in c('femur','tibia','lv6','lv7')){
  bone<-read.table(paste0('bvtv_trabec_',p,'_30oct2019.txt'),sep='\t',header=T); bone<-bone[order(bone$Bone),]; rownames(bone)<-bone$Bone
  bone<-merge(bone,info,by.x='Bone',by.y='Name')
  queen=subset(bone,breeder_status=="q")
  nonbreeder=subset(bone,breeder_status=="n")
    t1=t.test(queen$BV.TV.1,nonbreeder$BV.TV.1,paired=T)
    out[1,p]<-t1$p.value
} #close p
out #P-values from paired t-tests of trabecular bone volume/total volume for each bone type

###################################################
#Test effects of breeder status on cortical area/total area at the midshaft of femur and LV6 
###################################################
setwd('~/Desktop/github_dmr_bone/microct_analyses/')
#Femur:
df1<-read.table("femur_2D_midsection_paired_breeders_nonbreeders.txt",sep="\t",header=T)
df1$queen<-as.factor(df1$queen)
queen=subset(df1,queen=="q");nonbreeder=subset(df1,queen=="n")
t1=t.test(queen$CA.TA,nonbreeder$CA.TA,paired=T);t1

#LV6:
df2<-read.table("LV6_2D_midsection_paired_breeders_nonbreeders.txt",sep="\t",header=T)
df2$queen<-as.factor(df2$queen)
queen=subset(df2,queen=="q");nonbreeder=subset(df2,queen=="n")
t2=t.test(queen$CA.TA,nonbreeder$CA.TA,paired=T);t2

############################################################
#Test whether cortical bone density differs between queens and nonbreeders at femoral midshaft
############################################################
setwd('~/Desktop/github_dmr_bone/microct_analyses/')
df<-read.table('bvtv_fem100mid_data.txt',sep='\t',header=T)
head(df)
b<-subset(df,Queen==1)
n<-subset(df,Queen==0)

#Apparent density
t<-t.test(b$mean_apparent_density_tv,n$mean_apparent_density_tv,paired=T)
t

#Material density
t<-t.test(b$mean_material_density_bv,n$mean_material_density_bv,paired=T)
t

###################################################
#Test for effect of breeder status on cortical area/total area, marrow area, and periosteal area. Create corresponding plots for main figure 3e, 3f, 3g
###################################################
setwd('~/Desktop/github_dmr_bone/microct_analyses/')
this<-read.table('femur_2D_midsection_paired_breeders_nonbreeders.txt',sep='\t',header=T)
queen<-subset(this,queen=='q')
nonbreeder<-subset(this,queen=='n')

t1=t.test(queen$CA.TA,nonbreeder$CA.TA,paired=T);t1
t2=t.test(queen$marrow_area,nonbreeder$marrow_area,paired=T);t2
t3=t.test(queen$periosteal_area,nonbreeder$periosteal_area,paired=T);t3

tt1=paste('p =',as.character(round(t1$p.value,5)))
tt2=paste('p =',as.character(round(t2$p.value,5)))
tt3=paste('p =',as.character(round(t3$p.value,5)))

#Across all queens, test whether femoral cortical area/total area is explained by recent litter size (pups born within the past 30 days) or total number offspring
this<-read.table('femur_2D_midsection_all_animals.txt',sep='\t',header=T)
this$CA.TA<-this$CA/this$TA
queen<-subset(this,queen=='q')
l=lm(queen$CA.TA~queen$pups_30days_or_younger);summary(l)
l=lm(queen$CA.TA~queen$number_offspring_born);summary(l)

#Creat plots
library(coxme)
library(ggplot2)
library("gridExtra")
library(cowplot)
this<-read.table('femur_2D_midsection_paired_breeders_nonbreeders.txt',sep='\t',header=T)
rownames(this)<-this$Name
this$queen<-as.factor(this$queen)
this$color=c("#F8766D","#C77CFF")[this$queen]
p0 <- ggplot(this, aes(x = queen, y = CA.TA))+
  #  ggtitle(tt1)+
  ylab("CA/TA")+
  scale_x_discrete(labels = c("N","B"))+
  theme(axis.title.x=element_blank())+
  geom_line(aes(group = pair))+ 
  theme(axis.text=element_text(size=24,colour='black'),axis.title=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1 <- ggplot(this, aes(x = queen, y = marrow_area))+
  #  ggtitle(tt2)+
  ylab("Marrow area (mm2)")+
  scale_x_discrete(labels = c("N","B"))+
  theme(axis.title.x=element_blank())+
  geom_line(aes(group = pair))+
  theme(axis.text=element_text(size=24,colour='black'),axis.title=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2 <- ggplot(this, aes(x = queen, y = periosteal_area))+
  #  ggtitle(tt3)+
  ylab("Periosteal area (mm2)")+
  scale_x_discrete(labels = c("N","B"))+
  theme(axis.title.x=element_blank())+
  geom_line(aes(group = pair))+
  theme(axis.text=element_text(size=24,colour='black'),axis.title=element_text(size=24))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_grid(p0+theme(legend.position="none"),align='v',
          p1+theme(legend.position="none"),
          p2+theme(legend.position="none"),
          ncol=3,nrow=1)
#saved as 6x18 PDF

#######################################
#Evaluate effect of number of offspring on cortical thickness along decile sections of the femur. Create plot for Main Figure 4B.
#######################################
setwd('~/Desktop/github_dmr_bone/microct_analyses')
locs<-read.table('femur_deciles/GRF002_canonical_locs.txt',header=F,sep=",");colnames(locs)<-c('x','y','z')
head(locs)

l1=read.table('femur_deciles/DVF007_thick.txt');colnames(l1)<-'DVF007'
l2=read.table('femur_deciles/DVF008_thick.txt');colnames(l2)<-'DVF008'
l1b=read.table('femur_deciles/G10F023_thick.txt');colnames(l1b)<-'G10F023'
l2b=read.table('femur_deciles/G10F024_thick.txt');colnames(l2b)<-'G10F024'
l3=read.table('femur_deciles/GRF001_thick.txt');colnames(l3)<-'GRF001'
l4=read.table('femur_deciles/GRF002_thick.txt');colnames(l4)<-'GRF002'
l3b=read.table('femur_deciles/GRF003_thick.txt');colnames(l3b)<-'GRF003'
l4b=read.table('femur_deciles/GRF004_thick.txt');colnames(l4b)<-'GRF004'
l5=read.table('femur_deciles/LAF007_thick.txt');colnames(l5)<-'LAF007'
l6=read.table('femur_deciles/LAF008_thick.txt');colnames(l6)<-'LAF008'
l7=read.table('femur_deciles/MAF003_thick.txt');colnames(l7)<-'MAF003'
l8=read.table('femur_deciles/MAF006_thick.txt');colnames(l8)<-'MAF006'
l9=read.table('femur_deciles/NOF007_thick.txt');colnames(l9)<-'NOF007'
l10=read.table('femur_deciles/NOF009_thick.txt');colnames(l10)<-'NOF009'
l11=read.table('femur_deciles/TIF003_thick.txt');colnames(l11)<-'TIF003'
l12=read.table('femur_deciles/TIF005_thick.txt');colnames(l12)<-'TIF005'
l13=read.table('femur_deciles/WAF011_thick.txt');colnames(l13)<-'WAF011'
l14=read.table('femur_deciles/WAF014_thick.txt');colnames(l14)<-'WAF014'
l<-cbind(l1,l2,l1b,l2b,l3,l4,l3b,l4b,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14)
lall<-cbind(l,locs)
lall2<- lall[order(lall$z),] 
head(lall2)
l<-lall2[,c(1:18)]

#Break rows into 10 sections
a<-nrow(l)
a/10
dec=rep(0,10)
for (p in c(1:10)){
  dec[p]<-a*p/10
  dec[p]<-round(dec[p],0)
}
dec<-c(1,dec)

#Remove top and bottom of bone,then cut into 10 new slices (this is done to reduce thickness measurement error from inclusion of trabecular bone)
l<-l[dec[2]:dec[10],]
a<-nrow(l)
a/10
dec=rep(0,10)
for (p in c(1:10)){
  dec[p]<-a*p/10
  dec[p]<-round(dec[p],0)
}
dec<-c(1,dec)

#Univariate modeling
result<-as.data.frame(matrix(nrow=10,ncol=3))
colnames(result)<-c('beta','se','p')
info<-read.table('femur_2D_midsection_paired_breeders_nonbreeders.txt',header=T,sep="\t")
rownames(info)<-info$Name
info<-info[,c('queen','pair','number_offspring_born')]
head(info)
#Take average thickness from each decile for each bone, save in info
for (p in c(1:10)){
  t<-l[c(dec[p]:dec[p+1]),]
  t2<-colSums(t)/nrow(t)
  info[,paste0('dec',p)]<-t2
}

info2<-info[,c(4:13,2,3)]
head(info2)

#Subset to pairs in which queen had less than 6, or 6 or more offspring. Presented as supplemental table.
#info2<-subset(info2,pair != 'DV' & pair != 'G10' & pair != 'GR_b' & pair != 'LA' & pair != 'WA') #Subset to pairs in which queen had at least 6 total offspring
#info2<-subset(info2,pair == 'DV' & pair == 'G10' & pair == 'GR_b' & pair == 'LA' & pair == 'WA') #Subset to pairs in which queen had less than 6 total offspring

#Include pair as random effect
library(lme4)
library(lmerTest)
result<-as.data.frame(matrix(nrow=10,ncol=3))
colnames(result)<-c('beta','se','p')

for (p in c(1:10)){
d1 = lmer(info2[,p] ~ number_offspring_born + (1|pair),data=info2,REML=TRUE) #the 1 indicates that an intercept is to be fitted for each level of the random variable. Each row corresponds to pup
dd=summary(d1)
dd1<-round(dd$coefficients[2,1],4)
dd2<-dd$coefficients[2,2]
dd3<-round(dd$coefficients[2,5],5)
result[p,]<-c(dd1,dd2,dd3)
}
result2<-result[10:1,] #flip because "bottom" in data frame corresponds to top (proximal end) of femur
write.table(result2,file='lm_femur_deciles_mixed_numberoffspring_all27feb2020.txt',sep='\t')

#Create plot, Main Figure 4B, showing effect of number of offspring on cortical thickness along decile sections of the femur.
library(ggplot2)
library(cowplot)
a<-read.table('lm_femur_deciles_mixed_numberoffspring_all27feb2020.txt',sep='\t')
a$decile<-1:10
a$upper<-a$beta+a$se
a$lower<-a$beta-a$se
p <- ggplot(a, aes(x = decile, y = beta)) 
plot1<- p + 
  geom_hline(aes(yintercept = 0), size = 1, linetype = "dashed",col="white") +
  geom_errorbar(aes(ymax=upper, ymin=lower), data=a, width = .2, color = "white",size=1) +
  geom_point(data=a,size = 3, bg = "white",pch=21,col="white") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("Effect of # offspring\non cortical thickness") +
  theme(axis.text=element_text(size=30,color='white')) +
  theme(axis.title=element_text(size=30,color='white')) +
  xlab('Decile')+
  scale_x_continuous(breaks=c(1:10),labels=c(1:10))+
  theme(panel.background = element_rect(fill = '#323232', colour = 'white',size=2))+
  theme(plot.background = element_rect(fill = '#323232', colour = 'white'))
plot_grid(plot1+theme(legend.position="none"),
          nrow=1,ncol=1)
#save as 4x15 PDF 'decileplot2'

###################################################
#Test effects of body mass and number of offspring on cortical area (CA), Imin (2nd moment of inertia), and J (polar second moment of area)
###################################################
setwd('~/Desktop/github_dmr_bone/microct_analyses')
#Nonbreeders:
df<-read.table("femur_2D_midsection_all_animals.txt",sep="\t",header=T)
df<-subset(df,queen=='n')
ca_nonbreeders<-lm(df$CA~df$weight);summary(ca_nonbreeders)
iminnonbreeders<-lm(df$Imin~df$weight);summary(iminnonbreeders)
j_nonbreeders<-lm(df$J~df$weight);summary(j_nonbreeders)

#Queens:
df<-read.table("femur_2D_midsection_all_animals.txt",sep="\t",header=T)
df<-subset(df,queen=='q')
result<-lm(df$Imin~df$weight);summary(result)
result<-lm(df$J~df$weight);summary(result)
result<-lm(df$Imin~df$number_offspring_born);summary(result)
result<-lm(df$J~df$number_offspring_born);summary(result)
result<-lm(df$CA~df$number_offspring_born);summary(result)
result<-lm(df$CA~df$weight);summary(result)



###########################
#7may2020 Look animal max load vs cort data from Jepsen et al 2003
###########################
#Look at relationship, at femoral midsection, between cortical area and max load 
jj<-read.table('~/Desktop/tunglab/Kalahari/Molerats/molerat_bonebreaking/jepsenunadjusted.txt',sep='\t',header=T)
par(mfrow=c(1,1))
#Create plot, coloring points according to mouse strain
ccol<-c('blue','green','deepskyblue3','chartreuse','darkblue','orange','darkgreen','red')
plot(jj$Area,jj$Max_Ld,pch=19,col=ccol[as.factor(jj$strain)],xlab='Cortical area',ylab="Max load (Newtons)",cex.axis=1.4,cex.lab=1.4,main='R2=0.88, p<2.2e-16')
l<-lm(jj$Max_Ld~jj$Area);abline(l)
summary(l)




######################################################################################
#Perform Cox proportional hazards models, and corresponding main figure
######################################################################################
library(survival)
library(survminer)
library(Hmisc)
library(rms)

info<-read.table('~/Desktop/tunglab/Kalahari/Molerats/molerat_microct/full_3D_reconstructions/FEM/midshaft_analysis/second_moments_femur_midshaft_all_scaled30Nov2019.txt',header=T,sep="\t")
#info<-subset(info,Name!='Z3F022')
info<-subset(info,Name!='Z3F022' & Name!='L2F009' & Name!='L2F010')
#info<-read.table('~/Desktop/tunglab/Kalahari/Molerats/molerat_microct/full_3D_reconstructions/FEM/midshaft_analysis/second_moments_femur_midshaft_paired_scaled31Oct2019.txt',header=T,sep="\t")
j<-read.table('~/Desktop/tunglab/Kalahari/Molerats/molerat_bonebreaking/jepsenunadjusted.txt',sep='\t',header=T)
summary(info$CA)
mean(info$CA)
mean(j$Area)

#Obtain scaling factor of Damaraland mole-rat cortical area relative to mouse cortical area
mean(info$CA)/mean(j$Area)

#Scale down mole-rat femur areas to mouse by scaling factor of 3.6
info$Area<-info$CA/3.6 

#Get predicted mole-rat max loads according to linear relationship of max load vs cortical area in mouse
l<-lm(j$Max_Ld~j$Area)
#a provides the info we need to create the linear equation of max load vs cortical area.
a<-summary(l)
#b is the increase in max load for every 1 unit increase in cortical area (ie, the slope):
b<-a$coefficients[2,1]
#c is the y-intercept of the line:
c<-a$coefficients[1,1]
#Obtain predicted mole-rat max loads:
info$Max_Ld<-info$Area*b+c

info<-info[order(info$queen),]
#info$strain<-c(rep('y',23),rep('z',13)) #if using all midshafts
info$strain<-c(rep('y',21),rep('z',13)) #if removing L2F009 and L2F010 midshafts
info$spp<-'dmr'
info$pch1<-15

df<-info
head(df)
mean(df$Max_Ld)
nb<-subset(df,queen=='n')
df$Max_Ld2<-df$Max_Ld/median(nb$Max_Ld)
df$Area<-df$CA
df$censored<-df$Area!=100
censored<-df$censored
lifespan<-c(unlist(df$Max_Ld2))
df$lifespan<-lifespan
df$status<-2
#Making a 'survival object' which is the required response variable for a cox ph model
life.surv<-Surv(time= lifespan,event = censored) #specify that all measures are "events"/deaths
df$lifesurv<-life.surv

cph.model<-coxph(life.surv~df$queen)
summary(cph.model)
cox.zph(cph.model)

#Fit # offspring
cph.model<-coxph(life.surv~df$number_offspring_born)
summary(cph.model)
lifespan<-c(unlist(df$Max_Ld2))
df$lifespan<-lifespan
df$status<-2
#Making a 'survival object' which is the required response variable for a cox ph model
life.surv<-Surv(time= lifespan,event = censored) #specify that all measures are "events"/deaths
df$lifesurv<-life.surv
df_bla<-df

#
#Remove barren queens, and compare fertile queens to nonbreeders
df<-df_bla
#Create column specifying queens that had at least 6 offspring ('c') vs queens with less than 6 offspring ('b') vs nonbreeders 'a'.
df$fruitful2<-ifelse(df$number_offspring_born>=6,'c',ifelse(df$queen=='q','b','a'))
#Remove queens with less than 6 offspring
df<-subset(df,fruitful2!='b')
nb<-subset(df,queen=='n')
df$Max_Ld2<-df$Max_Ld/median(nb$Max_Ld)
df$Area<-df$CA
censored<-df$censored
lifespan<-c(unlist(df$Max_Ld2))
df$lifespan<-lifespan
df$status<-2
#Making a 'survival object' which is the required response variable for a cox ph model
life.surv<-Surv(time= lifespan,event = censored) #specify that all measures are "events"/deaths
df$lifesurv<-life.surv

#Fit fruitful queens vs nonbreeders
fit<-survfit(Surv(lifespan,censored)~queen,data=df)
#fit<-coxph(Surv(lifespan,censored)~queen,data=df)
fit

cph.model<-coxph(life.surv~df$queen)
summary(cph.model)
cox.zph(cph.model)

#
#


#Plot
#Remove barren, and compare fruitful to nonbreeders
df<-df_bla
df$fruitful2<-ifelse(df$number_offspring_born>=6,'c',ifelse(df$queen=='q','b','a'))
df2<-subset(df,fruitful2=='c') #all fruitful breeders, n=7
df2<-droplevels(df2)
df$fruitful2<-ifelse(df$queen=='q','b','a')
df3<-rbind(df,df2)

fit<-survfit(Surv(lifespan,censored)~fruitful2,data=df3)
surv_summary(fit,data=df3)
censored<-df3$censored
surv_diff <- survdiff(Surv(lifespan, censored) ~ fruitful2, data = df3)
surv_diff
pp<-ggsurvplot(fit,xlab='Relative force',conf.int = F, pval = F,palette=c('#009193','#41006A','#941100'))
fit
pp$plot
pp$plot + geom_vline(xintercept = 1,lty='dashed',col='#009193',lwd=1)+geom_vline(xintercept = 0.665,lty='dashed',col='#941100',lwd=1)+
  geom_vline(xintercept = 0.903,lty='dashed',col='#41006A',lwd=1)+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )+
  theme(axis.text.x = element_text(color="black",size=12),
        axis.text.y = element_text(color="black",size=12),
        axis.ticks = element_line(color="black"),
        axis.title.x = element_text(colour = "black",size=12),
        axis.title.y = element_text(colour = "black",size=12),
        axis.line = element_line(color="black"))

#saved as 3.5x6 'surv plot nonbreeder vs fruitful relative colored'
#saved as 3.5x6 'surv plot nonbreeder vs fruitful relative bright 10aug2020'

#saved as 3.5x6 'surv plot nonbreeder vs breeder vs fruitful' 9Mar2021
lifespan<-c(unlist(df2$Max_Ld2))
life.surv<-Surv(time= lifespan,event = censored) #specify that all measures are "events"/deaths
cph.model<-coxph(life.surv~df2$fruitful2)
summary(cph.model)
cox.zph(cph.model)

