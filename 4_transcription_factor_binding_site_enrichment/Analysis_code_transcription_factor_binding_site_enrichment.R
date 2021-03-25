###################################################
#This script starts with fastq files from ATAC-Seq libraries, calls open chromatin peaks, and tests whether breeder status-associated genes are enriched for trancription factor binding motifs.
#The script also creates plots for Main Figure 2a and 2b.
###################################################

###################################################
#trimmap.sh Trim and map OMNI-ATAC-Seq reads
for f in L31305 L31306 L31307 L31314 L31315 L31318 L31319; do cat trimmap.sh | sed -e s/FILEINFO/$f/g > trimmap_$f.sh; sbatch trimmap_$f.sh; done
###################################################

#!/bin/bash
#SBATCH -n 8 # Number of cores
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=32G 
#SBATCH -o trimmap_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e trimmap_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load TrimGalore
module load bwa
trim_galore -q 20 -o ../rawreads --fastqc -a AGATCGGAAGAGC --stringency 2 --length 25 --paired ../rawreads/FILEINFO_R1.fq ../rawreads/FILEINFO_R2.fq
bwa mem -t 8 /data/tunglab/raj/references/dmr1.0_bwa/dmr1.0 ../rawreads/FILEINFO_R1_val_1.fq ../rawreads/FILEINFO_R2_val_2.fq > ../mapped/FILEINFO.sam

###################################################
#unique.sh Filter reads to keep reads that mapped uniquely to genome.
for f in L31304 L31305 L31306 L31307 L31314 L31315 L31318 L31319; do cat unique.sh | sed -e s/FILEINFO/$f/g > unique_$f.sh; sbatch unique_$f.sh; done
###################################################
#!/bin/bash 
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-96:00 # Runtime in D-HH:MM 
#SBATCH --mem-per-cpu=1G # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o FILEINFOmerge.out # File to which STDOUT will be written 
#SBATCH -e FILEINFOmerge.err # File to which STDERR will be written 
#SBATCH --mail-type=BEGIN,END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load samtools
samtools view -bSq 1 ../mapped/FILEINFO.sam > ../mapped/FILEINFO.uniq.bam

###################################################
#merge.sh Merge bam files
###################################################

#!/bin/bash 
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-96:00 # Runtime in D-HH:MM 
#SBATCH --mem-per-cpu=1G # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o FILEINFOmerge.out # File to which STDOUT will be written 
#SBATCH -e FILEINFOmerge.err # File to which STDERR will be written 
#SBATCH --mail-type=BEGIN,END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load samtools
~/samtools-1.9/samtools merge ../mapped/merged.bam ../mapped/L31304.uniq.bam ../mapped/L31305.uniq.bam ../mapped/L31306.uniq.bam ../mapped/L31307.uniq.bam ../mapped/L31314.uniq.bam ../mapped/L31315.uniq.bam ../mapped/L31318.uniq.bam ../mapped/L31319.uniq.bam 

###################################################
#macs.sh Perform peak calling. Output file of interest is 'NA_peaks.narrowPeak', which we have posted.
###################################################
#!/bin/bash 
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 0-96:00 # Runtime in D-HH:MM 
#SBATCH --mem-per-cpu=64G # Memory pool for all cores (see also --mem-per-cpu) 
#SBATCH -o FILEINFOmacs2.out # File to which STDOUT will be written 
#SBATCH -e FILEINFOmacs2.err # File to which STDERR will be written 
#SBATCH --mail-type=BEGIN,END # Type of email notification- BEGIN,END,FAIL,ALL 
#SBATCH --mail-user=racheljohnston7@gmail.com # Email to which notifications will be sent
module load macs2
macs2 callpeak -t ../mapped/merged.bam --nomodel --outdir ../macsout --keep-dup all -q 0.05 -f BAMPE -n outputfileprefix

###################################################
#Footprinting analysis
###################################################
# The *.uniq.bam files and output from macs2, above, were used as input for footprinting analyses. 
# To perform the footprinting analyses, see the code 'tf_analysis.py' in the github subfolder 'footprinting'

###################################################
#Prepare input files for HOMER-- bed files specifying genomic coordinates upstream of breeder status associated genes that also lie within open chromatin regions. 
###################################################

#Read in gene expression linear modeling results (see Analysis_code_gene_expression_linear_modeling.R for how this file was generated)
real<-read.table("~/Desktop/github_dmr_bone/mediation_analysis/model_bone_bonequ_bonetreat_dayscult_EMMREML_regrlitterpropinfeatvoom.txt",sep="\t",header=T)
rownames(real)=real$Row.names

#Read in transcription start sites for all annotated transcripts of DMR genome. 
#I generated this file using Ensembl biomart, a very handy webtool:https://m.ensembl.org/biomart/martview/25ed49101f0a60e340dd8e334ffd2f0b
#Subset significant genes to 4 quadrants 
tss<-read.delim("DMR_v1.0.92_tss_transcript.txt")
model1<-merge(real,tss,by.x="Row.names",by.y="Gene")
results<-as.data.frame(matrix(nrow=nrow(real),ncol=ncol(model1)));rownames(results)<-rownames(real);colnames(results)<-colnames(model1)

#For each gene, take the 5' most TSS
for (i in c(1:nrow(real))){
  mygene=as.character(rownames(real[i,]))
  output=subset(model1,Row.names==mygene)
  output$Strand<-as.character(output$Strand)
  minval<-output[which(output$TSS==(min(output$TSS))),]
  maxval<-output[which(output$TSS==(max(output$TSS))),]
  if (output$Strand[1]=="+") {
    x<-minval
  } else {
    x<-maxval
  }
  #If multiple transcripts have same TSS, randomly sample 1
  results[i,]<-as.matrix(x[sample(nrow(x),1),]) 
} #close i

#Subset the genes to those that were significantly differentially expressed with breeder status in either the long bones or lumbar vertebrae.
#Divide these genes into four quadrants according to up- or down-regulation in long bones or vertebrae. We're interested in quadrant 2 (top right quadrant; genes up-regulated with breeder status in long bones and vert)
model<-subset(results,fdr_Queen_at_lb<=0.1|fdr_Queen_at_vert<=0.1)
aa<-subset(model,beta_bone0.queen1<=0&beta_bone1.queen1>=0)
bb<-subset(model,beta_bone0.queen1>=0&beta_bone1.queen1>=0)
cc<-subset(model,beta_bone0.queen1>=0&beta_bone1.queen1<=0)
dd<-subset(model,beta_bone0.queen1<=0&beta_bone1.queen1<=0)

#Write bed file for HOMER input:
quad1b<-aa[,c('Chromosome','TSS','TSS_1','ID','Strand','Strand')]
quad2b<-bb[,c('Chromosome','TSS','TSS_1','ID','Strand','Strand')]
quad3b<-cc[,c('Chromosome','TSS','TSS_1','ID','Strand','Strand')]
quad4b<-dd[,c('Chromosome','TSS','TSS_1','ID','Strand','Strand')]
bg<-results[,c('Chromosome','TSS','TSS_1','ID','Strand','Strand')] #Background set of all genes, for testing for enrichment.

write.table(quad1b,file="quad1_fdr0.1_11dec2019.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(quad2b,file="quad2_fdr0.1_11dec2019.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(quad3b,file="quad3_fdr0.1_11dec2019.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(quad4b,file="quad4_fdr0.1_11dec2019.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(bg,file="bg_11dec2019.txt",sep="\t",quote=F,row.names=F,col.names=F) #Background set of all genes, for testing for enrichment.

#Add 2000 bp to each end of TSS:
cd ~/Desktop/hiseq_allbmscs/
bedtools slop -b 2000 -i quad1_fdr0.1_11dec2019.txt -g Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.chrom.sizes.txt > quad1_fdr0.1_11dec2019_2000.txt
bedtools slop -b 2000 -i quad2_fdr0.1_11dec2019.txt -g Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.chrom.sizes.txt > quad2_fdr0.1_11dec2019_2000.txt
bedtools slop -b 2000 -i quad3_fdr0.1_11dec2019.txt -g Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.chrom.sizes.txt > quad3_fdr0.1_11dec2019_2000.txt
bedtools slop -b 2000 -i quad4_fdr0.1_11dec2019.txt -g Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.chrom.sizes.txt > quad4_fdr0.1_11dec2019_2000.txt
bedtools slop -b 2000 -i bg_11dec2019.txt -g ~/references/Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.chrom.sizes.txt > bg_11dec2019_2000.txt #Background set of all genes, for testing for enrichment.

#Only keep regions within 2000 bp of TSS that overlap open chromatin regions
cd ~/Desktop/hiseq_allbmscs/
bedtools intersect -a quad1_fdr0.1_11dec2019_2000.txt -b NA_peaks.narrowPeak > quad1_fdr0.1_11dec2019_2000peaks.txt
bedtools intersect -a quad2_fdr0.1_11dec2019_2000.txt -b NA_peaks.narrowPeak > quad2_fdr0.1_11dec2019_2000peaks.txt
bedtools intersect -a quad3_fdr0.1_11dec2019_2000.txt -b NA_peaks.narrowPeak > quad3_fdr0.1_11dec2019_2000peaks.txt
bedtools intersect -a quad4_fdr0.1_11dec2019_2000.txt -b NA_peaks.narrowPeak > quad4_fdr0.1_11dec2019_2000peaks.txt
bedtools intersect -a bg_11dec2019_2000.txt -b NA_peaks.narrowPeak > bg_11dec2019_2000peaks.txt #Background set of all genes, for testing for enrichment.

###################################################
#Run HOMER to perform enrichment analysis of transcription factor binding sites near breeder status-associated genes (specifically, up-regulated in both the long bones and lumbar vertebrae).
#At the time of writing this (Feb 2021), Homer is available at http://homer.ucsd.edu/homer/motif/motifDatabase.html
#The code runs pretty quickly (~20 minutes), so I typically run this interactively on our computer cluster
###################################################
srun -p interactive --pty /bin/bash
export PATH=$PATH:PATH_TO_HOMER
genome=dmrv1.fa
findMotifsGenome.pl quad2_fdr0.1_11dec2019_2000peaks.txt $genome quad2_fdr0.1_11dec2019_2000peaks -size given -bg bg_11dec2019_2000peaks.txt -mknown known.motifs -nomotif
#The file 'known.motifs' was obtained from Homer, and is a compilation of transcription factor binding motifs across vertebrates.
#-nomotif option tells Homer to not try to annotate/predict novel motifs from the genome.

###################################################
#Perform Fisher's exact tests on Homer output to identify significantly enriched transription factor binding sites. This creates the Supplementary Table of significant TFBS's, which is plotted in Main Figure 2.
###################################################
# To run FET, you need 4 numbers: outquad2_noTF	outquad2_TF	inquad2_noTF	inquad2_TF
# inquad2_TF is provided in knownResults.txt as 'numb_of_Target_Sequences_with_Motif_of_'.
# inquad2_noTF is calculated as constant number of query genes (provided in column 6 of knownResults.txt) minus inquad_TF.
# outquad2_TF is provided in knownResults.txt as '# of Background Sequences with Motif'
# outquad2_noTF is calculated as constant number of background query genes (provided in column 8 of knownResults.txt) minus outquad_TF.
# I used excel to add these values to the Homer output file, knownResults.txt, and saved the updated file as 'knownResults_significant2.txt'.
library(ggplot2)
library(cowplot)
require(gridExtra)
df<-read.table('~/Desktop/hiseq_allbmscs/motifs/3feb2020/quad2_fdr0.1_11dec2019_2000peaks/knownResults_significant2.txt',header=T,sep='\t')
head(df)
terms=df
dim(terms)
enrichment<-as.data.frame(matrix(nrow=nrow(terms),ncol=11))
colnames(enrichment)<-c('term','OR','query','color','Homer_p','lower_or','upper_or','log2or','log2lower','log2upper','FET_p')

#Run Fisher's exact tests: 
for (p in 1:nrow(terms)) {
  info<-terms[p,]
  aa<-info$outquad2_noTF
  bb<-info$inquad2_noTF
  cc<-info$outquad2_TF
  dd<-info$inquad2_TF
  values = matrix(c(aa,bb,cc,dd), nrow = 2) #create contingency table
  values
  this<-fisher.test(values,alternative="two.sided")
  this
  dircolor='purple'
  result<-1:11
  result[1]<-as.character(info$shortname);result[2]<-this$estimate; result[3]='Quad2';result[4]<-dircolor
  result[5]<-info$P.value;result[6]<-this$conf.int[1]
  result[7]<-this$conf.int[2];result[8]<-log2(this$estimate)
  result[9]<-log2(as.numeric(result[6]));result[10]<-log2(as.numeric(result[7]))
  result[11]<-this$p.value
  enrichment[p,]<-result
}

#Check p-value for enrichment of androgen response elements (ARE) and estrogen response elements (ERE), two motifs of interest: 
enrichment[which(enrichment$term=='ARE'),]
enrichment[which(enrichment$term=='ERE'),]

#Save significant results for plotting
enrichment$Homer_p<-as.numeric(enrichment$Homer_p)
enrichment$FET_p<-as.numeric(enrichment$FET_p)
e<-subset(enrichment,FET_p<=0.01)
e<-e[,c(1,2,8,9,10,11)]
colnames(e)<-c('term','OR','log2or','log2lower95','log2upper95','FET_p')
write.table(e,'significant_tfbs_p0.01.txt',sep='\t')

#################################
#Plot gene ontology results (Main Figure 2a):
#################################
library(ggplot2)
require(gridExtra)
library(cowplot)

setwd('~/Desktop/github_dmr_bone/transcription_factor_binding_site_enrichment/')

e<-read.table('goresultsforplotting5may2020.txt',header=T,sep='\t')

#Subset to highest-level terms 
e<-subset(e,shallowestterm==1) 
e<-subset(e,OR!="Inf")
enrichment2<-e
enrichment2$log2or<-as.numeric(enrichment2$log2or);enrichment2$log2lower95<-as.numeric(enrichment2$log2lower95);enrichment2$log2upper95<-as.numeric(enrichment2$log2upper95)
enrichment2$termname <- factor(enrichment2$termname, levels = enrichment2$termname[order(enrichment2$log2or)])

#Create plot
p <- ggplot(enrichment2, aes(x = log2or, y = termname)) 
plot1<- p + geom_vline(aes(xintercept = 0), size = .25, linetype = "solid") +
  geom_errorbarh(aes(xmax = log2upper95, xmin = log2lower95), data=enrichment2,size = 1, height = .2, color = "#41016A",position = position_nudge(y = 0.1)) +
  geom_point(data=enrichment2,size = 6, bg = "#41016A",pch=21,col="#41016A",position = position_nudge(y = 0.1)) +
  #scale_x_continuous(breaks=seq(0,10,by=2), labels=seq(0,10,by=2), limits=c(0,10)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("") +
  theme(plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=16,color='black')) +
  theme(axis.title=element_text(size=16,color='black')) +
  xlab(expression('log'[2]*'(OR)')) +
  theme(plot.title = element_text(hjust = 0.5,size=24))

#################################
#Obtain odds ratios from TFBS table (Main Figure 2b):
#################################
#Read in TFBS results
e<-read.table('significant_tfbs_p0.01.txt',header=T,sep='\t')
enrichment2<-e
enrichment2$log2or<-as.numeric(enrichment2$log2or);enrichment2$log2lower95<-as.numeric(enrichment2$log2lower95);enrichment2$log2upper95<-as.numeric(enrichment2$log2upper95)
enrichment2$termname<-enrichment2$term
enrichment2$termname <- factor(enrichment2$termname, levels = enrichment2$termname[order(enrichment2$log2or)])

#Create plot
p <- ggplot(enrichment2, aes(x = log2or, y = termname)) 
plot2<- p + geom_vline(aes(xintercept = 0), size = .25, linetype = "solid") +
  geom_errorbarh(aes(xmax = log2upper95, xmin = log2lower95), data=enrichment2,size = 1, height = .5, color = "#41016A",position = position_nudge(y = 0.1)) +
  geom_point(data=enrichment2,size = 6, bg = "#41016A",pch=21,col="#41016A",position = position_nudge(y = 0.1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("") +
  theme(plot.title = element_text(size=16)) +
  theme(axis.text=element_text(size=16,color='black')) +
  theme(axis.title=element_text(size=16,color='black')) +
  xlab(expression('log'[2]*'(OR)')) +
  theme(plot.title = element_text(hjust = 0.5,size=24))

#Render Figure 2 a,b plots
plot_grid(plot1+theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
          plot2+theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")),
          label_size = 32,
          nrow=2,ncol=1,align='v')
# Save as 9x8 PDF

