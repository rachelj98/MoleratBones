#####################################################
#This file starts with fastq files and ends with transcriptome-wide gene expression modeling results.
#The 47 fastq files used as input for this code are available in the NCBI Gene Expression Omnibus (series accession GSE152659).
######################################################

#####################################################
#trim.sh Trims reads. Create and submit multiple jobs (our computer cluster uses slurm) by replacing FILEINFO with sample id, using the command:
for f in `cat ids.txt`; do cat trim.sh | sed -e s/FILEINFO/${f}/g > trim_${f}.sh; sbatch trim_${f}.sh; done
#####################################################
#!/bin/bash
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G # Memory per core
#SBATCH -o trim_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e trim_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent
module load TrimGalore
module load cutadapt
module load samtools

inputdir=directory_containing_fastq
cutadapt -q 20 -e 0.2 --times 5 --overlap 2 -a AGATCGGAAGAGC -a "T{100}" --minimum-length=20 -o ${inputdir}/FILEINFO.cutadapt.fq.gz ${inputdir}/FILEINFO.fastq.gz

#####################################################
#build.sh Build genome indices for STAR, only need to do once. Make sure to specify GTF file:
#####################################################
#!/bin/bash
#SBATCH -n 24 # Number of cores
#SBATCH -N 2 # Ensure that all cores are on  machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G 
#SBATCH -o genomebuild.out # File to which STDOUT will be written
#SBATCH -e genomebuild.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01
STAR --runMode genomeGenerate --limitGenomeGenerateRAM 96000000000 --genomeDir DMR_v1.0_star --genomeFastaFiles Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.fa  --runThreadN 8 --sjdbGTFfile Fukomys_damarensis.DMR_v1.0.92.gtf 

#####################################################
#map.sh Map reads, single end, first pass.
for f in `cat ids.txt`; do cat map.sh | sed -e s/FILEINFO/$f/g > map_$f.sh; sbatch map_$f.sh; done
#####################################################

#!/bin/bash
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G # Memory per cpu
#SBATCH -o map_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e map_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01

genome_path=DMR_v1.0_star
path_out=mapped
path_in=trimmed
STAR --genomeDir $genome_path/ --outFileNamePrefix $path_out/mapped_FILEINFO --runThreadN 12 --readFilesCommand zcat --readFilesIn $path_in/FILEINFO.cutadapt.fq.gz

#####################################################
#spliceindex.sh From the first pass mapping, create an index of splice locations in preparation for 2nd pass mapping.
#####################################################
mkdir DMR_v1.0_star_secondpass-allbmsc

#!/bin/bash
#SBATCH -n 24 # Number of cores
#SBATCH -N 2 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=96G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o spliceindex.out # File to which STDOUT will be written
#SBATCH -e spliceindex.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01
path_out=DMR_v1.0_star_secondpass-allbmsc
genomefastapath=directory
splicelocationsdir=mapped
STAR --runMode genomeGenerate --outFileNamePrefix $path_out --runThreadN 12 --genomeDir $path_out --limitGenomeGenerateRAM 96000000000 --genomeFastaFiles $genomefastapath/Fukomys_damarensis.DMR_v1.0.dna_sm.toplevel.fa --sjdbFileChrStartEnd $splicelocationsdir/mapped_CRF002_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_CRF002_v_null_R1SJ.out.tab $splicelocationsdir/mapped_CRF007_v_or_lb_null_A_R1SJ.out.tab $splicelocationsdir/mapped_CRF007_v_orLlb_null_B_R1SJ.out.tab $splicelocationsdir/mapped_G10F026_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_G10F026_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_G10F026_v_null_R1SJ.out.tab $splicelocationsdir/mapped_G1F022_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_G1F022_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_G1F025_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_G1F025_v_null_R1SJ.out.tab $splicelocationsdir/mapped_G4F019_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_G4F019_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_G4F019_v_null_R1SJ.out.tab $splicelocationsdir/mapped_G4F020_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_G4F020_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_G4F020_v_null_R1SJ.out.tab $splicelocationsdir/mapped_GRF001_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_GRF001_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_GRF001_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_GRF001_v_null_R1SJ.out.tab $splicelocationsdir/mapped_GRF007_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_GRF007_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_GRF007_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_GRF007_v_null_R1SJ.out.tab $splicelocationsdir/mapped_HEF001_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_HEF001_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_HEF001_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_HEF001_v_null_R1SJ.out.tab $splicelocationsdir/mapped_L2F020_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_L2F020_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_L2F020_v_null_R1SJ.out.tab $splicelocationsdir/mapped_MAF003_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_MAF003_v_null_R1SJ.out.tab $splicelocationsdir/mapped_MAF006_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_MAF006_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_MAF006_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_MAF008_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_MAF008_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_MAF008_v_null_R1SJ.out.tab $splicelocationsdir/mapped_NOF007_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_NOF007_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_NOF007_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_NOF007_v_null_R1SJ.out.tab $splicelocationsdir/mapped_NOF009_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_NOF009_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_TIF008_V_A_R1SJ.out.tab $splicelocationsdir/mapped_TIF008_V_B_R1SJ.out.tab $splicelocationsdir/mapped_WEF001_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_WEF001_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_WEF001_v_null_R1SJ.out.tab $splicelocationsdir/mapped_WEF002_lb_E2_R1SJ.out.tab $splicelocationsdir/mapped_WEF002_lb_null_R1SJ.out.tab $splicelocationsdir/mapped_WEF002_v_null_R1SJ.out.tab $splicelocationsdir/mapped_Z3F022_v_E2_R1SJ.out.tab $splicelocationsdir/mapped_Z3F022_v_null_R1SJ.out.tab $splicelocationsdir2019/mapped_ARF002_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_ARF002_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_ARF002_v_E2SJ.out.tab $splicelocationsdir2019/mapped_ARF002_v_nullSJ.out.tab $splicelocationsdir2019/mapped_BEF001_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_BEF001_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_BEF001_v_E2SJ.out.tab $splicelocationsdir2019/mapped_BEF001_v_nullSJ.out.tab $splicelocationsdir2019/mapped_DVF007_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_DVF007_v_E2SJ.out.tab $splicelocationsdir2019/mapped_DVF007_v_nullSJ.out.tab $splicelocationsdir2019/mapped_DVF008_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_DVF008_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_DVF008_v_E2SJ.out.tab $splicelocationsdir2019/mapped_DVF008_v_nullSJ.out.tab $splicelocationsdir2019/mapped_LAF007_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_LAF007_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_LAF007_vert_nullSJ.out.tab $splicelocationsdir2019/mapped_LAF008_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_LAF008_vert_nullSJ.out.tab $splicelocationsdir2019/mapped_WAF011_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_WAF011_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_WAF011_v_E2SJ.out.tab $splicelocationsdir2019/mapped_WAF011_v_nullSJ.out.tab $splicelocationsdir2019/mapped_WAF014_lb_E2SJ.out.tab $splicelocationsdir2019/mapped_WAF014_lb_nullSJ.out.tab $splicelocationsdir2019/mapped_WAF014_v_E2SJ.out.tab $splicelocationsdir2019/mapped_WAF014_v_nullSJ.out.tab

#####################################################
#map2.sh Use the splice index to perform final alignment.
for f in `cat ids.txt`; do cat map2.sh | sed -e s/FILEINFO/$f/g > map2_$f.sh; sbatch map2_$f.sh; done
#####################################################

#!/bin/bash
#SBATCH -n 12 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem-per-cpu=4G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o map2_FILEINFO.out # File to which STDOUT will be written
#SBATCH -e map2_FILEINFO.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

module load STAR/2.5.0c-fasrc01
genome_path=DMR_v1.0_star_secondpass-allbmsc
path_in=trimmed
path_out=secondpass
STAR --genomeDir $genome_path --outFileNamePrefix $path_out/mapped2_FILEINFO --readFilesCommand zcat --readFilesIn $path_in/FILEINFO.cutadapt.fq.gz --runThreadN 12

############################################################
#htseq.sh Run HTSeq to count number of reads mapping to each gene. 
#Note that for running htseq, we extended the serpine gene 2kb in each direction in the gtf file to account for many reads mapping just adjacent to the Ensembl gene coordinates.
for f in `cat ids.txt`; do cat htseq.sh | sed -e s/FILEINFO/$f/g > htseq_$f.sh; sbatch htseq_$f.sh; done
############################################################
#!/bin/bash
#SBATCH -n 2 # Number of cores
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-96:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=myemail@gmail.com # Email to which notifications will be sent

f=FILEINFO
module load samtools
module load HTSeq
directory=secondpass #directory to sam files

cd $directory
samtools view -H mapped2_${f}Aligned.out.sam > head_${f}.sam #Creates file which just contains the header
awk '($5=="255") {print $0;}' mapped2_${f}Aligned.out.sam > select_${f}.sam #Filter to only retain reads that mapped to a single genomic location.
cat head_${f}.sam select_${f}.sam > ${f}_uniquelymapped.sam #Concatenate header and filtered sam file
samtools view -bS ${f}_uniquelymapped.sam | samtools sort -n - -o ${f}_uniquesorted.bam #Convert sam to bam and sort by read name for htseq
#Run htseq. Htseq requires sam format, so use samtools view to pipe the sam file to htseq.
samtools view ${f}_uniquesorted.bam | htseq-count --mode=union --stranded=yes --idattr=gene_id - Fukomys_damarensis.DMR_v1.0.92_serpineextended2kb.gtf > ${f}_HTSeq-counts_uniq_serpine2kb.txt

#Remove intermediate files:
rm head_${f}.sam
rm select_${f}.sam
rm ${f}_uniquelymapped.sam

############################################################
#In R, merge htseq output files to create read count matrix. 
############################################################
file1=read.delim("G1F022_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file2=read.delim("G1F022_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file3=read.delim("G1F025_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file4=read.delim("G1F025_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file5=read.delim("G4F019_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file6=read.delim("G4F019_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file7=read.delim("G4F019_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file8=read.delim("G4F020_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file9=read.delim("G4F020_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file10=read.delim("G4F020_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file11=read.delim("MAF003_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file12=read.delim("MAF003_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file13=read.delim("MAF006_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file14=read.delim("MAF006_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file15=read.delim("MAF006_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file16=read.delim("NOF007_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file17=read.delim("NOF007_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file18=read.delim("NOF007_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file19=read.delim("NOF007_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file20=read.delim("NOF009_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file21=read.delim("NOF009_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file22=read.delim("WEF001_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file23=read.delim("WEF001_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file24=read.delim("WEF001_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file25=read.delim("WEF002_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file26=read.delim("WEF002_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file27=read.delim("WEF002_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file28=read.delim("DVF007_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file29=read.delim("DVF007_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file30=read.delim("DVF007_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file31=read.delim("DVF008_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file32=read.delim("DVF008_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file33=read.delim("DVF008_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file34=read.delim("DVF008_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file35=read.delim("LAF007_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file36=read.delim("LAF007_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file37=read.delim("LAF007_vert_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file38=read.delim("LAF008_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file39=read.delim("LAF008_vert_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file40=read.delim("WAF011_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file41=read.delim("WAF011_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file42=read.delim("WAF011_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file43=read.delim("WAF011_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file44=read.delim("WAF014_lb_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file45=read.delim("WAF014_lb_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file46=read.delim("WAF014_v_E2_HTSeq-counts_uniq_serpine2kb.txt",header=F)
file47=read.delim("WAF014_v_null_HTSeq-counts_uniq_serpine2kb.txt",header=F)

this=cbind(file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12, file13, file14, file15, file16, file17, file18, file19, file20, file21, file22, file23, file24, file25, file26, file27, file28, file29, file30, file31, file32, file33, file34, file35, file36, file37, file38, file39, file40, file41, file42, file43, file44, file45, file46, file47)
this2=this[,seq(2,ncol(this),2)]

colnames(this2)=c("G1F022_lb_E2", "G1F022_lb_null", "G1F025_lb_null", "G1F025_v_null", "G4F019_lb_null", "G4F019_v_E2", "G4F019_v_null", "G4F020_lb_E2", "G4F020_lb_null", "G4F020_v_null", "MAF003_lb_null", "MAF003_v_null", "MAF006_lb_E2", "MAF006_lb_null", "MAF006_v_E2", "NOF007_lb_E2", "NOF007_lb_null", "NOF007_v_E2", "NOF007_v_null", "NOF009_lb_E2", "NOF009_lb_null", "WEF001_lb_null", "WEF001_v_E2", "WEF001_v_null", "WEF002_lb_E2", "WEF002_lb_null", "WEF002_v_null", "DVF007_lb_null", "DVF007_v_E2", "DVF007_v_null", "DVF008_lb_E2", "DVF008_lb_null", "DVF008_v_E2", "DVF008_v_null", "LAF007_lb_E2", "LAF007_lb_null", "LAF007_vert_null", "LAF008_lb_null", "LAF008_vert_null", "WAF011_lb_E2", "WAF011_lb_null", "WAF011_v_E2", "WAF011_v_null", "WAF014_lb_E2", "WAF014_lb_null", "WAF014_v_E2", "WAF014_v_null")
rownames(this2)<-file1[,1]
#Remove bottom 5 lines, which are summary statistics:
n<-nrow(this2)-4
this2<-this2[-c(n:nrow(this2)),]
write.table(this2,file="molerat_htseq_47bmscs.txt",sep="\t")

############################################################
#Perform gene-by-gene linear modeling 
############################################################
library(coxme)
library(edgeR)
library(limma)
library(spaa)
library(MASS)
library(ggplot2)
library(EMMREML)
library(stats)
library(cobs)
library(matrixcalc)
library(dplyr)

#Prepare traits file
traits<-read.table("traits_47bmscs.txt",sep="\t",header=T)
traits$cellcollectiondate <- as.Date(traits$cellcollectiondate,"%m/%d/%y")
traits$ageatsac_years<-as.numeric(traits$ageatsac_years)
traits$cellstrypsinized<-as.factor(traits$cellstrypsinized)
traits$bone<-as.factor(traits$bone)
traits$animal<-as.factor(traits$animal)
cols<-traits
cols$queen<-as.factor(cols$queen)
cols$bone<-as.factor(cols$bone)
cols$treatment<-as.factor(cols$treatment)
cols$treatment <- relevel(cols$treatment, ref = "null")

#Prepare counts file
rawcounts=read.table("molerat_htseq_47bmscs.txt",sep="\t",header=T)
rawcounts<-rawcounts[,as.character(traits$sample)]

#Calculate TPM, using each gene's average length (averaged across all transcripts annotated for that gene in the Ensembl gtf file)
gene_length<-read.table("Fukomys_damarensis.DMR_v1.0.92.gtf.metrics_transcriptLengths_meanGeneLength.txt",header=T)
gene_length<-gene_length[,c("gene_id","mean_gene_length")]
counts2<-merge(rawcounts,gene_length,by.x="row.names",by.y="gene_id")
rownames(counts2)<-counts2[,1];counts2<-counts2[,-c(1)]
total<-colSums(counts2[,1:nrow(traits)])

tpm<-counts2[,1:nrow(traits)]*100*10^6/counts2$mean_gene_length
#Set final row equal to T
tpm[nrow(tpm)+1,]<-apply(X = counts2[,1:nrow(traits)],2,function(x) sum(x*100/counts2$mean_gene_length))
for (i in 1:nrow(traits)){
  tpm[,i]<-tpm[,i]/tpm[nrow(tpm),i]
}
tpm<-tpm[-nrow(tpm),]

#Require 25% samples to have TPM>=2
counts2<-counts2[,-ncol(counts2)] #remove last column, which is mean_gene_length
rawcounts2<-counts2[which(apply(tpm,1,function(a){quantile(a,0.75)})>=2),] 

########################################
#Voom normalize read count data and prepare files for gene-by-gene linear analysis
########################################
#Voom normalize, removing litter effect and proportion reads mapping within feature
regressiondesign = model.matrix(~propwithinfeat_uniqmapped+natal_colony,data=traits)
dge <- DGEList(counts=rawcounts2)
dge <- calcNormFactors(dge)
v <- voom(dge,regressiondesign,plot=F)
fit <-lmFit(v,regressiondesign)
fit2 <- eBayes(fit)
exp<-residuals.MArrayLM(object=fit2, v)
#Save file to run mediation analyses later
write.table(exp,file='residuals_propinfeat_litter.txt',sep='\t')

#Read in relatedness matrix
K<-read.table("relatedness.txt",sep="\t",header=T)
K<-as.matrix(K[levels(traits$animal),levels(traits$animal)])
dim(K)

#Create Z, which is an incidence matrix mapping samples to individuals
random_effects=exp[,1:nrow(K)]
colnames(random_effects)=colnames(K)

Z=matrix(rep(0,nrow(cols)*ncol(K)),nrow=nrow(cols),ncol=ncol(K)) 
rownames(Z)=rownames(cols)
colnames(Z)=colnames(K)
for(i in 1:ncol(Z))
{
  set=which(cols$animal == colnames(Z)[i])
  Z[set,i]=1
}

#########################################################
#Run emmreml actual data:
#########################################################
#Prepare dataframe for output results
design = model.matrix(~bone+bone:queen+bone:treatment+daysinculture,data=cols)
res_full<-as.data.frame(matrix(nrow=nrow(exp),ncol=3*ncol(design))) #create new matrix
colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("se_",colnames(design))
colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_val_",colnames(design))
res_full<-data.frame(res_full)
rownames(res_full)<-rownames(exp)

#Run mixed model, gene by gene:
for(i in 1:nrow(exp))
{
  emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  random_effects[i,]=t(emma$uhat)
  res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
}

write.table(res_full,file="model_bone_bonequ_bonetreat_dayscult_EMMREML_regrlitterpropinfeatvoom.txt",sep="\t")

#########################################################
#To account for multiple hypothesis testing, we calculate here the empirically derived false discovery rate. 
#To do so, we first create a dataframe of all possible permutations of breeding status, while retaining in the permutations 
#the fact that only one animal in each pair of siblings can be a breeder:
#########################################################
library(combinat)
library(gtools)

#In the actual data, there were 8 animal pairs total, 5 pairs of which were breeder/nonbreeder and 3 pairs which were nonbreeder/nonbreeder. 
#We start with randomly assigning which 5 of the 8 pairs will contain 1 queen
x<-as.character(c(0,0,0,1,1,1,1,1)) #3 pairs without queen, 5 pairs with queen
this<-permn(x)
this2<-unique(this) #56 possible ways to assign which 5 of 8 pairs that contain a queen

#For each set of 5 pairs that contain a queen, we are essentially doing coinflips to determine if the first female is a queen or nonbreeder.
#There are 32 combinations of who are the 5 queens. This is specified in 'coinflips'.
xx=c("0","1") #possible values for coin flip
coinflips<-permutations(n=2,r=5,v=xx,repeats.allowed=T) #n=number of values to choose from (eg 2 for heads/tails); r=number of times you flip the coin)

#Declare flag variable (1 at the first entry per pair, 0 in the rest). We do this so that we can randomly assign queen statuses to the first *sample* of each animal. 
#Samples matching the same animal will then be filled in so that samples from same animal all have the same breeding status.
cols$flag=0
cols$natal_colony<-as.factor(cols$natal_colony)
for(i in 1:length(levels(cols$natal_colony)))
{
  cols$flag[which(cols$natal_colony == levels(cols$natal_colony)[i])[1]]=1
}

#Now perform permutations of breeding status;
cols_random<-cols
for (iters in 1:56) {
  cols_random$natal_colony2=0 #Create natal_colony2, which will assign, in the next line of code, 5 pairs to contain a queen
  cols_random$natal_colony2[which(cols_random$flag==1)]=this2[[iters]] 
  
  #For each pair containing a queen (specified in natal_colony2), if first animal is 0, then second animal is 1. Else second is 0. Create a separate column for this, repeating for all possible combinations of coinflips.
  for (f in 1:32) {
    #Assign the first animal of each of the 5 pairs according to 1 of 32 possibilities
    set0=which(cols_random$natal_colony2==1)
    colname<-paste0("perm",iters,"_",f,sep="")
    cols_random[,colname]=0
    cols_random[,colname][set0]=coinflips[f,]
    #Complete the assignment to the rest of the samples
    for(i in 1:length(levels(cols_random$natal_colony)))
    {
      set=which(cols_random$natal_colony==levels(cols_random$natal_colony)[i] & cols_random$firstinpair==0)
      set_ref=which(cols_random$natal_colony==levels(cols_random$natal_colony)[i] & cols_random$flag==1)
      if(cols_random[,colname][set_ref]==0 & cols_random$natal_colony2[set_ref]==1) {
        cols_random[,colname][set]=1
      } 
      set2=which(cols_random$natal_colony==levels(cols_random$natal_colony)[i] & cols_random$firstinpair==1)
      set_ref2=which(cols_random$natal_colony==levels(cols_random$natal_colony)[i] & cols_random$flag==1)
      cols_random[,colname][set2]=cols_random[,colname][set_ref2]
      
    }	#close i
  }	#close f
} #close iters		

cols_random2<-cols_random[,62:1853] #dim=41 rows x560 columns
cols_random3<-cols_random2
for (f in 1:1792) {
  cols_random3[,f]<-as.factor(cols_random3[,f])
}
cols_random2<-cols_random3
head(cols_random2[,c(1:6)]) #data frame in which columns represent all 1,792 possible permutations of breeding status
write.table(cols_random2,file='cols_random2.txt',sep='\t',header=T)
#Randomly sample 100 integers between 1 and 1792. This specifies which permutations from the permutations dataframe (cols_random2) we'll run.
randomsequence<-sample(1:1792,100,replace=F) 
#Save the permutations for posterity:
write.table(randomsequence,file='breeder_status_permutations100.txt',sep='\t')
#randomsequence vector from paper:
#randomsequence=c("1004","1009","1014","1016","1025","1082","1088","1095","1101","1104","1124","112","1143","1160","1175","1222","1226","1233","1237","1299","1321","1336","1378","1395","1407","142","1464","1487","1495","1517","1519","153","155","1561","1566","1575","158","15","1616","1621","1631","1658","1659","1685","1715","1721","1764","1766","1767","1790","230","253","255","260","261","264","265","269","317","319","339","342","385","430","433","447","449","464","501","504","505","506","508","53","565","582","585","605","651","656","659","65","697","6","732","757","761","785","794","798","834","851","873","893","915","936","954","957","962","970")

#########################################################
#Now run some permutations!
#########################################################
cols_randoms2<-read.table('cols_random2.txt',sep='\t',header=T)
randomsequence<-read.table('breeder_status_permutations100.txt',sep='\t',header=T)

for (j in 1:100) {
  p<-randomsequence[j] 
  cols$queen<-cols_random2[,p] #permute queen
  design = model.matrix(~bone+bone:queen+bone:treatment+daysinculture,data=cols)
  res_full<-as.data.frame(matrix(nrow=nrow(exp),ncol=3*ncol(design))) #create new matrix
  colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
  colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("se_",colnames(design))
  colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_val_",colnames(design))
  res_full<-data.frame(res_full)
  
  for(i in 1:nrow(exp))
  {
    emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
    this<-paste(p,"Row",i,sep=" ")
    print(this)
  } #close i
  write.table(res_full,file=paste0("model_bone_bonequ_queentreat_dayscult_EMMREML_regrlitterpropinfeatvoom_permbs",j,".txt",sep=""),sep="\t")
  print(j)
} #close j

############################################################
## Create dataframe *shuffled_queen_pvals_qinlb* with pvals from permutations (one column=1 permutation)
############################################################
for (j in 1:100) {
  perm<-paste0("model_bone_bonequ_queentreat_dayscult_EMMREML_regrlitterpropinfeatvoom_permbs",j,".txt",sep="")
  permfile<-read.table(perm,header=T,sep="\t")
  
  if(p=="1")
  {
    shuffled_queen_pvals_qinlb <-data.frame(x=permfile[,"p_val_bone0.queen1"])
    rownames(shuffled_queen_pvals_qinlb)=rownames(permfile)
  } else {
    shuffled_queen_pvals_qinlb <- cbind(shuffled_queen_pvals_qinlb,x=permfile[,"p_val_bone0.queen1"])
  }
}

############################################################
## Create dataframe *shuffled_queen_pvals_qinvert* with pvals from permutations (one column=1 permutation)
############################################################
for (j in 1:100) {
  perm<-paste0("model_bone_bonequ_queentreat_dayscult_EMMREML_regrlitterpropinfeatvoom_permbs",j,".txt",sep="")
  permfile<-read.table(perm,header=T,sep="\t")
  
  if(j=="1")
  {
    shuffled_queen_pvals_qinvert <-data.frame(x=permfile[,"p_val_bone1.queen1"])
    rownames(shuffled_queen_pvals_qinvert)=rownames(permfile)
  } else {
    shuffled_queen_pvals_qinvert <- cbind(shuffled_queen_pvals_qinvert,x=permfile[,"p_val_bone1.queen1"])
  }
}

############################################################
#Correct for Multiple testing and filter results to report
#perm.fdr is the function we use to calculate the empirical FDR by comparing the observed p-values to the permuted p-values.
############################################################
library(cobs)
perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}

#perm.fdr=function(input_df,permuted_df,input_column_name_observed_pvalues,output_column_name_prefix)

#Read in observed data
real<-read.table("model_bone_bonequ_bonetreat_dayscult_EMMREML_regrlitterpropinfeatvoom.txt",header=T,sep="\t")

#Apply perm.fdr to calculate FDR for effect of breeder status in the long bones:
res_full=perm.fdr(data.frame(real),shuffled_queen_pvals_qinlb,"p_val_bone0.queen1","Queen_at_lb")
#Order the dataframe by rownames (important if you're going to run perm.fdr on another column)
res_full=res_full[order(rownames(res_full)),]

#Apply perm.fdr to calculate FDR for effect of breeder status in the lumber vertebrae:
res_full=perm.fdr(data.frame(res_full),shuffled_queen_pvals_qinvert,"p_val_bone1.queen1","Queen_at_vert")
res_full=res_full[order(rownames(res_full)),]

write.table(res_full,file="model_bone_bonequ_bonetreat_dayscult_permpc1_EMMREML_regrlitterpropfeat_voom_FDR.txt",sep="\t")
