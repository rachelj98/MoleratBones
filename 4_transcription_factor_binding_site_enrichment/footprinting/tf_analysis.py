'''
Author : Saideep Gona

This PYTHON3 script is intended for running TF enrichment and
activity analysis on paired condition samples. 

Motif_Enrichment_Matching.png shows the whole workflow

PREREQS (Midway Cluster):

module load python
conda env create <python environment>
conda activate <python environment>
pip install luigi
OR
conda install -c conda-forge luigi
pip install RGT
OR
conda install -c bioconda rgt


STEPS

1.) Make sure a python environment is created as above, and that the prerequisites are
installed

2.) Create a "metadata" file for the ATACseq alignments. 

    This should be a whitespace delimited table file with the following columns:
    individual_id	condition	bam

    Individual-condition combinations should be unique. If not, please merge files
    ahead of time.

    Where the individual_id is the source of the sample, the condition describes the 
    experimental status of the sample, and bam is the path to the alignment file.

    For activity analysis, the comparison will be done per individual. This means 
    that you should have samples for each condition paired by individual to run the
    analysis.

    templates are provided : example_metadata.(xlsx,txt)

3.) Create a "luigi.cfg" config file (template is provided, example_luigi.cfg) and fill in the
    fields. Place it in a directory of interest (preferablly the output directory)

    This pipeline uses the regulatory genomics toolbox (http://www.regulatory-genomics.org/rgt/basic-introduction/).
    Some config parameters (GENOME,FILTER_PARAMS,MOTIF_DB) are related to RGT.

    REGION_FILES should be a directory of unique bed files containing regions of 
    interest to be analyzed individually.

    NOTE: The config file must specifically be called "luigi.cfg"

5.) Copy the luigi_wrapper.sh file and place it next to the "luigi.cfg" file

6.) To run (with rough runtime estimates):

(long)    bash luigi_wrapper.sh Condition_Footprints --condition1 *Condition1* --condition2 *Condition2*

(short)   bash luigi_wrapper.sh Meta_Footprint --condition1 *Condition1* --condition2 *Condition2*

(short)   bash luigi_wrapper.sh Meta_Motif_Match --condition1 *Condition1* --condition2 *Condition2*
    
    bash luigi_wrapper.sh Enrichment --condition1 *Condition1* --condition2 *Condition2*

        --- AND/OR ---
    
(long)  bash luigi_wrapper.sh Activity --condition1 *Condition1* --condition2 *Condition2*

(short) bash luigi_wrapper.sh Collect_Activity --condition1 *Condition1* --condition2 *Condition2*

    NOTE: Condition1 and Condition2 are conditions from the metadata file to be compared
    in a particular run sequencesqueu

'''

import sys, os, argparse
from glob import glob
from datetime import datetime
from os.path import join

import linecache

import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi.util import requires
# from bioluigi.scheduled_external_program import ScheduledExternalProgramTask
import shlex

import threading
import subprocess
import time

import scipy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Boilerplate
cwd = os.getcwd()           # Store current working directory
print("Current Working Directory: ", cwd)
dtime = str(datetime.now())
print("Time of Running ", dtime)
print("Python Interpreter: ", sys.executable)

class CustomWrapperTask(luigi.WrapperTask):
    CACHED_REQUIRES = []

    def cached_requires(self):
        return self.CACHED_REQUIRES or self.requires()

    def complete(self):
        return all(r.complete() for r in self.cached_requires())

    def deps(self):
        return self.cached_requires()

#####################################################################################
# CONFIG AND PREP INPUT FILES
#####################################################################################

class TF_ANALYSIS(luigi.Config):

    OUTPUT_DIR = luigi.Parameter()

    ALIGNMENT_METADATA = luigi.Parameter()

    REGION_FILES = luigi.Parameter()

    # RGT-Specific parameters
    # rgt-motifanalysis matching --help
    # NA will result in default
    GENOME = luigi.Parameter()
    FILTER_PARAMS = luigi.Parameter()
    MOTIF_DB = luigi.Parameter()

    PYTHON_ENVIRONMENT = luigi.Parameter()

    sbatch_templates = luigi.Parameter()

cfg = TF_ANALYSIS()

os.system("mkdir -p " + cfg.OUTPUT_DIR + "/logs")

class ProduceBam(luigi.Task):
    """
    Produce the bam files that relate to a sample.
    """

    sample_id = luigi.Parameter()

    bams = []
    def output(self):

        with open(cfg.ALIGNMENT_METADATA, "r") as am:
            c = 0
            for line in am:
                if c == 0:
                    c += 1
                    continue
                p_line = line.strip().split()
                if p_line[0] == self.sample_id:
                    f = p_line[2]
                    break

        return luigi.LocalTarget(f)

#####################################################################################
# Footprinting
#####################################################################################

class MergeRegions(luigi.Task):
    """
    Merge region files to create region superset
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):
        self.output().makedirs()

        os.system("module load bedtools")
        regions = glob(cfg.REGION_FILES+"/*")

        if len(regions) == 1:
            os.system("cp "+regions[0]+" "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","merged_regions.bed"))
        else:
            os.system("cat "+" ".join(regions)+" > "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all.bed"))
            os.system("sort -k1,1 -k2,2n "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all.bed")+" > "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all.bed"))
            os.system("bedtools merge -i "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all.bed") + " > " + join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","merged_regions.bed"))

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","merged_regions.bed"))


@requires(MergeRegions)
class Condition_Footprints(luigi.Task):
    """
    Merge bam files for each condition and then call footprints
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):

        if self.condition1 == self.condition2:
            print("Condition1 and Condition2 are the same")
            sys.exit()

        with open(cfg.sbatch_templates+"/basic.sbatch", "r") as b:
            script = str(b.read())
            script+=("\nconda activate "+cfg.PYTHON_ENVIRONMENT+"\n")
            script+="module load samtools\n"

        merge1 = [
            "samtools",
            "merge",
            "-f",
        ]
        merge1.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition1+"_merged.bam"))
        with open(cfg.ALIGNMENT_METADATA, "r") as am:
            c = 0
            for line in am:
                if c == 0:
                    c += 1
                    continue
                p_line = line.strip().split()
                if p_line[1] == self.condition1:
                    merge1.append(p_line[2])
        sort1 = [
            "samtools",
            "sort",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition1+"_merged.bam"),
            "-o", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition1+"_merged.bam").replace(".bam",".sorted.bam")
        ]
        footprint1 = [
            "rgt-hint",
            "footprinting",
            "--atac-seq",
            "--organism", cfg.GENOME,
            "--output-prefix",self.condition1+"_footprints",
            "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting"),
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition1+"_merged.bam").replace(".bam",".sorted.bam"),
            self.input().path
        ]

        with open(join(cfg.OUTPUT_DIR,"logs",self.condition1+"_footprints.sbatch"), "w") as first:
            first.write(script)
            first.write(" ".join(merge1) + " \n\n")
            first.write(" ".join(sort1) + " \n\n")
            first.write(" ".join(footprint1))
        os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition1+"_footprints.sbatch"))

        merge2 = [
            "samtools",
            "merge",
            "-f",
        ]
        merge2.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition2+"_merged.bam"))
        with open(cfg.ALIGNMENT_METADATA, "r") as am:
            c = 0
            for line in am:
                if c == 0:
                    c += 1
                    continue
                p_line = line.strip().split()
                if p_line[1] == self.condition2:
                    merge2.append(p_line[2])
        sort2 = [
            "samtools",
            "sort",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition2+"_merged.bam"),
            "-o", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition2+"_merged.bam").replace(".bam",".sorted.bam")
        ]
        footprint2 = [
            "rgt-hint",
            "footprinting",
            "--atac-seq",
            "--organism", cfg.GENOME,
            "--output-prefix",self.condition2+"_footprints",
            "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting"),
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition2+"_merged.bam").replace(".bam",".sorted.bam"),
            self.input().path
        ]

        with open(join(cfg.OUTPUT_DIR,"logs",self.condition2+"_footprints.sbatch"), "w") as second:
            second.write(script)
            second.write(" ".join(merge2) + " \n\n")
            second.write(" ".join(sort2) + " \n\n")
            second.write(" ".join(footprint2))
        os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition2+"_footprints.sbatch"))

class Meta_Footprint(luigi.Task):
    """
    Merge footprints into meta footprint file
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):

        with open(cfg.sbatch_templates+"/motif_match.sbatch", "r") as b:
            script = str(b.read())
            script+=("\nconda activate "+cfg.PYTHON_ENVIRONMENT+"\n")
            script+="module load samtools\n"
            script+="module load bedtools\n"

        cat = [
            "cat",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition1+"_footprints.bed"),
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",self.condition2+"_footprints.bed"),
            ">",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed")
        ]
        sort = [
            "sort -k1,1 -k2,2n",
            "-o", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed"),
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed")
        ]
        remove_trailing_white = [
            r"sed 's/^[ \t]*//;s/[ \t]*$//'",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed"),
            ">",
            join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed.clean")
        ]
        merge = [
            "bedtools",
            "merge",
            "-i", 
        ]

        merge.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","all_footprints.bed.clean"))
        merge.append(">")
        merge.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","meta_footprints.bed"))

        os.system("mkdir "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match"))

        match = [
            "rgt-motifanalysis matching",
            "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match"),
            "--organism", cfg.GENOME
        ]
        if cfg.FILTER_PARAMS == "NA":
            pass
        else:
            match.append("--filter")
            match.append(cfg.FILTER_PARAMS)
        match.append("--input-files")
        match.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","meta_footprints.bed"))

        with open(join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_meta_footprint.sbatch"), "w") as out:
            
            out.write(script)
            out.write("\n")
            out.write(" ".join(cat) + " \n\n")
            out.write(" ".join(sort) + " \n\n")
            out.write(" ".join(remove_trailing_white) + " \n\n")
            out.write(" ".join(merge) + " \n\n")

        os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_meta_footprint.sbatch"))


#####################################################################################
# Run Enrichment Analysis
#####################################################################################

class Enrichment(luigi.Task):
    """
    Run Enrichments on a set of conditions
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):

        regions = glob.glob(cfg.REGION_FILES+"/*")
        for r in regions:
            r_tag = regions.split("/")[-1].split(".")[0]


            # CONDITION 1
            with open(cfg.sbatch_templates+"/enrichment.sbatch", "r") as b:
                script = str(b.read())
                script+=("\nconda activate "+cfg.PYTHON_ENVIRONMENT+"\n")
                script+="module load samtools\n"
                script+="module load bedtools\n"
                script+="module load rstudio\n"

            os.system("mkdir -p "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,"motif_match"))

            # Intersect condition-specific footprints with region file
            intersect_footprints2_1 = [
                "bedtools",
                "intersect",
                "-a", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting", self.condition1+"FOOTPRINTS"),
                "-b", r,
                "-wa", ">",
                join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,rtag,"motif_match",self.condition1+"_regional_footprints.bed")
            ]

            # Motif Match on the intersected bed file
            match = [
                "rgt-motifanalysis matching",
                "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,rtag,"motif_match"),
                "--organism", cfg.GENOME
            ]
            if cfg.FILTER_PARAMS == "NA":
                pass
            else:
                match.append("--filter")
                match.append(cfg.FILTER_PARAMS)
            match.append("--input-files")
            match.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","meta_footprints.bed"))

            # Run Enrichment
            os.system("mkdir -p "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,"enrichment"))

            enrichment1 = [
                "rgt-motifanalysis enrichment",
                "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,"enrichment"),
                "--organism", cfg.GENOME,
            ]
            if cfg.FILTER_PARAMS == "NA":
                pass
            else:
                enrichment1.append("--filter")
                enrichment1.append(cfg.FILTER_PARAMS)
            enrichment1.append([
                "--matching-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,rtag,"motif_match"),
                join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match")
            ])


# rgt-motifanalysis enrichment --output-location /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/open \
# --organism hg19 --filter "species:sapiens" \
# --matching-location /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/open \
# footprint_background/All_merged.bed \
# /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/open/open_footprints.bed



            with open(join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_motif_match.sbatch"), "w") as out:
                
                out.write(script)
                out.write(" ".join(cat) + " \n\n")
                out.write(" ".join(sort) + " \n\n")
                out.write(" ".join(merge) + " \n\n")

            os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition1+"_footprints.sbatch"))




# rgt-motifanalysis matching --output-location /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/change \
# --organism hg19 --rand-proportion 10 --filter "species:sapiens" \
# --input-files /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/change/change_footprints.bed


#####################################################################################
# Run Meta Motif Matching
#####################################################################################

class Meta_Motif_Match(luigi.Task):
    """
    Compute Meta Motif Matches
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):

        regions = glob(cfg.REGION_FILES+"/*")
        print(cfg.REGION_FILES)
        print(regions)

        for r in regions:
            print(r)
            r_tag = r.split("/")[-1].split(".")[0]


            # Open sbatch
            with open(cfg.sbatch_templates+"/motif_match.sbatch", "r") as b:
                script = str(b.read())
                script+=("\nconda activate "+cfg.PYTHON_ENVIRONMENT+"\n")
                script+="module load bedtools\n"


            # Intersect meta footprints with region file
            intersect_footprints = [
                "bedtools",
                "intersect",
                "-a", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting","meta_footprints.bed"),
                "-b", r,
                "-wa", ">",
                join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",r_tag+"_footprints.bed")
            ]

            os.system("mkdir -p "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match",r_tag))
            # Motif Match on the intersected bed file
            match = [
                "rgt-motifanalysis matching",
                "--output-location", join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match",r_tag),
                "--organism", cfg.GENOME
            ]
            if cfg.FILTER_PARAMS == "NA":
                pass
            else:
                match.append("--filter")
                match.append(cfg.FILTER_PARAMS)
            match.append("--input-files")
            match.append(join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"footprinting",r_tag+"_footprints.bed"))

            script += ("\n" + " ".join(intersect_footprints) + "\n")
            script += ("\n" + " ".join(match) + "\n")

            with open(join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_"+r_tag+"_meta_match.sbatch"), "w") as out:
            
                out.write(script)
                
            os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_"+r_tag+"_meta_match.sbatch"))


#####################################################################################
# Run Activity Analysis
#####################################################################################

class Activity(luigi.Task):
    """
    Use the meta-footprints and meta-motif matches to compute sample-specific
    activity scores
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):

        regions = glob(cfg.REGION_FILES+"/*")
        for r in regions:
            r_tag = r.split("/")[-1].split(".")[0]
            with open(cfg.ALIGNMENT_METADATA, "r") as am:
                c = 0
                individuals = {}
                individuals_bams = {}
                current_individuals = set()
                for line in am:
                    if c == 0:
                        c += 1
                        continue

                    p_line = line.strip().split()
                    if p_line[0] not in individuals:
                        individuals[p_line[0]] = [False,False]
                        individuals_bams[p_line[0]] = [None,None]
                    if p_line[1] == self.condition1:
                        individuals[p_line[0]][0] = True
                        individuals_bams[p_line[0]][0] = p_line[2]
                    if p_line[1] == self.condition2:
                        individuals[p_line[0]][1] = True
                        individuals_bams[p_line[0]][1] = p_line[2]
                    if individuals[p_line[0]] == [True,True]:
                        current_individuals.add(p_line[0])
            print(current_individuals)
            for indiv in current_individuals:   
                print(indiv)        

                with open(cfg.sbatch_templates+"/activity.sbatch", "r") as b:
                    script = str(b.read())
                    script+=("\nconda activate "+cfg.PYTHON_ENVIRONMENT+"\n")

                # Run Activity
                os.system("mkdir -p "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,indiv,"activity"))

                activity = [
                    "rgt-hint differential",
                    "--organism", cfg.GENOME,
                    "--bc",
                    "--nc", "28",
                    "--mpbs-files="+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match",r_tag,r_tag+"_footprints_mpbs.bed")+","+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,"motif_match",r_tag,r_tag+"_footprints_mpbs.bed"),
                    "--conditions="+self.condition1+","+self.condition2,
                    "--colors=red,blue",
                    "--output-profiles",
                    "--reads-files="+individuals_bams[indiv][0]+","+individuals_bams[indiv][1],
                    "--output-location="+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,indiv,"activity")
                ]
            # if cfg.FILTER_PARAMS == "NA":
            #     pass
            # else:
            #     match.append("--filter")
            #     match.append(cfg.FILTER_PARAMS)

                script += ("\n" + " ".join(activity) + "\n")

                with open(join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_"+r_tag+"_"+indiv+"_activity.sbatch"), "w") as out:
                
                    out.write(script)
                
                os.system("sbatch "+join(cfg.OUTPUT_DIR,"logs",self.condition1+"_"+self.condition2+"_"+r_tag+"_"+indiv+"_activity.sbatch"))

class Collect_Activity(luigi.Task):
    """
    Use the meta-footprints and meta-motif matches to compute sample-specific
    activity scores
    """

    condition1 = luigi.Parameter()
    condition2 = luigi.Parameter()

    def run(self):


        regions = glob(cfg.REGION_FILES+"/*")
        for r in regions:
            r_tag = r.split("/")[-1].split(".")[0]
            with open(cfg.ALIGNMENT_METADATA, "r") as am:
                c = 0
                individuals = {}
                individuals_bams = {}
                current_individuals = set()
                for line in am:
                    if c == 0:
                        c += 1
                        continue

                    p_line = line.strip().split()
                    if p_line[0] not in individuals:
                        individuals[p_line[0]] = [False,False]
                        individuals_bams[p_line[0]] = [None,None]
                    if p_line[1] == self.condition1:
                        individuals[p_line[0]][0] = True
                        individuals_bams[p_line[0]][0] = p_line[2]
                    if p_line[1] == self.condition2:
                        individuals[p_line[0]][1] = True
                        individuals_bams[p_line[0]][1] = p_line[2]
                    if individuals[p_line[0]] == [True,True]:
                        current_individuals.add(p_line[0])   

            os.system("mkdir "+join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,"multisample_activity"))
            out_root = join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,"multisample_activity")
            exclude = []
            inds = set()
            for ind in current_individuals:

                if ind in inds:
                    continue
                else:
                    inds.add(ind)
                outdir_cur = join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,ind,"activity")

                if not os.path.exists(join(outdir_cur,"differential_statistics.txt")):
                    print(ind, " does not have diff statistics")
                    exclude.append(ind)
            tfs = []
            inds = set()
            for ind in current_individuals:

                if ind in exclude:
                    continue
                if ind in inds:
                    continue
                else:
                    inds.add(ind)
                outdir_cur = join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,ind,"activity")

                with open(join(outdir_cur,"differential_statistics.txt"), "r") as stats:
                    c=0
                    for line in stats:
                        if c == 0:
                            c+=1
                            continue
                        p_line = line.strip().split("\t")
                        tfs.append(p_line[0])

                break

            cols = 0
            inds = set()
            for ind in current_individuals:

                if ind in exclude:
                    continue
                if ind in inds:
                    continue
                else:
                    inds.add(ind)

                cols += 1

            all_p = np.zeros((len(tfs),cols))
            all_a = np.zeros((len(tfs),cols))

            t = 0
            tf_inds = {}
            for tf in tfs:
                tf_inds[tf] = t
                t+=1

            e=0
            inds = set()
            with open(join(out_root,"activity_samples.txt"), "w") as osamp:

                for ind in current_individuals:

                    if ind in exclude:
                        continue
                    if ind in inds:
                        continue
                    else:
                        inds.add(ind)
                    osamp.write(ind + "\n")

                    outdir_cur = join(cfg.OUTPUT_DIR, self.condition1+"_"+self.condition2,r_tag,ind,"activity")

                    with open(join(outdir_cur,"differential_statistics.txt"), "r") as stats:
                        k=0
                        for line in stats:
                            if k == 0:
                                k+=1
                                continue
                            p_line = line.strip().split("\t")
                            tf = p_line[0]
                            c = tf_inds[tf]
                            all_p[c,e] = float(p_line[8])
                            all_a[c,e] = p_line[6]
                            c+=1

                    e+=1


            outs = np.zeros((len(tfs),3))
            with open(join(out_root,"tfs.txt"), "w") as otfs:
                for x in range(len(tfs)):
                    # print(scipy.stats.combine_pvalues(all_p[x,:]))
                    outs[x,1] = scipy.stats.combine_pvalues(all_p[x,:])[1]
                    outs[x,0] = np.mean(all_a[x,:])
                    outs[x,2] = np.var(all_a[x,:])
                    otfs.write(tfs[x] + "\n")

            np.savetxt(join(out_root,"multisample_act.txt"), outs)
            np.savetxt(join(out_root,"multisample_a.txt"), all_a)
            np.savetxt(join(out_root,"multisample_p.txt"), all_p)

