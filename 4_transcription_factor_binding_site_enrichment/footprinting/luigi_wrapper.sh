#!/bin/sh

# PYTHONPATH to location of bioluigi and tf_analysis.py
PYTHONPATH="/project2/lbarreiro/programs/TF_Activity_Analysis/"
module load python
module load bedtools
/project2/lbarreiro/programs/luigi/bin/luigi --module luigi_pipeline.tf_analysis --local-scheduler $@
