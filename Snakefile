import os
import sys
import pandas as pd  
import numpy as np   

# Snakefile for running the CNV analysis pipeline

# Define the number of cores and other resources
cores = 8
walltime = "47:59:00"
memory = "64g"
vmemory = "64g"

# Rule to run Python script for CNV analysis
rule run_cnv_python:
    input:
        "src/MAIN_CNV.py"
    output:
        temp("output_cnv.tsv")  # Temporary output
    resources:
        mem=memory, vmem=vmemory
    shell:
        "module load python/3.9.2_torch_gpu && "
        "python {input} > {output}"

# Rule to run R script for CNV analysis
rule run_r_script:
    input:
        "src/R_test.r"
    output:
        temp("output_r.tsv")  # Temporary output
    resources:
        mem=memory, vmem=vmemory
    shell:
        "module load R/4.1.0 && "
        "chmod +x {input} && "
        "Rscript {input} > {output}"

# Rule to run the Python wrapper for CNV
rule run_wrapper:
    input:
        "src/Wrapper.py"
    output:
        "p-values_lncRNA_CNVs.tsv"  # Final output file
    resources:
        mem=memory, vmem=vmemory
    params:
        arg1="lncRNA", arg2="CNV"
    shell:
        "python {input} {params.arg1} {params.arg2} > {output}"

