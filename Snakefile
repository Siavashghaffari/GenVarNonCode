"""
Snakefile for Genetic Impact Pipeline for Non-coding Regions (GenVarNonCode)
Handles both CNV and SV analysis workflows
"""

import os

# Configuration
configfile: "config.yaml"

# Default values if config file doesn't exist
CORES = config.get("cores", 8)
MEMORY = config.get("memory", "64g")
PYTHON = config.get("python", "python")
RSCRIPT = config.get("rscript", "Rscript")

# Define output directories
OUTPUT_DIR = config.get("output_dir", "results")
CNV_DIR = f"{OUTPUT_DIR}/CNV"
SV_DIR = f"{OUTPUT_DIR}/SV"

# Default rule - runs both CNV and SV analysis
rule all:
    input:
        f"{CNV_DIR}/p-values_lncRNA_CNVs.tsv",
        f"{SV_DIR}/p-values_lncRNA_SVs.tsv"

# Rule to run only CNV analysis
rule cnv:
    input:
        f"{CNV_DIR}/p-values_lncRNA_CNVs.tsv"

# Rule to run only SV analysis
rule sv:
    input:
        f"{SV_DIR}/p-values_lncRNA_SVs.tsv"

# ========================================
# CNV Analysis Workflow
# ========================================

# Step 1: Run Python script for CNV analysis
rule run_cnv_python:
    input:
        script="src/MAIN_CNV.py"
    output:
        temp(f"{CNV_DIR}/cnv_python.done")
    log:
        f"{CNV_DIR}/logs/cnv_python.log"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {CNV_DIR}/logs
        {PYTHON} {input.script} > {log} 2>&1 && touch {output}
        """

# Step 2: Run R script for CNV analysis
rule run_cnv_r:
    input:
        script="src/R_test.r",
        depends=f"{CNV_DIR}/cnv_python.done"
    output:
        temp(f"{CNV_DIR}/cnv_r.done")
    log:
        f"{CNV_DIR}/logs/cnv_r.log"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {CNV_DIR}/logs
        chmod +x {input.script}
        {RSCRIPT} {input.script} > {log} 2>&1 && touch {output}
        """

# Step 3: Run Python wrapper for CNV
rule run_cnv_wrapper:
    input:
        script="src/Wrapper.py",
        depends=f"{CNV_DIR}/cnv_r.done"
    output:
        f"{CNV_DIR}/p-values_lncRNA_CNVs.tsv"
    log:
        f"{CNV_DIR}/logs/cnv_wrapper.log"
    params:
        arg1="lncRNA",
        arg2="CNV"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {CNV_DIR}/logs
        {PYTHON} {input.script} {params.arg1} {params.arg2} > {output} 2> {log}
        """

# ========================================
# SV Analysis Workflow
# ========================================

# Step 1: Run Python script for SV analysis
rule run_sv_python:
    input:
        script="src/MAIN_SV.py"
    output:
        temp(f"{SV_DIR}/sv_python.done")
    log:
        f"{SV_DIR}/logs/sv_python.log"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {SV_DIR}/logs
        {PYTHON} {input.script} > {log} 2>&1 && touch {output}
        """

# Step 2: Run R script for SV analysis
rule run_sv_r:
    input:
        script="src/R_SV.r",
        depends=f"{SV_DIR}/sv_python.done"
    output:
        temp(f"{SV_DIR}/sv_r.done")
    log:
        f"{SV_DIR}/logs/sv_r.log"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {SV_DIR}/logs
        chmod +x {input.script}
        {RSCRIPT} {input.script} > {log} 2>&1 && touch {output}
        """

# Step 3: Run Python wrapper for SV
rule run_sv_wrapper:
    input:
        script="src/Wrapper.py",
        depends=f"{SV_DIR}/sv_r.done"
    output:
        f"{SV_DIR}/p-values_lncRNA_SVs.tsv"
    log:
        f"{SV_DIR}/logs/sv_wrapper.log"
    params:
        arg1="lncRNA",
        arg2="SV"
    threads: CORES
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 64000
    shell:
        """
        mkdir -p {SV_DIR}/logs
        {PYTHON} {input.script} {params.arg1} {params.arg2} > {output} 2> {log}
        """

# ========================================
# Utility Rules
# ========================================

# Clean up temporary files
rule clean:
    shell:
        """
        rm -rf {OUTPUT_DIR}
        echo "Cleaned up output directory"
        """

# Generate workflow diagram
rule dag:
    output:
        "workflow_dag.png"
    shell:
        """
        snakemake --dag | dot -Tpng > {output}
        echo "Workflow DAG saved to {output}"
        """
