# Genetic Impact Pipeline for Non-coding Regions


## Overview
**GenVarNonCode** is a pipeline designed to evaluate the impact of genetic variations, specifically **Copy Number Variations (CNV)** and **Structural Variants (SV)**, on non-coding regions of the genome. This tool helps identify the effects of these variations on gene regulation, expression, and other biological functions that are often overlooked when focusing on coding regions alone.

## Features
- **Copy Number Variation (CNV) Analysis**: Detects and analyzes CNVs in non-coding regions, assessing potential impacts on gene regulation and expression.
- **Structural Variant (SV) Analysis**: Identifies and evaluates SVs, including insertions, deletions, and rearrangements, in non-coding areas.

## Requirements
- **Python 3.9+**
- **Dependencies**: (e.g., `pandas`, `numpy`, `scipy`)
- **Other Tools**: R 4.1+ (for statistical analysis and visualization)

## Installation
Clone the repository and install the required dependencies:
```bash
git clone https://github.com/Siavashghaffari/GenVarNonCode
cd GenVarNonCode
pip install -r requirements.txt


## Usage

To run the pipeline, you can execute either:

- `run_CNV.sh` or `run_SV.sh` shell scripts for CNV and SV analysis, respectively.
- Alternatively, use the `Snakefile` for parallelized execution.

### Running with Snakemake
Run the pipeline with Snakemake:
```bash
snakemake -j 8

Adjust -j based on the number of available CPUs. Ensure to set paths correctly for your environment.

## Authors and acknowledgment

This work was authored by **Siavash Ghaffari**. For any queries or further information, please feel free to contact Siavash directly. Feedback is always appreciated to help improve and refine the pipeline.
