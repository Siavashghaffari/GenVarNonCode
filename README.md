# Genetic Impact Pipeline for Non-coding Regions

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.1+-blue.svg)](https://www.r-project.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-7.0+-green.svg)](https://snakemake.readthedocs.io/)

## Overview

**GenVarNonCode** is a computational pipeline designed to evaluate the impact of genetic variations, specifically **Copy Number Variations (CNV)** and **Structural Variants (SV)**, on non-coding regions of the genome. This tool helps identify the effects of these variations on gene regulation, expression, and other biological functions that are often overlooked when focusing on coding regions alone.

## Key Features

- **Copy Number Variation (CNV) Analysis**: Detects and analyzes CNVs in non-coding regions, assessing potential impacts on gene regulation and expression
- **Structural Variant (SV) Analysis**: Identifies and evaluates SVs, including insertions, deletions, and rearrangements, in non-coding areas
- **Snakemake Workflow**: Parallelized execution with automatic dependency management
- **Flexible Execution**: Run both analyses together or separately
- **Organized Output**: Structured results with detailed logging
- **Environment Validation**: Built-in scripts to verify dependencies and configuration

## Requirements

### Software Requirements
- **Python 3.9+**
- **R 4.1+** (for statistical analysis)
- **Snakemake 7.0+** (for workflow management)

### Python Dependencies
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- seaborn
- biopython
- rpy2
- pyyaml
- snakemake

All Python dependencies are listed in `requirements.txt` and can be installed automatically.

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/Siavashghaffari/GenVarNonCode.git
cd GenVarNonCode
```

### 2. Install Python Dependencies
```bash
pip install -r requirements.txt
```

### 3. Verify Your Environment
Run the environment check script to ensure all dependencies are installed:
```bash
python check_environment.py
```

This script will check:
- Python version (3.9+)
- Required Python packages
- R installation
- Snakemake installation

### 4. Validate Workflow Configuration
Verify that the Snakefile and configuration are set up correctly:
```bash
python validate_snakefile.py
```

## Quick Start

For a quick start guide, see [QUICKSTART.md](QUICKSTART.md).

### Run Both CNV and SV Analysis
```bash
snakemake --cores 8
```

### Run Only CNV Analysis
```bash
snakemake --cores 8 cnv
```

### Run Only SV Analysis
```bash
snakemake --cores 8 sv
```

## Configuration

The pipeline is configured through `config.yaml`. You can customize:

```yaml
# Computational resources
cores: 8              # Number of CPU cores to use
memory: "64g"         # Memory allocation

# Executables
python: "python"      # Python executable (or path to specific version)
rscript: "Rscript"    # R executable (or path to specific version)

# Output
output_dir: "results" # Directory for all output files
```

Edit `config.yaml` before running the pipeline to match your system configuration.

## Usage

### Using Snakemake (Recommended)

Snakemake provides parallelized execution, automatic dependency management, and checkpoint/restart capabilities.

#### Basic Commands

**Run the complete pipeline (CNV + SV):**
```bash
snakemake --cores 8
```

**Run only CNV analysis:**
```bash
snakemake --cores 8 cnv
```

**Run only SV analysis:**
```bash
snakemake --cores 8 sv
```

#### Advanced Options

**Dry run (see what will be executed without running):**
```bash
snakemake -n
```

**Dry run with detailed output:**
```bash
snakemake -n -p
```

**Generate workflow visualization:**
```bash
snakemake --dag | dot -Tpng > workflow_dag.png
```
*Requires Graphviz to be installed*

**Run with custom number of cores:**
```bash
snakemake --cores 16
```

**Clean up results:**
```bash
snakemake clean
```
or manually:
```bash
rm -rf results/
```

**Force re-run of specific rule:**
```bash
snakemake --cores 8 --forcerun run_cnv_python
```

**Unlock working directory (if previous run was interrupted):**
```bash
snakemake --unlock
```

### Using Shell Scripts (HPC with PBS)

For HPC clusters with PBS scheduling system, you can use the provided shell scripts:

**For CNV analysis:**
```bash
qsub run_CNV.sh
```

**For SV analysis:**
```bash
qsub run_SV.sh
```

Or run directly:
```bash
bash run_CNV.sh
bash run_SV.sh
```

## Workflow Details

### CNV Analysis Pipeline

1. **Python Analysis** (`src/MAIN_CNV.py`)
   - Processes CNV data
   - Identifies variants in non-coding regions
   - Generates intermediate results

2. **Statistical Analysis** (`src/R_test.r`)
   - Performs statistical tests
   - Calculates significance values

3. **Final Processing** (`src/Wrapper.py`)
   - Combines results
   - Generates p-values for lncRNA-CNV associations
   - Output: `results/CNV/p-values_lncRNA_CNVs.tsv`

### SV Analysis Pipeline

1. **Python Analysis** (`src/MAIN_SV.py`)
   - Processes SV data
   - Identifies variants in non-coding regions
   - Generates intermediate results

2. **Statistical Analysis** (`src/R_SV.r`)
   - Performs statistical tests
   - Calculates significance values

3. **Final Processing** (`src/Wrapper.py`)
   - Combines results
   - Generates p-values for lncRNA-SV associations
   - Output: `results/SV/p-values_lncRNA_SVs.tsv`

## Output Structure

Results are organized in the `results/` directory:

```
results/
├── CNV/
│   ├── p-values_lncRNA_CNVs.tsv    # Final CNV analysis results
│   └── logs/                        # Execution logs for CNV workflow
│       ├── cnv_python.log           # Python analysis log
│       ├── cnv_r.log                # R analysis log
│       └── cnv_wrapper.log          # Wrapper script log
└── SV/
    ├── p-values_lncRNA_SVs.tsv     # Final SV analysis results
    └── logs/                        # Execution logs for SV workflow
        ├── sv_python.log            # Python analysis log
        ├── sv_r.log                 # R analysis log
        └── sv_wrapper.log           # Wrapper script log
```

### Output Files

- **p-values_lncRNA_CNVs.tsv**: Contains p-values and statistics for CNV impacts on lncRNA
- **p-values_lncRNA_SVs.tsv**: Contains p-values and statistics for SV impacts on lncRNA
- **logs/*.log**: Detailed execution logs for debugging and verification

## Helper Scripts

### check_environment.py

Validates your environment setup:
```bash
python check_environment.py
```

Checks:
- Python version compatibility
- Required Python packages
- R installation
- Snakemake availability

### validate_snakefile.py

Validates the workflow configuration:
```bash
python validate_snakefile.py
```

Checks:
- Snakefile presence
- config.yaml validity
- Source script availability
- Directory structure

## Troubleshooting

### Common Issues

**Issue: Snakemake not found**
```bash
pip install snakemake
```

**Issue: R not found**
- Install R from https://www.r-project.org/
- Ensure R is in your PATH
- Update `config.yaml` with full path to Rscript if needed

**Issue: Permission denied for R scripts**
- The workflow automatically makes R scripts executable
- If issues persist, manually run: `chmod +x src/*.r`

**Issue: Workflow locked**
```bash
snakemake --unlock
```

**Issue: Missing dependencies**
```bash
pip install -r requirements.txt --upgrade
```

**Issue: Out of memory**
- Reduce the number of cores: `snakemake --cores 4`
- Edit `config.yaml` to reduce memory allocation
- Run CNV and SV separately instead of together

### Getting Help

1. Run environment check: `python check_environment.py`
2. Validate configuration: `python validate_snakefile.py`
3. Check logs in `results/*/logs/`
4. Use dry run to debug: `snakemake -n -p`

## Project Structure

```
GenVarNonCode/
├── Snakefile                  # Main workflow definition
├── config.yaml                # Configuration file
├── requirements.txt           # Python dependencies
├── README.md                  # This file
├── QUICKSTART.md             # Quick start guide
├── check_environment.py      # Environment validation script
├── validate_snakefile.py     # Workflow validation script
├── run_CNV.sh                # PBS script for CNV analysis
├── run_SV.sh                 # PBS script for SV analysis
└── src/                      # Source code directory
    ├── MAIN_CNV.py           # Main CNV analysis script
    ├── MAIN_SV.py            # Main SV analysis script
    ├── R_test.r              # R statistical tests for CNV
    ├── R_SV.r                # R statistical tests for SV
    ├── Wrapper.py            # Results wrapper and formatter
    └── ...                   # Additional analysis modules
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Citation

If you use GenVarNonCode in your research, please cite:

```
Ghaffari, S. (2024). GenVarNonCode: Genetic Impact Pipeline for Non-coding Regions.
GitHub repository: https://github.com/Siavashghaffari/GenVarNonCode
```

## Authors and Acknowledgment

This work was authored by **Siavash Ghaffari**. For any queries or further information, please feel free to contact Siavash directly. Feedback is always appreciated to help improve and refine the pipeline.



## Version History

- **v2.0** (2024): Added Snakemake workflow, configuration system, and validation tools
- **v1.0** (2024): Initial release with shell scripts

---

For a quick start guide, see [QUICKSTART.md](QUICKSTART.md).
