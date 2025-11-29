# GenVarNonCode - Quick Start Guide

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Siavashghaffari/GenVarNonCode
   cd GenVarNonCode
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Verify your environment:**
   ```bash
   python check_environment.py
   ```

4. **Validate the Snakefile configuration:**
   ```bash
   python validate_snakefile.py
   ```

## Running the Pipeline

### Option 1: Using Snakemake (Recommended)

**Run both CNV and SV analysis:**
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

**Dry run (see what will be executed without running):**
```bash
snakemake -n
```

**View workflow diagram:**
```bash
snakemake --dag | dot -Tpng > workflow_dag.png
```

### Option 2: Using Shell Scripts (for HPC with PBS)

**For CNV analysis:**
```bash
bash run_CNV.sh
```

**For SV analysis:**
```bash
bash run_SV.sh
```

## Configuration

Edit `config.yaml` to customize:
- `cores`: Number of CPU cores to use (default: 8)
- `memory`: Memory allocation (default: "64g")
- `python`: Python executable path (default: "python")
- `rscript`: Rscript executable path (default: "Rscript")
- `output_dir`: Output directory location (default: "results")

## Output

Results are saved in the `results/` directory:
```
results/
├── CNV/
│   ├── p-values_lncRNA_CNVs.tsv    # Final CNV results
│   └── logs/                        # CNV processing logs
└── SV/
    ├── p-values_lncRNA_SVs.tsv     # Final SV results
    └── logs/                        # SV processing logs
```

## Troubleshooting

**If Snakemake is not installed:**
```bash
pip install snakemake
```

**If R is not installed:**
Visit: https://www.r-project.org/

**To clean up results:**
```bash
snakemake clean
```
or
```bash
rm -rf results/
```

## Workflow Steps

### CNV Analysis:
1. Run Python analysis (MAIN_CNV.py)
2. Run R statistical tests (R_test.r)
3. Generate final results (Wrapper.py with "lncRNA" "CNV")

### SV Analysis:
1. Run Python analysis (MAIN_SV.py)
2. Run R statistical tests (R_SV.r)
3. Generate final results (Wrapper.py with "lncRNA" "SV")

## Support

For questions or issues, please contact Siavash Ghaffari or open an issue on GitHub.
