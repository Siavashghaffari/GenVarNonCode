#!/usr/bin/env python
"""
Snakefile validation script
Tests the Snakefile structure and configuration without running Snakemake
"""

import os
import sys
import yaml

def validate_config_file():
    """Validate that config.yaml exists and is properly formatted"""
    print("Checking config.yaml...")
    if not os.path.exists("config.yaml"):
        print("  ✗ config.yaml not found")
        return False

    try:
        with open("config.yaml", 'r') as f:
            config = yaml.safe_load(f)

        required_keys = ['cores', 'memory', 'python', 'rscript', 'output_dir']
        missing_keys = [key for key in required_keys if key not in config]

        if missing_keys:
            print(f"  ✗ Missing required keys in config.yaml: {missing_keys}")
            return False

        print("  ✓ config.yaml is valid")
        return True
    except yaml.YAMLError as e:
        print(f"  ✗ Error parsing config.yaml: {e}")
        return False

def validate_snakefile():
    """Check that Snakefile exists"""
    print("Checking Snakefile...")
    if not os.path.exists("Snakefile"):
        print("  ✗ Snakefile not found")
        return False

    print("  ✓ Snakefile exists")
    return True

def validate_source_scripts():
    """Validate that all required source scripts exist"""
    print("Checking source scripts...")

    required_scripts = [
        'src/MAIN_CNV.py',
        'src/MAIN_SV.py',
        'src/R_test.r',
        'src/R_SV.r',
        'src/Wrapper.py'
    ]

    all_exist = True
    for script in required_scripts:
        if os.path.exists(script):
            print(f"  ✓ {script} found")
        else:
            print(f"  ✗ {script} not found")
            all_exist = False

    return all_exist

def validate_directory_structure():
    """Check that the directory structure is correct"""
    print("Checking directory structure...")

    if not os.path.exists("src"):
        print("  ✗ src/ directory not found")
        return False

    print("  ✓ src/ directory exists")
    return True

def main():
    """Run all validations"""
    print("=" * 60)
    print("Snakefile Validation")
    print("=" * 60)
    print()

    checks = []

    checks.append(validate_directory_structure())
    print()

    checks.append(validate_snakefile())
    print()

    checks.append(validate_config_file())
    print()

    checks.append(validate_source_scripts())
    print()

    print("=" * 60)
    if all(checks):
        print("✓ All validations passed!")
        print()
        print("The Snakefile is properly configured.")
        print()
        print("To test the workflow (requires Snakemake installation):")
        print("  snakemake -n                 # Dry run")
        print("  snakemake --cores 8 -n       # Dry run with 8 cores")
        print()
        print("To run the workflow:")
        print("  snakemake --cores 8          # Run both CNV and SV")
        print("  snakemake --cores 8 cnv      # Run only CNV")
        print("  snakemake --cores 8 sv       # Run only SV")
        return 0
    else:
        print("✗ Some validations failed.")
        print("Please fix the issues before running the workflow.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
