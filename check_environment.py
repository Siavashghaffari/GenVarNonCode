#!/usr/bin/env python
"""
Environment check script for GenVarNonCode
Verifies that all required dependencies are installed
"""

import sys
import subprocess

def check_python_version():
    """Check if Python version is 3.9+"""
    version = sys.version_info
    if version.major >= 3 and version.minor >= 9:
        print(f"✓ Python {version.major}.{version.minor}.{version.micro} found")
        return True
    else:
        print(f"✗ Python 3.9+ required, found {version.major}.{version.minor}.{version.micro}")
        return False

def check_python_package(package_name, import_name=None):
    """Check if a Python package is installed"""
    if import_name is None:
        import_name = package_name

    try:
        __import__(import_name)
        print(f"✓ {package_name} installed")
        return True
    except ImportError:
        print(f"✗ {package_name} not installed")
        return False

def check_r_installation():
    """Check if R is installed"""
    try:
        result = subprocess.run(['Rscript', '--version'],
                              capture_output=True,
                              text=True,
                              timeout=5)
        if result.returncode == 0:
            version = result.stderr.split('\n')[0] if result.stderr else "version unknown"
            print(f"✓ R installed ({version})")
            return True
        else:
            print("✗ R not found")
            return False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("✗ R not found or not in PATH")
        return False

def check_snakemake():
    """Check if Snakemake is installed"""
    try:
        result = subprocess.run(['snakemake', '--version'],
                              capture_output=True,
                              text=True,
                              timeout=5)
        if result.returncode == 0:
            version = result.stdout.strip()
            print(f"✓ Snakemake {version} installed")
            return True
        else:
            print("✗ Snakemake not found")
            return False
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("✗ Snakemake not found or not in PATH")
        return False

def main():
    """Run all checks"""
    print("=" * 60)
    print("GenVarNonCode Environment Check")
    print("=" * 60)
    print()

    checks = []

    print("Checking Python version...")
    checks.append(check_python_version())
    print()

    print("Checking required Python packages...")
    packages = [
        ('numpy', 'numpy'),
        ('pandas', 'pandas'),
        ('scipy', 'scipy'),
        ('scikit-learn', 'sklearn'),
        ('matplotlib', 'matplotlib'),
        ('seaborn', 'seaborn'),
        ('biopython', 'Bio'),
        ('rpy2', 'rpy2'),
        ('pyyaml', 'yaml')
    ]

    for package_name, import_name in packages:
        checks.append(check_python_package(package_name, import_name))
    print()

    print("Checking R installation...")
    checks.append(check_r_installation())
    print()

    print("Checking Snakemake...")
    checks.append(check_snakemake())
    print()

    print("=" * 60)
    if all(checks):
        print("✓ All checks passed! Your environment is ready.")
        print()
        print("To run the pipeline:")
        print("  snakemake --cores 8        # Run both CNV and SV")
        print("  snakemake --cores 8 cnv    # Run only CNV")
        print("  snakemake --cores 8 sv     # Run only SV")
        return 0
    else:
        print("✗ Some checks failed. Please install missing dependencies.")
        print()
        print("To install Python packages:")
        print("  pip install -r requirements.txt")
        print()
        print("To install R (if needed):")
        print("  Visit: https://www.r-project.org/")
        return 1

if __name__ == "__main__":
    sys.exit(main())
