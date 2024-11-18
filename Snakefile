import os
import sys
import pandas as pd  
import numpy as np   

rule all:
    input:
        "Wrapper_output.txt"  # Final expected output file, adjust based on actual results

rule main_sv:
    output:
        "main_sv_output.txt"  # Placeholder for actual output, modify as needed
    shell:
        """
        cd PIPE
        module load python/3.9.2_torch_gpu
        python /home/sghaffari/PIPE/MAIN_SV.py > {output}
        """

rule r_sv:
    input:
        "main_sv_output.txt"
    output:
        "r_sv_output.txt"  # Placeholder for actual output, modify as needed
    shell:
        """
        module load R/4.1.0
        chmod +x ./R_SV.r
        Rscript ./R_SV.r > {output}
        """

rule wrapper:
    input:
        "r_sv_output.txt"
    output:
        "Wrapper_output.txt"  # Placeholder for actual output, modify as needed
    shell:
        """
        python /home/sghaffa