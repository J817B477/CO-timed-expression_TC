"""
Docstring for scripts.infer_strand:
Finds consensus for strand direction based on output files of infer_experiment.py. Value in output file is used for strand direction flag in FeatureCounts calls.
"""

import re

forward_count = 0
reverse_count = 0

for file in snakemake.input:
    with open(file) as f:
        text = f.read()
        
        # gets proportion of reads found to be forward for .bam
        forward = float(re.search(r'"\+\+,--":\s+([0-9\.]+)', text).group().split(":")[1])
        # gets proportion of reads found to be reversed for .bam
        reverse = float(re.search(r'"\+-,-\+":\s+([0-9\.]+)', text).group().split(":")[1])

    # scores if .bam direction consensus is forward or reverse
    if reverse + forward > 0.66:
        if forward > reverse:
            forward_count += 1
        else:
            reverse_count += 1

# establishes direction consensus for all bam files
if forward_count > reverse_count:
    fc_strand = 1
elif forward_count < reverse_count:
    fc_strand = 2
else:
    fc_strand = 0

with open(snakemake.output[0], "w") as out:
    out.write(f"{str(fc_strand)}\n")
    
