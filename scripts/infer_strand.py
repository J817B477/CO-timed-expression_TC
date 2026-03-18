import re

forward = 0
reverse = 0

for file in snakemake.input:
    with open(file) as f:
        text = f.read()
        if "1++,1--" in text:
            forward += 1
        if "1+-,1-+" in text:
            reverse += 1

if reverse > forward:
    strand = 2
elif forward > reverse:
    strand = 1
else:
    strand = 0

with open(snakemake.output[0], "w") as out:
    out.write(str(strand))