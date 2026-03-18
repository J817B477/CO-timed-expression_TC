import os
import re
import pandas as pd
from io import StringIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time


# Snakemake variables
sampleQC_list = snakemake.input.qc_txts


overrep_df = pd.DataFrame()
overrep_map = {}
fa_sequences = []


# creates sample level fastqc dict
for sample in sampleQC_list:
  # modules = {}

  sample_name = re.split("/",sample)[1].split(".")[0]
  # sample_names.append(sample_name)
  with open(sample) as f:
    fastqc_content = f.read()
  

  #---- generates data frame of overrepresented sequences ----#
  # for sample_name in sample_names:

  overrep_module = fastqc_content.split(">>Overrepresented sequences")[1].split(">>END_MODULE")[0].strip()
  clean_module = overrep_module.splitlines()[1:]
  overrep_module = "\n".join(clean_module)
  
  if len(overrep_module) > 0: 
    
    temp_df = pd.read_csv(StringIO(overrep_module), 
                          sep='\t')
    temp_df["sample"] = sample_name
    
    overrep_df = pd.concat([overrep_df,temp_df], axis = 0, ignore_index=True)

#---- defines overep fastA strings ----#
if not overrep_df.empty and 'Possible Source' in overrep_df.columns:
  if not os.path.isfile("overrepresented_sequences.csv"):
    overrep_seqs = overrep_df[overrep_df['Possible Source'] == 'No Hit']['#Sequence'].unique()

    for seq in overrep_seqs:
      result_handle = NCBIWWW.qblast("blastn",
                                    "nt", 
                                    seq, 
                                    entrez_query="txid9606[ORGN]")

      blast_record = NCBIXML.read(result_handle)

      if blast_record.alignments:
        human_source = True if "Homo sapien" in blast_record.alignments[0].hit_def else False
      else:
        human_source = None

      overrep_map[seq] = human_source

      time.sleep(2)

    overrep_df["human_source"] = overrep_df['#Sequence'].map(overrep_map)

    overrep_df.to_csv("overrepresented_sequences.csv",
                      index= False)

  else:
    overrep_df = pd.read_csv("overrepresented_sequences.csv")

else:
    # If there's no data, ensure we have a column to look for later to avoid crashes
    overrep_df["human_source"] = None


if "human_source" in overrep_df.columns:
    fa_sequences = overrep_df[overrep_df["human_source"] == False]['#Sequence'].unique()
else:
    fa_sequences = []

with open(snakemake.output.custom_fa, "w") as f:
  if len(fa_sequences) > 0:
    for i, seq in enumerate(fa_sequences, start=1):
        f.write(f">contaminant_{i}\n{seq}\n")

  else:
    # writing a comment prevents some parsers from complaining.
    f.write("# No non-human overrepresented sequences found\n")