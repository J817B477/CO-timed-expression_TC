import os
import re
import pandas as pd
from io import StringIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
from statistics import median
from ast import literal_eval
import yaml



# Snakemake variables
sampleQC_list = snakemake.input.qc_txts
perBaseQual_df = pd.read_csv(snakemake.input.multiqc_perbase_seq_qual, sep = "\t")
seq_len_df = pd.read_csv(snakemake.input.multiqc_seq_len, sep = "\t")

# script global variables:
sample_names = []
samples_dict = {}
overrep_df = pd.DataFrame()
read_type = snakemake.params.read_type
adapter_dir = os.path.join(os.environ["CONDA_PREFIX"],
                           "share",
                           "trimmomatic",
                           "adapters")
kit_fa = ''
combined_fa_str = ''
custom_fa = ''
adapterkit_call = ''
polyA_call = ''
polyG_call = ''
overrep_map = {}
overrep_fa_string = ""
headcrop_call = ""
crop_call = ""
sliding_window_call = ""
minlen_call = ""

# creates sample level fastqc dict
for sample in sampleQC_list:
  modules = {}
  sample_name = re.split("/",sample)[1].split(".")[0]
  sample_names.append(sample_name)
  with open(sample) as f:
    fastqc_content = f.read()
  
  txt_chunks = fastqc_content.split(">>END_MODULE")

  #---- creates map of fastqc modules ----#
  for chunk in txt_chunks:
    # removes any leading/trailing empty space 
    chunk = chunk.strip() 
    if not chunk:
      continue
    
    parts = chunk.split('\n',1)
    
    if len(parts) <2:
      continue
    
    module_header, chunk = parts
    module_name = module_header.replace('>>','').split('\t')[0].strip()

    modules[module_name] = chunk

  samples_dict[sample_name] = modules

print(samples_dict[sample_name].keys())


#---- generates data frame of overrepresented sequences ----#
for sample_name in sample_names:

  if "Overrepresented sequences" in (sample := samples_dict[sample_name]): 
    overrep_module = sample["Overrepresented sequences"]
    temp_df = pd.read_csv(StringIO(overrep_module), 
                          sep='\t')
    
    temp_df["sample"] = sample_name
    
    overrep_df = pd.concat([overrep_df,temp_df], axis = 0, ignore_index=True)

#---- defines call params based on adapter presence ----#
sources = overrep_df["Possible Source"]
if "Illumina Universal Adapter" in sources:
  kit_fa = f"TruSeq3-{read_type}.fa"
elif "Nextera Transposase Sequence" in sources:
  kit_fa = "NexteraPE-PE.fa"

if kit_fa:
  with open(kit_fa) as f:
    combined_fa_str = f.read()


#---- defines overep fastA strings ----#

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

  overrep_df["human_source"] = overrep_df['#Sequence'].map(overrep_map
  )

  overrep_df.to_csv("overrepresented_sequences.csv",
                    index= False)

else:
  overrep_df = pd.read_csv("overrepresented_sequences.csv")



fa_sequences = overrep_df[overrep_df["human_source"] == False]['#Sequence']

if not fa_sequences.empty:
    overrep_fa_string = "".join(
        f">seq{i}\n{seq}\n"
        for i, seq in enumerate(fa_sequences, start=1))

    combined_fa_str += overrep_fa_string

if combined_fa_str:
  custom_fa = "combined_adapters.fa"
  with open(custom_fa, "w") as f:
    f.write(combined_fa_str)


#---- defines call strings ----# 
if "PolyA" in sources:
  polyA_call = 'cutadapt -a "A{10}" -o out.fq in.fq'

if "PolyG" in sources:
  polyG_call = "fastp -i in.fq -o out.fq --trim_poly_g"
    
if custom_fa:
  adapterkit_call = f"ILLUMINACLIP:{custom_fa}:2:30:10"


#---- defines postional crop parameters ----#

headcrop_positions = []
crop_positions = []

# unwraps per base sequence content to id candidate positions
for sample_name in sample_names:
  last_position = 0
  if "Per base sequence content" in (sample := samples_dict[sample_name]):

    per_base_module = sample["Per base sequence content"]
    temp_df = pd.read_csv(StringIO(per_base_module), 
                          sep='\t')
    
    bases = ['G','A','T','C'] 

    # data for determining headcrop call
    temp_df['stable'] = ((temp_df[bases] >= 20) & (temp_df[bases] <= 30)).all(axis=1)
    headcrop_position = int(temp_df[temp_df['stable'] == True].iloc[0,0].split('-')[0])
    headcrop_positions.append(headcrop_position)

    # data for determining cropping call
    reversed_df = temp_df.iloc[::-1].reset_index(drop=True)
    crop_position = max([int(x) for x in reversed_df[reversed_df['stable'] == True].iloc[0,0].split('-')])
    crop_positions.append(crop_position)
    last_position = max([int(x) for x in reversed_df.iloc[0,0].split('-')]+ [last_position])


# assigns positional cropping parameters if positions are inside full range 
if any(x > 1 for x in set(headcrop_positions)):
  headcrop_call = f"HEADCROP:{median(headcrop_positions)}"

if any(x < last_position for x in set(crop_positions)):
  crop_call = f"CROP:{median(crop_positions)}"



#---- defines sliding window parameters ----#

window_size = 4
min_quality_threshold = 20

# Convert string tuples to actual tuples
for col in perBaseQual_df.columns[1:]:
    perBaseQual_df[col] = perBaseQual_df[col].apply(literal_eval)
    
# Extract just the mean quality
for col in perBaseQual_df.columns[1:]:
    perBaseQual_df[col] = perBaseQual_df[col].apply(lambda x: x[1])

quality_df = perBaseQual_df.iloc[:, 1:].astype(float)

rolling_mean = (
    quality_df.T
    .rolling(window=window_size)
    .mean()
    .T)

observed_min = rolling_mean.min().min()  # lowest rolling mean across samples

if observed_min < min_quality_threshold:
    # only call sliding window if quality drops below threshold
    sliding_window_call = f"SLIDINGWINDOW:{window_size}:{min_quality_threshold}"

#---- determines min length parameter ----#

for col in seq_len_df.columns[1:]:
    seq_len_df[col] = seq_len_df[col].apply(
        lambda x: literal_eval(x) if pd.notna(x) else None
    )



lengths = []
counts = []

for row in seq_len_df.iloc[:,1:].itertuples(index=False):
  for cell in row:
    if cell is not None:
      length, count = cell
      lengths.append(length)
      counts.append(count)

length_df = pd.DataFrame({"length": lengths,
                          "count": counts})

agg = (
        length_df
        .groupby("length", as_index=False)
        .sum()
        .sort_values("length")
    )

# cumulative retention from longest downward
agg["cumulative_prop"] = (
    agg["count"][::-1].cumsum()[::-1] /
    agg["count"].sum()
)

# determine MINLEN
minlen = agg.loc[
    agg["cumulative_prop"] >= .95,
        "length"].min()

minlen = max(minlen, 36)

minlen_call = f"MINLEN:{minlen}"

trim_params = {
   "adapter": adapterkit_call,
    "headcrop": headcrop_call,
    "crop": crop_call,
    "sliding_window": sliding_window_call,
    "minlen": minlen_call,
    "polyG": polyG_call,
    "ployA": polyA_call
}

with open("trim_params.yml", "w") as f:
    yaml.dump(trim_params, f)