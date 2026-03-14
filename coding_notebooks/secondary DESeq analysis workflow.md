You are now working **downstream of differential expression**, with a dataset representing the **union of DEGs across contrasts**. The useful structure is:

- rows: genes
    
- columns: contrasts
    
- objects already created
    
    - **significance matrix** (`padj < 0.05`)
        
    - **direction matrix** (`-1, 0, +1` or categorical up/down)
        

The goal is to extract **biological patterns across contrasts**, not compute DE again. The workflow below assumes you still retain the original **DESeq2 results tables** (log2FC, padj, baseMean, etc.), which should remain accessible because they are needed later.

---

# Recommended analysis workflow

## 1. Organize the core matrices

Maintain three matrices derived from the DESeq results.

**1. Significance matrix**

```
S[g, c] ∈ {0,1}
```

gene g significant in contrast c.

**2. Direction matrix**

```
D[g, c] ∈ {-1,0,1}
```

-1 = downregulated  
0 = not significant  
1 = upregulated

**3. Effect size matrix**

```
L[g, c] = log2FoldChange
```

These serve different purposes:

|Matrix|Use|
|---|---|
|S|DEG overlap patterns|
|D|regulatory direction patterns|
|L|clustering and trend analysis|

---

# Step 2 — Basic DEG distribution analysis

Before clustering or enrichment, quantify **how many genes change per contrast**.

Objects required:

- significance matrix `S`
    

Analyses:

### DEG counts per contrast

```
deg_counts = S.sum(axis=0)
```

Interpretation:

- identifies strongest perturbations
    
- detects contrasts with weak signal
    

### DEG frequency per gene

```
deg_freq = S.sum(axis=1)
```

Interpretation:

- genes responding to many conditions
    
- genes specific to one treatment
    

Useful visualizations:

- histogram of `deg_freq`
    
- ranked DEG counts per contrast
    

This step identifies **global response structure**.

---

# Step 3 — DEG overlap analysis

Use the significance matrix.

Goal: determine **which contrasts share transcriptional responses**.

### Contrast similarity

Compute similarity between columns of `S`.

Examples:

**Jaccard similarity**

```
intersection / union
```

**Hamming distance**

```
mean(S[:,i] != S[:,j])
```

Produces:

```
contrast_similarity_matrix
```

Use:

- clustering
    
- heatmap
    

Interpretation:

- which treatments/timepoints produce similar DEG sets
    

---

# Step 4 — Gene response patterns

Use the **direction matrix**.

Goal: identify **genes with shared regulation patterns** across contrasts.

Example gene pattern:

```
[0,0,1,1,0,-1]
```

Meaning:

- induced in some contrasts
    
- suppressed in others
    

Useful operations:

### Unique response signatures

```
group genes by row pattern of D
```

Example clusters:

|Pattern|Interpretation|
|---|---|
|up only in APAP contrasts|APAP-specific genes|
|up in CO and CO+APAP|CO driven genes|
|opposite regulation|antagonistic pathways|

This step extracts **logical regulatory rules**.

---

# Step 5 — clustering genes by expression response

Use **log2FC matrix (L)**.

Recommended approach:

### K-means clustering

Input:

```
L[g,c]
```

Purpose:

- identify **co-regulated gene groups**
    

Typical preprocessing:

```
standardize rows
```

Clustering output:

```
cluster_labels[g]
```

Then compute:

```
cluster_centroids
```

Interpretation:

Example clusters:

|Cluster|Pattern|
|---|---|
|early APAP induction||
|CO suppression||
|late response genes||
|sustained response||

This captures **temporal trends** across contrasts.

---

# Step 6 — treatment-specific gene sets

Use metadata from your design:

```
co_time
apap_time
```

Define logical gene groups such as:

### CO-specific genes

```
significant in CO contrasts
not significant in APAP contrasts
```

### APAP-specific genes

```
significant in APAP contrasts
not significant in CO contrasts
```

### interaction genes

```
significant only in combined CO+APAP
```

These groups directly address **treatment interactions**.

---

# Step 7 — time-dependent response analysis

Your data contains timepoints:

```
4hr
8hr
24hr
```

Use the log2FC matrix to track:

### temporal trends

Example:

```
CO_4hr
CO_24hr
```

Check genes where:

```
|L_CO24| > |L_CO4|
```

Interpretation:

- progressive induction
    
- transient responses
    

You can cluster genes based on **time-ordered contrasts**.

---

# Step 8 — gene enrichment analysis

Perform enrichment on gene sets derived above.

Common inputs:

- clusters from k-means
    
- treatment-specific gene sets
    
- high DEG frequency genes
    

Typical analyses:

### Gene Ontology enrichment

categories:

- Biological Process
    
- Molecular Function
    
- Cellular Component
    

### pathway enrichment

examples:

- KEGG
    
- Reactome
    

Interpretation:

|cluster|enriched pathway|
|---|---|
|early APAP genes|oxidative stress|
|CO genes|hypoxia response|
|interaction genes|inflammatory signaling|

---

# Step 9 — regulatory network inference (optional)

Although you are not building WGCNA networks, you can still infer regulation.

Examples:

### transcription factor enrichment

Check if gene clusters are enriched for targets of specific TFs.

Common databases:

- ChEA
    
- ENCODE
    
- TRANSFAC
    

This identifies **candidate regulators** of DEG programs.

---

# Step 10 — integrate DESeq statistics

Do not ignore the full DESeq output.

Useful attributes include:

|column|use|
|---|---|
|baseMean|filter low expression genes|
|log2FoldChange|effect size clustering|
|lfcSE|reliability|
|stat|ranking genes|
|padj|significance|

These allow additional analyses such as:

### weighted clustering

cluster using **log2FC weighted by significance**

### gene ranking

identify **most extreme responders**

---

# Data objects you should keep available

Essential objects for downstream analysis:

```
S  = significance matrix
D  = direction matrix
L  = log2FoldChange matrix
P  = padj matrix
C  = contrast metadata
G  = gene annotation table
```

Metadata should include:

```
contrast
co_time
apap_time
treatment_type
```

---

# Summary workflow

1. Build matrices
    
    - significance
        
    - direction
        
    - log2FC
        
2. Quantify DEG distributions
    
    - genes per contrast
        
    - contrasts per gene
        
3. Compare contrasts
    
    - DEG overlap
        
    - similarity clustering
        
4. Identify gene response signatures
    
    - pattern grouping
        
5. Cluster genes by log2FC patterns
    
6. Define treatment-specific gene sets
    
7. Analyze temporal trends
    
8. Perform enrichment on gene groups
    
9. Identify transcriptional regulators
    
10. Integrate full DESeq statistics
    

---

If useful, the next step can be constructing a **structured analysis plan specific to your contrasts**, because the **CO / APAP factorial time design** creates several biologically distinct comparisons that can be exploited systematically rather than treated as unrelated contrasts.