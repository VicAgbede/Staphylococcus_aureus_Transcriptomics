# Transcriptomic Profiling of a cvfE Virulence Regulator Knockout in Staphylococcus aureus

## Project Overview
This repository contains an end-to-end computational biology pipeline profiling differential gene expression (DGE) in *Staphylococcus aureus* following the knockout of the key virulence regulator gene, *cvfE*. By analyzing downstream transcriptomic shifts, this project maps how *S. aureus* alters its metabolic machinery and virulence expression, offering critical insights into novel therapeutic targets to combat antimicrobial resistance (AMR).

---

## Biological Context & Dataset
* **Reference Genome:** *Staphylococcus aureus* NCTC 8325 (Alignment Reference)
* **Experimental Strain Background:** *Staphylococcus aureus* RN4220
* **Dataset Source:** NCBI Sequence Read Archive (SRA) / BioProject
* **BioProject ID:** PRJDB5479
* **Reference Publication:** Panthee et al. (PMC9660545)
* **Experimental Design:** Transcriptional profiling across 9 replicates evaluating a three-condition genetic matrix: Wild-Type baseline, ∆cvfE knockout mutant, and the complemented rescue strain.

---

## Technical Pipeline & Workflow
The analysis transitioned from raw sequencing reads to biological interpretation using a structured bioinformatics infrastructure:

1. **Quality Control & Trimming:** Raw FastQ sequencing data processing and adapter removal.
2. **Alignment & Quantification:** Mapping sequence reads against the *S. aureus* reference genome using a Linux environment via **WSL Ubuntu**.
3. **Statistical Modeling:** Differential expression analysis, dispersion estimation, and log-fold change shrinkage executed using **DESeq2** in R.
4. **Data Transformations:** Calculating gene-by-gene Z-scores to standardize highly variable expression scales for downstream clustering.

---

## Key Visualizations & Biological Insights

### 1. The Global Transcriptomic Landscape (Volcano Plot)
To assess the overall magnitude and statistical significance of the genetic shift between the Mutant and Wild-Type strains, a customized volcano plot was generated using `ggplot2` and `ggrepel`.

<img width="543" height="479" alt="Mutant_vs_WildType_Vol_top10_gene_plot" src="https://github.com/user-attachments/assets/1e7a5a82-b016-4d84-9284-6906c23b2ba3" />

* **Key Insight:** The plot demonstrates a massive vertical separation of high-significance responders. A notable heavy-hitter on the upregulated side includes **`asnC`** (associated with amino acid stress response modulation), which climbs to an astronomical statistical ceiling near $-\log_{10}(P_{\text{adj}}) = 78$. Conversely, **`pyrE`** (associated with pyrimidine down-regulation) displays deep, undeniable down-regulation on the left wing, marking it as a critical downstream component suppressed by the *cvfE* knockout. Key components of the staphyloxanthin biosynthesis cluster were also strongly implicated as drivers of the transcriptional shift.

### 2. Hierarchical Clustering (Expression Signatures Heatmap)
To visualize clean visual contrast and directional expression profiles across all biological replicates, a diverging Red-Yellow-Blue (`RdYlBu`) heatmap was constructed using `ComplexHeatmap`. The matrix displays a balanced panel of the top 5 most significant upregulated and top 5 most significant downregulated genes.

<img width="543" height="479" alt="Mutant_vs_WildType_HeatMap" src="https://github.com/user-attachments/assets/3dd0ffb1-c146-4145-94d2-f38a408c973f" />

* **Key Insight:** The dual dendrograms reveal replicate consistency. The algorithm splits the 9 samples cleanly by biological condition without prior grouping inputs, which explicitly map to:
  - **Wild-Type Baseline** (DRR084259–61)
  - **∆cvfE Knockout Mutant** (DRR084262–64)
  - **Complemented Strain** (DRR084265–67)
  
  The heatmap beautifully shows that when the regulator is knocked out, the gene expression patterns completely flip—and when the gene is complemented, the expression profile reverts right back to the wild-type baseline, proving that reintroducing the gene successfully reversed the mutant defects.

---

## Repository Structure
* `/data` : Contains final processed data tables and top gene lists.
* `/scripts` : Production R scripts for DESeq2 pipelines and figure rendering.
* `README.md` : Project abstract, pipeline description, and visual findings.

---

## Tools & Libraries Used
* **Environment:** Linux (WSL Ubuntu), RStudio
* **Packages:** `DESeq2`, `ComplexHeatmap`, `ggplot2`, `ggrepel`, `dplyr`, `circlize`
