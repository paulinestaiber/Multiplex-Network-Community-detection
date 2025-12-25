# Disease-Focused Community Detection in a Multiplex Network 

This repository implements a reproducible workflow for identifying and
functionally characterizing disease-enriched gene modules in multiplex biological networks using Infomap community detection.

The pipeline was developed in the context of a Master's internship project
at the Ludwig Boltzmann Institute for Network Medicine.

## Workflow Overview

The analysis consists of four main steps:

1. LCC-based layer selection  
2. Construction of a disease-focused multilayer network for Infomap  
3. Identification of disease-gene-enriched communities  
4. Functional enrichment analysis of disease-gene-enriched communities 

## Dependencies
### Python
- Python 3.9
- pandas
- numpy
- scipy
- gseapy
- matplotlib
- plotly

### R
- R ≥ 4.0
- InfoWalkR
- igraph
- dplyr
- tidyr
- ggplot2

### External tools
- Infomap (command-line version)  
https://www.mapequation.org/infomap/#Install

## Workflow 
## 1. LCC-based Layer Selection

This step identifies informative layers in which disease-associated genes
form larger connected components compared to random expectation.

### Input
- Network data provided as graph files for multiple layers  
- A list of disease-associated genes 

### Method
For each layer, the size of the largest connected component (LCC) formed by
the disease genes is computed.  
The observed LCC size is compared against a null model generated from random
gene sets of identical size.  

Based on the null distribution, a z-score is calculated for each layer,
quantifying how strongly the observed LCC size deviates from random
expectation.  
Layers in which disease genes form a significantly larger LCC than expected
by chance are considered disease-informative.

### Output
- A set of selected layers and corresponding z-scores in which disease genes
  show statistically significant connectivity



## 2. Construction of Infomap Input for Multilayer Community Detection

This step constructs an Infomap-compatible multilayer network using the
disease-informative layers identified in step 1.

### Input
- Network graph files  
- `.csv` file containing selected layers and corresponding z-scores  

### Method
- All nodes are assigned unique IDs  
- Intra-layer edges are included with weight *w = 1*  
- Inter-layer edges connect node replicas across layers and are weighted
  using layer-specific z-scores derived from the LCC analysis 

Optional index mapping files are generated for reproducibility:
- Layer index file  
- Vertex (node) index file  

### Output
- A `.net` file in Infomap’s multilayer format containing:
  - Node definitions  
  - Intra-layer edges  
  - Inter-layer edges with z-score-based coupling  

Infomap can then be executed from the command line, for example:
infomap input_file.net output_folder --clu --tree -N 100



## 3. Identification of Modules Enriched in Disease Genes

This step identifies multilayer communities that are significantly enriched
for disease-associated genes.

### Input
- Infomap community assignments (`.clu`)  
- Vertex index file  
- List of disease-associated genes  

### Method
Nodes from the Infomap output are mapped back to gene names using the
vertex index file.  
For each community, the total number of genes and the number of disease genes
are counted.  
A hypergeometric overrepresentation test is applied to assess whether the
observed number of disease genes in a community is higher than expected by
chance.

### Output
- A CSV file listing significantly disease-enriched modules, including:
  - Module ID  
  - Number of genes  
  - Number of disease genes  
  - p-values  
  - adjusted p-values
  - enrichment  
- Additionally, a merged dataframe containing node ID, gene name, and module
assignment is generated.



## 4. Functional Enrichment Analysis of Disease-Enriched Modules

This step performs functional enrichment analysis on genes from
disease-enriched communities.

### Input
- CSV file of significant modules  
- Merged node–module dataframe  
- Functional annotation libraries supported by **gseapy**
  (e.g. GO, Reactome)  

### Method
Gene sets corresponding to disease-enriched modules are extracted and
submitted to overrepresentation analysis (ORA) using **gseapy**.  
Significantly enriched functional terms are identified.

### Output
- Bar plot showing the five most significantly enriched functional terms
  for the disease-associated module, using all genes in the
  disease-focused network as background  
- Bar plot showing the five most significantly enriched functional terms
  for the disease-associated genes within the disease-enriched module,
  using all genes in the disease-enriched module as background 







