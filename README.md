# PIP-eco
Pathotype Identification Pipeline for *Escherichia coli* (PIP-eco)

Overview:
PIP-eco is a systematic framework designed to identify the pathotype of *Escherichia coli* (*E. coli*) using Whole Genome Sequencing (WGS) data, especially for those with hybrid pathotypes. PIP-eco pipeline consists of three infrastructures: **Marker gene alignment** process, **Pan-phylogenetic analysis** process, and **Pathogenicity Islands (PAIs) analysis** process. It accepts assembled bacterial strain collections as input, which can be either NCBI RefSeq records or user's own data in fasta format.
In the PIP-eco pipeline, genome annotation on the input WGS data is performed. Follwing this, the pathotype is determined based on marker genes. Additionally, by conducting phylogenetic analysis based on pan-genome analysis, the genetic distances are investigated, thus effectively discriminating hybrid pathotypes. Through these processes, the PIP-eco pipeline is utilized not only for pathotype assignment but also for tracing the trajectories of pathogenic factors. 
The Processing within the PIP-eco pipeline uses publicly available tools: PROKKA, USEARCH, MUSCLE, and MAFFT. 

## Table of contents
  * [Overview of dependencies](#overview-of-dependencies)
  * [Pipeline overview](#pipeline-overview)
  * [Quick start and installing dependencies](#quick-start-and-installing-dependencies)
  * [Usage](#usage)  


## Pipeline overview
![PIPeco](/PIPeco.png)

## Quick start and installing dependencies

#### 1. Download RADAR pipeline on Github and syncronized environment.
```
conda create -y pathotype.yaml
git clone https://github.com/SBL-Kimlab/PIP-eco.git
```
#### 2. Installing dependencies 
##### Overview of dependencies:
  * Genome annotation: [Prokka](https://github.com/tseemann/prokka)
  * Local alignment tool: [Usearch](https://www.drive5.com/usearch/)

## Usage

In the RADAR pipeline, there are eight different modules in detail. Each process is performed according to defined modules. Users can directly use the individual modules as shown below, so all processes can be executed at once.


```
#Before executing the PIP-eco pipeline, it needs to declare /include/include.ipynb.

import os; import os.path as path
from time import sleep

path_root = path.abspath( path.join( os.getcwd(), "..", ".." ) )
path_local = path_root +  "/pipeco"; path_include = path_local +  "/include"
file_include = path_include +  "/include.ipynb"
%run $file_include
```

```
#PIP-eco pipeline excution 
os.chdir( path_local )
pipeco = pipeline()

pipeco.method.prokka() 
pipeco.method.marker_alignment() 
pipeco.method.draft_assignment() 
pipeco.method.pan_phylo() 
pipeco.method.pai_analysis() 
```

## Reference
1. *Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069.*
2. *Edgar, R. C. (2010). **Search and clustering orders of magnitude faster than BLAST.** _Bioinformatics_, _26_(19), 2460-2461.*
3. *Edgar, R. C. (2004). MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic acids research, 32(5), 1792-1797.*
4. *Katoh, K., Misawa, K., Kuma, K. I., & Miyata, T. (2002). MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic acids research, 30(14), 3059-3066.*
