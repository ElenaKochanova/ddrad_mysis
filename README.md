# ddRAD Mysis spp Bioinformatic pipeline

This repository contains the code for the manuscript:"Phylogeography and population genomics of three freshwater Mysis species at their post-glacial contact in Eastern Fennoscandia", investigating the structure and colonization history of Mysis relicta, M. salemaai and M. segerstralei using ddRAD sequencing data.

Raw sequence data is available on NCBI SRA (BioProject PRJNA1439230).

The filtered genomic datasets are available on DRYAD

## Overview

This pipeline includes the scripts 
- filtering+assembly_stacks.sh, which covers Clone filtering, Adapter trimming with cutadapt, and De novo assembly using Stacks
- Principal_component_analysis.R doing PCAs with adegenet
- Dxy_distances.R script calculating Dxy distances from ddRAD data using vcfR and dplyr packages
- snmf.sh describing snmf admixture analysis 

