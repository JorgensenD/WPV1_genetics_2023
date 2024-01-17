# Code for the paper "Genetic epidemiology of wild type 1 poliovirus in the last endemic region"

This repository provides the code used to analyse and plot the date presented in the above paper. This code covers both bayesian phylogeographic analysis in BEAST and the downstream analysis and plotting of the outputs of this work in R.

A brief description of the functionality of each of the included scripts is provided below.

## [1_BEAST_templates](./1_BEAST_templates)
Scripts are provided to handle sequence data and prepare for Markov Jump Discrete Trait Analysis (MJ-DTA) in BEAST. These analyses are carried out using R (.r scripts) and [BEAST version 1.10.4](https://beast.community/) (.xml scripts). The BEAST templates have also been tested in the 1.10.5pre version. R scripts include code to install and load required packages. [IQ-TREE (v1)](http://www.iqtree.org/) is used to build the starting trees for beast analysis but this could be substituted for any other maximum likelihood tree building software. If installed to the system PATH, IQ-TREE can be run from within R using the included pipeline.

## [2_BF_supported_movements](./2_BF_supported_movements)
Following analysis with BEAST, Bayes Factor (BF) support for the returned numbers of transitions between regions is assessed with [spreaD3 (v.0.9.7)](https://rega.kuleuven.be/cev/ecv/software/SpreaD3). Only transitions with BF support above 10 are included in plots of the numbers of movements between regions shown in the manuscript (Fig. 1 E&F).

The numbers of transitions and bounds on these numbers are matched manually from the output of the MJ-DTA visualised in Tracer 1.7.2 to the BF supported transitions returned from spreaD3. The matched data are provided in 'XXXX.csv'

## [3_Local_Transmsision_Lineages](./3_Local_Transmsision_Lineages)
A [consensus phylogeny](./3_Local_Transmsision_Lineages/) is generated from the MJ-DTA analysis using the treeannotator software bundled with BEAST and the burn-in as specified in the manuscript. Local Transmission Lineages (LTLs) are then produced in R by splitting the phylogeny wherever a change in state (location) is inferred on the tree. This can be used to produce figures showing all of the LTLs or further summarised. Code to carry out the tree splitting, plotting and summary plotting are provided. 

## [4_Orphan_Lineages](./4_Orphan_Lineages)
Orphan lineages are defined in the Global Polio Eradication Initiative (GPEI) based on the pairwise genetic distance between sequences rather than on phylogenetic relationships. The R code used to generate these distances and produce the plots in the manuscript is provided. This code extracts the location of genetic sequences from the dataframe produced for the eariler analyses, although these could be extracted directly from the initial metadata without running the preceeding steps with modifications to the code provided.

## Other supplementary information
Background WPV1 sequecnes were extracted from GenBank. These are not reproduced here as they are publicly available on that platform. Accession numbers for the included sequences are provided here.
Non-polio AFP rates and numbers of environmental surveillance sites provided in the supplementary material are extracted from the World Health Organization Polio Information System. [A public version of this database is available](https://extranet.who.int/polis/public/CaseCount.aspx). More complex data is available on request from the WHO. 
As these figures cannot be recreated without the required data access the code is not reporduced here. Code can be provided on request.



