# Code for the paper "Genetic epidemiology of wild type 1 poliovirus in the last endemic region"

This repository provides the code used to analyse and plot the date presented in the above paper. This code covers both bayesian phylogeographic analysis in BEAST and the downstream analysis and plotting of the outputs of this work in R. 

A brief description of the functionality of each of the included scripts is provided below.

The code has been tested on R version 4.4.1, BEAST version 1.10.4 and 1.10.5 and SpreaD3 version 0.9.7. Package versions for R are included in the respective code files, BEAST does not have any additional dependencies. The code is independent of operating system but is tested primarily on MAC OS 15. The Included scripts in [Figures](./Figures) and minimal working data in [Plot_data](./Plot_data) can be used to recreate the figures in the main text of the paper. This is possible by cloning this repository locally or by manually downloading the relevant files. Plotting from these data should be almost instantaneuos on a standard desktop computer. To run the full analysis takes around 5 days from start to finish, with each iteration of the tree generating process (8 used here) and the split tree plotting (here we use a subset of 300 posterior phylogenies) split across multiple CPU cores. Here we use 8 cores each with 4gb of memory for the tree generation process and 12 cores with a pooled 16gb of ram for the tree plotting. Both of these can be carried out with more restricted computational resources but the analysis will take significantly longer to run.


## [1_BEAST_templates](./1_BEAST_templates)
**Sequence data must be provided with names in the format "SEQUENCENAME_SAMPLETYPE_ADMIN0NAME_SAMPLEDECIMALDATE_CLUSTER_REGION".** A map of the regions used in this analysis is shown in Fig.1.

Scripts are provided to handle sequence data and prepare for Markov Jump Discrete Trait Analysis (MJ-DTA) in BEAST. These analyses are carried out using R (.r scripts) and [BEAST version 1.10.4](https://beast.community/) (.xml scripts). The BEAST templates have also been tested in the 1.10.5pre version. R scripts include code to install and load required packages. [IQ-TREE (v1)](http://www.iqtree.org/) is used to build the starting trees for beast analysis but this could be substituted for any other maximum likelihood tree building software. If installed to the system PATH, IQ-TREE can be run from within R using the included pipeline.

In this analysis tree generation and inference of movement events are split into two independent processes with BEAST run twice. BEAST is first run in parallel to generate a set of phylogenies over which movements can be inferred. These are then combined with logcombiner to create a set of 1000 trees which are iteratively investigated with MJ-DTA. R files `1_meta.R`, `2_treegen.R` and `3_template_load.R` take the data from cleaned metadata tagged fasta data to starting trees for an analysis in BEAST and then feed this into a BEAST template file `TREEGEN_TEMPLATE.xml`. `4_template_MJ.R` feeds the phylogenies generated from the first BEAST run into an MJ-DTA analysis. After generating a consensus phylogeny with treeannotator, this can be visualised with `5_MCC_tree_viz.R`.

## [2_BF_supported_movements](./2_BF_supported_movements)
Following analysis with BEAST, Bayes Factor (BF) support for the returned numbers of transitions between regions is assessed with [spreaD3 (v.0.9.7)](https://rega.kuleuven.be/cev/ecv/software/SpreaD3). Only transitions with BF support above 10 are included in plots of the numbers of movements between regions shown in the manuscript (Fig. 1 E&F).

The numbers of transitions and bounds on these numbers are matched manually from the output of the MJ-DTA visualised in Tracer 1.7.2 to the BF supported transitions returned from spreaD3. The matched data are provided in 'supported_out.csv'. Total numbers of movements into and  out of each region are estimated separately in the MJ-DTA analysis. These estimates and their bounds are given in`movement_bounds.csv`. Code to produce circos and bar plots from these data is provided in `movement_plots.R`.

## [3_Local_Transmsision_Lineages](./3_Local_Transmsision_Lineages)
A [consensus phylogeny](./MJ_MCC_CA_DS.trees/) is generated from the MJ-DTA analysis using the treeannotator software bundled with BEAST and the burn-in as specified in the manuscript. Local Transmission Lineages (LTLs) are then produced in R by splitting the phylogeny wherever a change in state (location) is inferred on the tree. This can be used to produce figures showing all of the LTLs or can be further summarised. Code to carry out the tree splitting, plotting and summary plotting are provided. The code in `bezier_map_anim.R` can be used to animate the LTLs over time on a map.

## Tipswap analysis 
The tipswap analysis is carried out using the template provided: [MJ_TIPSWAP_TEMPLATE.xml](./1_BEAST_templates/MJ_TIPSWAP_TEMPLATE.xml). This is the same as the template for a standard markov jump analysis with an additional tipswap operator. To ensure at least 95% of the tips are swapped at each iteration the weight of the tipswap operator is calculated as follows:

Where threshold = 0.95, we calculate X based on the number of tips used in the analysis

`x = log(1-threshold)/log((nb.tips-2)/nb.tips))`

We can then use this to calculate the total number of randomisations using the number of samples to be collected during the run (logevery/chain length)

`nb.randomistations = x*nb.samples`

The operator weight is then calculated with this number of randomisations. Total operators is the sum of the weight of all other operators.

`weight = nb.randomisations/MCMC.length*total.operators`

The posterior probabilities for the tipswap analysis are calculated as outlined in section 2 and can be used to correct the Bayes Factors estimated in the main analysis:

`BF = [PP_dta/(1-PP_dta)] / [PP_tsw/(1-PP_tsw)]`

Where either PP is equal to 1 it is approximated with `pp = (nb.samples - 1) / nb.samples` and where both are equal to 1 the BF is set to 1.

## Other supplementary information
Non-polio AFP rates and numbers of environmental surveillance sites provided in the supplementary material are extracted from the World Health Organization Polio Information System. [A public version of this database is available](https://extranet.who.int/polis/public/CaseCount.aspx). More complex data is available on request from the WHO. 
As these figures cannot be recreated without the required data access the code is not reporduced here. This can be provided on request. Code used to produce animated and static maps of viral movements can also be provided on request.



