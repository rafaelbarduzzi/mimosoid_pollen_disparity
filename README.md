# Analyses of pollen morphological disparity in the Mimosa and Stryphnodendron clades (Mimoseae, Leguminosae)

## Reference publication

For details on methodology and results, see the published paper:  
https://doi.org/10.1093/aob/mcaf213

## Description

This repository contains data and R scripts for analyzing pollen morphological disparity across the Mimosa and Stryphnodendron clades. The data include pollen morphological traits and biome occurrence for the Mimosa clade (genera *Adenopodia*, *Mimosa*, and *Piptadenia*), the Stryphnodendron clade (genera *Gwilymia*, *Marlimorimia*, *Microlobius*, *Naiadendron*, *Parapiptadenia*, *Pityrocarpa*, and *Stryphnodendron*), and the outgroups *Inga*, *Lachesiodendron* and *Senegalia*. The analyses include ancestral biome estimation, ordination of morphological data, calculation of morphospace occupancy metrics, and statistical tests to assess differences between groups/distributions (Wilcoxon rank-sum test and Bhattacharyya coefficient). See the referenced publication for more details.

## Data folder

- **eco_data.csv**: Biome occurrence data. The coded discrete states are provided in **eco_data_states.csv**.

- **eco_data_states.csv**: States for the coded data matrix **eco_data.csv**.

- **mimo_stryph_phylo.tre**: Complete phylogenetic hypothesis synthesized from Vasconcelos et al. 2020, Borges et al. 2022, and Lima et al. 2022.

- **pollen_data.csv**: Pollen morphological data, including discrete and continuous traits. The coded discrete states are provided in the main text. Individual matrices are also available in Morphobank as .nexus and .tnt files (http://dx.doi.org/10.7934/P5877).

## Scripts folder

- **1.anc_biomes_estimation.R**: Estimates ancestral biomes states for the study groups.

- **2.disparity_analyses.R**: Performs all data preparation (e.g., discretization of continuous traits) and downstream analyses (ordination of morphological data, pollen morphospace plotting, calculation of morphospace occupancy metrics, and statistical tests) presented in the main results.

- **3.supp_sum_of_ranges.R**: Calculates an additional morphospace occupancy metric (Sum of Ranges).

- **4.disparity_analyses_preoase.R**: Performs the same analyses as **disparity_analyses.R**, but with pre-ordination ancestral state estimation for the generation of the phylomorphospace.

### Usage

1. Install R (version >= 4.0) and required packages (see individual scripts for required library calls).
2. Run scripts in the numbered order.

### Citation

If you use the data or scripts available here, please cite both the referenced publication and this repository: 
Barduzzi, R. 2025. Data from "Pollen morphological disparity analyses for Mimosa and Stryphnodendron clades". GitHub. Available at https://github.com/rafaelbarduzzi/mimosoid_pollen_disparity. Deposited September 17, 2025.

### Contact

For questions, please feel free to contact me at rfbarduzzi@gmail.com.