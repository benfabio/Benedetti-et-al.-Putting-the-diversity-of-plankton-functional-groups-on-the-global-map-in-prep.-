# 19/05/2021

# Benedetti et al. Putting the diversity of plankton functional groups on the global map (in prep.)
This repository contains parts of the R scripts developed during my postdoc at ETH Zürich, D-USYS, IBP, UP group.
The present R scripts were developed to produce the results of the manuscript in preparation where we simulatenously model and map the global distribution of species diversity and composition for various key plankton functional groups (Coccolithophores, Diatoms, Krill, Copepods, Jellyfish...).

- Scripts labelled as 'RSCRIPTBATCH' are those that were run on a local cluster to perform those functions defined within the script in parallel using the 'parallel' R package (R Core Team, 2020) within a mclapply().
- Scripts that were given a number correspond to those R scripts were I formated and mined the main datasets involved. They usually contain organized sequences of code where data are being read, examined, reformatted, analyzed and plotted.
- "OVERSEE" just corresponds to name I gave to the project overarching my main research activities at ETH Zürich (modelling and forecast of global plankton biodiversity). Therefore, other repositories corresponding to other studies/manuscripts also fall within the 'OVERSEE' umbrella.
- The end of the scripts' names usually indicate their purposes. A more detailed list of the scripts' goals and content is usually given in the beginning of the scripts.
