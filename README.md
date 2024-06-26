# 07/10/2023

# Benedetti et al. (2023) Global gradients in species richness of marine plankton functional groups (Journal of Plankton Research)
This repository contains parts of the R scripts developed during my postdoc at ETH Zürich, D-USYS, IBP, UP group.
The present R scripts were developed to produce the results of Benedetti et at. (2023) where we simultaneously model and map the global distribution of species diversity and composition for various key plankton functional groups (Coccolithophores, Diatoms, Krill, Copepods, Jellyfish...).
See: https://academic.oup.com/plankt/article/45/6/832/7308688?searchresult=1 

- Scripts labelled as 'RSCRIPTBATCH' are those that were run on a local cluster to perform those functions defined within the script in parallel using the 'parallel' R package (R Core Team, 2020) within a mclapply().
- Scripts that were given a number correspond to those R scripts where I formatted and mined the main datasets involved. They usually contain organized sequences of code where data are being read, examined, reformatted, analyzed and plotted. These scripts often correspond to the preparation. The number given to a script mainly corresponds to a temporal marker, and not necessarily a direct sequence from one numbered script to another. Considering that multiple studies fall under the OVERSEE project, some of the numbered scripts are deposited in another repository. This is why there might be gaps between the numbers of the present R scripts.
- "OVERSEE" just corresponds to the name I gave to the project overarching my main research activities at ETH Zürich (modelling and forecast of global plankton biodiversity). Therefore, other repositories corresponding to other studies/manuscripts also fall within the 'OVERSEE' umbrella.
- The end of the scripts' names usually indicate their main purposes. A more detailed list of the scripts' goals and content is usually given in the beginning of the numbered scripts.
