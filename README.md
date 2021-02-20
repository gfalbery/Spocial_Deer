![banner](https://github.com/gfalbery/Spocial_Deer/blob/master/Banner.jpeg)

## Fit similarity matrices and spatial locations in linear models!

This repository includes code to help with fitting Animal Models to Spatial/Social data in `INLA`.

The response variables in this analysis are social network position traits in a social network of wild red deer.

The variance components (random effects) are: 
- Genetic similarity matrix
- Home range overlap matrix (representing space sharing; constructed in the package `AdeHabitatHR`) 
- INLA spatial SPDE effects (representing point location effects)  

The fixed explanatory variables are a set of phenotypic and environmental traits for the deer, all of which were compared with the spatiotemporal effects, alongside two spatially varying components:
- Home range size (representing individual ranging capacity)  
- Local population density (representing the local spatial availability of social partners)  

All analyses are carried out in `INLA`, using animal models (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3737164/).

The component manipulations should be easily transferrable to your data. If they are not, shoot me an email at gfalbery@gmail.com and I'll help troubleshoot.

Accompanies this paper: https://onlinelibrary.wiley.com/doi/10.1111/ele.13684
