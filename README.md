![banner](https://github.com/gfalbery/INLA_N_Out/blob/master/INLA_N_Out.jpg)

## Fit similarity matrices and spatial locations in linear models!

This repository includes code to help with fitting Animal Models to Spatial/Social data in `INLA`.

The response variables in this toy analysis are social centrality traits in a real social network of wild red deer.

The variance components (random effects) are: 
- genetic similarity matrix
- social association matrix (constructed with the package `igraph`)
- home range overlap matrix (constructed in the package `AdeHabitatHR`) 
- INLA spatial SPDE effects

The fixed explanatory variables are a set of phenotypic and environmental traits for the deer (that we're not all that interested in here).

All analyses are carried out in `INLA`

No data are available to share, but the component manipulations should be easily transferrable to your data. If they are not, shoot me an email at gfalbery@gmail.com and I'll help troubleshoot.

Accompanies this preprint: 

And this review: 
