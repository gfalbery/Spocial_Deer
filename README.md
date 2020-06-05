![banner](https://github.com/gfalbery/Spocial_Deer/blob/master/Banner.jpg)

## Fit similarity matrices and spatial locations in linear models!

This repository includes code to help with fitting Animal Models to Spatial/Social data in `INLA`.

The response variables in this analysis are social network position traits in a social network of wild red deer.

The variance components (random effects) are: 
- genetic similarity matrix
- home range overlap matrix (constructed in the package `AdeHabitatHR`) 
- INLA spatial SPDE effects (lifetime centroids, annual centroids, and spatiotemporally varying)

The fixed explanatory variables are a set of phenotypic and environmental traits for the deer, all of which were compared with the spatiotemporal effects.

All analyses are carried out in `INLA`.

The component manipulations should be easily transferrable to your data. If they are not, shoot me an email at gfalbery@gmail.com and I'll help troubleshoot.

Accompanies this preprint: 
