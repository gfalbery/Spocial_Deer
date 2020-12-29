
# X_Making SocialHinds.rds ####

library(tidyverse)

Covar <- c("Age", "ReprodStatus", "Year", "PopN", "NObs")

SocResps <- Resps <- c("GroupSize", "Degree",
                       "Strength", "Strength_Mean",
                       "Eigenvector", "Eigenvector_Weighted",
                       "Betweenness", "Clustering")

SpatialVars <- c(
  "AnnualDensity", "HRA", 
  "f(IndexPhylo, model='generic0', Cmatrix = GRMatrix, constr = TRUE, param = c(0.5, 0.5))", 
  "f(IndexSpace, model='generic0', Cmatrix = SpaceContacts, constr = TRUE, param = c(0.5, 0.5))"
  
)

Cols <- c("Name", 
          Resps, 
          Covar, SpatialVars[1:2],
          "E", "N", 
          "LifetimeE", "LifetimeN", 
          "MeshInclude")

SocialHinds %>% dplyr::select(Cols) %>% 
  mutate_at("Name", ~as.numeric(as.factor(.x))) %>% 
  saveRDS("SocialHinds.rds")
