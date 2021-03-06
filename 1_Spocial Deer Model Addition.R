
# 1_Untangling Model Addition ####

library(INLA); library(ggregplot); library(tidyverse); library(fs); library(magrittr)

detach(package:raster); 
library(dplyr)

# Setting up ####

dir_create("Output Files")

CentralityList <- MeshList <- SPDEList <- TestDFList <- list()

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

# Making a dummy data frame ####

r <- 1

SocTestHinds <- SocialHinds %>% 
  
  dplyr::select(Cols) %>%
  mutate_at(c("Degree", "Strength"), sqrt) %>%
  mutate_at(c("Strength_Mean", "GroupSize"), ~kader:::cuberoot(.x)) %>%
  mutate_at(c("Betweenness"), ~log(.x + 1)) %>%
  mutate_at(c("Clustering"), ~.x^2) %>%
  mutate(fYear = as.factor(Year)) %>% 
  
  filter(Name%in%DeerNames, MeshInclude == 1) %>% 
  
  arrange(Name) %>%
  na.omit() %>%
  as.data.frame() # removing NA's and individuals that died halfway through the year

N <- nrow(SocTestHinds); print(N)

Classes <- SocTestHinds %>% sapply(class)

ToScale <- names(Classes[Classes %in% c("integer", "numeric")]) %>% 
  setdiff(c("Eigenvector","Eigenvector2","E","N"))

SocTestHinds %<>% mutate_at(ToScale, ~c(scale(.x)))

# Adding Covariates ####

SpaceContacts <- as(solve(HRO[DeerNames, DeerNames]),"dgCMatrix")
GRMatrix <- as(solve(GRM[DeerNames, DeerNames]),"dgCMatrix")

SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))

DeerNames <- reduce(list(as.character(SocTestHinds$Name), 
                         colnames(HRO), 
                         colnames(GRM)), intersect)

DeerNames <- sort(DeerNames)

FamilyList <- rep("gaussian", length(Resps))

# Removing outliers

OutlierLimits <- 
  list(c(-Inf, 50),
       c(-Inf, Inf),
       c(-Inf, 15),
       c(-Inf, 0.375),
       c(-Inf, Inf),
       c(-Inf, Inf),
       c(-Inf, 1000),
       c(0.375, Inf))

SpocialList <- list()

# Looping through response variables ####
# Use a loop, not a map/apply: these models take time.
# After this, run the script `X_Appending Centrality Lists.r` 
# to combine model files into a single list.

r <- 1

for(r in r:length(Resps)){
  
  print(SocResps[r])
  
  SocialHinds %>% pull(SocResps[r]) %>% 
    between(OutlierLimits[[r]][1], OutlierLimits[[r]][2]) %>%
    which() %>%
    slice(SocialHinds, .) ->
    
    SocTestHinds # Test data frame
  
  # Cleaning test data frame ####
  
  SocTestHinds %>% 
    
    dplyr::select(Cols) %>%
    mutate_at(c("Degree", "Strength"), sqrt) %>%
    mutate_at(c("Strength_Mean", "GroupSize"), ~kader:::cuberoot(.x)) %>%
    mutate_at(c("Betweenness"), ~log(.x + 1)) %>%
    mutate_at(c("Clustering"), ~.x^2) %>%
    mutate(fYear = as.factor(Year)) %>% 
    
    filter(Name%in%DeerNames, MeshInclude == 1) %>% 
    
    arrange(Name) %>%
    na.omit() %>%
    as.data.frame() %>%
    mutate_at(c(SocResps, "NObs"), ~c(scale(.x))) ->
    
    SocTestHinds 
  
  SpaceContacts <- as(solve(HRO[DeerNames, DeerNames]),"dgCMatrix")
  GRMatrix <- as(solve(GRM[DeerNames, DeerNames]),"dgCMatrix")
  
  SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
  SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
  
  if(1){ # Log-transforming HRA
    
    SocTestHinds %<>% 
      mutate_at("HRA", ~log(.x))
    
  }
  
  Classes <- SocTestHinds %>% sapply(class)
  
  ToScale <- names(Classes[Classes %in% c("integer", "numeric")]) %>% 
    setdiff(c("Eigenvector","Eigenvector2","E","N", "IndexSpace", "IndexPhylo"))
  
  SocTestHinds %<>% mutate_at(ToScale, ~c(scale(.x))) # Scaling
  
  # Conducting model addition (function found in the `ggregplot` package) ####
  # Adds each variable, keeping the best-fitting one, until eventually all of them are fitted
  # Then adds SPDE effect
  
  IM1 <- INLAModelAdd(
    
    Response = Resps[r],
    Data = SocTestHinds,
    Explanatory = Covar,
    Add = SpatialVars,
    # AllModels = T,
    Random = c("Name", "fYear"), 
    RandomModel = rep("iid", 2),
    Family = "gaussian", 
    Delta = -Inf,
    AddSpatial = T, Coordinates = c("E", "N"), #Boundary = RumBoundary,
    Groups = T, GroupVar = "Year"
    
  )
  
  SpocialList[[Resps[r]]] <- IM1
  
  IM1 %>% saveRDS(paste0("Output Files/", Resps[r], "Models.rds")) # Saving
  
}
