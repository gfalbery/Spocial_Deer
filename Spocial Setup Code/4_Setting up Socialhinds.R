
# Setting up test data frames ####

library(tidyverse); library(igraph); library(adehabitatHR)

# Making genomic relatedness matrix ####

if(file.exists("Output Files/GRM.Rdata")) load("Output Files/GRM.Rdata") else{
  makeGRM <- function(grm.object, ids.object, id.vector = NULL) {
    
    require(MASS)
    
    elements <- nrow(ids.object)
    X <- diag(elements)
    
    X[upper.tri(X, diag = TRUE)] <- grm.object[,ncol(grm.object)]
    X.grm <- X + t(X) - diag(diag(X))
    nam <- as.factor(ids.object[, 2])
    
    # THE FOLLWING PART IS ESSENTIAL TO MAKE THE MATRIX CONTAIN THE EXACT SAME IDS AS THE PHENOTYPE FILE
    rownames(X.grm) <- nam
    colnames(X.grm) <- nam
    
    if(!is.null(id.vector)) X.grm <- X.grm[which(rownames(X.grm) %in% id.vector),which(colnames(X.grm) %in% id.vector)]
    
    X.grm
    
  }
  
  grm.auto <- read.table("D:/Scripy/Deer Analyses/Genomic relatedness matrix code/Deer33.grm.gz")
  grm.id <- read.table("D:/Scripy/Deer Analyses/Genomic relatedness matrix code/Deer33.grm.id")
  
  GRM <- makeGRM(grm.auto, grm.id)
  save(GRM, file = "GRM.Rdata")
}

# making home range overlap ####

if(file.exists("Output Files/HRO.Rdata")) load("Output Files/HRO.Rdata") else{
  
  memory.limit(56000)
  
  Censuses2 <- droplevels(Censuses[Censuses$Code%in%names(table(Censuses$Code)[table(Censuses$Code)>5])&!Censuses$Code=="",])
  
  Censuses2 %>% mutate(
    
    EastingJ = Easting + runif(n(), -10, 10),
    NorthingJ = Northing + runif(n(), -10, 10)
    
  ) %>% filter(Hind == "Y") %>% dplyr::select(EastingJ, NorthingJ, Code) %>% na.omit ->
   
  Censuses2
  
  Censuses2 <- droplevels(Censuses2[Censuses2$Code%in%names(table(Censuses2$Code))[table(Censuses2$Code)>5],])
  
  spdf <- SpatialPointsDataFrame(data = Censuses2, coords = Censuses2[,c("EastingJ","NorthingJ")])
  
  spdf <- spdf[,"Code"]
  
  kudl <- kernelUD(spdf, same4all = TRUE, grid = 500)
  
  HRO <- kerneloverlaphr(kudl, percent = 95, method = "BA")
  
  save(HRO, file = "HRO.Rdata")
  
}

# Making final dataset ####

DeerNames <- reduce(list(as.character(SocialHinds$Name), 
                         colnames(HRO), 
                         colnames(GRM)), intersect)

DeerNames <- sort(DeerNames)
