
# Setting up test data frames ####

library(tidyverse)

# Making genomic relatedness matrix ####

if(file.exists("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/GRM.Rdata")) load("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/GRM.Rdata") else{
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

if(file.exists("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/HRO.Rdata")) load("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/HRO.Rdata") else{
  
  memory.limit(56000)
  
  Censuses2<-droplevels(Censuses[Censuses$Code%in%names(table(Censuses$Code)[table(Censuses$Code)>5])&!Censuses$Code=="",])
  
  Censuses2<-droplevels(na.omit(Censuses2[Censuses2$Hind == "Y",c("EastingJ","NorthingJ","Code")]))
  
  Censuses2 <- droplevels(Censuses2[Censuses2$Code%in%names(table(Censuses2$Code))[table(Censuses2$Code)>5],])
  
  spdf <- SpatialPointsDataFrame(data = Censuses2, coords = Censuses2[,c("EastingJ","NorthingJ")])
  
  spdf<-spdf[,"Code"]
  
  kudl <- kernelUD(spdf, same4all=TRUE, grid=500)
  
  HRO <- kerneloverlaphr(kudl, percent = 95, method = "BA")
  
  save(HRO, file = "HRO.Rdata")
  
}

# Making social overlap ####

if(file.exists("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/AM.Rdata")) load("D:/Scripy/Deer Analyses/Extrication Social Spatial Code/AM.Rdata") else{
  
  Censuses3 <- droplevels(Censuses[Censuses$Season == "Rut"&Censuses$Hind == "Y",])
  
  M <- with(Censuses3, table(Code, GroupDate)) %>% as.matrix
  SocGraph <- graph.incidence(M, weighted = T)
  DeerProj <- bipartite.projection(SocGraph)$proj1
  AssMat <- DeerProj %>% get.adjacency(attr = "weight") %>% as.matrix
  diag(AssMat) <- table(Censuses3$Code)
  
  N <- nrow(AssMat)
  
  A <- matrix(rep(table(Censuses3$Code), N), N)
  B <- matrix(rep(table(Censuses3$Code), each = N), N)
  
  AM <- AssMat/(A + B - AssMat)
  
  AM[AM>0.25] <- 0.25
  
  save(AM, file = "AM.Rdata")
}

# Making final dataset ####

DeerNames <- reduce(list(as.character(SocialHinds$Name), 
                         colnames(HRO), 
                         colnames(GRM),
                         colnames(AM)), intersect)

DeerNames <- sort(DeerNames)

Resps <- c("Degree", "Strength", "Eigenvector", "Eigenvector2", "GroupSize")

Cols <- c("Name", Resps, "Year", "fYear", "PopN", "ReprodStatus", "Age", "Age2", "Survived1", "E", "N")

TestHinds <- na.omit(SocialHinds[, Cols]) # removing NA's and individuals that died halfway through the year

WITCH <- which(sapply(TestHinds,is.numeric))
WITCH <- WITCH[!names(WITCH)%in%c("Eigenvector","E","N")]
TestHinds[,WITCH] <- apply(TestHinds[,WITCH], 2, scale)

TestHinds <- TestHinds %>% filter(Name%in%DeerNames)



Cordf <- data.frame(
  
  Space = c(HRO[DeerNames, DeerNames]),
  Social = c(AM[DeerNames, DeerNames]),
  GRM = c(GRM[DeerNames, DeerNames])
  
)

Cordf[,c("Deer1","Deer2")] <- expand.grid(DeerNames,DeerNames)

Cordf <- Cordf[lower.tri(HRO[DeerNames, DeerNames]),]
Cordf[,c("Space", "Social", "GRM")] %>% GGally::ggpairs(lower = list(continuous = "smooth"))
