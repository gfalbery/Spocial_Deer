
# Analysing social network ####

library(INLA); library(ggregplot); library(tidyverse); library(fs)

dir_create("Output Files")

CentralityList <- MeshList <- SPDEList <- TestDFList <- list()

Covar <- c("Age", "ReprodStatus", "Year", "PopN", "NObs")

SocResps <- Resps <- c("GroupSize", "Degree",
                       "Strength", "Strength_Mean",
                       "Eigenvector", "Eigenvector_Weighted",
                       "Betweenness", "Clustering")

Cols <- c("Name", 
          Resps, 
          Covar,
          "E", "N", 
          "LifetimeE", "LifetimeN", 
          "MeshInclude")

# Dummy df ####

r <- 1

SocTestHinds <- SocialHinds %>% 
  
  select(Cols) %>%
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

AnnualLocations = cbind(SocTestHinds$E, SocTestHinds$N) #cbind(SocTestHinds$E, SocTestHinds$N)

AnnualMesh <- inla.mesh.2d(loc = AnnualLocations, 
                           boundary = RumBoundary,
                           max.edge = c(10), cutoff = 1)

AnnualA3 <- inla.spde.make.A(AnnualMesh, loc = AnnualLocations) # Making A matrix

AnnualSPDE = inla.spde2.pcmatern(mesh = AnnualMesh, 
                                 prior.range = c(10, 0.5), 
                                 prior.sigma = c(.5, .5)) # Making SPDE

w.index.annual <- inla.spde.make.index('wAnnual', n.spde = AnnualSPDE$n.spde)

LifetimeLocations = cbind(SocTestHinds$LifetimeE, SocTestHinds$LifetimeN)

LifetimeMesh <- inla.mesh.2d(loc = LifetimeLocations, 
                             boundary = RumBoundary,
                             max.edge = c(10), cutoff = 1)

LifetimeA3 <- inla.spde.make.A(LifetimeMesh, loc = LifetimeLocations) # Making A matrix

LifetimeSPDE = inla.spde2.pcmatern(mesh = LifetimeMesh, 
                                   prior.range = c(10, 0.5), 
                                   prior.sigma = c(.5, .5)) # Making SPDE

w.index.lifetime <- inla.spde.make.index('wLifetime', n.spde = LifetimeSPDE$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(Covar, collapse = " + "))), 
                   data = SocTestHinds)

X <- as.data.frame(Xm)

X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("Survived10",
                                                  "ReprodStatusTrue.Yeld"))]) # Model Matrix


# Establishing model formulae ####

IID <- c("Name", "fYear")

f1 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ",
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))'))

f2 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ", 
                       'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + ", 
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))'))

f3 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ",
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + f(wAnnual, model = AnnualSPDE)"))

f4 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ", 
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + f(wLifetime, model = LifetimeSPDE)"))

f5 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ",
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))', 
                       " + ",
                       'f(IndexSpace, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + f(wAnnual, model = AnnualSPDE)"))

f6 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ",
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))', 
                       " + ", 
                       'f(IndexSpace, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + f(wLifetime, model = LifetimeSPDE)"))

FormulaList <- list(f1, f3, f4, f2, f5, f6)

names(FormulaList) <- c("Base", "SPDE", "LifetimeSPDE", "HRO", "HROSPDE", "HROLifetimeSPDE")

# Establishing distance matrices ####

SpaceContacts <- as(solve(HRO[DeerNames, DeerNames]),"dgCMatrix")
GRMatrix <- as(solve(GRM[DeerNames, DeerNames]),"dgCMatrix")

SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))

DeerNames <- reduce(list(as.character(SocTestHinds$Name), 
                         colnames(HRO), 
                         colnames(GRM)), intersect)

DeerNames <- sort(DeerNames)

FamilyList <- rep("gaussian", length(Resps))

OutlierLimits <- 
  list(c(-Inf, 50),
       c(-Inf, Inf),
       c(-Inf, 15),
       c(-Inf, 0.375),
       c(-Inf, Inf),
       c(-Inf, Inf),
       c(-Inf, 1000),
       c(0.375, Inf))

r = 1

for(r in r:length(SocResps)){
  
  print(SocResps[r])
  
  SocialHinds %>% pull(SocResps[r]) %>% 
    between(OutlierLimits[[r]][1], OutlierLimits[[r]][2]) %>%
    which() %>%
    slice(SocialHinds, .) ->
    
    SocTestHinds
  
  SocTestHinds %>% 
    
    select(Cols) %>%
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
  
  N <- nrow(SocTestHinds); print(N)
  
  WITCH <- which(sapply(SocTestHinds,is.numeric))
  
  WITCH <- WITCH[!names(WITCH)%in%c("Eigenvector","Eigenvector_Weighted",
                                    "E", "N", 
                                    "LifetimeE", "LifetimeN")]
  
  SocTestHinds %>% mutate_at(WITCH, ~c(scale(.x))) -> SocTestHinds
  
  SubDeerNames <- intersect(DeerNames, SocTestHinds$Name)
  
  SpaceContacts <- as(solve(HRO[SubDeerNames, SubDeerNames]),"dgCMatrix")
  GRMatrix <- as(solve(GRM[SubDeerNames, SubDeerNames]),"dgCMatrix")
  
  SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(SubDeerNames==a)))
  SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(SubDeerNames==a)))
  
  AnnualLocations = cbind(SocTestHinds$E, SocTestHinds$N) #cbind(SocTestHinds$E, SocTestHinds$N)
  
  AnnualMesh <- inla.mesh.2d(loc = AnnualLocations, 
                             boundary = RumBoundary,
                             max.edge = c(10), cutoff = 1)
  
  AnnualA3 <- inla.spde.make.A(AnnualMesh, loc = AnnualLocations) # Making A matrix
  
  AnnualSPDE = inla.spde2.pcmatern(mesh = AnnualMesh, 
                                   prior.range = c(10, 0.5), 
                                   prior.sigma = c(.5, .5)) # Making SPDE
  
  w.index.annual <- inla.spde.make.index('wAnnual', n.spde = AnnualSPDE$n.spde)
  
  LifetimeLocations = cbind(SocTestHinds$LifetimeE, SocTestHinds$LifetimeN)
  
  LifetimeMesh <- inla.mesh.2d(loc = LifetimeLocations, 
                               boundary = RumBoundary,
                               max.edge = c(10), cutoff = 1)
  
  LifetimeA3 <- inla.spde.make.A(LifetimeMesh, loc = LifetimeLocations) # Making A matrix
  
  LifetimeSPDE = inla.spde2.pcmatern(mesh = LifetimeMesh, 
                                     prior.range = c(10, 0.5), 
                                     prior.sigma = c(.5, .5)) # Making SPDE
  
  w.index.lifetime <- inla.spde.make.index('wLifetime', n.spde = LifetimeSPDE$n.spde)
  
  # Making the models ####
  
  Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(Covar, collapse = " + "))), 
                     data = SocTestHinds)
  
  X <- as.data.frame(Xm)
  
  X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("Survived10",
                                                    "ReprodStatusTrue.Yeld"))]) # Model Matrix
  
  EffectList <- list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    IndexSpace = SocTestHinds$IndexSpace,
    IndexPhylo = SocTestHinds$IndexPhylo)
  
  EffectList[c(IID)] <- lapply(c(IID), function(a){
    
    SocTestHinds[,a]
    
  })
  
  EffectList$wAnnual <- w.index.annual
  EffectList$wLifetime <- w.index.lifetime
  
  AList <- as.list(rep(1, length(EffectList) - 2)) %>% 
    append(list(AnnualA3)) %>% 
    append(list(LifetimeA3))
  
  # names(AList) <- names(EffectList)
  
  CentStack <- inla.stack(
    
    data = list(y = SocTestHinds[,SocResps[r]]),  
    
    A = AList,
    
    effects = EffectList) # Leave
  
  i <- 1
  
  for(i in i:length(FormulaList)){
    
    print(names(FormulaList)[i])
    
    IM <-
      
      inla(
        FormulaList[[i]], 
        family = FamilyList[r],
        data = inla.stack.data(CentStack),
        control.compute = list(dic = TRUE),
        control.predictor = list(A = inla.stack.A(CentStack))
      )
    
    CentralityList[[SocResps[r]]][[i]] <- IM
    
    beepr::beep()
    
  }
  
  names(CentralityList[[SocResps[r]]]) <- names(FormulaList)
  
  TestDFList[[SocResps[r]]] <- SocTestHinds
  
  MeshList[[SocResps[r]]] <- SPDEList[[SocResps[r]]] <- list()
  
  MeshList[[SocResps[r]]]$AnnualMesh <- AnnualMesh
  MeshList[[SocResps[r]]]$LifetimeMesh <- LifetimeMesh
  
  SPDEList[[SocResps[r]]]$AnnualSPDE <- AnnualSPDE
  SPDEList[[SocResps[r]]]$LifetimeSPDE <- LifetimeSPDE
  
  
  saveRDS(CentralityList, file = "Output Files/CentralityList.rds")
  
  saveRDS(SPDEList, file = "Output Files/SPDEList.rds")
  saveRDS(MeshList, file = "Output Files/MeshList.rds")
  saveRDS(TestDFList, file = "Output Files/TestDFList.rds")
  
}

# The Temporal part ####

f7 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste0(colnames(X), collapse = " + "),
                       " + ",
                       paste0(paste0("f(", IID[1:2], ", model = 'iid')"), collapse = " + "), 
                       " + ",
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))', 
                       " + ",
                       'f(IndexSpace, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE, param = c(0.5, 0.5))',
                       " + f(wTemporal, model = TemporalSPDE, replicate = wTemporal.repl)"))

r = 1

inla.setOption(num.threads = 8)

for(r in r:length(SocResps)){ # Takes a while I bet
  
  SocialHinds %>% pull(SocResps[r]) %>% 
    between(OutlierLimits[[r]][1], OutlierLimits[[r]][2]) %>%
    which() %>%
    slice(SocialHinds, .) ->
    
    SocTestHinds
  
  SocTestHinds %>% 
    
    select(Cols) %>%
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
  
  SocTestHinds %>% mutate(GroupVar = as.numeric(as.factor(Year))) -> SocTestHinds
  
  NGroup <- nunique(SocTestHinds$GroupVar)
  
  N <- nrow(SocTestHinds); print(N)
  
  WITCH <- which(sapply(SocTestHinds,is.numeric))
  
  WITCH <- WITCH[!names(WITCH)%in%c("Eigenvector","Eigenvector_Weighted", "GroupVar",
                                    "E", "N", 
                                    "LifetimeE", "LifetimeN")]
  
  SocTestHinds %>% mutate_at(WITCH, ~c(scale(.x))) -> SocTestHinds
  
  SubDeerNames <- intersect(DeerNames, SocTestHinds$Name)
  
  SpaceContacts <- as(solve(HRO[SubDeerNames, SubDeerNames]),"dgCMatrix")
  GRMatrix <- as(solve(GRM[SubDeerNames, SubDeerNames]),"dgCMatrix")
  
  SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(SubDeerNames==a)))
  SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(SubDeerNames==a)))
  
  TemporalLocations = cbind(SocTestHinds$E, SocTestHinds$N) #cbind(SocTestHinds$E, SocTestHinds$N)
  
  TemporalMesh <- inla.mesh.2d(loc = TemporalLocations, 
                               boundary = RumBoundary,
                               max.edge = c(10), cutoff = 1)
  
  TemporalA3 <- inla.spde.make.A(TemporalMesh, 
                                 loc = TemporalLocations,
                                 repl = SocTestHinds$GroupVar,
                                 n.repl = NGroup) # Making A matrix
  
  TemporalSPDE = inla.spde2.pcmatern(mesh = TemporalMesh, 
                                     prior.range = c(10, 0.5), 
                                     prior.sigma = c(.5, .5)) # Making SPDE
  
  w.index.temporal <- inla.spde.make.index('wTemporal', 
                                           n.spde = TemporalSPDE$n.spde,
                                           n.repl = NGroup)
  
  # Making the models ####
  
  Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(Covar, collapse = " + "))), 
                     data = SocTestHinds)
  
  X <- as.data.frame(Xm)
  
  X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("Survived10",
                                                    "ReprodStatusTrue.Yeld"))]) # Model Matrix
  
  print(SocResps[r])
  
  EffectList <- list(
    Intercept = rep(1, N), # Leave
    X = X, # Leave
    IndexSpace = SocTestHinds$IndexSpace,
    IndexPhylo = SocTestHinds$IndexPhylo)
  
  EffectList[c(IID)] <- lapply(c(IID), function(a){
    
    SocTestHinds[,a]
    
  })
  
  EffectList$wTemporal <- w.index.temporal
  
  AList <- as.list(rep(1, length(EffectList) - 1)) %>% 
    append(list(TemporalA3))
  
  # names(AList) <- names(EffectList)
  
  CentStack <- inla.stack(
    
    data = list(y = SocTestHinds[,SocResps[r]]),  
    
    A = AList,
    
    effects = EffectList) # Leave
  
  IM <-
    
    inla(
      f7, 
      family = FamilyList[r],
      data = inla.stack.data(CentStack),
      control.compute = list(dic = TRUE),
      control.predictor = list(A = inla.stack.A(CentStack))
    )
  
  CentralityList[[SocResps[r]]]$Temporal <- IM
  
  beepr::beep()
  
  TestDFList[[SocResps[r]]] <- SocTestHinds
  
  MeshList[[SocResps[r]]] <- SPDEList[[SocResps[r]]] <- list()
  
  MeshList[[SocResps[r]]]$Temporal <- TemporalMesh
  SPDEList[[SocResps[r]]]$Temporal <- TemporalSPDE
  
  saveRDS(CentralityList, file = "Output Files/Temporal/CentralityList.rds")
  
  saveRDS(SPDEList, file = "Output Files/Temporal/SPDEList.rds")
  saveRDS(MeshList, file = "Output Files/Temporal/MeshList.rds")
  saveRDS(TestDFList, file = "Output Files/Temporal/TestDFList.rds")
  
}
