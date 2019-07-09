
# Analysing social network ####

library(INLA); library(ggregplot)

covar <- c("Age", "Survived1", "ReprodStatus", "Year", "PopN", "NObs")

SocResps <- Resps <- c("GroupSize","Degree", "Strength","Strength2")
CentCovar <- covar

Cols <- c("Name", SocResps, 
          "Year", "PopN", covar,
          "Age2", "E", "N")

SocTestHinds <- SocialHinds[!NARows(SocialHinds,Cols),] %>% select(Cols) %>%
  mutate(
    Degree = as.vector(scale(sqrt(Degree))),
    Strength = as.vector(scale(sqrt(Strength))),
    GroupSize = as.vector(scale(kader:::cuberoot(GroupSize))),
    NObs = scale(NObs) %>% c,
    fYear = as.factor(Year)
  ) %>% filter(Name%in%DeerNames) %>% na.omit() %>%
  as.data.frame()# removing NA's and individuals that died halfway through the year

WITCH <- which(sapply(SocTestHinds,is.numeric))
WITCH <- WITCH[!names(WITCH)%in%c("Eigenvector","Eigenvector2","E","N")]
SocTestHinds[,WITCH] <- apply(SocTestHinds[,WITCH], 2, scale)

HostLocations = cbind(SocTestHinds$E, SocTestHinds$N)
SocMesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
A3 <- inla.spde.make.A(SocMesh, loc = HostLocations) # Making A matrix
spde = inla.spde2.pcmatern(mesh = SocMesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

# Making the models ####

Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(CentCovar, collapse = " + "))), 
                   data = SocTestHinds)
N <- nrow(SocTestHinds)
X <- as.data.frame(Xm[, -which(colnames(Xm) %in%c("Survived10"))]) # Model Matrix

# Establishing distance matrices ####
SpaceContacts <- as(solve(HRO[DeerNames, DeerNames]),"dgCMatrix")
SocialContacts <- as(solve(AM[DeerNames, DeerNames]),"dgCMatrix")
GRMatrix <- as(solve(GRM[DeerNames, DeerNames]),"dgCMatrix")

SocTestHinds$IndexSpace = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
SocTestHinds$IndexSocial = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))
SocTestHinds$IndexPhylo = unlist(sapply(SocTestHinds$Name, function(a) which(DeerNames==a)))

# Establishing model formulae ####
f1 <- as.formula(paste("y ~ -1 + Intercept + ", 
                       paste(names(X), collapse = " + "),
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f2 <- as.formula(paste0("y ~ -1 + Intercept + ",
                        paste(names(X), collapse = " + "),
                        " + f(w, model = spde)",
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f3 <-  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr=TRUE,param = c(0.5, 0.5))',
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f4 <-  as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                       constr=TRUE,param = c(0.5, 0.5))',
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f5 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr=TRUE,param = c(0.5, 0.5))',
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f6 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 
                       'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       " + f(w, model = spde)",
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f7 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 
                       'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       " + f(w, model = spde)",
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f8 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 
                       'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                       constr = TRUE,param = c(0.5, 0.5))',
                       " + f(w, model = spde)",
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f9 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                       " + ", 
                       'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       " + ", 
                       'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                       " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f10 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f11 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))


f12 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f13 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(w, model = spde)",
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))


f14 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                       constr = TRUE,param = c(0.5, 0.5))',
                        " + f(w, model = spde)",
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f15 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSpace, model="generic0", Cmatrix = SpaceContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(w, model = spde)",
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))

f16 <- as.formula(paste("y ~ -1 + Intercept + ", paste(names(X), collapse = " + "), 
                        " + ", 
                        'f(IndexSocial, model="generic0", Cmatrix = SocialContacts,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + ", 
                        'f(IndexPhylo, model="generic0", Cmatrix = GRMatrix,
                        constr = TRUE,param = c(0.5, 0.5))',
                        " + f(w, model = spde)",
                        " + f(Name, model = 'iid') + f(fYear, model = 'iid')"))


FormulaList <- list(f1, f2, 
                    f3, f4, f5, 
                    f6, f7, f8, 
                    f9, f10, f11, 
                    f12, f13,
                    f14, f15, f16)

ModelNames <- c("Base", 
                "SPDE", "SpaceMat", "SocMat", "GRMat",
                "SpaceMat:w","SocMat:w","GRMat:w",
                "SpaceMat:SocialMat","SpaceMat:GRMat","SocialMat:GRMat",
                "MatMatMat","MatMatMatw",
                "SpaceMat:SocialMat:w","SpaceMat:GRMat:w","SocialMat:GRMat:w")

FamilyList <- rep("gaussian", 4)

CentralityList <- list()

r = 3

for(r in r:length(SocResps)){ # Takes a while I bet
  
  print(SocResps[r])
  
  CentralityList[[SocResps[r]]] <- list()
  
  CentStack <- inla.stack(
    data = list(y = SocTestHinds[,SocResps[r]]),  
    A = list(1, 1, 1, 1, 1, 1, 1, A3), # Vector of Multiplication factors              
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      IndexSpace = SocTestHinds$IndexSpace,
      IndexSocial = SocTestHinds$IndexSocial,
      IndexPhylo = SocTestHinds$IndexPhylo,
      Name = SocTestHinds$Name,
      fYear = SocTestHinds$fYear,
      w = w.index)) # Leave
  
  q = 1
  
  for(q in q:length(ModelNames)){
    print(ModelNames[q])
    CentralityList[[SocResps[r]]][[ModelNames[[q]]]] <-
      inla(
        FormulaList[[q]], 
        family = FamilyList[r],
        data = inla.stack.data(CentStack),
        control.compute = list(dic = TRUE),
        control.predictor = list(A = inla.stack.A(CentStack))
      )
  }
}

save(CentralityList, file = "Output Files/CentralityList.Rdata")

a = 2

INLADICFig(CentralityList[[a]], ModelNames = ModelNames) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Efxplot(CentralityList[[a]], ModelNames = ModelNames)

# Plotting out ####

lapply(CentralityList[[1]][c(6:8,13:16)], function(a) a %>%
         ggField(SocMesh)) %>% arrange_ggplot2(nrow = 3)

lapply(1:length(CentralityList), function(a) INLADICFig(CentralityList[[a]], ModelNames = ModelNames)) %>% 
  arrange_ggplot2(nrow = 3)

lapply(CentralityList, function(a) Efxplot(a, ModelNames = ModelNames)) %>% 
  arrange_ggplot2(ncol = 3)

x=1

lapply(1:length(CentralityList[[x]]), function(a){ 
  qplot(SocTestHinds[,SocResps[x]],
        CentralityList[[x]][[a]]$summary.fitted.values$mean[1:dim(SocTestHinds)[1]]) + 
    ggtitle(ModelNames[a]) + labs(x = paste("Data", SocResps[x]), y = paste("Fitted", SocResps[x]))
}) %>%
  arrange_ggplot2

lapply(CentralityList, function(a) INLARep(a[[7]]))

# Plotting field ####

ggField(CentralityList[[3]][[2]], SocMesh) + 
  geom_path(data = WorldMap, inherit.aes = F, aes(long/50000, lat/50000, group = group)) +
  geom_point(data = SocTestHinds[,c("LongMean", "LatMean")]/50000, aes(LongMean, LatMean), inherit.aes = F) + 
  scale_fill_brewer(palette = AlberPalettes[2])

