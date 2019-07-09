
# Analysing social network ####

library(INLA); library(ggregplot); library(tidyverse)

FitCovar <- c("Age", "Year", "PopN", "ReprodStatus")

FitResps <- c("cBirthDate.Calf","BirthWt.Calf","Pregnant", "Survived1")

Cols <- c("Name", #FitResps[3],
          "Year", "PopN", FitCovar,
          "Age2", "E", "N")

ModelNames <- c("Base", "SPDE",
                "SpaceMat", "SocMat", "GRMat",
                "SpaceMat:w","SocMat:w","GRMat:w",
                "SpaceMat:SocialMat","SpaceMat:GRMat","SocialMat:GRMat",
                "MatMatMat","MatMatMatw",
                "SpaceMat:SocialMat:w","SpaceMat:GRMat:w","SocialMat:GRMat:w")

FamilyList <- c("gaussian", "gaussian", "binomial","binomial")

FitdfList <- FitMeshList <- list()
FitnessList <- list()

r = 1

for(r in r:length(FitResps)){ # Takes a while I bet

  print(FitResps[r])

  TestHinds <- SocialHinds %>% filter(Shot == 0) %>% 
    slice(-NARows(SocialHinds,c(Cols,FitResps[r]))) %>%
    select(Cols,FitResps[r]) %>%
    mutate(
      #Survived1 = as.numeric(Survived1)-1,
      fYear = as.factor(Year)
    ) %>% filter(Name%in%DeerNames)#,
  #cBirthDate.Calf>-50)# removing NA's and individuals that died halfway through the year

  if(r==1){
    TestHinds <- TestHinds %>% filter(cBirthDate.Calf>-50&cBirthDate.Calf<50)
  }

  if(r==2){
    TestHinds <- TestHinds %>% filter(BirthWt.Calf>2.5&BirthWt.Calf<10)
  }

  if(r==4){
    #TestHinds <- TestHinds %>% filter(!is.na(Survived1)) %>%
    #  mutate(Survived1 = as.numeric(as.character(Survived1)))
    
    TestHinds <- TestHinds %>% #filter(!is.na(Survived1)) %>%
      mutate(Survived1 = as.numeric(as.character(Survived1)))
    
    for(i in (TestHinds %>% filter(is.na(Survived1)))$Name){
      Subdf <- TestHinds[is.na(TestHinds$Survived1)&TestHinds$Name==i,]
      
    }
    
  }

  FitDeerNames <- intersect(DeerNames, TestHinds$Name)

  WITCH <- which(sapply(TestHinds,is.numeric))
  WITCH <- WITCH[!names(WITCH)%in%c("Eigenvector","Eigenvector2","Survived1","Pregnant","E","N")]
  TestHinds[,WITCH] <- apply(TestHinds[,WITCH], 2, scale)

  HostLocations = cbind(TestHinds$E, TestHinds$N)
  FitMeshList[[r]] <- Mesh <- inla.mesh.2d(loc = HostLocations, max.edge = c(10, 25), cutoff = 10)
  A3 <- inla.spde.make.A(Mesh, loc = HostLocations) # Making A matrix
  spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
  w.index <- inla.spde.make.index('w', n.spde = spde$n.spde)

  # Making the models ####

  Xm <- model.matrix(as.formula(paste0("~ -1 + ", paste(FitCovar, collapse = " + "))),
                     data = TestHinds)
  N <- nrow(TestHinds)
  X <- as.data.frame(Xm)[, -which(colnames(Xm) %in%c("Survived10", "ReprodStatusTrue.Yeld"))] # Model Matrix

  if(r%in%1:2){
    X <- X %>% select(-starts_with("ReprodStatus"))
  }

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

  # Establishing distance matrices ####
  SpaceContacts <- as(solve(HRO[FitDeerNames, FitDeerNames]),"dgCMatrix")
  SocialContacts <- as(solve(AM[FitDeerNames, FitDeerNames]),"dgCMatrix")
  GRMatrix <- as(solve(GRM[FitDeerNames, FitDeerNames]),"dgCMatrix")

  TestHinds$IndexSpace = unlist(sapply(TestHinds$Name, function(a) which(FitDeerNames==a)))
  TestHinds$IndexSocial = unlist(sapply(TestHinds$Name, function(a) which(FitDeerNames==a)))
  TestHinds$IndexPhylo = unlist(sapply(TestHinds$Name, function(a) which(FitDeerNames==a)))

  FitdfList[[r]] <- TestHinds

  FitnessList[[FitResps[r]]] <- list()

  FitStack <- inla.stack(
    data = list(y = TestHinds[,FitResps[r]]),
    A = list(1, 1, 1, 1, 1, 1, 1, A3), # Vector of Multiplication factors
    effects = list(
      Intercept = rep(1, N), # Leave
      X = X, # Leave
      IndexSpace = TestHinds$IndexSpace,
      IndexSocial = TestHinds$IndexSocial,
      IndexPhylo = TestHinds$IndexPhylo,
      Name = TestHinds$Name,
      fYear = TestHinds$fYear,
      w = w.index)) # Leave

  q = 12

  for(q in q:length(ModelNames)){
    print(ModelNames[q])
    FitnessList[[FitResps[r]]][[ModelNames[[q]]]] <-
      inla(
        FormulaList[[q]],
        family = FamilyList[r],
        data = inla.stack.data(FitStack),
        control.compute = list(dic = TRUE),
        control.predictor = list(A = inla.stack.A(FitStack))
      )
  }
}

save(FitnessList, file = "Output Files/FitnessList.Rdata")
