
# Predictions etc ####

library(INLA); library(ggregplot); library(tidyverse); library(cowplot); library(ggpointdensity);
library(colorspace)

theme_set(theme_cowplot())

CentralityList <- readRDS("Output Files/FullCentralityList.rds")

SPDEList <- readRDS("Output Files/FullSPDEList.rds")
TestDFList <- readRDS("Output Files/FullTestDFList.rds")
MeshList <- readRDS("Output Files/FullMeshList.rds")

names(SPDEList) <- names(TestDFList) <- names(MeshList) <- Resps[1:length(SPDEList)]

FitList <- list()

Covar <- c("Age", "ReprodStatus", "Year", "PopN", "NObs")

for(cl in Resps){
  
  print(cl)
  
  TestDF <- TestDFList[[cl]][[1]]
  
  Fits <- INLAFit(Model = CentralityList[[cl]]$HROSPDE, 
                  TestDF = TestDF,
                  FixedCovar = Covar, 
                  Locations = TestDF[,c("E", "N")],
                  Mesh = MeshList[[cl]]$AnnualMesh, 
                  SPDE = SPDEList[[cl]]$AnnualSPDE)
  
  TestDF$Fits <- Fits
  
  FitsMatrix <- INLAFit(Model = CentralityList[[cl]]$HROSPDE, 
                        TestDF = TestDF,
                        FixedCovar = Covar, 
                        Locations = TestDF[,c("E", "N")],
                        Mesh = MeshList[[cl]]$AnnualMesh, 
                        SPDE = SPDEList[[cl]]$AnnualSPDE,
                        Return = "Matrix") %>% 
    
    as.data.frame
  
  names(FitsMatrix)
  
  AlterVar <- c("Intercept", "Age", "ReprodStatus", "Year", "PopN", "NObs", 
                "Name", "fYear", 
                "IndexSpace", "IndexPhylo", "W")
  
  #AlterVar <- AlterVar[-length(AlterVar)]
  
  PlotList <- list()
  
  for(r in AlterVar){
    
    print(r)
    
    FitsMatrix %>% select(starts_with(r)) %>% unlist -> Values
    
    Values[!Values == 0] -> Values
    
    Mean <- mean(Values, na.rm = T)
    
    FitsMatrix %>% mutate_at(vars(starts_with(r)), ~ifelse(.x == 0, 0, Mean)) %>% rowSums -> HoldFits
    
    TestDF[[paste0("Fits.", r)]] <- HoldFits
    
    PlotList[[r]] <- ggplot(TestDF, aes(Fits, HoldFits)) + 
      geom_abline() +
      geom_pointdensity() + 
      # scale_color_continuous_sequential(palette = AlberPalettes[[4]]) +
      scale_color_continuous_sequential(palette = "Terrain") +
      ggpubr::stat_cor() + 
      theme(legend.position = "none") +
      ggtitle(r)
    
  }
  
  TestDF$Observed <- TestDF[,cl]
  
  FitList[[cl]] <- TestDF
  
}

Resps[1:length(FitList)] %>% 
  
  map_dfr(~FitList[[.x]] %>% 
            dplyr::select(starts_with("Fits"), Resps, Observed) %>% 
            gather("Key", "Value", -c(Fits, Resps, Observed)) %>% 
            group_by(Key) %>% 
            summarise(Cor = cor(Value, Observed, use = "complete.obs")) %>% 
            t %>% as.data.frame %>% GregHeader, 
          
          .id = "Response") %>%
  
  mutate(Response = Resps[1:length(FitList)]) -> R2DF

R2DF %>% 
  gather("Variable", "R2", -Response) %>% 
  mutate_at("Variable", ~str_remove(.x, "Fits.")) %>%
  group_by(Response) %>%
  mutate_at("R2", as.numeric) %>%
  mutate_at("R2", ~ ((max(.x) - .x)/sum(max(.x)  - .x))*max(.x)) %>%
  ungroup %>%
  mutate_at("Response", ~factor(.x, levels = Resps)) -> 
  
  RelativeR2

saveRDS(RelativeR2, file = "Output Files/RelativeR2DF.rds")

RelativeR2 %>% 
  ungroup %>%
  mutate_at("Response", ~factor(.x, levels = Resps)) %>%
  ggplot(aes(Response, R2, fill = Variable)) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.6) +
  geom_col(position = "stack", colour = "black") + 
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]])

CentralityList %>% map("HROSPDE") %>% Efxplot(ModelNames = Resps)

CentralityList %>% map("Temporal") %>% Efxplot(ModelNames = Resps)

# Just HRO Model ####

FitList <- list()

Covar <- c("Age", "ReprodStatus", "Year", "PopN", "NObs")

for(cl in Resps){
  
  print(cl)
  
  TestDF <- TestDFList[[cl]][[1]]
  
  Fits <- INLAFit(Model = CentralityList[[cl]]$HRO, 
                  TestDF = TestDF,
                  FixedCovar = Covar, 
                  Locations = TestDF[,c("E", "N")],
                  Mesh = MeshList[[cl]]$AnnualMesh, 
                  SPDE = SPDEList[[cl]]$AnnualSPDE)
  
  TestDF$Fits <- Fits
  
  FitsMatrix <- INLAFit(Model = CentralityList[[cl]]$HRO, 
                        TestDF = TestDF,
                        FixedCovar = Covar, 
                        Locations = TestDF[,c("E", "N")],
                        Mesh = MeshList[[cl]]$AnnualMesh, 
                        SPDE = SPDEList[[cl]]$AnnualSPDE,
                        Return = "Matrix") %>% 
    
    as.data.frame
  
  names(FitsMatrix)
  
  AlterVar <- c("Intercept", "Age", "ReprodStatus", "Year", "PopN", "NObs", 
                "Name", "fYear", 
                "IndexSpace", "IndexPhylo")
  
  #AlterVar <- AlterVar[-length(AlterVar)]
  
  PlotList <- list()
  
  for(r in AlterVar){
    
    print(r)
    
    FitsMatrix %>% select(starts_with(r)) %>% unlist -> Values
    
    Values[!Values == 0] -> Values
    
    Mean <- mean(Values, na.rm = T)
    
    FitsMatrix %>% mutate_at(vars(starts_with(r)), ~ifelse(.x == 0, 0, Mean)) %>% rowSums -> HoldFits
    
    TestDF[[paste0("Fits.", r)]] <- HoldFits
    
    PlotList[[r]] <- ggplot(TestDF, aes(Fits, HoldFits)) + 
      geom_abline() +
      geom_pointdensity() + 
      # scale_color_continuous_sequential(palette = AlberPalettes[[4]]) +
      scale_color_continuous_sequential(palette = "Terrain") +
      ggpubr::stat_cor() + 
      theme(legend.position = "none") +
      ggtitle(r)
    
  }
  
  TestDF$Observed <- TestDF[,cl]
  
  FitList[[cl]] <- TestDF
  
}

Resps[1:length(FitList)] %>% 
  
  map_dfr(~FitList[[.x]] %>% 
            dplyr::select(starts_with("Fits"), Resps, Observed) %>% 
            gather("Key", "Value", -c(Fits, Resps, Observed)) %>% 
            group_by(Key) %>% 
            summarise(Cor = cor(Value, Observed, use = "complete.obs")) %>% 
            t %>% as.data.frame %>% GregHeader, 
          
          .id = "Response") %>%
  
  mutate(Response = Resps[1:length(FitList)]) -> R2DF

R2DF %>% 
  gather("Variable", "R2", -Response) %>% 
  mutate_at("Variable", ~str_remove(.x, "Fits.")) %>%
  group_by(Response) %>%
  mutate_at("R2", as.numeric) %>%
  mutate_at("R2", ~ ((max(.x) - .x)/sum(max(.x)  - .x))*max(.x)) %>%
  ungroup %>%
  mutate_at("Response", ~factor(.x, levels = Resps)) -> 
  
  RelativeR2

saveRDS(RelativeR2, file = "Output Files/HRORelativeR2DF.rds")

RelativeR2 %>% 
  ungroup %>%
  mutate_at("Response", ~factor(.x, levels = Resps)) %>%
  ggplot(aes(Response, R2, fill = Variable)) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.6) +
  geom_col(position = "stack", colour = "black") + 
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]])

