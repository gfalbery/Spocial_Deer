# Fitness Extricating Figures ####


# 4_Extricating Figures ####

library(INLA); library(ggregplot); library(tidyverse); library(cowplot); library(gganimate); library(colorspace)

ModelOrder <- c(1:11,14:16,12:13)

# Strength DIC changes ####

INLADICFig(CentralityList$Strength[ModelOrder], ModelNames[ModelOrder], Just = T) +
  labs(x = NULL) + theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  transition_reveal(Model)

anim_save("Strength_DIC.gif")

INLADICFig(CentralityList$Strength[ModelOrder], ModelNames[ModelOrder], Just = T) +
  labs(x = NULL) + theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave("Strength_DIC.jpeg", units = "mm", height = 150, width = 150, dpi = 600)

# Strength variance components ####

INLARepPlot(CentralityList$Strength, ModelNames, Just = 1, DICOrder = T, Residual = F) + 
  #ggtitle(names(CentralityList)[a]) +
  scale_fill_brewer(palette = AlberPalettes[[1]]) +
  coord_fixed(ratio = 15) +
  ggsave("Strength_Var.jpeg", units = "mm", height = 150, width = 150, dpi = 600)

# DIC Output ####

c(1:4) %>% lapply(function(a){
  
  INLADICFig(CentralityList[[a]], ModelNames, Just = T) + ggtitle(names(CentralityList)[a]) +
    labs(x = NULL) + theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
}) %>% arrange_ggplot2(nrow = 2)

# Fields

c(1:4) %>% lapply(function(a){
  
  ggField(CentralityList[[a]][[2]], 
          SocMesh, 
          Fill = "Continuous") + ggtitle(names(CentralityList)[a]) + theme_cowplot() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_continuous_sequential(AlberPalettes[[2]])
  
}) %>% arrange_ggplot2(nrow = 4)

# Variance Component Extraction ####

c(1:4) %>% lapply(function(a){
  
  INLARepPlot(CentralityList[[a]], ModelNames, Just = 1, DICOrder = T, Family = FamilyList[[a]], Residual = F) + 
    ggtitle(names(CentralityList)[a]) +
    scale_fill_brewer(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(nrow = 2)

# Checking no perfect fit ####

a = 1

CentralityList[[a]] %>% lapply(function(b) qplot(SocTestHinds[,Resps[a]], 
                                                 b$summary.fitted.values$mean[1:nrow(SocTestHinds)],
                                                 geom = c("point", "smooth")) + labs(x = "Observed", y = "Fitted")
) %>% 
  arrange_ggplot2

a = 4

FitnessList[[a]] %>% lapply(function(b) qplot(FitdfList[[a]][,FitResps[a]], 
                                              b$summary.fitted.values$mean[1:nrow(FitdfList[[a]])],
                                              geom = c("point", "smooth")) + labs(x = "Observed", y = "Fitted")
) %>% 
  arrange_ggplot2

# Fitness DIC Plots ####

c(1:4) %>% lapply(function(a){
  
  INLADICFig(FitnessList[[a]], ModelNames, Just = 1) + theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(names(FitnessList)[a])
  
}) %>% arrange_ggplot2(nrow = 2, ncol = 2)



# Fitness Extricating #####

FitnessFamilyList <- c("gaussian", "gaussian", "binomial", "binomial")

c(1,2,3,4) %>% lapply(function(a){
  
  INLARepPlot(FitnessList[[a]][ModelOrder], ModelNames[ModelOrder], Just = 1, DICOrder = T, Residual = F, Family = FitnessFamilyList[a]) + 
    ggtitle(names(FitnessList)[a]) +
    scale_fill_brewer(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(nrow = 2, ncol = 2)


# Figures for Presentations ####

names(CentralityList[[2]])[c(2,6,7,13)] %>% 
  lapply(function(a) ggField(CentralityList[[2]][[a]], SocMesh) + ggtitle(a) + scale_fill_discrete_sequential(palette = AlberPalettes[[2]])) %>% 
  arrange_ggplot2(nrow = 1)


ggplot(Cordf, aes(Space, Social)) + geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + geom_smooth(colour = "black") +
  coord_fixed() + THEME + ggpubr::stat_cor() + labs(x = "Home Range Overlap", y = "Association Index") +
  ggsave(file = "~/SocSpat.tiff", units = "mm", width = 100, height = 100, dpi = 600)

ggplot(Cordf, aes(Space, GRM)) + geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + geom_smooth(colour = "black") +
  coord_fixed() + THEME + ggpubr::stat_cor() + labs(x = "Home Range Overlap", y = "Genetic Relatedness") +
  ggsave(file = "~/SpaceGRM.tiff", units = "mm", width = 100, height = 100, dpi = 600)

ggplot(Cordf, aes(GRM, Social)) + geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + geom_smooth(colour = "black") +
  coord_fixed() + THEME + ggpubr::stat_cor() + labs(x = "Genetic Relatedness", y = "Association Index") +
  ggsave(file = "~/GRMSocial.tiff", units = "mm", width = 100, height = 100, dpi = 600)

Efxplot(CentralityList[[2]], ModelNames = ModelNames) + ggsave("ModelOutputs.jpeg", units = "mm", height = 150, width = 200)


lapply(1:4, function(a){
  
  ggField(FitnessList[[a]]$SPDE, Mesh = FitMeshList[[a]], Fill = "continuous") + 
    scale_fill_continuous_sequential(palette = AlberPalettes[[2]]) +
    ggtitle(names(FitnessList)[a])
  
}) %>% arrange_ggplot2(nrow = 1)





