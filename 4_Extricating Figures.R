
# 4_Extricating Figures ####

library(INLA); library(ggregplot); library(tidyverse); 
library(cowplot); library(gganimate); library(colorspace); library(RColorBrewer); library(GGally)

ModelLabels = c("Base", "SPDE", "HRO", "Social", "GRM", "HRO+SPDE", "Social+SPDE","GRM+SPDE","HRO+Social","HRO+GRM","Social+GRM","HRO+Social+GRM","HRO+Social+GRM+SPDE","HRO+Social+SPDE","HRO+GRM+SPDE","Social+GRM+SPDE")
RespLabels = c("Group Size", "Degree", "Strength", "Mean Strength")

ModelOrder <- c(1:11,14:16,12:13)

SPDEModels <- c(2,6:8,13:16)

# Checking no perfect fit ####

a = 1

CentralityList[[a]] %>% lapply(function(b) qplot(SocTestHinds[,Resps[a]], 
                                                 b$summary.fitted.values$mean[1:nrow(SocTestHinds)],
                                                 geom = c("point", "smooth")) + labs(x = "Observed", y = "Fitted")
) %>% 
  arrange_ggplot2

# Main text Figures ####

# Figure 1: SPDE effect ####

ggField(CentralityList$GroupSize$`MatMatMatw`, SocMesh, Fill = "Continuous") +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]]) +
  theme_cowplot() + 
  labs(fill = "Group Size") +
  ggsave("Figure1.jpeg", units = "mm", height = 100, width = 100, dpi = 600)

# Figure 2: Variance Component Extraction ####

INLARepPlot(CentralityList[[a]], ModelNames, Just = 1, DICOrder = T, CutOff = 10, Family = "gaussian", Residual = F) + 
  ggtitle(RespLabels[a]) +
  scale_fill_brewer(palette = AlberPalettes[[1]], 
                    labels = c("Name", "Year", "SPDE", "HRO", "Social", "GRM")) + 
  labs(x = NULL, fill = "Variable") + 
  scale_x_discrete(labels = ModelLabels[rev(order(sapply(CentralityList[[a]],MDIC)))]) +
  ggsave("Figure2.jpeg", units = "mm", height = 120, width = 150, dpi = 600)

jpeg("FullFigure2.jpeg", units = "mm", height = 150, width = 200, res = 600)

c(1:4) %>% lapply(function(a){
  
  INLARepPlot(CentralityList[[a]], ModelNames, Just = 1, DICOrder = T, CutOff = 10, Family = "gaussian", Residual = F) + 
    ggtitle(RespLabels[a]) +
    scale_fill_brewer(palette = AlberPalettes[[1]], 
                      labels = c("Name", "Year", "SPDE", "HRO", "Social", "GRM")) + 
    labs(x = NULL, fill = "Variable") + 
    scale_x_discrete(labels = ModelLabels[rev(order(sapply(CentralityList[[a]],MDIC)))]) +
    theme(axis.text.x = element_text(size = 6))
  
}) %>% arrange_ggplot2(nrow = 2)

dev.off()

# Supplementary Information Figures ####

# Figure SI1 ####

Cor1 = vegan::mantel(HRO[DeerNames,DeerNames],AM[DeerNames,DeerNames])
Cor2 = vegan::mantel(HRO[DeerNames,DeerNames],GRM[DeerNames,DeerNames])
Cor3 = vegan::mantel(GRM[DeerNames,DeerNames],AM[DeerNames,DeerNames])

list(
  ggplot(Cordf, aes(Space, Social)) + 
    geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + 
    geom_text(data = Cordf[1,], aes(label = paste0("R=",round(Cor1$statistic,2)), x = 0.15, y = 0.7)) +
    geom_smooth(colour = "black") +
    labs(x = "Home Range Overlap", y = "Association Index"),
  
  ggplot(Cordf, aes(Space, GRM)) + 
    geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + 
    geom_text(data = Cordf[1,], aes(label = paste0("R=",round(Cor2$statistic,2)), x = 0.15, y = 0.7)) +
    geom_smooth(colour = "black") +
    labs(x = "Home Range Overlap", y = "Genetic Relatedness"),
  
  ggplot(Cordf, aes(GRM, Social)) + 
    geom_point(alpha = 0.1, colour = brewer.pal(9,AlberPalettes[[3]])[7]) + 
    geom_text(data = Cordf[1,], aes(label = paste0("R=",round(Cor1$statistic,2)), x = 0, y = 0.7)) +
    geom_smooth(colour = "black") +
    labs(x = "Genetic Relatedness", y = "Association Index")
) %>% 
  plot_grid(plotlist = ., labels = "AUTO", ncol = 3) %>%
  save_plot(., file = "SIFigures/MatrixCors.jpeg", base_width = 4)


# Figure SI2 ####

jpeg("SIFigures/ResponseCorrelations.jpeg", units = "mm", width = 200, height = 200, res = 300)

ggpairs(SocialHinds[,SocResps], 
        lower = list(continuous = wrap("smooth", 
                                       colour = AlberColours[2], 
                                       alpha = 0.3, 
                                       method = "loess")), 
        columnLabels = RespLabels,
        upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "white"))

dev.off()

# Figure SI3: INLA DIC Values ####

jpeg("FullDICFig.jpeg", units = "mm", height = 150, width = 200, res = 600)

c(1:4) %>% lapply(function(a){
  
  INLADICFig(CentralityList[[a]][ModelOrder], ModelLabels[ModelOrder], Just = T, CutOff = 10) + 
    ggtitle(RespLabels[a]) +
    labs(x = NULL) + theme_cowplot() + 
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
  
}) %>% arrange_ggplot2(nrow = 2)

dev.off()

# Figure SI4: Model Ouptuts ####

jpeg("FullModelOutputs.jpeg", units = "mm", height = 250, width = 200, res = 600)

c(1:4) %>% lapply(function(a){
  
  Efxplot(CentralityList[[a]][ModelOrder], ModelNames = ModelLabels[ModelOrder])  + 
    ggtitle(RespLabels[a]) +
    guides(colour = guide_colourbar(cols = 2))
  
}) %>% arrange_ggplot2(nrow = 4)

dev.off()

# Figure SI5: Spatial Fields ####

c(1:4) %>% lapply(function(a){
  
  lapply(1:length(SPDEModels), function(b){
    ggField(CentralityList[[a]][[SPDEModels[b]]], SocMesh, Fill = "Continuous") + 
      theme_cowplot() + theme(legend.position = "top") +
      scale_fill_continuous_sequential(palette = AlberPalettes[[2]])
  }) %>% plot_grid(plotlist = ., ncol = length(SPDEModels))
  
}) %>% plot_grid(plotlist = ., nrow = 4) %>% save_plot(base_aspect_ratio = 5)

# Figure SI6: INLA Ranges ####

INLARange(CentralityList[[a]][intersect(ModelOrder,SPDEModels)], Mesh = SocMesh, maxrange = 50, ModelLabels[intersect(ModelOrder,SPDEModels)]) + theme_cowplot() +
  ggsave("SIFigures/FigureSI6.jpeg", units = "mm", width = 150, height = 150, dpi = 300)

jpeg("INLARanges.jpeg", units = "mm", height = 200, width = 250, res = 600)

c(1:4) %>% lapply(function(a){
  INLARange(CentralityList[[a]][intersect(ModelOrder,SPDEModels)], 
            Mesh = SocMesh, maxrange = 50, ModelLabels[intersect(ModelOrder,SPDEModels)]) + 
    theme_cowplot() + ggtitle(RespLabels[a])
}) %>% plot_grid(plotlist = ., nrow = 2)

dev.off()
