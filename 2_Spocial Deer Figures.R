
# 4_Extricating Figures ####

{
  
  library(INLA); library(ggregplot); library(tidyverse); library(GGally); library(patchwork)
  library(cowplot); library(gganimate); library(colorspace); library(RColorBrewer); library(MCMCglmm)
  library(ggrepel)
  
  theme_set(theme_cowplot() + 
              theme(strip.background = element_rect(fill = "white")))
  
  AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
  AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
  
  AlberColours[length(AlberColours)+1:2] <- 
    
    RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]
  
  AlberColours <- append(AlberColours, list(Pink = "#FD6396", Blue = "#3C78D8")) %>% unlist
  
  #Palette <- 
  #  get_colors_from_image(
  #    "https://raw.githubusercontent.com/HughSt/mappalettes/master/images/nathan-lindahl-1j18807_ul0-unsplash.jpg",
  #    5)
  
  # Palette2 <- colorRampPalette(Palette)
  
  CentralityList <- readRDS("Output Files/FullCentralityList.rds")
  
  SPDEList <- readRDS("Output Files/FullSPDEList.rds")
  TestDFList <- readRDS("Output Files/FullTestDFList.rds")
  MeshList <- readRDS("Output Files/FullMeshList.rds")
  
  names(SPDEList) <- names(TestDFList) <- names(MeshList) <- Resps
  
  ModelLabels = c("Base", "SPDE", "HRO", "Social", "GRM", "HRO+SPDE", "Social+SPDE","GRM+SPDE","HRO+Social","HRO+GRM","Social+GRM","HRO+Social+SPDE","HRO+GRM+SPDE","Social+GRM+SPDE","HRO+Social+GRM","HRO+Social+GRM+SPDE")
  
  RespLabels = c("Group Size", "Degree", "Strength", 
                 "Mean Strength", "Eigenvector", "Weighted Eigenvector", 
                 "Betweenness", #"Closeness", 
                 "Clustering")
  
  SmallRespLabels = c("Group\nSize", "Degree", "Strength", 
                      "Mean\nStrength", "Eigenvector", "Weighted\nEigenvector", 
                      "Between-\nness", #"Closeness", 
                      "Clustering")
  
  VarNames <- c("Residual", 
                "Age", "Year", "Population", "Observations", 
                "ID", "Reproduction", 
                "Home Range", 
                "GRM", "SPDE")
  
  VarOrder <- c("Residual", 
                "Year", "Population", "Observations", 
                "ID", "Age", "Reproduction", 
                "GRM", 
                "Home Range", "SPDE")
  
}

# Main text Figures ####

# Figure 2: SPDE effect and spatial distribution of the network ####

x = 2016
Censuses2 <- droplevels(Censuses[Censuses$Year == x&Censuses$Season == "Rut"&Censuses$Hind == "Y",])

AssMatrix<-matrix(unlist(lapply(levels(Censuses2$Code),Assoce)),nrow=nlevels(Censuses2$Code))
dimnames(AssMatrix)<-list(levels(Censuses2$Code),levels(Censuses2$Code))
NObs<-diag(as.matrix(AssMatrix))
AssTotals<-matrix(rep(apply(AssMatrix,2,max),length(AssMatrix[1,]))+
                    rep(apply(AssMatrix,1,max),length(AssMatrix[1,])),nrow=nrow(AssMatrix))
UnionAss<-AssTotals-as.matrix(AssMatrix)
Am2<-as.matrix(AssMatrix)/UnionAss
AM <- Am2

NAMES <- intersect(colnames(AM), DeerNames)

AM <- AM[NAMES, NAMES]
HRO <- HRO[NAMES, NAMES]

# Social Network Links
links2<-reshape2::melt(Am2)
links2$identifier<-rownames(links2)
colnames(links2)[1:3] <-c("from","to","weight")
links2<-merge(links2,longlist3[longlist3$Year==2016,c("E","N","Name")],by.x="from",by.y="Name",all.x=TRUE,sort=FALSE)
links2<-merge(links2,longlist3[longlist3$Year==2016,c("E","N","Name")],by.x="to",by.y="Name",all.x=TRUE,sort=FALSE)

linksfrom<-reshape2::melt(links2[links2$weight>0,],value.name = "E",variable.name = "Coord",id.vars=c("to","from","weight","identifier"))

Cutlinks<-dim(linksfrom)[1]/4
Newlinks<-rbind(cbind(linksfrom[1:Cutlinks,],N=linksfrom[(Cutlinks+1):(2*Cutlinks),"E"]),
                cbind(linksfrom[(2*Cutlinks+1):(3*Cutlinks),],N=linksfrom[(3*Cutlinks+1):(4*Cutlinks),"E"])) %>%
  mutate(MeshInclude = MeshInclude(E, N, 1355,1385,7997.5,8050))

list(
  ggplot(Newlinks %>% filter(MeshInclude == 1) , aes(E, N)) +
    geom_polygon(data = boundary %>% mutate(E = Easting, N = Northing), 
                 fill = NA, colour = "black", size = 1) + 
    geom_line(aes(group = identifier, alpha = weight), colour = AlberColours[[3]]) +
    coord_fixed() + scale_alpha_continuous(range = c(0,0.25)) + 
    theme_cowplot() +
    labs(x = "Easting", y = "Northing") +
    theme(legend.position = "none"),
  
  ggField(CentralityList$Strength$`MatMatMatw`, 
          SocMesh, Fill = "Discrete", Boundary = boundary[,c("Easting", "Northing")]) +
    scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) +
    theme_cowplot() + 
    labs(fill = "Strength", x = "Easting", y = "Northing")
) %>% 
  plot_grid(plotlist = ., labels= "AUTO", rel_widths = c(1, 1.35)) %>%
  save_plot(filename = "Untangling Code/Figures/Figure1.jpeg")

# Figure 3: Model output ####

CentralityList %>% 
  map_dfr(~map_df(.x, c("dic", "dic")), .id = "Response") -> 
  
  DICDF

DICDF %>% mutate_at("Response", ~{
  
  .x %>% str_split("_") %>% 
    map_chr(function(a) sapply(a, function(b) substr(b, 1, 1)) %>% 
              unlist %>% paste0(collapse = ""))
  
}) %>% pull(Response) %>% 
  str_replace_all(c("EW" = "WE", "SM" = "MS")) -> DICDF$ResponseLabel

DICDF %>%
  gather("Model", "DIC", -c(Response, ResponseLabel)) %>%
  group_by(Response) %>%
  mutate_at("DIC", ~.x - min(.x)) -> LongDICDF

LongDICDF %>% arrange(desc(DIC), Model) %>%
  mutate_at("Model", ~.x %>% 
              str_replace_all(c("HROSPDE" = "HRO.AnnualSPDE",
                              "HROLifetimeSPDE" = "HRO.LifetimeSPDE",
                              "^SPDE$" = "AnnualSPDE",
                              "^Temporal$" = "SpatiotemporalSPDE"))) %>% 
  mutate_at("Model", ~factor(.x, levels = unique(.x))) -> 
  
  LongDICDF

LongDICDF %>% 
  ungroup %>%
  mutate(Response = factor(Response, levels = names(CentralityList) %>% 
                             rev)) ->
  
  LongDICDF

LongDICDF %>%
  ggplot(aes(as.numeric(as.factor(Model)), DIC)) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_line(aes(colour = Response)) +
  geom_point(aes(colour = Response)) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = NULL, y = "DeltaDIC") +
  scale_x_continuous(limits = c(0.25, 7), breaks = 1:7, 
                     labels = levels(LongDICDF$Model)) +
  ggrepel::geom_label_repel(data = LongDICDF %>% filter(Model == "Base"), 
                            aes(x = 1, hjust = 1,
                                label = ResponseLabel, 
                                colour = Response),
                            xlim = c(0.1, 1), 
                            force = 10) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]],
                                   nmax = 12, 
                                   order = 12:4) ->
  
  Plot1

RelativeR2 <- readRDS("Output Files/RelativeR2DF.rds")

RelativeR2 %>% 
  filter(Response == "Clustering") %>% 
  mutate(Label = Variable %>% 
           str_replace_all(c("IndexPhylo" = "GRM",
                             "ReprodStatus" = "Status",
                             "IndexSpace" = "HRO",
                             "W" = "SPDE"))) %>%
  mutate(Colour = as.factor(Label %in% c("SPDE", "Year", "Status"))) %>%
  arrange(desc(Variable)) %>%
  mutate(Y = c(0, cumsum(R2[1:(n() - 1)]))) %>%
  mutate(Y2 = R2/2) %>%
  mutate(Y3 = Y + Y2) -> 
  
  LabelsRight

Plot2 <- 
  RelativeR2 %>%
  ggplot(aes(as.numeric(Response), R2, fill = Variable)) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.6) +
  geom_col(position = "stack", colour = "black") + 
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = 1:8,
                     labels = RespLabels) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
  geom_label_repel(data = LabelsRight %>% filter(R2>0.0), 
                   direction = "y",
                   #fill = "white",
                   #label.colour = "black",
                   force = 5,
                   colour = c(rep("white", 3), rep("black", 7)),
                   segment.color = "black",
                   segment.alpha = c(0.6),
                   aes(label = Label, 
                       #colour = Variable,
                       x = 8, y = Y3), nudge_y = 0.2,
                   xlim = c(9.0, 9.0) + 1) +
  coord_cartesian(xlim = c(0.5, 10)) +
  theme(legend.position = "none") +
  labs(y = expression(R^{2}), x = NULL)

ModelEffects <- 
  CentralityList[8:1] %>% 
  map("HROSPDE") %>%
  Efxplot(Size = 3,
    # SigAlpha = T, 
    Alpha1 = 1, Alpha2 = 0.4,
    ModelNames = names(CentralityList) %>% 
            str_replace_all(c("GroupSize" = "Group size",
                              "Eigenvector_Weighted" = "Weighted Eigenvector",
                              "Strength_Mean" = "Mean strength")) %>% 
      rev, 
          Intercept = F, 
          VarNames = rev(c("Intercept", "Age", "Summer_Yeld", "Winter_Yeld", "Milk", 
                           "Year", "Population", "Observations"))) +
  theme(legend.position = c(0.6, 0.75), 
        #legend.box.margin = c(1),
        legend.background = element_rect(fill = "white", 
                                         colour = "white", size = 8-1),
        legend.box.background = element_rect(colour = "white",
                                             size = 10-1)) +
  labs(colour = NULL) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]],
                                   nmax = 12, 
                                   order = 12:4) +
  guides(colour = guide_legend(reverse = T)) +
  geom_point(colour = "black", aes(group = Model), 
             #alpha = graph$Alpha,
             position = position_dodge(w = 0.5), size = 4) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper, group = Model),
                width = 0.1, 
                #alpha = graph$Alpha,
                position = position_dodge(w = 0.5),
                colour = "black") +
  geom_point(position = position_dodge(w = 0.5), size = 3)

((Plot1|Plot2)/ModelEffects + plot_annotation(tag_levels = "A")) + 
  plot_layout(heights = c(1, 1.5)) +
  ggsave("Untangling Code/Figures/Figure3.jpeg", units = "mm", 
         width = 250, height = 280, dpi = 600)


# Getting errors ####

CentralityList %>% map(~INLARep(.x,
                                SPDEModel = DeerSPDE,
                                Return = "Summary",
                                Draw = T, NDraw = 1000
)) ->
  
  VarianceEstimates

CentralityList %>% map(~INLARep(.x,
                                SPDEModel = DeerSPDE,
                                Return = "Summary",
                                SPDEComponent = "Tau",
                                Draw = T, NDraw = 1000
)) ->
  
  VarianceEstimatesTau

VarianceEstimates %>% bind_rows(.id = "Model") %>% 
  bind_cols(VarianceEstimatesTau %>% bind_rows(.id = "Model") %>% 
              rename_all(~paste0(.x, ".Tau"))) -> CompareVariances

CompareVariances %>%
  ggplot(aes(Mean, Mean.Tau)) + geom_point(aes(colour = Component))

CentralityList %>% map(c("summary.hyperpar", "mean")) %>% map_dbl(6)

CompareVariances %>% mutate(Ratio = (Mean - Mean.Tau)/Mean.Tau) %>%
  filter(Component == "SPDE") -> SPDECompareVariances

CentralityList %>% map(c("summary.hyperpar", "mean")) %>% map_dbl(6) ->
  
  SPDECompareVariances$Range

SPDECompareVariances %>% ggplot(aes(Range, Ratio)) + geom_point()

CentralityList %>% map(~INLARep(.x,
                                SPDEModel = DeerSPDE
)) ->
  
  VarianceEstimatesMean

CentralityList %>% map(~INLARep(.x,
                                SPDEModel = DeerSPDE,
                                Draw = T, NDraw = 1000
)) ->
  
  VarianceEstimatesRaw

VarianceEstimatesRaw[[1]] %>% 
  gather("Component", "Value") %>%
  ggplot(aes(Component, Value)) + 
  ggforce::geom_sina(scale = "width") +
  geom_point(data = VarianceEstimatesMean[[1]] %>% 
               mutate(Component = rownames(VarianceEstimatesMean[[1]])), aes(y = Mean),
             size = 5, colour = AlberColours[[2]]) +
  geom_point(data = VarianceEstimates[[1]], aes(y = Mean),
             size = 5, colour = AlberColours[[1]])

VarianceEstimatesRaw %>% bind_rows(.id = "Model") %>% 
  mutate(Ratio = SPDE/IndexSpace) %>% 
  ggplot(aes(Model, Ratio)) + 
  geom_boxplot() +
  ggforce::geom_sina()

VarianceEstimatesRaw %>% bind_rows(.id = "Model") %>% 
  mutate(Ratio = SPDE/(SPDE + IndexSpace)) %>% 
  mutate(Model = fct_reorder(Model, Ratio)) %>%
  ggplot(aes(Model, Ratio)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggforce::geom_sina()

VarianceEstimates[[6]] %>% #
  
  ggplot(aes(Component, Mean)) + geom_point() +
  
  geom_errorbar(aes(ymin = Lower, ymax = Upper))


VarianceEstimates %>% lapply(function(a) a %>% gather %>% ggplot(aes(key, value)) + 
                               ggforce::geom_sina(scale = "width")
) %>% arrange_ggplot2

# Figure 4: Spatial fields ####

jpeg(filename = "Untangling Code/Figures/Figure4.jpeg", 
     units = "mm", height = 300, width = 150, res = 600)

c(1:length(CentralityList)) %>% lapply(function(a){
  
  ggField(CentralityList[[a]]$HROSPDE, 
          MeshList[[a]]$AnnualMesh, 
          #Fill = "Continuous", 
          Boundary = boundary[,c("Easting", "Northing")]) + 
    theme_cowplot() +
    labs(fill = SmallRespLabels[a]) +
    theme(legend.title = element_text(size = 9)) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) +
    labs(x = NULL, y = NULL) +
    scale_fill_discrete_sequential(palette = AlberPalettes[[3]])
}) %>%
  plot_grid(plotlist = ., nrow = 2, labels = "AUTO") 

dev.off()

jpeg(filename = "Untangling Code/Figures/Figure4B.jpeg", 
     units = "mm", height = 150, width = 300, res = 600)

# Supplementary Information Figures ####

# Figure SI1 ####

SocialHinds %>% select(Name, LifetimeE, LifetimeN) %>%
  unique %>% arrange(Name) %>% 
  as.data.frame %>% na.omit-> DistDF

vegan::vegdist(DistDF[,2:3]) %>% as.matrix(ncol = nrow(DistDF)) -> 
  
  DistMatrix

dimnames(DistMatrix) <- list(DistDF[,1], DistDF[,1])

SubDeerNames <- intersect(rownames(DistMatrix), DeerNames)

DistMatrix <- DistMatrix[SubDeerNames, SubDeerNames]
HRO <- HRO[SubDeerNames, SubDeerNames]

mantel(DistMatrix, HRO)
mantel(DistMatrix, GRM[SubDeerNames, SubDeerNames])
mantel(GRM[SubDeerNames, SubDeerNames], HRO)
HRO %>% vegan::mantel()
GRM %>% vegan::mantel()

DistMatrix %>% 
  reshape2::melt() %>%
  as.data.frame() %>% rename(Distance = value) %>%
  slice(which(lower.tri(DistMatrix))) %>%
  left_join(HRO %>% 
              reshape2::melt() %>%
              as.data.frame() %>% rename(HRO = value) %>%
              slice(which(lower.tri(HRO)))) %>%
  left_join(GRM %>% 
              reshape2::melt() %>%
              as.data.frame() %>% rename(GRM = value) %>%
              slice(which(lower.tri(GRM)))) -> CompDF

CompDF %>% ggplot(aes(Distance, HRO)) + 
  geom_point(alpha = 0.01, colour = AlberColours[[3]]) +
  geom_smooth(colour = "black") +
  ggpubr::stat_cor(label.y = 1.1) ->
  
  Plot1

CompDF %>% ggplot(aes(Distance, GRM)) +
  geom_point(alpha = 0.01, colour = AlberColours[[3]]) +
  geom_smooth(colour = "black") +
  ggpubr::stat_cor(label.y = 0.7) ->
  
  Plot2

CompDF %>% ggplot(aes(HRO, GRM)) + 
  geom_point(alpha = 0.01, colour = AlberColours[[3]]) +
  geom_smooth(colour = "black") +
  ggpubr::stat_cor(label.y = 0.7) ->
  
  Plot3

(Plot1|Plot2|Plot3) + 
  ggsave("Untangling Code/Figures/SIFigure1.jpeg",
         units = "mm", width = 250, height = 100,
         dpi = 600)

# Figure SI2 ####

jpeg("Untangling Code/SIFigures/FigureSI2.jpeg", units = "mm", width = 250, height = 250, res = 300)

ggpairs(SocialHinds[,SocResps], 
        lower = list(continuous = wrap("smooth", 
                                       colour = AlberColours[2], 
                                       alpha = 0.3, 
                                       method = "loess")), 
        columnLabels = RespLabels,
        upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 8))

dev.off()

# Figure SI3: INLA Ranges ####

ggsave("Untangling Code/SIFigures/FigureSI7.jpeg", units = "mm", height = 300, width = 250, res = 600)

MeshList <- list()
MeshList[1:8] <- list(SocMesh)

INLARange(ModelList = CentralityList, 
          MeshList = MeshList, 
          MaxRange = 50, 
          ModelNames = RespLabels,
          Return = "Figure") + 
  scale_x_continuous(breaks = c(0:5*10), labels = c(0:5)) +
  labs(x = "Distance (km)") + 
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]]) +
  ggsave("Untangling Code/Figures/FigureSI7.jpeg", units = "mm", height = 200, width = 250, dpi = 600)


# Figure SI4: HRO R2 contributions ####

HRORelativeR2 <- readRDS("Output Files/HRORelativeR2DF.rds")

HRORelativeR2 %>% 
  filter(Response == "Clustering") %>% 
  mutate(Label = Variable %>% 
           str_replace_all(c("IndexPhylo" = "GRM",
                             "ReprodStatus" = "Status",
                             "IndexSpace" = "HRO",
                             "W" = "SPDE"))) %>%
  mutate(Colour = as.factor(Label %in% c("SPDE", "Year", "Status"))) %>%
  arrange(desc(Variable)) %>%
  mutate(Y = c(0, cumsum(R2[1:(n() - 1)]))) %>%
  mutate(Y2 = R2/2) %>%
  mutate(Y3 = Y + Y2) -> 
  
  HROLabelsRight

HROPlot <- 
  HRORelativeR2 %>%
  ggplot(aes(as.numeric(Response), R2, fill = Variable)) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.6) +
  geom_col(position = "stack", colour = "black") + 
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = 1:8,
                     labels = RespLabels) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
  scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
  geom_label_repel(data = HROLabelsRight %>% filter(R2>0.0), 
                   direction = "y",
                   #fill = "white",
                   #label.colour = "black",
                   force = 5,
                   colour = c(rep("white", 2), rep("black", 7)),
                   segment.color = "black",
                   segment.alpha = c(0.6),
                   aes(label = Label, 
                       #colour = Variable,
                       x = 8, y = Y3), nudge_y = 0.2,
                   xlim = c(9.0, 9.0) + 0.5 ) +
  coord_cartesian(xlim = c(0.5, 10)) +
  theme(legend.position = "none") +
  labs(y = expression(R^{2}), x = "Response")

HROPlot + ggsave(file = "SIFigure4.jpeg", units = "mm", width = 200, height = 200, dpi = 600)

# Figures SI5-12: Spatiotemporal fields ####

a <- 1

1:length(CentralityList) %>% map(function(a){
  
  ggField(Model = CentralityList[[a]]$Temporal, 
          Mesh = MeshList[[a]]$Temporal,
          Groups = nunique(TestDFList[[a]][[2]]$GroupVar),
          #Fill = "Continuous", 
          Boundary = boundary[,c("Easting", "Northing")]) + 
    labs(fill = SmallRespLabels[a]) +
    theme(legend.title = element_text(size = 9)) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) +
    labs(x = NULL, y = NULL) +
    facet_wrap(~Group, ncol = 8) +
    theme(legend.position = "top") +
    guides(fill = guide_legend(reverse = F,
                               direction = "horizontal",
                               title.position = "left",
                               title.vjust = 0.25, title.hjust = 1,
                               label.position = "top",
                               label.hjust = 0.5,
                               label.vjust = 1.5,
                               label.theme = element_text(angle = 0), nrow = 1)) +
    scale_fill_discrete_sequential(palette = AlberPalettes[[3]]) +
    ggsave(file = paste0("Untangling Code/Figures/SIFigure", a + 4, ".jpeg"),
           height = 9, width = 7)
  
})

# Little figures for schematic ####

Censuses %>% 
  ggplot(aes(Easting, Northing)) + 
  geom_point(alpha = 0.01) + theme_void() + coord_fixed() -> ObservationPlot

SocTestHinds %>% 
  ggplot(aes(E, N)) + 
  geom_point(alpha = 0.5) + theme_void() + coord_fixed() -> CentroidPlot

ObservationPlot +
  ggsave("Untangling Code/SIFigures/Observations.jpeg", units = "mm", width = 120, height = 150)

CentroidPlot +
  ggsave("Untangling Code/SIFigures/Centroids.jpeg", units = "mm", width = 120, height = 150)



