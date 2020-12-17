
# 3_Untangling Figures ####

{
  
  library(INLA); library(ggregplot); library(tidyverse); library(GGally); library(patchwork)
  library(cowplot); library(gganimate); library(colorspace); library(RColorBrewer); library(MCMCglmm)
  library(ggrepel); library(magrittr)
  
  theme_set(theme_cowplot() + 
              theme(strip.background = element_rect(fill = "white")))
  
  AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
  AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
  
  AlberColours[length(AlberColours)+1:2] <- 
    
    RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]
  
  AlberColours <- append(AlberColours, list(Pink = "#FD6396", Blue = "#3C78D8")) %>% unlist
  
  CentralityList <- readRDS("Untangling Code/Output Files/FullCentralityList.rds")
  
  SPDEList <- CentralityList %>% map(c("Spatial", "SPDE"))
  TestDFList <- CentralityList %>% map(c("Data"))
  MeshList <- CentralityList %>% map(c("Spatial", "Mesh"))
  
  names(SPDEList) <- names(TestDFList) <- names(MeshList) <- Resps
  
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
  map("DIC") %>% 
  map(~.x[2:length(.x)]) %>% 
  map(c(DICTableGet)) %>% 
  bind_rows(.id = "Response") %>% 
  mutate_at("Variable", ~.x %>% str_split("f[()]") %>% 
              map_chr(last) %>% str_split(", ") %>% 
              map_chr(first) %>% str_remove_all("[*]")
  ) -> 
  
  DICDF

CentralityList %>% 
  map("dDIC") %>% 
  map(c(DICTableGet)) %>% 
  bind_rows(.id = "Response") %>% 
  mutate_at("Variable", ~.x %>% str_split("f[()]") %>% 
              map_chr(last) %>% str_split(", ") %>% 
              map_chr(first) %>% str_remove_all("[*]")
  ) -> 
  
  dDICDF

dDICDF %>% filter(Kept == 1) %>% 
  bind_rows(expand.grid(Response = Resps, Round = "0", Variable = "Base", Delta = 0, Kept = 1)) %>% 
  arrange(Response, Round) -> DICDF

CentralityList %>% map_dbl(c("FinalModel", "dic", "dic")) -> FinalDIC

CentralityList %>% map_dbl(c("Spatial", "Model", "dic", "dic")) %>% 
  subtract(FinalDIC) -> SpatialDIC

CentralityList %>% map_dbl(c("Spatial", "SpatiotemporalModel", "dic", "dic")) %>% 
  subtract(CentralityList %>% map_dbl(c("Spatial", "Model", "dic", "dic"))) -> 
  SpatiotemporallDIC

list(SpatialDIC, SpatiotemporallDIC) %>% map(as.data.frame) %>% 
  map(~rownames_to_column(.x, "Response") %>% rename(Delta = ".x[[i]]")) %>% 
  bind_rows(.id = "Variable") %>% 
  mutate(Round = Variable) %>% 
  mutate_at("Variable", ~c("Spatial", "Spatiotemporal")[as.numeric(.x)]) %>% 
  bind_rows(DICDF, .) %>% 
  group_by(Response) %>% 
  mutate(DIC = cumsum(Delta)) %>% 
  mutate(Round = 1:n()) %>% ungroup %>% 
  mutate_at("Response", ~factor(.x, levels = rev(names(CentralityList)))) %>% 
  arrange(Response) -> CombinedDICDF

CombinedDICDF %>% group_by(Variable) %>% 
  summarise_at("Delta", mean) %>% arrange(Delta) %>% 
  pull(Variable) %>% rev -> VarOrder

VarLabels <- c("Base", "GRM", "HRA", "HRO", "Density", "SPDE", "tSPDE")

names(VarLabels) <- VarOrder

CombinedDICDF %>% 
  
  mutate_at("Variable", ~factor(.x, levels = VarOrder)) ->
  CombinedDICDF

CombinedDICDF %>% mutate_at("Response", ~{
  
  .x %>% str_split("_") %>% 
    map_chr(function(a) sapply(a, function(b) substr(b, 1, 1)) %>% 
              unlist %>% paste0(collapse = ""))
  
}) %>% pull(Response) %>% 
  str_replace_all(c("EW" = "WE", "SM" = "MS")) %>% 
  factor(., levels = rev(unique(.))) -> CombinedDICDF$ResponseLabel

CombinedDICDF %>%
  ggplot(aes(as.numeric(as.factor(Variable)), GregCube(Delta))) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_line(aes(colour = Response)) +
  geom_point(colour = "black", size = 2) + 
  geom_point(aes(colour = Response)) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.title.y = element_text(vjust = -10)) + 
  labs(x = NULL, y = "DeltaDIC") +
  scale_x_continuous(limits = c(1, 8), breaks = 1:7, 
                     labels = VarLabels[levels(CombinedDICDF$Variable)]) +
  #scale_colour_discrete_sequential(palette = AlberPalettes[[2]],
  #                                 nmax = 12, 
  #                                 order = 12:4) +
  scale_colour_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(limits = c(-15, 0), breaks = -c(0:3*5), labels = -c(0:3*5)^3) +
  ggrepel::geom_label_repel(data = CombinedDICDF %>% filter(Variable == "Spatiotemporal"), 
                            aes(x = 7, hjust = 0,
                                label = ResponseLabel, 
                                fill = Response),
                            alpha = 1,
                            colour = "black",
                            xlim = c(7.1, 8), 
                            force = 10) ->
  
  Plot1

RelativeR2 <- readRDS("Output Files/RelativeR2DF.rds")

RelativeR2 %>% 
  filter(Response == "Clustering") %>% 
  mutate(Label = Variable %>% 
           str_replace_all(c("IndexPhylo" = "GRM",
                             "ReprodStatus" = "Status",
                             "IndexSpace" = "HRO",
                             "W" = "SPDE",
                             "AnnualDensity" = "Density"))) %>%
  mutate(Colour = as.factor(Label %in% c("SPDE", "Year", "Status"))) %>%
  arrange(desc(Variable)) %>%
  mutate(Y = c(0, cumsum(na.omit(R2[1:(n() - 1)])))) %>%
  mutate(Y2 = R2/2) %>%
  mutate(Y3 = Y + Y2) -> 
  
  LabelsRight

RelativeR2 %>% 
  filter(Response == "GroupSize") %>% 
  mutate(Label = Variable %>% 
           str_replace_all(c("IndexPhylo" = "GRM",
                             "ReprodStatus" = "Status",
                             "IndexSpace" = "HRO",
                             "W" = "SPDE",
                             "AnnualDensity" = "Density"))) %>%
  mutate(Colour = as.factor(Label %in% c("SPDE", "Year", "Status"))) %>%
  arrange(desc(Variable)) %>%
  mutate(Y = c(0, cumsum(na.omit(R2[1:(n() - 1)])))) %>%
  mutate(Y2 = R2/2) %>%
  mutate(Y3 = Y + Y2) -> 
  
  LabelsLeft

RelativeR2 %>% 
  group_by(Response) %>% 
  mutate(Label = Variable %>% 
           str_replace_all(c("IndexPhylo" = "GRM",
                             "ReprodStatus" = "Status",
                             "IndexSpace" = "HRO",
                             "W" = "SPDE",
                             "AnnualDensity" = "Density"))) %>%
  mutate(Colour = as.factor(Label %in% c("SPDE", "Year", "Status"))) %>%
  arrange(desc(Variable)) %>% na.omit %>% 
  mutate(Y = c(0, cumsum(R2[1:(n() - 1)]))) %>%
  mutate(Y2 = R2/2) %>%
  mutate(Y3 = Y + Y2) -> 
  
  LabelsCentre

Plot2 <- 
  
  RelativeR2 %>%
  ggplot(aes(as.numeric(Response), R2, fill = Variable)) + 
  geom_hline(yintercept = 1, lty = 2, alpha = 0.6) +
  geom_col(position = "stack", colour = "black") + 
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = 1:8,
                     labels = c(RespLabels)) +
    scale_fill_discrete_diverging(palette = "Tropic") +
    scale_colour_discrete_diverging(palette = "Tropic") +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) +
    scale_colour_discrete_sequential(palette = AlberPalettes[[1]]) +
  #scale_fill_discrete_qualitative(palette = "warm") +
  #scale_colour_discrete_qualitative(palette = "warm") +
  geom_label_repel(data = LabelsLeft %>% filter(!Variable == "Intercept"), 
                   direction = "y",
                   force = 5,
                   #colour = "black", #c(rep("white", 3), rep("black", 9)),
                   colour = c(rep("white", 3), rep("black", 9)),
                   segment.color = "black",
                   segment.alpha = c(0.6),
                   aes(label = Label, 
                       x = 0.7, y = Y3), nudge_y = 0.2,
                   xlim = c(0, 0.5)-1) +
  coord_cartesian(xlim = c(0, 9)-1) +
  theme(legend.position = "none") +
  labs(y = expression(R^{2}), x = NULL) +
  geom_text(data = LabelsCentre %>% filter(R2>0.05), 
            aes(y = Y3, label = Label %>% substr(1, 2)), colour = "black")

ModelEffects <- 
  CentralityList[8:1] %>% 
  map(c("Spatial", "Model")) %>%
  Efxplot(Size = 3,
          # SigAlpha = T, 
          Alpha1 = 1, Alpha2 = 0.4,
          ModelNames = (names(CentralityList)) %>% 
            str_replace_all(c("GroupSize" = "Group size",
                              "Eigenvector_Weighted" = "Weighted Eigenvector",
                              "Strength_Mean" = "Mean strength")) %>% 
            rev, 
          Intercept = F, 
          VarNames = rev(c("Intercept", "Age", "Summer_Yeld", "Winter_Yeld", "Milk", 
                           "Year", "Population", "Observations", "HRA", "AnnualDensity"))) +
  theme(legend.position = c(0.6, 0.75), 
        legend.background = element_rect(fill = "white", 
                                         colour = "white", size = 8-1),
        legend.box.background = element_rect(colour = "white",
                                             size = 10-1)) +
  labs(colour = NULL) +
  scale_colour_brewer(palette = "Spectral") +
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

# Figure 4: Spatial fields ####

jpeg(filename = "Untangling Code/Figures/Figure4.jpeg", 
     units = "mm", height = 150, width = 300, res = 600)

c(1:length(CentralityList)) %>% lapply(function(a){
  
  ggField(CentralityList[[a]]$Spatial$Model, 
          MeshList[[a]], 
          #Fill = "Continuous", 
          Boundary = Boundary[,c("Easting", "Northing")]) + 
    theme_cowplot() +
    labs(fill = SmallRespLabels[a]) +
    theme(legend.title = element_text(size = 9)) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) +
    labs(x = NULL, y = NULL) +
    scale_fill_discrete_sequential(palette = AlberPalettes[[3]])
}) %>%
  plot_grid(plotlist = ., nrow = 2, ncol = 4, labels = "AUTO") 

dev.off()

# Figure 5: Density Effects ####

SpocialList <- CentralityList

Covar <- c("ReprodStatus", "Age", "AnnualDensity",
           "Year", "PopN", "NObs")

a <- Resps[1]

Figure1aList <- LabelYList <- LabelList <- LabelXList <- list()

for(a in Resps){
  
  print(a)
  
  TestHinds <- TestDFList[[a]]
  SocMesh <- MeshList[[a]]
  
  X = seq(min(TestDFList[[a]]$AnnualDensity), 
          max(TestDFList[[a]]$AnnualDensity), 
          length.out = nrow(TestHinds)) %>% c
  
  B = summary(SpocialList[[a]]$Spatial$Model)$fixed["AnnualDensity","mean"]
  A = summary(SpocialList[[a]]$Spatial$Model)$fixed[,"mean"]%*%t(SpocialList[[a]]$Spatial$Model$model.matrix) %>% mean
  Y = A + B*X
  
  SpocialList[[a]]$Spatial$Model %>% 
    INLAFit(TestDFList[[a]], 
            Locations = TestDFList[[a]][,c("E", "N")],
            SPDEModel = SPDEList[[a]],
            Mesh = MeshList[[a]],
            FixedCovar = Covar, 
            NDraw = 100, Draw = T) %>% map_dbl(mean) -> Intercepts
  
  SpocialList[[a]]$Spatial$Model %>% 
    GetEstimates("AnnualDensity", NDraw = 100, Draws = T) -> Slopes
  
  1:length(Slopes) %>% map(~data.frame(X = X,
                                       Y = X*Slopes[[.x]] + Intercepts[[.x]]) %>% 
                             slice(1, n())) -> SlopeDF
  
  SlopeDF %<>% bind_rows(.id = "Rep")
  
  FitLine <- data.frame(
    AnnualDensity = rep(X, 1),
    Y = c(Y),
    Model = c("AnnualDensity")
  )
  
  TestHinds$Y <- TestHinds[,a]
  
  SpocialList[[a]]$Spatial$Model %>% 
    
    GetEstimates(Variable = "AnnualDensity") %>%
    paste0("; P = ", SpocialList[[a]]$Spatial$Model  %>% INLAPValue("AnnualDensity")) -> 
    LabelList[[a]]
  
  LabelYList[[a]] <- max(TestHinds$Y + diff(range(TestHinds$Y))/10)
  
  LabelXList[[a]] <- median(X)
  
  LabelDf <- data.frame(
    
    AnnualDensity = LabelXList[[a]],
    Y = LabelYList[[a]],
    Label = LabelList[[a]]
    
  )
  
  Figure5List[[a]] <- 
    ggplot(data = TestHinds, aes(AnnualDensity, Y)) + 
    geom_point(alpha = 0.1, colour = AlberColours[[1]]) + 
    geom_line(alpha = 0.1, data = SlopeDF, aes(X, Y, group = Rep)) +
    #geom_line(alpha = 0.1, data = SlopeDF2, aes(X, Y, group = Rep), 
    #          colour = "grey") +
    geom_line(data = FitLine, size = 1)  +
    labs(y = a, x = "Annual density") +
    scale_linetype_manual(values = 1:2,
                          breaks = c("AnnualDensity", "HROSPDE.AnnualDensity"),
                          labels = c("Base", "SPDE")) +
    geom_text(data = LabelDf, aes(label = Label), hjust = 0.5,
              colour = AlberColours[[1]], size = 3)# %>% plot
  
}

Figure1aList %>% ArrangeCowplot(

)

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

# Model Effects ####

CentralityList %>% map("Base") %>% Efxplot

1:length(CentralityList) %>% 
  map(~CentralityList[[.x]][c("Base","SPDE")] %>% 
        map("summary.fixed")) %>% 
  unlist(recursive = F) %>% 
  map(c(as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Model") -> Output

Output$Model %>% table() %>% extract2(1) %>% 
  divide_by(length(CentralityList)) %>% 
  rep(names(CentralityList), each = .) %>% rep(each = 2) ->
  Output$Response

Output %>% 
  filter(!rowname == "Intercept") %>% 
  ggplot(aes(rowname, mean, colour = Model)) + 
  geom_hline(yintercept = 0, lty = 2, size = 1, alpha = 0.3) +
  geom_errorbar(position = position_dodge(w = 0.5),
                aes(ymin = `0.025quant`, ymax = `0.975quant`)) +
  geom_point(position = position_dodge(w = 0.5)) + 
  facet_wrap(~Response, ncol = 2) + coord_flip() +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  ggsave("Untangling Code/Figures/FigureSI13.jpeg", units = "mm", width = 200, height = 300)

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

# Figure SI5: Spatiotemporal density etc. ####

SpocialList %>% map("dDIC") %>% 
  map(DICTableGet) %>% bind_rows(.id = "Response") %>% 
  mutate_at("Variable", ~str_remove_all(.x, "[*]")) -> DensityVarDIC

DensityVarDIC %>% 
  group_by(Variable) %>% 
  summarise_at("Delta", mean) %>% arrange(Variable) %>% 
  pull(Variable) %>% rev -> DensityVarOrder

DensityVarDIC %>% 
  
  mutate_at("Variable", ~factor(.x, levels = DensityVarOrder)) ->
  DensityVarDIC

DensityVarDIC %>% mutate_at("Response", ~{
  
  .x %>% str_split("_") %>% 
    map_chr(function(a) sapply(a, function(b) substr(b, 1, 1)) %>% 
              unlist %>% paste0(collapse = ""))
  
}) %>% pull(Response) %>% 
  str_replace_all(c("EW" = "WE", "SM" = "MS")) %>% 
  factor(., levels = rev(unique(.))) -> DensityVarDIC$ResponseLabel

DensityVarDIC %>%
  ggplot(aes(as.numeric(as.factor(Variable)), GregCube(Delta))) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_line(aes(colour = Response)) +
  geom_point(colour = "black", size = 2) + 
  geom_point(aes(colour = Response)) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(axis.title.y = element_text(vjust = -10)) + 
  labs(x = NULL, y = "DeltaDIC") +
  scale_x_continuous(limits = c(1, 8), breaks = 1:7, 
                     labels = VarLabels[levels(CombinedDICDF$Variable)]) +
  #scale_colour_discrete_sequential(palette = AlberPalettes[[2]],
  #                                 nmax = 12, 
  #                                 order = 12:4) +
  scale_colour_brewer(palette = "Spectral") +
  scale_fill_brewer(palette = "Spectral") +
  scale_y_continuous(limits = c(-15, 0), breaks = -c(0:3*5), labels = -c(0:3*5)^3) +
  ggrepel::geom_label_repel(data = CombinedDICDF %>% filter(Variable == "Spatiotemporal"), 
                            aes(x = 7, hjust = 0,
                                label = ResponseLabel, 
                                fill = Response),
                            alpha = 1,
                            colour = "black",
                            xlim = c(7.1, 8), 
                            force = 10) ->
  
  Plot1


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

# SISomething: Model effect comparisons ####

CentralityList %>% 
  map(~.x[c("Base", "HROSPDEDensity")] %>% 
        Efxplot(ModelNames = c("Base", "Spatial"))) -> EffectComparisonList

EffectComparisonList[[1]] + 
  EffectComparisonList[[2]] + 
  EffectComparisonList[[3]] + 
  EffectComparisonList[[4]] + 
  EffectComparisonList[[5]] + 
  EffectComparisonList[[6]] + 
  EffectComparisonList[[7]] + 
  EffectComparisonList[[8]] + 
  plot_layout(nrow = 4, guides = "collect") +
  plot_annotation(tag_levels = "A")

# Little figures for schematic ####

Censuses %>% 
  ggplot(aes(Easting, Northing)) + 
  geom_point(alpha = 0.01) + theme_void() + coord_fixed() -> ObservationPlot

SocTestHinds %>% 
  ggplot(aes(E, N)) + 
  geom_point(alpha = 0.5) + theme_void() + coord_fixed() -> CentroidPlot

LifetimeKUDL@coords %>% as.data.frame() %>% 
  rename(Y = Var1, X = Var2) %>% 
  mutate(Density = LifetimeKUDL$ud) %>% 
  filter(Density>0.000005) %>% 
  ggplot(aes(X, Y, fill = Density)) + 
  geom_tile() + coord_fixed() +
  scale_fill_gradient2(low = "grey", high = "black") + 
  theme_void() + 
  theme(legend.position = "none") +
  geom_contour(aes(z = Density), 
               colour = "white", size = 1, alpha = 0.8) ->
  DensityPlot

ObservationPlot +
  ggsave("Untangling Code/SIFigures/Observations.jpeg", units = "mm", width = 120, height = 150)

CentroidPlot +
  ggsave("Untangling Code/SIFigures/Centroids.jpeg", units = "mm", width = 120, height = 150)

DensityPlot +
  ggsave("Untangling Code/SIFigures/Density.jpeg", units = "mm", width = 120, height = 150)
