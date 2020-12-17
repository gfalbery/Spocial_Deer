
# X_Affixing the centrality list ####

library(tidyverse); library(cowplot)

theme_set(theme_cowplot())

i <- 1

CentralityList <- list()

for(i in 1:length(Resps)){
  
  print(i)
  
  Models <- readRDS(paste0("Output Files/", Resps[i], "Models.rds"))
  
  CentralityList[[Resps[i]]] <- Models
  
}

saveRDS(CentralityList, file = "Output Files/FullCentralityList.rds")
