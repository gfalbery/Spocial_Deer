# Subsampling Social and Spatial behvaiour Metrics

rm(list = ls())

Root <- "Data"

# Initial data import and cleaning ####

library(tidyverse);library(reshape2);library(igraph); library(ggregplot)

# Cleaning the phenotypic data ####

Censuses <- read.csv(paste0(Root, "/FullCensuses.csv"))
Dates <- read.csv(paste0(Root, "/Dates.csv"))
Dates <- Dates[order(Dates$Ndate),]
Censuses <- merge(Censuses,Dates[,c("Ndate","DateUK")],by.x="Date",by.y="DateUK",all.x=T)

YearStarts <- Dates[substr(Dates$DateUK,1,5)=="01/01","Ndate"]
SpringStarts <- Dates[substr(Dates$DateUK,1,5)=="01/05","Ndate"][2:52]

Timing <- data.frame(Year=1969:2019, YearStarts, SpringStarts)

Censuses$Year <- as.numeric(as.character(cut(Censuses$Ndate,breaks=YearStarts,labels=1969:2018)))
Censuses$DeerYear <- as.numeric(as.character(cut(Censuses$Ndate,breaks=SpringStarts,labels=1969:2018)))
Censuses$Season <- ifelse(Censuses$Year == Censuses$DeerYear, "Rut", "Spring")
Censuses$GroupDate<-paste(Censuses$Ndate,Censuses$Group,sep=",")

Censuses<-Censuses[-which(is.na(Censuses$Code)|Censuses$Code==""),
                   c("Date","Code","Easting","Northing","GroupSize","Ndate","Year","DeerYear","Season","GroupDate")]

Individuals <- read.csv(paste0(Root, "/Individuals.csv"))
Names <- read.csv(paste0(Root, "/Names.csv"))

Individuals<-merge(Individuals,Names,by="Code",all.x=TRUE)

Individuals$GivenName<-as.character(Individuals$GivenName)
Individuals$FamilyName<-as.character(Individuals$FamilyName)

Individuals[is.na(Individuals$GivenName),"GivenName"]<-Individuals[is.na(Individuals$GivenName),"FamilyName"]
Individuals$Animal<-Individuals$GivenName
colnames(Individuals)[colnames(Individuals)=="Birth.Date"]<-"BirthYear"

Individuals$Sex<-cut(Individuals$Sex,breaks=c(0,1.5,2.5,3.5),labels=c("F","M","3"))

for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
  Individuals[,x]<-factor(as.factor(Individuals[,x]),levels=c(0:10000,paste(0,1:10,sep="")))
}

for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
  Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x]<-paste(0,Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x],sep="")
}

Individuals$Birth.Date<-as.character(factor(with(Individuals,paste(BirthDay,BirthMonth,BirthYear,sep="/"))))
Individuals$Death.Date<-as.character(factor(with(Individuals,paste(DeathDay,DeathMonth,DeathYear,sep="/"))))

Individuals[Individuals$Birth.Date=="NA/NA/NA","Birth.Date"]<-NA
Individuals[Individuals$Death.Date=="NA/NA/NA","Death.Date"]<-NA

Individuals<-merge(Individuals,Dates[,c("DateUK","Ndate")],by.x="Birth.Date",by.y="DateUK",all.x=T)
Individuals<-merge(Individuals,Dates[,c("DateUK","Ndate")],by.x="Death.Date",by.y="DateUK",all.x=T)

colnames(Individuals)[c(length(colnames(Individuals))-1,length(colnames(Individuals)))]<-c("BirthDate","DeathDate")

Individuals$BirthYear<-as.numeric(as.character(Individuals$BirthYear))
Individuals$DeathYear<-as.numeric(as.character(Individuals$DeathYear))
Individuals$Name=Individuals$Code

Censuses <- merge(Censuses, Individuals[,c("Sex", "Name", "BirthYear")], 
                  by.x = c("Code"), by.y = c("Name"))

Censuses$Age <- with(Censuses, Year - BirthYear)
Censuses$Hind <- with(Censuses, ifelse(Age>2&Sex == "F", "Y", "N"))

# to make annual records for all sightings ####

SEASON = "Rut"

Records = 5

RutAMList <- longlist <- list()

for(x in (min(Censuses$Year,na.rm=TRUE):2017)){
  
  print(x)
  
  Censuses2 <- droplevels(Censuses[Censuses$DeerYear == x&
                                     Censuses$Season == SEASON&
                                     Censuses$Hind == "Y",])
  
  if(nrow(Censuses2)>0){
    M <- with(Censuses2, table(Code, GroupDate))
    SocGraph <- graph.incidence(M, weighted = T)
    DeerProj <- bipartite.projection(SocGraph)$proj1
    AssMat <- DeerProj %>% get.adjacency(attr = "weight") %>% as.matrix 
    
    NObs <- diag(AssMat) <- table(Censuses2$Code)
    
    N <- nrow(AssMat)
    
    # Making proportional ####
    
    A <- matrix(rep(table(Censuses2$Code), N), N)
    B <- matrix(rep(table(Censuses2$Code), each = N), N)
    
    AM <- AssMat/(A + B - AssMat)
    
    diag(AssMat) <- diag(AM) <- 0
    
    RutAMList[[x-(min(Censuses$Year,na.rm=TRUE)-1)]] <- AM
    
    # Making graph and deriving node traits ####
    
    DeerGraph <- graph_from_adjacency_matrix(AM, weighted = TRUE, mode = "undirected")
    
    Eigens <- data.frame(Name = rownames(AM),
                         Degree = colSums(AM>0), 
                         Strength = colSums(AM), 
                         Eigenvector = eigen_centrality(DeerGraph, scale=TRUE)$vector,
                         Eigenvector_Weighted = eigen_centrality(DeerGraph, weights = NA, scale=TRUE)$vector,
                         Betweenness = betweenness(DeerGraph),
                         Closeness = closeness(DeerGraph),
                         Clustering = transitivity(DeerGraph, type = "local")
    ) %>%
      mutate(Strength_Mean = Strength/Degree) %>%
      mutate(Strength_Mean = ifelse(is.na(Strength_Mean), 0, Strength_Mean))
    
    Eigens <- Eigens[order(Eigens$Name),]
    
    Censuses2$EastingJ <- as.numeric (Censuses2$Easting + as.integer (runif (nrow (Censuses2),-1,1)))
    Censuses2$NorthingJ <- as.numeric (Censuses2$Northing + as.integer (runif (nrow (Censuses2),-1,1)))
    
    E <- with(Censuses2,tapply(Easting,Code,mean))
    N <- with(Censuses2,tapply(Northing,Code,mean))
    
    Subdf <- data.frame(Year = x,
                        Name = levels(Censuses2$Code),
                        NObs = c(NObs),
                        GroupSize = as.numeric(with(Censuses2,tapply(GroupSize,Code,function(x) mean(x, na.rm = T)))),
                        E = as.numeric(E),N = as.numeric(N),
                        RiverDistance = as.numeric(abs(1363-E))
    )
    
    rownames(Subdf) <- Subdf$Name
    
    longlist[[x-(min(Censuses$Year,na.rm=TRUE)-1)]] <- Subdf
    
    longlist[[x-(min(Censuses$Year,na.rm=TRUE)-1)]]$Reg6 <- with(longlist[[x-(min(Censuses$Year,na.rm=TRUE)-1)]], LocToReg6(E, N))
    longlist[[x-(min(Censuses$Year,na.rm=TRUE)-1)]] <- merge(longlist[[x-(min(Censuses$Year,na.rm=TRUE)-1)]], Eigens, by="Name",all.X=TRUE)
    
  }
}

FullSociality <- bind_rows(longlist)

FullSociality <- FullSociality %>% filter(!Year == 1992)

Resps <- c("GroupSize", "Degree",
           "Strength", "Strength_Mean",
           "Eigenvector", "Eigenvector_Weighted",
           "Betweenness", "Clustering")

