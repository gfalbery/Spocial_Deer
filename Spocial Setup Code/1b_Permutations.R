
# Making network permutations for Josh's tests ####
  
library(tidyverse)

if(file.exists("Output Files/Perm Orders.Rdata")) load("Output Files/Perm Orders.Rdata") else{
  
  make.network<-function(gbi){
    yb<-xab<-matrix(NA,ncol(gbi),ncol(gbi))
    colsums<-colSums(gbi)
    for(i in 1:ncol(gbi)){
      sums<-gbi[which(gbi[,i]>0),]
      if(is.null(dim(sums))==T) {xab[,i]<-sums} else {xab[,i]<-colSums(sums)}}
    yb<-matrix(colsums,ncol(gbi),ncol(gbi)); ya<-t(yb)
    am<-xab/(xab+(ya-xab)+(yb-xab))
    diag(am)<-0;am[is.nan(am)]<-0;rownames(am)<-colnames(am)<-colnames(gbi)
    am}
  
  mround <- function(x,base){
    base*round(x/base)
  }
  
  for(SEASON in c("Rut","Spring")){
    
    #Make the Censuses object comparable to the network objects (i.e. only used groups)
    c.c<- droplevels(Censuses[Censuses$Season == SEASON&Censuses$Hind == "Y",])
    # c.c<- droplevels(Censuses[Censuses$Season == "Spring"&Censuses$Hind == "Y",])
    c.c<-c.c[! (is.na(c.c$GroupDate) | is.na(c.c$Code) | is.na(c.c$Easting)),]
    c.c$Row<-1:nrow(c.c)
    c.c$Code<-as.character(c.c$Code);c.c$GroupDate<-as.character(c.c$GroupDate)
    #Number of focal inds in the group:
    c.c$FocInds<-table(c.c$GroupDate)[c.c$GroupDate]
    
    #Change into a GBI
    gbi<-as.matrix(table(c.c$GroupDate,c.c$Code))
    gbi<-gbi[order(rownames(gbi)),order(colnames(gbi))]
    
    #Create group info
    gd<-unique(c.c[,c("Date","Ndate","Year","DeerYear","Season","GroupDate")])
    gd<-gd[order(gd$GroupDate),]
    gd$Row<-1:nrow(gd)
    
    #Get group mid-points
    mean.x<-tapply(c.c$Easting,c.c$GroupDate,mean,na.rm=T)
    mean.y<-tapply(c.c$Northing,c.c$GroupDate,mean,na.rm=T)
    gd$X<-mean.x[match(gd$GroupDate,names(mean.x))]
    gd$Y<-mean.x[match(gd$GroupDate,names(mean.x))]
    
    #Make list objects:
    l.gbi<-lapply(tapply(gd$Row,gd$DeerYear,unique),function(a){b<-gbi[a,];b[,colSums(b)>0]})
    l.gd<-lapply(tapply(gd$Row,gd$DeerYear,unique),function(a){b<-gd[a,]})
    l.c.c<-lapply(tapply(c.c$Row,c.c$DeerYear,unique),function(a){b<-c.c[a,]})
    l.am<-lapply(l.gbi,make.network)
    #DeerYearly individual list for needed info
    l.ind<-lapply(l.c.c,function(a)data.frame(Code=sort(unique(a$Code)),
                                              NObs=as.numeric(table(a$Code)),
                                              Obs=tapply(a$GroupDate,a$Code,as.character),
                                              Groups=tapply(a$GroupDate[a$FocInds>1],as.factor(a$Code)[a$FocInds>1],as.character),
                                              Xs=tapply(a$Easting,a$Code,as.numeric),Ys=tapply(a$Northing,a$Code,as.numeric),
                                              Pos=tapply(paste(a$Easting,a$Northing),a$Code,as.character),
                                              Grid.5=tapply(paste(mround(a$Easting,5),mround(a$Northing,5)),a$Code,as.character), 
                                              X.mean=tapply(a$Easting,a$Code,mean,na.rm=T),
                                              Y.mean=tapply(a$Northing,a$Code,mean,na.rm=T),stringsAsFactors=F))
    
    #Chosen Permutation procedure begins here:
    nperms<-1000 #How many we are going to run
    #As we aren't using gbi permutations, the overall structure of the network is maintained, so we just swap identities
    perm.orders<-sapply(l.am,function(a){b<-rownames(a);names(b)<-b;b}) #Vector for the new orders for individual IDs
    perm.orders<-rep(list(perm.orders),nperms)
    
    
    #Spatial Grid Permutation Procedure
    #Assign each individual to a grid space (where grid positions for each observation are the rounded x and y , and then allow swaps between individuals in the same grid space
    for(p in 1:nperms){
      for(y in 1:length(l.am)){
        #In each permutation, this permutation assigns each individual to a grid square each DeerYear, 
        # and then allows swaps of identity between those observed in the same assigned grid square for that DeerYear.
        
        l.ind.r<-l.ind[[y]] 
        
        #First, place some restrictions on the swaps
        #-In DeerYears when some individuals are observed once and others arent, only swap those observed more than once 
        # (because we won't include these 'observed once' as focal individuals in the age analysis, but need to be kept here for the network structure):
        # -note that in DeerYears when all individuals are observed once (i.e. 19), this doesn't need to be done as they'll be cut from analysis anyway
        
        if(sum(l.ind.r$NObs>1)>1){
          l.ind.r<-l.ind.r[l.ind.r$NObs>1,]}
        
        #Now, randomly assign each individual to grid position
        l.ind.r$Grid.5.p<-sapply(l.ind.r$Grid.5,function(a)sample(as.character(a),1))
        
        #Now permute so individuals swap with those in same grid space
        l.ind.r<-l.ind.r[order(l.ind.r$Grid.5.p),]
        l.ind.r$Code.p<-l.ind.r$Code[unlist(tapply(1:nrow(l.ind.r),l.ind.r$Grid.5.p,function(a)as.numeric(sample(as.character(a)))))]
        #Now save the newly ordered identities to the perm.orders vector
        perm.orders[[p]][[y]][l.ind.r$Code]<-l.ind.r$Code.p
      }}
    
    # The perm.orders object now holds all the info needed for carrying out same models but on permuted data. 
    # The format is the following, there is a vector within each list for each DeerYear within each list for each run of the permutations: 
    # e.g. perm.orders[[100]][[10]] would be the 100th permutation for the 10th DeerYear. 
    # This vector has the names of the observed individuals, with the characters set as the individual they are to be swapped identities with. 
    # Using this, the identity swap can take place at the network level (i.e. reassigning the response variables within each DeerYear into the new orders), 
    # or at the level of each predictor factor (just reassigning the order of one predictor at a time)
    
    PermOrderList[[SEASON]] <- perm.orders
    
  }
  
  save(PermOrderList,file="Output Files/Perm Orders.Rdata")
  
}

PermList <- YearPermList <- list()

perm.orders <- PermOrderList[[SEASON]]

i = 1

for(i in 1:1000){
  
  if(i %% 100 == 0) 
    print(i)
  
  Order = perm.orders[[i]]
  
  for(j in 1:length(Order)){
    
    if(length(longlist[[j]])>0){
      
      YearPermList[[j]] <- longlist[[j]]
      
      rownames(YearPermList[[j]]) <- YearPermList[[j]]$Name
      
      SubDF <- YearPermList[[j]][Order[[j]], Resps] %>%
        mutate(Name = names(Order[[j]]))
      
      YearPermList[[j]] <- YearPermList[[j]] %>% left_join(SubDF, by = "Name", suffix = c("_A","")) %>%
        suppressWarnings()
      
      names(YearPermList[[j]])
      
    }
  }
  
  PermList[[i]] <- bind_rows(YearPermList) %>% mutate(Perm_GroupSize = kader:::cuberoot(GroupSize),
                                                      Perm_Degree = sqrt(Degree),
                                                      Perm_Strength = sqrt(Strength),
                                                      Perm_Strength_Mean = as.vector(scale(kader:::cuberoot(Strength_Mean))),
                                                      Perm_GroupSize = as.vector(scale(kader:::cuberoot(GroupSize))),
                                                      Perm_Modularity = as.vector(scale(Modularity^2)),
                                                      Perm_Betweenness = as.vector(log(Betweenness + 1)),
                                                      Perm_Eigenvector = Eigenvector,
                                                      Perm_Eigenvector_Weighted = Eigenvector_Weighted)
  
  PermList[[i]]$Permutation <- i
  
  YearPermList <- list()
  
}

FinalPermList <- PermList %>% bind_cols() %>% select(Name, Year, starts_with("Perm"))

saveRDS(FinalPermList, file = paste0("Output Files/",SEASON,"PermList.rds"))
