#Attaching Social Network

{ 
  
  longlist3 <- FullSociality
  
  LBS <- read.csv(paste0(Root, "/LBS.csv"))
  PopSize <- read.csv(paste0(Root, "/PopSize.csv"))
  HindStatus <- read.csv(paste0(Root, "/HindStatus.csv"), header=TRUE)
  BirthWts <- read.csv(paste0(Root, "/BirthWts.csv"), header=TRUE)
  
  levels(HindStatus$ReprodStatus) <- c("Milk","Naive","Summer.Yeld","True.Yeld","Winter.Yeld")
  HindStatus$ReprodStatus <- factor(HindStatus$ReprodStatus,levels=c("Naive","True.Yeld","Summer.Yeld","Winter.Yeld","Milk"))
  
  Names <- read.csv(paste0(Root, "/Names.csv"), header=T)
  
  Individuals <- read.csv(paste0(Root, "/Individuals.csv"))
  
  Individuals <- merge(Individuals,Names,by="Code",all.x=TRUE)
  
  Individuals$GivenName <- as.character(Individuals$GivenName)
  Individuals$FamilyName <- as.character(Individuals$FamilyName)
  
  Individuals[is.na(Individuals$GivenName),"GivenName"] <- Individuals[is.na(Individuals$GivenName),"FamilyName"]
  Individuals$Animal <- Individuals$GivenName
  colnames(Individuals)[colnames(Individuals)=="Birth.Date"] <- "BirthYear"
  
  Individuals$Sex <- cut(Individuals$Sex,breaks=c(0,1.5,2.5,3.5),labels=c("F","M","3"))
  
  for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
    Individuals[,x] <- factor(as.factor(Individuals[,x]),levels=c(0:10000,paste(0,1:10,sep="")))
  }
  
  for(x in which(colnames(Individuals)=="BirthDay"):which(colnames(Individuals)=="DeathYear")){
    Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x] <- paste(0,Individuals[(sapply(Individuals[,x],function(y) nchar(substr(y,1,2)))==1)&!is.na(Individuals[,x]),x],sep="")
  }
  
  Individuals$Birth.Date <- as.character(factor(with(Individuals,paste(BirthDay,BirthMonth,BirthYear,sep="/"))))
  Individuals$Death.Date <- as.character(factor(with(Individuals,paste(DeathDay,DeathMonth,DeathYear,sep="/"))))
  
  Individuals[Individuals$Birth.Date=="NA/NA/NA","Birth.Date"] <- NA
  Individuals[Individuals$Death.Date=="NA/NA/NA","Death.Date"] <- NA
  
  Individuals <- merge(Individuals,Dates[,c("DateUK","Ndate")],by.x="Birth.Date",by.y="DateUK",all.x=T)
  Individuals <- merge(Individuals,Dates[,c("DateUK","Ndate")],by.x="Death.Date",by.y="DateUK",all.x=T)
  
  colnames(Individuals)[c(length(colnames(Individuals))-1,length(colnames(Individuals)))] <- c("BirthDate","DeathDate")
  
  Individuals$BirthYear <- as.numeric(as.character(Individuals$BirthYear))
  Individuals$DeathYear <- as.numeric(as.character(Individuals$DeathYear))
  Individuals$Name=Individuals$Code
  
  longlist3 <- merge(longlist3,PopSize,by.x="Year",by.y="DeerYear",sort=FALSE,all.x=T)
  longlist3 <- merge(longlist3,LBS[,c("Name","LRS","LBS")],by="Name",sort=FALSE,all.x=T)
  
  longlist3 <- merge(longlist3,Individuals[,c("Code","BirthDate",
                                              "DeathDate","DeathMonth","Sex","Animal","MumCode","Name","BirthYear","DeathYear","DeathType")],by.x="Name",by.y="Name",sort=FALSE,all.x=T)
  
  longlist3 <- longlist3 %>% 
    mutate(DeathYear = ifelse(DeathMonth%in%c("01","02","03","04"), DeathYear - 1, DeathYear)) %>%
    mutate(Age = Year-BirthYear,
           Longevity = DeathYear - BirthYear) %>%
    mutate(Age2 = Age^2,
           Hind = ifelse(Age>2&Sex == "F", "Y", "N"),
           Y2d = DeathYear - Year) %>%
    mutate(Survived = as.numeric(DeathYear == Year),
           Survived1 = as.numeric(DeathYear == Year + 1))
  
  longlist3 <- merge(longlist3,Individuals[,c("Code","MumCode","Sex","BirthDate","BirthYear","DeathDate","DeathMonth","DeathYear","DeathType")],by.x=c("Name","Year"),by.y=c("MumCode","BirthYear"),sort=FALSE,all.x=T,suffixes=c("",".Calf"))
  longlist3$Sex.Calf <- factor(longlist3$Sex.Calf,levels=c("F","M","3",""))
  longlist3$Sex.Calf[is.na(longlist3$Sex.Calf)] <- ""
  longlist3$DeathYear.Calf <- as.numeric(as.character(longlist3$DeathYear.Calf))
  longlist3$CalfSurvived <- factor(with(longlist3,ifelse(DeathYear.Calf>Year|is.na(DeathYear.Calf),"Y","N")),ordered=TRUE)
  longlist3$CalfMatured <- with(longlist3,ifelse((DeathYear.Calf>Year+2)|is.na(DeathYear.Calf),"Y","N"))
  longlist3$DeathType <- as.character(longlist3$DeathType)
  longlist3[is.na(longlist3$DeathType),"DeathType"] <- "Unknown"
  
  longlist3 <- longlist3[longlist3$Y2d>-1|is.na(longlist3$Y2d),]
  longlist3 <- longlist3[order(longlist3$Year),]
  
  longlist3 <- merge(longlist3,HindStatus[,c("Female","ReprodStatus","DeerYear")],by.x=c("Name","Year"),by.y=c("Female","DeerYear"),all.x=TRUE,sort=FALSE)
  
  longlist3$Year2 <- longlist3$Year+1
  longlist3$Year3 <- longlist3$Year-1
  
  longlist3 <- merge(longlist3,HindStatus[,c("Female","ReprodStatus","DeerYear")],by.x=c("Name","Year2"),by.y=c("Female","DeerYear"),all.x=TRUE,sort=FALSE,suffixes=c("",".t2"))
  longlist3 <- merge(longlist3,HindStatus[,c("Female","ReprodStatus","DeerYear")],by.x=c("Name","Year3"),by.y=c("Female","DeerYear"),all.x=TRUE,sort=FALSE,suffixes=c("",".t0"))
  
  longlist3$CalfSurvived <- with(longlist3,ifelse(Hind=="N",NA,
                                                  ifelse(ReprodStatus=="Milk","Y","N")))
  
  Individuals$BirthYear2 <- Individuals$BirthYear-1
  longlist3 <- merge(longlist3,Individuals[,c("Code","MumCode","Sex","BirthDate","BirthYear2","DeathDate","DeathMonth","DeathYear")],by.x=c("Name","Year"),by.y=c("MumCode","BirthYear2"),sort=FALSE,all.x=T,suffixes=c("",".Calft2"))
  longlist3$Sex.Calft2 <- factor(longlist3$Sex.Calft2,levels=c("F","M","3",""))
  longlist3$Sex.Calft2[is.na(longlist3$Sex.Calft2)] <- ""
  longlist3$DeathYear.Calft2 <- as.numeric(as.character(longlist3$DeathYear.Calft2))
  longlist3$Calft2Survived <- factor(with(longlist3,ifelse(DeathYear.Calft2>Year+1|is.na(DeathYear.Calft2),"Y","N")),ordered=TRUE)
  longlist3$Calft2Matured <- with(longlist3,ifelse((DeathYear.Calft2>Year+3)|is.na(DeathYear.Calft2),"Y","N"))
  
  longlist3$Pregnant <- with(longlist3,ifelse(Hind=="N",NA,
                                              ifelse(ReprodStatus.t2%in%c("Naive","True.Yeld"),0,1)))
  
  longlist3$cBirthDate.Calft2 = longlist3$BirthDate.Calft2 - with(longlist3,tapply(BirthDate.Calft2, Year, function(a) mean(a, na.rm = T)))[as.character(longlist3$Year)]
  longlist3$cBirthDate.Calf = longlist3$BirthDate.Calf - with(longlist3,tapply(BirthDate.Calf, Year, function(a) mean(a, na.rm = T)))[as.character(longlist3$Year)]
  
  longlist3 <- merge(longlist3, BirthWts[,c("Code","BirthWt")], by.x = "Name", by.y = "Code", all.x = T)
  longlist3 <- merge(longlist3, BirthWts[,c("Code","BirthWt")], by.x = "Code.Calf", by.y = "Code", all.x = T, suffixes = c("",".Calf"))
  
  SocialFactors <- c("Eigenvector","Degree","GroupSize","HRA")
  longlist3 <- merge(longlist3,Individuals[,c("Code","Status")], by.x = "Name", by.y="Code",all.x=TRUE,sort=FALSE,suffixes=c("",".Death"))
  
  longlist3 <- longlist3[!is.na(longlist3$Name),]
  
  # Setting up test data frames ####
  
  SocialHinds <- longlist3
  all(which(is.na(SocialHinds$Degree)) == which(is.na(SocialHinds$Eigenvector)))
  SocialHinds <- SocialHinds[-which(SocialHinds$ReprodStatus == "Naive"), ]
  SocialHinds$ReprodStatus <- factor(SocialHinds$ReprodStatus, 
                                     levels = levels(SocialHinds$ReprodStatus)[2:5])
  
  SocialHinds %<>%
    arrange(Name, Year) %>% 
    mutate(Reprod = as.numeric(!ReprodStatus == "True.Yeld")) %>% 
    group_by(Name) %>%  
    mutate(CumReprod = cumsum(Reprod)) %>% 
    ungroup()
  
  SocialHinds$AgeCat <- with(SocialHinds, cut(Age, breaks = c(0,4,10,25), 
                                              label = c("Young", "Med", "Old")))
  
  SocialHinds$Age <- scale(SocialHinds$Age)[,1]
  SocialHinds$Age2 <- scale(poly(SocialHinds$Age,2)[,2])[,1]
  
  SocialHinds <- SocialHinds[-(which(SocialHinds$Y2d==0)),]
  
  SocialHinds <- SocialHinds %>% 
    mutate(fY2d = cut(Y2d, breaks = c(-1,0,1,2,30), labels = c(0,1,2,"3Plus")),
           Terminal = as.factor(as.numeric(Y2d==0)),
           Strength2 = Strength/Degree) %>%
    mutate(Strength2 = ifelse(is.na(Strength2), 0, Strength2))
  
  SocialHinds <- SocialHinds %>% 
    mutate(RelStrength = Strength - NoRelStrength,
           RelDegree = Degree - NoRelDegree)
  
  SocialHinds <- SocialHinds %>% mutate(Shot = as.numeric(DeathType == "S"),
                                        fYear = as.factor(Year))
  
  SocialHinds <- SocialHinds %>% 
    filter(!Year %in% c(2003, 1992)) %>%
    mutate(MeshInclude = MeshInclude(E, N, 1355, 1384.777, 7997.5, 8050)) %>% 
    group_by(Name) %>% 
    mutate(LifetimeE = mean(E, na.rm = T), LifetimeN = mean(N, na.rm = T)) %>%
    ungroup()
  
  SaveSocialHinds <- SocialHinds
  
}

