
# Sorting Rum boundary ####

if(file.exists("Output Files/RumBoundary.rds")){
  
  RumBoundary <- readRDS("Output Files/RumBoundary.rds")
  
}else{
  
  library(sf);library(spdep); library(ggplot2)
  
  Root <- "Data"
  
  myroute <- st_read(dsn = paste0(Root,"/myroute.kml"))
  
  myroute_line <- myroute[st_is(myroute, type = "LINESTRING"),]
  
  myroute_points <- myroute[st_is(myroute, type = "POINT"),]
  
  myroute_line_xy <- st_zm(myroute_line)
  
  myroute_rd_string <- st_transform(myroute_line_xy, crs = 28992)
  
  boundary <- as.data.frame(
    spdep::Rotation(
      data.frame(E=unlist(myroute_rd_string$geometry)[c(1:7,15:180)],N=unlist(myroute_rd_string$geometry)[c(8:14,181:346)]),
      6*pi/180))
  
  names(boundary)<-c("E","N")
  
  c1<-boundary[67,"E"]-1360
  c2<-boundary[67,"N"]-8040
  
  boundary$Easting<-boundary$E-c1
  boundary$Northing<-boundary$N-c2
  
  #  ggplot(boundary,aes(Easting,Northing))+geom_path()+geom_point()+coord_fixed()
  # ggplot(boundary,aes(E,N))+geom_path()+geom_point()+coord_fixed()+geom_point(data=boundary[151,],colour="red")
  
  c3<-(boundary[151,"E"]-boundary[65,"E"])/20.5
  
  boundary$Easting<-boundary$E/c3
  boundary$Northing<-boundary$N/c3
  
  
  c1<-boundary[67,"Easting"]-1360
  c2<-boundary[67,"Northing"]-8040
  
  boundary$Easting<-boundary$Easting-c1
  boundary$Northing<-boundary$Northing-c2
  
  #boundary2<-as.data.frame(Rotation(boundary[,c("Easting","Northing")], 355*pi/180))
  
  boundary[64:83,"Easting"]<-seq(boundary[c(64,83),"Easting"][1],boundary[c(64,83),"Easting"][2],length.out=length(64:83))
  boundary[64:83,"Northing"]<-boundary[64,"Northing"]+(boundary[64:83,"Easting"]-boundary[64,"Easting"])^2/16.15
  
  boundary<-boundary[-c(20:45),]
  
  boundary$Easting[146]<-1385
  boundary$Easting[13]<-1355
  
  boundary<-boundary[-which(boundary$Easting>1385|boundary$Easting<1355),]
  boundary[c(dim(boundary)[1]+1,dim(boundary)[1]+2),c("Easting","Northing")]<-cbind(c(1385,1355),c(7997.5,7997.5))
  
  Boundary <- boundary[,3:4]
  
  write.csv(Boundary, file = "RumBoundary.csv", row.names = F)
  
  poly1<-Polygon(boundary[,c("Easting","Northing")])
  
  RumBoundary <- SpatialPolygons(list(Polygons(list(poly1), ID = '0')))
  
  saveRDS(RumBoundary, file = "Output Files/RumBoundary.rds")
  
}
