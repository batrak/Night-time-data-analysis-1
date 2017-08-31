# Nighttime-Data-Analysis


##India's Nighttime Spatial Analysis (State-wise)



VIIRS Nighttime data obtained from: https://ngdc.noaa.gov/eog/viirs/download_dnb_composites.html (2017 monthly data, June month, Tile3_75N060E)

India Shapefile: http://www.gadm.org/country (.rds file type)

India Population Data: http://www.census2011.co.in/states.php (used as a .csv file)

```{r}
library(doParallel)
library(foreach)
library(raster)
library(sp)
library(rgdal)
library(ggmap)
library(plotly)
library(ggplot2)
library(devtools)
getwd()
setwd("C:/Users/admin1/Desktop/imagery")

imagery = "[path]/imagery"

##Obtain a list of TIF files, load in the first file in list. Stored as 'imagery' folder
tifs = list.files(imagery,pattern = "\\.tif")
rast <- raster("image.tif")
rast

wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
projection(rast) <- CRS(wgs84)
 
      
bound<-readRDS("IND_adm1.rds")

projection(bound) <- CRS(wgs84)
names(bound)

require(rgdal)
# Read SHAPEFILE.shp from the current working directory (".")
shape <- readOGR(dsn = "C:/Users/admin1/Desktop/imagery", layer = "IND_adm1")
projection(shape) <- CRS(wgs84)
s<-shapefile("IND_adm1")

pop<- read.csv("pca.csv")
pop$NAME <- as.character(pop$State) 




```

Now, after selecting the states to be compared, we geocode the states and create their spatial bounding box. This is then extended via longitude and latitide to include a greater statespatial box and cropped from the raster tile. This cropped raster formes the sample size within which we first identify clusters of 15 and create a separate  data frame which we can then plot. 
```{r}
states<- c("Maharashtra,Mh","Bihar,Bi","Gujarat,Gu","West Bengal,WB","Kerala,Kr","Madhya Pradesh,MP")
par(mai=c(0,0,0,0),mfrow = c(1,1),bg='#001a4d', bty='n')

coords <- data.frame() ##place holder

for(i in 1:length(states)){
  
  ##Coords
  temp_coord <- geocode(states[i], source = "google")
  coords <- rbind(coords,temp_coord)
  
  e <- extent(temp_coord$lon - 4, temp_coord$lon + 4,
              temp_coord$lat - 2, temp_coord$lat + 2)
  rc <- crop(rast, e)    
  
  ##Rescale brackets
  sampled <- as.vector(rc)
  clusters <- 15
  clust <- kmeans(sampled,clusters)$cluster
  combined <- as.data.frame(cbind(sampled,clust))
  brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = T, asp=1.5); plot(bound, add=T, border="white")
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(states[i],1,regexpr(",",states[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}

```
For inter-state comparisons, now the entire raster is the sample size from which clusters are identified. Now, these clusters are sorted by radiance. Post this, individual spatial bound boxes for each state are made and plotted. 

```{r}

par(mai=c(0,0,0,0),mfrow = c(1,1),bg='#001a4d', bty='n')

#Run clustering
set.seed(123) #set seed for reproducibility
sampled <- sample(rast, 20000) #sample 20,000 pixels
clusters <- 15 ##15 clusters
clust <- kmeans(sampled,clusters)$cluster
combined <- as.data.frame(cbind(sampled,clust))
brk <- sort(aggregate(combined[,1], list(combined[,2]), max)[,2])

##Loop through each city
for(i in 1:length(states)){
  
  temp_coord <- coords[i,] ##re-use the coordinates 
  e <- extent(temp_coord$lon - 4, temp_coord$lon + 4,
              temp_coord$lat - 2, temp_coord$lat + 2)
  rc <- crop(rast, e)    
  
  #Plots
  plot(rc, breaks=brk, col=colorRampPalette(c("#001a4d","#0066FF", "yellow"))(clusters), 
       legend=F,yaxt='n',xaxt='n',frame = F, asp=1.5)
  plot(bound, add=T, border="white")
  text(temp_coord$lon ,temp_coord$lat + 0.15,
       substr(states[i],1,regexpr(",",states[i])-1), 
       col="white", cex=1.25)
  
  rm(combined)
}

```

For statistical comparison, we extract radiance values of all states' geoshapes through a new function, masq.
```{r}


masq <- function(s,rast,i){
  #Extract one polygon based on index value i
  polygon <- s[i,] #extract one polygon
  extent <- extent(polygon) #extract the polygon extent 
  
  #Raster extract
  outer <- crop(rast, extent) #extract raster by polygon extent
  inner <- mask(outer,polygon) #keeps values from raster extract that are within polygon
  
  #Convert cropped raster into a vector
  #Specify coordinates
  coords <- expand.grid(seq(extent@xmin,extent@xmax,(extent@xmax-extent@xmin)/(ncol(inner)-1)),
                        seq(extent@ymin,extent@ymax,(extent@ymax-extent@ymin)/(nrow(inner)-1)))
  #Convert raster into vector
  rastdata <- as.vector(inner)
  
  #package data in neat dataframe
  rastdata <- cbind(as.character(s@data$NAME_1[i]),coords, rastdata) 
  colnames(rastdata)<-c("GEOID","lon","lat","avg_rad") #note that 
  rastdata <- rastdata[!is.na(rastdata$avg_rad),] #keep non-NA values only
  
  return(rastdata)
}




```

To create a histogram for the 6 state comparison:

Average radiances are extracted from the raster and shape files for the states and plotted. However, theses plots do not take into account the area of each state and analysis based on the histograms should be made keeping this fact in mind.

```{r}
 

skt<-c(5,11,17,19,20,36)

  radiances <- data.frame() 
 
for(i in skt){
  
    print(skt)
    
    #Extract i polygon
      shp_temp <- shape[shape@data$ID_1==i,]
    
   #if(regexpr(" ",as.character(shp_temp@data$NAME_1)[1])[1]==-1){
        loc = as.character(shp_temp@data$NAME_1)[1]
      #} else{
        #loc = (as.character(shp_temp@data$NAME_1)[0])
      #}
    
    #Extract the radiances, append to radiances placeholder
      rad <- masq(shp_temp,rast,1)$avg_rad 
      temp <- data.frame(loc = as.character(paste(loc,"(TNL = ",round(sum(rad),0),")",sep="")), avg_rad = rad) 
      radiances <- rbind(radiances,temp)
      #print(radiances)
}

#Use ggplot to create histograms by States.
  ggplot(radiances, aes(x=log(avg_rad))) +
    geom_histogram(position="identity", alpha=0.6) +
    facet_grid(. ~ loc)

#Remove all axes labels for style
    x <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    )
    y <- list(
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
    ) 
    
#Initiate a plotly graph without axes
  ggplotly()  %>% layout(xaxis=x, yaxis=y)
  
```
Scatter plot of Total Nighttime Light and Population to see the underlying correlation between the two. 

For this, extracted GEOIDs and average radiance data frame extracted is merged with the population data frame and then plotted. 

There exists a positive correlation between the two parameters. 

```{r}
library(sp)
  registerDoParallel(cores=2)
    extract <- foreach(i=1:nrow(s@data),.packages='raster', .combine=rbind) %dopar% {
        rastdata <- masq(shape,rast,i)
       data.frame(GEOID = rastdata$GEOID[1],sum = sum(rastdata$avg_rad))
    }
   #extract$GEOID <- as.numeric(as.character(extract$GEOID))
    
   
   
  # D <- as.numeric(as.character(extract$GEOID))
    
  ##Join in data
    
  joined<-merge(extract, pop[,c("Density","NAME","Population")],by.x="GEOID",by.y="NAME")
   
    colnames(joined) <- c("State","TNL","Density","Population")
    
    aTNL<-log(joined$TNL)
    bPop<- log(joined$Population)
    cState<-joined$State
    
    joinedfd<-cbind(joined,aTNL,bPop,cState)
    attach(joinedfd)
    plot_ly(joinedfd, 
            x = aTNL,
            y = bPop, 
            text = paste("State: ",cState),
            mode = "markers", 
            color = TNL,colors="PuOr")  %>% 
            layout(title="Total Nighttime Light vs. Population", showlegend = F)


```




```{r}
