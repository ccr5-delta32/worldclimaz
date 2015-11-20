library(raster)
library(rgdal)

# Load the below image to prevent long raster processing times
# it contains data for the Azores plus Cologne and Lisbon
load('../Azores_Cologne_Lisbon.Rdata')

lpath <- '~/server/biodata/common/Bjorn/WorldClim_data/'
spath <- '/biodata/dep_tsiantis/common/Bjorn/WorldClim_data/'

location <- 'local' # either 'local' or 'cluster'

if (location == 'local') { 
  path <- lpath
} else if (location == 'cluster') {
  path <- spath
} else {
  stop(paste("location should be either 'local' or 'cluster' while ", location, " was specified", sep=''))
}

####################################################################################
#          Open raster objects and extract data with spatial coordinates           #
####################################################################################

EastAz <- raster(paste(path, 'current_EastAzores/tmax_15_tif/tmax1_15.tif', sep=''))
WestAz <- raster(paste(path, 'current_WestAzores/tmax_14_tif/tmax1_14.tif', sep=''))
Az <- merge(EastAz, WestAz)
plot(Az, col=grey(1:99/100), main="tmax_1", xlim=c(-31.6, -24.7), ylim=c(36, 40))

# Extract cell data with latitude and longitude
Az.pts <- rasterToPoints(Az, spatial=TRUE)
proj4string(Az.pts)

geo.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
Az.pts <- spTransform(Az.pts, CRS(geo.proj))
proj4string(Az.pts)
Az.pts@data <- data.frame(Az.pts@data, long=coordinates(Az.pts)[,1], lat=coordinates(Az.pts)[,2])

####################################################################################
#                  Determine locations to extract layer data from                  #
####################################################################################

# The script ./find_locations_for_comparison.R was used to prepare ./locations_for_extraction.tbl
select <- read.table('./locations_for_extraction.tbl', header=TRUE, stringsAsFactors=FALSE)

####################################################################################
#                    Extract layer data at selected locations                      #
####################################################################################

Az.pts.pos <- Az.pts[which(Az.pts$layer >= 0),]
select$tmax1 <- vector(length=length(select[,1]), mode='numeric')

for (x in 1:length(select[,1])) {
  select$tmax1[x] <- Az.pts$layer[which(abs(Az.pts$long - select$longitude[x]) == min(abs(Az.pts$long - select$longitude[x])) & abs(Az.pts$lat - select$latitude[x]) == min(abs(Az.pts$lat - select$latitude[x])))]
}

####################################################################################
#       Now extract layer data at selected locations for all variables/times       #
####################################################################################

# tmax
tmax <- matrix(nrow=length(select[,1]), ncol=12)
for (x in 1:12) {
  EastAz <- raster(paste(path, 'current_EastAzores/tmax_15_tif/tmax', x, '_15.tif', sep=''))
  WestAz <- raster(paste(path, 'current_WestAzores/tmax_14_tif/tmax', x, '_14.tif', sep=''))
  cologn <- raster(paste(path, 'current_Cologne/tmax_16_tif/tmax', x, '_16.tif', sep=''))
  Az <- merge(EastAz, WestAz, cologn)
  Az.pts <- rasterToPoints(Az, spatial=TRUE)
  proj4string(Az.pts)
  geo.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  Az.pts <- spTransform(Az.pts, CRS(geo.proj))
  proj4string(Az.pts)
  Az.pts@data <- data.frame(Az.pts@data, long=coordinates(Az.pts)[,1], lat=coordinates(Az.pts)[,2])
  for (y in 1:length(select[,1])) {
    tmax[y, x] <- Az.pts$layer[which(abs(Az.pts$long - select$longitude[y]) == min(abs(Az.pts$long - select$longitude[y])) & abs(Az.pts$lat - select$latitude[y]) == min(abs(Az.pts$lat - select$latitude[y])))]
  }
}
tmax <- tmax/10

# tmean
tmean <- matrix(nrow=length(select[,1]), ncol=12)
for (x in 1:12) {
  EastAz <- raster(paste(path, 'current_EastAzores/tmean_15_tif/tmean', x, '_15.tif', sep=''))
  WestAz <- raster(paste(path, 'current_WestAzores/tmean_14_tif/tmean', x, '_14.tif', sep=''))
  cologn <- raster(paste(path, 'current_Cologne/tmean_16_tif/tmean', x, '_16.tif', sep=''))
  Az <- merge(EastAz, WestAz, cologn)
  Az.pts <- rasterToPoints(Az, spatial=TRUE)
  proj4string(Az.pts)
  geo.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  Az.pts <- spTransform(Az.pts, CRS(geo.proj))
  proj4string(Az.pts)
  Az.pts@data <- data.frame(Az.pts@data, long=coordinates(Az.pts)[,1], lat=coordinates(Az.pts)[,2])
  for (y in 1:length(select[,1])) {
    tmean[y, x] <- Az.pts$layer[which(abs(Az.pts$long - select$longitude[y]) == min(abs(Az.pts$long - select$longitude[y])) & abs(Az.pts$lat - select$latitude[y]) == min(abs(Az.pts$lat - select$latitude[y])))]
  }
}
tmean <- tmean/10

# tmin
tmin <- matrix(nrow=length(select[,1]), ncol=12)
for (x in 1:12) {
  EastAz <- raster(paste(path, 'current_EastAzores/tmin_15_tif/tmin', x, '_15.tif', sep=''))
  WestAz <- raster(paste(path, 'current_WestAzores/tmin_14_tif/tmin', x, '_14.tif', sep=''))
  cologn <- raster(paste(path, 'current_Cologne/tmin_16_tif/tmin', x, '_16.tif', sep=''))
  Az <- merge(EastAz, WestAz, cologn)
  Az.pts <- rasterToPoints(Az, spatial=TRUE)
  proj4string(Az.pts)
  geo.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  Az.pts <- spTransform(Az.pts, CRS(geo.proj))
  proj4string(Az.pts)
  Az.pts@data <- data.frame(Az.pts@data, long=coordinates(Az.pts)[,1], lat=coordinates(Az.pts)[,2])
  for (y in 1:length(select[,1])) {
    tmin[y, x] <- Az.pts$layer[which(abs(Az.pts$long - select$longitude[y]) == min(abs(Az.pts$long - select$longitude[y])) & abs(Az.pts$lat - select$latitude[y]) == min(abs(Az.pts$lat - select$latitude[y])))]
  }
}
tmin <- tmin/10

# precipitation 
prec <- matrix(nrow=length(select[,1]), ncol=12)
for (x in 1:12) {
  EastAz <- raster(paste(path, 'current_EastAzores/prec_15_tif/prec', x, '_15.tif', sep=''))
  WestAz <- raster(paste(path, 'current_WestAzores/prec_14_tif/prec', x, '_14.tif', sep=''))
  cologn <- raster(paste(path, 'current_Cologne/prec_16_tif/prec', x, '_16.tif', sep=''))
  Az <- merge(EastAz, WestAz, cologn)
  Az.pts <- rasterToPoints(Az, spatial=TRUE)
  proj4string(Az.pts)
  geo.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  Az.pts <- spTransform(Az.pts, CRS(geo.proj))
  proj4string(Az.pts)
  Az.pts@data <- data.frame(Az.pts@data, long=coordinates(Az.pts)[,1], lat=coordinates(Az.pts)[,2])
  for (y in 1:length(select[,1])) {
    prec[y, x] <- Az.pts$layer[which(abs(Az.pts$long - select$longitude[y]) == min(abs(Az.pts$long - select$longitude[y])) & abs(Az.pts$lat - select$latitude[y]) == min(abs(Az.pts$lat - select$latitude[y])))]
  }
}

####################################################################################
#                           Plot together with daylength                           #
####################################################################################

library(BP.Rfunlib)

midmonth <- c(16.0, 45.5, 75.0, 105.5, 136.0, 166.5, 197.0, 228.0, 258.5, 289.0, 319.5, 350.0)
intermonth <- c(31.5, 59.5, 90.5, 120.5, 151.5, 181.5, 212.5, 243.5, 273.5, 304.5, 334.5)
#cols <- c(552, 452, 578, 267, 148, 173)
cols <- c("red", "darkorchid2", "darkorange2", "darkgreen", "cyan3", "burlywood4", "blue", "magenta")

pdf("Monthly_dayl_temp_prec_AzColLis.pdf", width=13, height=9)
  close.screen(all.screens=TRUE)
  split.screen(c(3,1))
  screen(1)
    par(mar=c(2,3.5,0.1,0.1))
    plot(1:365, type='n', axes=FALSE, xaxs='i', yaxs='i', xlim=c(1,365), ylim=c(0,24.4), xlab='', ylab='')
    box()
    abline(v=intermonth, col='grey82')
    for (y in 1:length(select[,1])) {
      for (x in 1:365) {
        points(x, dayl(select$latitude[y], x, 0.8333), col=cols[y], pch=15, cex=0.3)
      }
    }
    axis(side=2, at=seq(0, 24, by=2), tck=-0.012, las=1)
    mtext('daylength (hours)', side=2, line=2.6) 
    legend('top', horiz=TRUE, legend=select$location, lty=1, col=cols, box.lwd=0, box.col='white', bg='white', lwd=1.5)
    box()
  screen(2)
    par(mar=c(2,3.5,0.1,0.1))
    for (y in 1:length(select[,1])) {
      tmaxmat <- matrix(ncol=2, nrow=14, data=c(c(midmonth[1]-((365-midmonth[12])+midmonth[1]),midmonth, midmonth[12]+((365-midmonth[12])+midmonth[1]),c(tmax[y,12], tmax[y,], tmax[y,1]))))
      tmeanmat <- matrix(ncol=2, nrow=14, data=c(c(midmonth[1]-((365-midmonth[12])+midmonth[1]),midmonth, midmonth[12]+((365-midmonth[12])+midmonth[1]),c(tmean[y,12], tmean[y,], tmean[y,1]))))
      tminmat <- matrix(ncol=2, nrow=14, data=c(c(midmonth[1]-((365-midmonth[12])+midmonth[1]),midmonth, midmonth[12]+((365-midmonth[12])+midmonth[1]),c(tmin[y,12], tmin[y,], tmin[y,1]))))
      if (y == 1) {
        plot(midmonth, tmean[y,], pch=15, cex=0.5, col=cols[y], ylim=c(min(tmin), max(tmax)), ylab='', xlim=c(1,365), xaxs='i', axes=FALSE)
        abline(v=intermonth, col='grey82')
        lines(tmaxmat, col=cols[y], lwd=1, type='l', lty=2)
        lines(tmeanmat, col=cols[y], lwd=1, type='l')
        lines(tminmat, col=cols[y], lwd=1, type='l', lty=2)
      } else {
        points(midmonth, tmean[y,], pch=15, cex=0.5, col=cols[y])
        lines(tmaxmat, col=cols[y], lwd=1, type='l', lty=2)
        lines(tmeanmat, col=cols[y], lwd=1, type='l')
        lines(tminmat, col=cols[y], lwd=1, type='l', lty=2)
      }
    }
    axis(side=2, at=seq(-2, 26, by=2), tck=-0.012, las=1)
    mtext(~degree~C, side=2, line=2.6)
    legend("bottom", legend=c('mean temperature', 'min/max temperature'), lty=c(1,2), lwd=1.5, box.lwd=0, box.col='white', bg='white')
    box()
  screen(3)
    par(mar=c(2,3.5,0.1,0.1))
    for (y in 1:length(select[,1])) {
      precmat <- matrix(ncol=2, nrow=14, data=c(c(midmonth[1]-((365-midmonth[12])+midmonth[1]),midmonth, midmonth[12]+((365-midmonth[12])+midmonth[1]),c(prec[y,12], prec[y,], prec[y,1]))))
      if (y == 1) {
        plot(midmonth, prec[y,], pch=15, cex=0.5, col=cols[y], ylim=c(min(prec), max(prec)), ylab='', xlim=c(1,365), xaxs='i', axes=FALSE)
        abline(v=intermonth, col='grey82')
        lines(precmat, col=cols[y], lwd=1, type='l')
      } else {
        points(midmonth, prec[y,], pch=15, cex=0.5, col=cols[y])
        lines(precmat, col=cols[y], lwd=1, type='l')
      }
    }
    axis(side=1, at=midmonth, labels=c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), tck=-0.012)
    axis(side=2, at=seq(0, 220, by=50), tck=-0.012, las=1)
    mtext('precipitation (mm)', side=2, line=2.6)
    box()
dev.off()
