pdata <- read.table('./photoperiods', header=TRUE, stringsAsFactors=FALSE)
pdata$timevector <- 24 * as.numeric(times(pdata$dayLength))
library(BP.Rfunlib)

locations <- list()
locations[['Ponta Delgada']] <- 37.8244
locations[['Horta']] <- 38.5794
locations[['Santa Cruz das Flores']] <- 39.45151
locations[['Cologne']] <- 50.9364 
locations[['Oxford']] <- 51.7519 

plot(1:365, type='n', axes=FALSE, xaxs='i', yaxs='i', xlim=c(1,365), ylim=c(0,24), xlab='', ylab='')
midmonth <- vector(length=12, mode="numeric")
intermonth <- vector(length=12, mode="numeric")
for (x in 1:12) {
  abline(v=max(pdata$days[pdata$Month == unique(pdata$Month)[x]]), col='grey82')
  midmonth[x] <- min(pdata$days[pdata$Month == unique(pdata$Month)[x]]) + (max(pdata$days[pdata$Month == unique(pdata$Month)[x]])-min(pdata$days[pdata$Month == unique(pdata$Month)[x]]))/2
  intermonth[x] <- max(pdata$days[which(pdata$Month == unique(pdata$Month)[x])]) + 0.5
}
axis(side=1, at=midmonth, labels=unique(pdata$Month))
axis(side=2, at=seq(0,24, by=1), las=1)
box()

custom_cols <- list()
custom_cols[['col1']] <- rgb(255,100,100,150, maxColorValue=255)
custom_cols[['col2']] <- rgb(200,150,100,150, maxColorValue=255)
custom_cols[['col3']] <- rgb(200,100,150,150, maxColorValue=255)
custom_cols[['col4']] <- rgb(150,100,255,150, maxColorValue=255)
custom_cols[['col5']] <- rgb(100,200,55,150, maxColorValue=255)

for (y in 1:length(locations)) {
  for (x in 1:365) {
    points(x, dayl(locations[[names(locations)[y]]], x, 0.8333), col=custom_cols[y][[1]], pch=15, cex=0.4)
  }
} 
