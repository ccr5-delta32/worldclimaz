####################################################################################
#         Get interesting spatial coordinates of Az haplotypes per island          #
####################################################################################

snpdata <- read.table('./AzSPL9_SNP_genotype_all.tbl', header=TRUE, stringsAsFactors=FALSE)

floresAz1 <- snpdata[which(snpdata$minorarea == 'Flores' & snpdata$SNP_20580327 == 'C' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
floresAz1.avglat <- mean(floresAz1$latitude)
floresAz1.avglon <- mean(floresAz1$longitude)
# calculate distance between collection site and centroid using Pythagoras Theorem 
floresAz1$cendist <- sqrt((floresAz1$latitude-floresAz1.avglat)^2 + (floresAz1$longitude-floresAz1.avglon)^2)
floresAz1.select <- floresAz1[which(floresAz1$longitude == min(floresAz1$longitude)),][1,] 

faialAz1 <- snpdata[which(snpdata$minorarea == 'Faial' & snpdata$SNP_20580327 == 'C' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
faialAz1.avglat <- mean(faialAz1$latitude)
faialAz1.avglon <- mean(faialAz1$longitude)
# calculate distance between collection site and centroid using Pythagoras Theorem
faialAz1$cendist <- sqrt((faialAz1$latitude-faialAz1.avglat)^2 + (faialAz1$longitude-faialAz1.avglon)^2)
faialAz1.select <- faialAz1[which(faialAz1$cendist == min(faialAz1$cendist)),][1,]

# for the Az2 allele on Faial find the most Western one
faialAz2 <- snpdata[which(snpdata$minorarea == 'Faial' & snpdata$SNP_20580327 == 'G' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
faialAz2.select <- faialAz2[which(faialAz2$longitude == min(faialAz2$longitude)),][1,]

picoAz1 <- snpdata[which(snpdata$minorarea == 'Pico' & snpdata$SNP_20580327 == 'C' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
picoAz1.avglat <- mean(picoAz1$latitude)
picoAz1.avglon <- mean(picoAz1$longitude)
# calculate distance between collection site and centroid using Pythagoras Theorem 
picoAz1$cendist <- sqrt((picoAz1$latitude-picoAz1.avglat)^2 + (picoAz1$longitude-picoAz1.avglon)^2)
picoAz1.select <- picoAz1[which(picoAz1$longitude == min(picoAz1$longitude)),][1,] 

# for the Az2 allele on Faial find the most Eastern one
picoAz2 <- snpdata[which(snpdata$minorarea == 'Faial' & snpdata$SNP_20580327 == 'G' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
picoAz2.select <- picoAz2[which(picoAz2$longitude == max(picoAz2$longitude)),][1,]

# Sao Miguel has no Az1 so take the Az2 closest to the centroid
sao_miguelAz2 <- snpdata[which(snpdata$minorarea == 'Sao_Miguel' & snpdata$SNP_20580327 == 'G' & !is.na(snpdata$latitude) &!is.na(snpdata$longitude)),]
sao_miguelAz2.avglat <- mean(sao_miguelAz2$latitude)
sao_miguelAz2.avglon <- mean(sao_miguelAz2$longitude)
# calculate distance between collection site and centroid using Pythagoras Theorem 
sao_miguelAz2$cendist <- sqrt((sao_miguelAz2$latitude-sao_miguelAz2.avglat)^2 + (sao_miguelAz2$longitude-sao_miguelAz2.avglon)^2)
sao_miguelAz2.select <- sao_miguelAz2[which(sao_miguelAz2$longitude == min(sao_miguelAz2$longitude)),][1,] 
