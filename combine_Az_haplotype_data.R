lpath <- '~/server/home/x-perimentz/Azores/Az_SPL9_sequencing/analysis/'
spath <- '/home/pieper/MPIPZ/x-perimentz/Azores/Az_SPL9_sequencing/analysis/'

location <- 'local' # either 'local' or 'cluster'

if (location == 'local') { 
  path <- lpath
} else if (location == 'cluster') {
  path <- spath
} else {
  print(paste("location should be either 'local' or 'cluster' while ", location, " was specified", sep=''))
}

snps <- read.table(paste(path, "SNP_genotypes", sep=''), header=TRUE, stringsAsFactors=FALSE)
coll <- read.table(paste(path, "collection_data", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
wayp <- read.table(paste(path, "waypoint_data", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
sq96 <- read.table(paste(path, "Seq96_snps", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
sali <- read.table(paste(path, "alias_wayp", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
alls <- read.table(paste(path, "alias_all_seqd", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')
gps  <- read.table(paste(path, "Az_coll1_waypoint_LizzieNew", sep=''), header=TRUE, stringsAsFactors=FALSE, sep='\t')

for (x in 1:length(snps[,1])) {
  snps$country[x] = coll$country[which(coll$waypoint == as.numeric(snps$line[x]))][1]
  snps$majorarea[x] = coll$majorarea[which(coll$waypoint == snps$line[x])][1]
  if (!is.na(snps$add[x]) & snps$add[x] == 'Ibe') {
    snps$minorarea[x] = 'Iberian'
  } else { 
    snps$minorarea[x] = coll$minorarea[which(coll$waypoint == snps$line[x])][1]
  }
  snps$latitude[x] = wayp$latitude_dec[which(wayp$waypoint == snps$line[x])][1]
  snps$longitude[x] = wayp$longitude_dec[which(wayp$waypoint == snps$line[x])][1]
}

for (x in 1:length(alls$alias)) {
  alls$country[x] = coll$country[which(coll$waypoint == alls$number[x])]
  alls$majorarea[x] = coll$majorarea[which(coll$waypoint == alls$number[x])]
  alls$minorarea[x] = coll$minorarea[which(coll$waypoint == alls$number[x])]
  alls$latitude[x] = gps$latitude[which(gps$waypoint == alls$number[x])][1]
  alls$longitude[x] = gps$longitude[which(gps$waypoint == alls$number[x])][1]
  if (alls$alias[x] %in% sq96$allele[which(sq96$position == 20580327)]) {
    alls$SNP_20580327[x] = sq96$consensus[which(sq96$allele == alls$alias[x] & sq96$position == 20580327)]
  } else {
    alls$SNP_20580327[x] = sq96$Oxford[which(sq96$position == 20580327)][1] 
  }
  if (alls$alias[x] %in% sq96$allele[which(sq96$position == 20580657)]) {
    alls$SNP_20580657[x] = sq96$consensus[which(sq96$allele == alls$alias[x] & sq96$position == 20580657)]
  } else {
    alls$SNP_20580657[x] = sq96$Oxford[which(sq96$position == 20580657)][1] 
  }
  if (alls$alias[x] %in% sq96$allele[which(sq96$position == 20581044)]) {
    alls$SNP_20581044[x] = sq96$consensus[which(sq96$allele == alls$alias[x] & sq96$position == 20581044)]
  } else {
    alls$SNP_20581044[x] = sq96$Oxford[which(sq96$position == 20581044)][1] 
  }
}
alls <- alls[1:61,]

both <- data.frame(c(snps$line, alls$number), c(snps$suffix, alls$suffix), c(snps$SNP_20580327, alls$SNP_20580327), c(snps$majorarea, alls$majorarea), c(snps$minorarea, alls$minorarea), c(snps$latitude, alls$latitude), c(snps$longitude, alls$longitude))
write.table(both, 'AzSPL9_SNP_genotype_all.tbl', quote=FALSE, sep='\t', row.names=FALSE) # colnames have been edited after saving
