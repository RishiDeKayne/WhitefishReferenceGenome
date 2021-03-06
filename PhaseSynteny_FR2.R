#PhaseSynteny_FR2.R
#analysis of Falcon genome assembly for De-Kayne, Zoller and Feulner 2019

#load in chromosome lengths and linkage map lengths
setwd("Dropbox/RishiMAC/Genome/PostPhaseScaffoldingFastas/")
lengthsdf <- read.csv(file = "FalconR2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
lm <- read.csv(file = "LM_stats.txt", header = FALSE, sep = " ")
lengthsdf$Scaffold <- 0:39

#load in samfile already separated columns in excel...
samfile <- read.csv(file = "SA_linkagemap_FalconR2PhaseMapped_Filtered.csv")
samfile <- subset(samfile, as.character(samfile$PN.bwa) == "PGA")

#make synteny map structure file with each linkage group/chromosome and start and end
#find equivalent chromosomes
for (i in 1:40){
  linkagegroup <- levels(lm$V1)[i]
  new_linkagegroup <- gsub("Calb", "SA", linkagegroup)
  sam_subset <- subset(samfile, as.character(samfile$X) == as.character(new_linkagegroup))
  scaffoldvect <- as.vector(sam_subset$X.2)
  abundant <- names(sort(summary(as.factor(scaffoldvect)), decreasing=T)[1:1])
  lm$scaffold[i] <- abundant
  lm$tot_markernumber[i] <- length(sam_subset$X.PG)
  abundant_number <- subset(sam_subset, as.character(sam_subset$X.2) == as.character(abundant))
  lm$scaf_markernumber[i] <- length(abundant_number$X.PG)
  lm$proportion_abundant <- ((lm$scaf_markernumber/lm$tot_markernumber)*100)
}

orderedbyscaffolds <- lm[order(lm$scaffold), ]

#now get lg data for plotting
lg_dat_1 <- as.data.frame(lm$V1)
lg_dat_1$cM <- lm$V3

lg_dat_2 <- as.data.frame(lm$V1) 
lg_dat_2$cM <- 0

lg_dat <- rbind(lg_dat_2, lg_dat_1)

lg_dat <- lg_dat[order(lg_dat$`lm$V1`),]

#now get scaffold order for plotting:
chrom_order_df <- unique(lm$scaffold) 
chrom_order_df

###••••••••#####
missing <- c("scaffold34", "scaffold37","scaffold38", "scaffold39", "scaffold6", "scaffold7")
fullchroms <- c(chrom_order_df, missing)
chromosomes <- as.data.frame(fullchroms)

for (i in 1:40){
  searchscaffold <- as.character(chromosomes$fullchroms[i])
  new_searchscaffold <- gsub("scaffold", "", searchscaffold)
  lengthsubset <- subset(lengthsdf, as.character(lengthsdf$Scaffold) == as.character(new_searchscaffold))
  chromosomes$bp[i] <- lengthsubset$V1
}

#newrename chromosomes
chromosomes$newchrom <- 1:40
chromosomes$newchrom <- paste("W", as.character(chromosomes$newchrom), sep = "")
for (i in 1:length(chromosomes$newchrom)){
  if (nchar(as.character(chromosomes$newchrom))[i] == 2){
    chromosomes$newchrom[i] <- as.character(gsub("W", "W0", as.character(chromosomes$newchrom)[i]))
  }
}

chromosomes$cm <- chromosomes$bp
renamedchromosomes <- chromosomes[,3:4]

totalbp <- sum(as.numeric(chromosomes$bp))
totalcm <- sum(lg_dat$cM)


#work out bp in cm
conversion <- totalcm/totalbp

chrom_dat_1 <- renamedchromosomes
chrom_dat_1$cm <- as.numeric(chromosomes$bp)*conversion

chrom_dat_2 <- chrom_dat_1
chrom_dat_2$cm <- 0

chrom_dat <- rbind(chrom_dat_2, chrom_dat_1)

chrom_dat <- chrom_dat[order(chrom_dat$newchrom),]

#now merge these two files:
#first make column the same
colnames(lg_dat) <- c("segment", "cM")
lg_dat$segment <- gsub("Calb", "c", lg_dat$segment)

colnames(chrom_dat) <- c("segment", "cM")
full_dat <- rbind (lg_dat, chrom_dat)

#now try circos plot, track outlines are given using 'complete' - dots can be added from this - then want to add links
library(circlize)
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

#now rename LG and scaffolds in samfile
linkfile <- samfile[,1:7]
#change bp to cM for plot
linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188*conversion

#find chromosomes conversion in df: chromosomes
for (i in 1:length(linkfile$X.PG)){
  oldchrom <- as.character(linkfile$X.2[i])
  oldlg <- as.character(linkfile$X)[i]
  newchrom_name_df <- subset(chromosomes, as.character(chromosomes$fullchroms) == as.character(oldchrom))
  linkfile$newchrom[i] <- newchrom_name_df$newchrom
  linkfile$newLG[i] <- gsub("SA", "c", oldlg)
}

circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.5, col = "black")
}
##################################################################
######################################
#load in chromosome lengths and linkage map lengths
lengthsdf <- read.csv(file = "FalconR2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
lm <- read.csv(file = "LM_stats.txt", header = FALSE, sep = " ")
lengthsdf$Scaffold <- 0:39

#load in samfile already separated columns in excel...
samfile <- read.csv(file = "SA_linkagemap_FalconR2PhaseMapped_Filtered.csv")
samfile <- subset(samfile, as.character(samfile$PN.bwa) == "PGA")

#make synteny map structure file with each linkage group/chromosome and start and end
#find equivalent chromosomes
for (i in 1:40){
  linkagegroup <- levels(lm$V1)[i]
  new_linkagegroup <- gsub("Calb", "SA", linkagegroup)
  sam_subset <- subset(samfile, as.character(samfile$X) == as.character(new_linkagegroup))
  scaffoldvect <- as.vector(sam_subset$X.2)
  abundant <- names(sort(summary(as.factor(scaffoldvect)), decreasing=T)[1:1])
  lm$scaffold[i] <- abundant
  lm$tot_markernumber[i] <- length(sam_subset$X.PG)
  abundant_number <- subset(sam_subset, as.character(sam_subset$X.2) == as.character(abundant))
  lm$scaf_markernumber[i] <- length(abundant_number$X.PG)
  lm$proportion_abundant <- ((lm$scaf_markernumber/lm$tot_markernumber)*100)
}

orderedbyscaffolds <- lm[order(lm$scaffold), ]

#now get lg data for plotting
lg_dat_1 <- as.data.frame(lm$V1)
lg_dat_1$cM <- lm$V3

lg_dat_2 <- as.data.frame(lm$V1) 
lg_dat_2$cM <- 0

lg_dat <- rbind(lg_dat_2, lg_dat_1)

lg_dat <- lg_dat[order(lg_dat$`lm$V1`),]

#now get scaffold order for plotting:
chrom_order_df <- unique(lm$scaffold) 
chrom_order_df
missing <- c("scaffold34", "scaffold37","scaffold38", "scaffold39", "scaffold6", "scaffold7")
fullchroms <- c(chrom_order_df, missing)
chromosomes <- as.data.frame(fullchroms)

for (i in 1:40){
  searchscaffold <- as.character(chromosomes$fullchroms[i])
  new_searchscaffold <- gsub("scaffold", "", searchscaffold)
  lengthsubset <- subset(lengthsdf, as.character(lengthsdf$Scaffold) == as.character(new_searchscaffold))
  chromosomes$bp[i] <- lengthsubset$V1
}

#newrename chromosomes
chromosomes$newchrom <- 40:1
chromosomes$newchrom <- paste("Z", as.character(chromosomes$newchrom), sep = "")
for (i in 1:length(chromosomes$newchrom)){
  if (nchar(as.character(chromosomes$newchrom))[i] == 2){
    chromosomes$newchrom[i] <- as.character(gsub("Z", "Z0", as.character(chromosomes$newchrom)[i]))
  }
}

chromosomes$cm <- chromosomes$bp
renamedchromosomes <- chromosomes[,3:4]

totalbp <- sum(as.numeric(chromosomes$bp))
totalcm <- sum(lg_dat$cM)


#work out bp in cm
conversion <- totalcm/totalbp

chrom_dat_1 <- renamedchromosomes
chrom_dat_1$cm <- as.numeric(chromosomes$bp)*conversion

chrom_dat_2 <- chrom_dat_1
chrom_dat_2$cm <- 0

chrom_dat <- rbind(chrom_dat_2, chrom_dat_1)

chrom_dat <- chrom_dat[order(chrom_dat$newchrom),]

#now merge these two files:
#first make column the same
colnames(lg_dat) <- c("segment", "cM")
lg_dat$segment <- gsub("Calb", "c", lg_dat$segment)


colnames(chrom_dat) <- c("segment", "cM")
chrom_dat$cM <- (chrom_dat$cM)-(2*(chrom_dat$cM))


full_dat <- rbind (lg_dat, chrom_dat)

#now try circos plot, track outlines are given using 'complete' - dots can be added from this - then want to add links
library(circlize)
circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })

#now rename LG and scaffolds in samfile
linkfile <- samfile[,1:8]
#change bp to cM for plot
linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188*conversion

linkfile$VN.0.7.17.r1188 <- linkfile$VN.0.7.17.r1188-(2*as.numeric(linkfile$VN.0.7.17.r1188))


#find chromosomes conversion in df: chromosomes
for (i in 1:length(linkfile$X.PG)){
  oldchrom <- as.character(linkfile$X.2[i])
  oldlg <- as.character(linkfile$X)[i]
  newchrom_name_df <- subset(chromosomes, as.character(chromosomes$fullchroms) == as.character(oldchrom))
  linkfile$newchrom[i] <- newchrom_name_df$newchrom
  linkfile$newLG[i] <- gsub("SA", "c", oldlg)
}

circos.clear()
#then set circos plot parameters
circos.par("track.height" = 0.05, start.degree=90, cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors = full_dat$segment, x = full_dat$cM)
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(5, "mm"), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })
#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h=0.7, col = add_transparency("darkgrey", transparency = 0.35))
  circos.lines(c(as.integer(b)), 1, as.character(a), col = "darkgrey", type = 'h')
  circos.lines(c(as.integer(d)), 1, as.character(e), col = "darkgrey", type = 'h')
}

#####


library(circlize)
circos.clear()
#name of file/parameters
tiff("FalconR2Hi_C.vs.LinkageMap.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# Customize plot parameters
par(fig=c(0,1,0,1),mar=c(1,1,1,1))
#create gaps vector
gapsnorm <- c(1)
gapswide <- c(1)
wfgaps <- rep(gapsnorm, 39)
salmgaps <- rep(gapswide, 39)
startgap <- c(9)
endgap <- c(9)
gapsall <- c(wfgaps, startgap, salmgaps, endgap)

#create names vector
Wnamevector <- c("LG1", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", 
                 "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", 
                 "LG21", "LG22", "LG23", "LG24", "LG25", "LG26", "LG27", "LG28", "LG29", "LG30", 
                 "LG31", "LG32", "LG33", "LG34", "LG35", "LG36", "LG37", "LG38", "LG39", "LG40")
Snamevector <- c("Chr40", "Chr39", "Chr38", "Chr37", "Chr36", "Chr35", "Chr34", "Chr33", "Chr32",
                 "Chr31", "Chr30", "Chr29", "Chr28", "Chr27", "Chr26", "Chr25", "Chr24", "Chr23", "Chr22", "Chr21", "Chr20", 
                 "Chr19", "Chr18", "Chr17", "Chr16", "Chr15", "Chr14", "Chr13", "Chr12", "Chr11", "Chr10", 
                 "Chr09", "Chr08", "Chr07", "Chr06", "Chr05", "Chr04", "Chr03", "Chr02", "Chr01")

namevector <- c(Wnamevector, Snamevector)

#then set circos plot parameters including the rotation, gaps between segments and gaps between tracks
circos.par("track.height" = 0.08 , start.degree=87, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall, track.margin = c(0.0045, 0.0045))
#initialize
circos.initialize(factors = full_dat$segment, x = full_dat$cM)

Wlet <- "W"
vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
#initialize first track - chromosomes/arms
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               if(grepl(Wlet, namevector[CELL_META$sector.numeric.index])){
                 circos.text(CELL_META$xcenter, 0.9 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               } else {
                 circos.text(CELL_META$xcenter, 1.45 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               }
             })

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  fact <- as.factor(linkfile$newLG)[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = colvec[fact])
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "red", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  }
  # if ((linkfile$CL.bwa)[row] < 60){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  # }
  # if ((linkfile$CL.bwa)[row] > 59){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  # }
}

legend("topright", c("30-39 - n=112", "40-49 - n=179", "50-59 - n=189", "60 - n=3168"), col = c("red", "orchid1", "lightblue", "royalblue1"), lwd = 1.5, cex = 0.55, title = "MapQ of 3648 links")

dev.off()

HQ <- subset(linkfile, linkfile$CL.bwa > 59)
HLQ <- subset(linkfile, linkfile$CL.bwa < 60)
LHQ <- subset(HLQ, HLQ$CL.bwa < 50)
LQ <- subset(LHQ, LHQ$CL.bwa < 40)
nrow(HQ)
nrow(HLQ)
nrow(LHQ)
nrow(LQ)

sum(nrow(HQ),nrow(HLQ),nrow(LHQ),nrow(LQ))
nrow(linkfile)



#--------------------------------
# plot synteny in loop
#--------------------------------
for (i in 1:1){
  circos.clear()
  #name of file/parameters
  filename <- paste("FalconHi_C.vs.LinkageMap_LG", as.character(i), ".tiff", sep = "")
  if (nchar(as.character(i)) < 2){
    LGNAME <- paste("c", "0", as.character(i), sep = "")
  }
  if (nchar(as.character(i)) > 1){
    LGNAME <- paste("c", as.character(i), sep = "")
  }
  #tiff(as.character(filename),height=8,width=8,units="in",res=300,compression="lzw")
  # Customize plot parameters
  par(fig=c(0,1,0,1),mar=c(1,1,1,1))
  #create gaps vector
  gapsnorm <- c(1)
  gapswide <- c(1)
  wfgaps <- rep(gapsnorm, 39)
  salmgaps <- rep(gapswide, 39)
  startgap <- c(9)
  endgap <- c(9)
  gapsall <- c(wfgaps, startgap, salmgaps, endgap)
  
  #create names vector
  Wnamevector <- c("LG1", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", 
                   "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", 
                   "LG21", "LG22", "LG23", "LG24", "LG25", "LG26", "LG27", "LG28", "LG29", "LG30", 
                   "LG31", "LG32", "LG33", "LG34", "LG35", "LG36", "LG37", "LG38", "LG39", "LG40")
  Snamevector <- c("Chr40", "Chr39", "Chr38", "Chr37", "Chr36", "Chr35", "Chr34", "Chr33", "Chr32",
                   "Chr31", "Chr30", "Chr29", "Chr28", "Chr27", "Chr26", "Chr25", "Chr24", "Chr23", "Chr22", "Chr21", "Chr20", 
                   "Chr19", "Chr18", "Chr17", "Chr16", "Chr15", "Chr14", "Chr13", "Chr12", "Chr11", "Chr10", 
                   "Chr09", "Chr08", "Chr07", "Chr06", "Chr05", "Chr04", "Chr03", "Chr02", "Chr01")
  
  namevector <- c(Wnamevector, Snamevector)
  
  #then set circos plot parameters including the rotation, gaps between segments and gaps between tracks
  circos.par("track.height" = 0.08 , start.degree=87, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall, track.margin = c(0.0045, 0.0045))
  #initialize
  circos.initialize(factors = full_dat$segment, x = full_dat$cM)
  
  Wlet <- "W"
  vec <- c(add_transparency("red", transparency = 0.8), 
           add_transparency("orange", transparency = 0.8), 
           add_transparency("green", transparency = 0.8), 
           add_transparency("blue", transparency = 0.8), 
           add_transparency("purple", transparency = 0.8))
  colvec <- rep(vec, 8)
  #initialize first track - chromosomes/arms
  circos.track(factors = full_dat$segment, ylim = c(0,1),
               panel.fun = function(x, y) {
                 if(grepl(Wlet, namevector[CELL_META$sector.numeric.index])){
                   circos.text(CELL_META$xcenter, 0.9 + uy(2, "mm"), 
                               namevector[CELL_META$sector.numeric.index], cex = 0.3)
                 } else {
                   circos.text(CELL_META$xcenter, 1.45 + uy(2, "mm"), 
                               namevector[CELL_META$sector.numeric.index], cex = 0.3)
                 }
               })
  
  #loop through linkage data frame and extract info to make links
  for (row in 1:(nrow(linkfile))){
    a <- linkfile$newLG[row]
    b <- linkfile$X.1[row]
    e <- linkfile$newchrom[row]
    d <- linkfile$VN.0.7.17.r1188[row]
    fact <- as.factor(linkfile$newLG)[row]
    #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
    #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "lightgrey")
    if ((linkfile$CL.bwa)[row] > 59){
      circos.lines(c(as.integer(b)), 1, as.character(a), col = "royalblue1", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 60){
      circos.lines(c(as.integer(b)), 1, as.character(a), col = "lightblue", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 50){
      circos.lines(c(as.integer(b)), 1, as.character(a), col = "orchid1", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 40){
      circos.lines(c(as.integer(b)), 1, as.character(a), col = "red", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] > 59){
      circos.lines(c(as.integer(d)), 1, as.character(e), col = "royalblue1", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 60){
      circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 50){
      circos.lines(c(as.integer(d)), 1, as.character(e), col = "orchid1", type = 'h')
    }
    if ((linkfile$CL.bwa)[row] < 40){
      circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
    }
    # if ((linkfile$CL.bwa)[row] < 60){
    #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
    # }
    # if ((linkfile$CL.bwa)[row] > 59){
    #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
    # }
    if (as.character(a) == as.character(LGNAME)){
      circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "red")
    }
  }  
  for (row in 1:(nrow(linkfile))){
    a <- linkfile$newLG[row]
    b <- linkfile$X.1[row]
    e <- linkfile$newchrom[row]
    d <- linkfile$VN.0.7.17.r1188[row]
    fact <- as.factor(linkfile$newLG)[row]
    #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
    #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))
    if (as.character(a) == as.character(LGNAME)){
      circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "red")
    }
  }
  legendtitle <- paste(as.character(gsub("c", "LG", LGNAME)), "MapQ of 3648 links", sep = " ")
  legend("topright", c("30-39 - n=112", "40-49 - n=179", "50-59 - n=189", "60 - n=3168"), col = c("red", "orchid1", "lightblue", "royalblue1"), lwd = 1.5, cex = 0.55, title = as.character(legendtitle))
  #dev.off()
}

#--------------------------------
#synteny plot no loop
#--------------------------------
circos.clear()
#name of file/parameters
#tiff("CanuHi_C.vs.LinkageMap.tiff",height=8,width=8,units="in",res=300,compression="lzw")
# Customize plot parameters
par(fig=c(0,1,0,1),mar=c(1,1,1,1))
#create gaps vector
gapsnorm <- c(1)
gapswide <- c(1)
wfgaps <- rep(gapsnorm, 39)
salmgaps <- rep(gapswide, 39)
startgap <- c(9)
endgap <- c(9)
gapsall <- c(wfgaps, startgap, salmgaps, endgap)

#create names vector
Wnamevector <- c("LG1", "LG02", "LG03", "LG04", "LG05", "LG06", "LG07", "LG08", "LG09", "LG10", 
                 "LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", 
                 "LG21", "LG22", "LG23", "LG24", "LG25", "LG26", "LG27", "LG28", "LG29", "LG30", 
                 "LG31", "LG32", "LG33", "LG34", "LG35", "LG36", "LG37", "LG38", "LG39", "LG40")
Snamevector <- c("Chr40", "Chr39", "Chr38", "Chr37", "Chr36", "Chr35", "Chr34", "Chr33", "Chr32",
                 "Chr31", "Chr30", "Chr29", "Chr28", "Chr27", "Chr26", "Chr25", "Chr24", "Chr23", "Chr22", "Chr21", "Chr20", 
                 "Chr19", "Chr18", "Chr17", "Chr16", "Chr15", "Chr14", "Chr13", "Chr12", "Chr11", "Chr10", 
                 "Chr09", "Chr08", "Chr07", "Chr06", "Chr05", "Chr04", "Chr03", "Chr02", "Chr01")

namevector <- c(Wnamevector, Snamevector)

#then set circos plot parameters including the rotation, gaps between segments and gaps between tracks
circos.par("track.height" = 0.08 , start.degree=87, cell.padding = c(0.005, 0, 0.005, 0), gap.degree = gapsall, track.margin = c(0.0045, 0.0045))
#initialize
circos.initialize(factors = full_dat$segment, x = full_dat$cM)

Wlet <- "W"
vec <- c(add_transparency("red", transparency = 0.8), 
         add_transparency("orange", transparency = 0.8), 
         add_transparency("green", transparency = 0.8), 
         add_transparency("blue", transparency = 0.8), 
         add_transparency("purple", transparency = 0.8))
colvec <- rep(vec, 8)
#initialize first track - chromosomes/arms
circos.track(factors = full_dat$segment, ylim = c(0,1),
             panel.fun = function(x, y) {
               if(grepl(Wlet, namevector[CELL_META$sector.numeric.index])){
                 circos.text(CELL_META$xcenter, 0.9 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               } else {
                 circos.text(CELL_META$xcenter, 1.45 + uy(2, "mm"), 
                             namevector[CELL_META$sector.numeric.index], cex = 0.3)
               }
             })

#loop through linkage data frame and extract info to make links
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  fact <- as.factor(linkfile$newLG)[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))
  circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "lightgrey")
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(b)), 1, as.character(a), col = "red", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] > 59){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "royalblue1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 60){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 50){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "orchid1", type = 'h')
  }
  if ((linkfile$CL.bwa)[row] < 40){
    circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  }
  # if ((linkfile$CL.bwa)[row] < 60){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "red", type = 'h')
  # }
  # if ((linkfile$CL.bwa)[row] > 59){
  #   circos.lines(c(as.integer(d)), 1, as.character(e), col = "lightblue", type = 'h')
  # }
  if (as.character(a) == "c01"){
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "red")
  }
}  
for (row in 1:(nrow(linkfile))){
  a <- linkfile$newLG[row]
  b <- linkfile$X.1[row]
  e <- linkfile$newchrom[row]
  d <- linkfile$VN.0.7.17.r1188[row]
  fact <- as.factor(linkfile$newLG)[row]
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c(as.integer(d)), h=0.9, col = "black")
  #circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = add_transparency("darkgrey", transparency = 0.35))
  if (as.character(a) == "c01"){
    circos.link(as.character(a), c(as.integer(b)), as.character(e), c((as.integer(d))), h.ratio = 0.7, col = "red")
  }
}
legend("topright", c("30-39 - n=112", "40-49 - n=179", "50-59 - n=189", "60 - n=3168"), col = c("red", "orchid1", "lightblue", "royalblue1"), lwd = 1.5, cex = 0.55, title = "MapQ of 3648 links")
# loop to write names and make LG name for plotting
for (i in 1:40){
  filename <- paste("FalconHi_C.vs.LinkageMap.tiff_LG", as.character(i), ".tiff", sep = "")
  if (nchar(as.character(i)) < 2){
    LGNAME <- paste("c", "0", as.character(i), sep = "")
  }
  if (nchar(as.character(i)) > 1){
    LGNAME <- paste("c", as.character(i), sep = "")
  }
}

#--------------------------------
#now plot recombination rate
#--------------------------------
dev.off()
old.par <- par(mar = c(1, 1, 1, 1))
par(old.par)

for (i in 1:40){
  LG_name <- levels(samfile$X)[i]
  LG_subset <- subset(samfile, as.character(samfile$X) == as.character(LG_name))
  abundant_chrom <- names(sort(summary(as.factor(LG_subset$X.2)), decreasing=T)[1:1])
  new_LG_subset <- subset(LG_subset, as.character(LG_subset$X.2) == as.character(abundant_chrom))
  title = paste(as.character(LG_name), "vs", as.character(abundant_chrom))
  lo <- loess(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1)
  plot(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1, main = title)
  abline(lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1), col = "green")
  lines(predict(lo), col='red', lwd=2)
}

abline(lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1), col = "green")
lm_stats <- lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1)
lm_stats$coefficients[2] < 0

?plot
LG01 <- subset(samfile, as.character(samfile$X) == as.character("SA01"))
abundant_chrom <- names(sort(summary(as.factor(LG01$X.2)), decreasing=T)[1:1])
new_LG01 <- subset(LG01, as.character(LG01$X.2) == as.character(abundant_chrom))
plot(new_LG01$X.1 ~ new_LG01$VN.0.7.17.r1188)

LG22 <- subset(samfile, as.character(samfile$X) == as.character("SA22"))
abundant_chrom <- names(sort(summary(as.factor(LG22$X.2)), decreasing=T)[1:1])
new_LG22 <- subset(LG22, as.character(LG22$X.2) == as.character(abundant_chrom))
plot(new_LG22$VN.0.7.17.r1188 ~ new_LG22$X.1)
new_LG22$X.1 <- (new_LG22$X.1)-(2*(new_LG22$X.1))

hist(samfile$CL.bwa)

#1/3 plot each linkage group, want bp on x axis and cm on y axis - plot ab line and determine if need to swap cM around
recomb_samfile <- samfile 
#make empty vector to contain r-estimated spar values for smooth.spline
sparvector <- c()

for (i in 1:40){
  LG_name <- levels(recomb_samfile$X)[i]
  LG_subset <- subset(recomb_samfile, as.character(recomb_samfile$X) == as.character(LG_name))
  abundant_chrom <- names(sort(summary(as.factor(LG_subset$X.2)), decreasing=T)[1:1])
  new_LG_subset <- subset(LG_subset, as.character(LG_subset$X.2) == as.character(abundant_chrom))
  title = paste(as.character(LG_name), "vs", as.character(abundant_chrom))
  #plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
  #abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
  lm_stats <- lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1)
  if (lm_stats$coefficients[2] < 0){
    new_LG_subset$X.1 <- ((max(new_LG_subset$X.1)+1)-new_LG_subset$X.1)
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1)
    sparvector[i] <- d$spar
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="p",col="blue",ylab="rec [cM/Mb]", main = title)
  }
  if (lm_stats$coefficients[2] >= 0 ){
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1)
    sparvector[i] <- d$spar
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="p",col="blue",ylab="rec [cM/Mb]", main = title)
  }
}

#2/3 new graph with predict smooth spline with mean value of spar determined from sparvector (==0.8)
for (i in 1:40){
  LG_name <- levels(recomb_samfile$X)[i]
  LG_subset <- subset(recomb_samfile, as.character(recomb_samfile$X) == as.character(LG_name))
  abundant_chrom <- names(sort(summary(as.factor(LG_subset$X.2)), decreasing=T)[1:1])
  new_LG_subset <- subset(LG_subset, as.character(LG_subset$X.2) == as.character(abundant_chrom))
  title = paste(as.character(LG_name), "vs", as.character(abundant_chrom))
  #plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
  #abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
  lm_stats <- lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1)
  if (lm_stats$coefficients[2] < 0){
    new_LG_subset$X.1 <- ((max(new_LG_subset$X.1)+1)-new_LG_subset$X.1)
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1, spar = mean(sparvector))
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="p",col="blue",ylab="rec [cM/Mb]", main = title)
  }
  if (lm_stats$coefficients[2] >= 0 ){
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1, spar = mean(sparvector))
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="p",col="blue",ylab="rec [cM/Mb]", main = title)
  }
}

#3/3and now plotting again with lines - need to first order each new_LG_subset by bp
#new graph with predict smooth spline with set value of spar determined from sparvector
for (i in 1:22){
  LG_name <- levels(recomb_samfile$X)[i]
  LG_subset <- subset(recomb_samfile, as.character(recomb_samfile$X) == as.character(LG_name))
  abundant_chrom <- names(sort(summary(as.factor(LG_subset$X.2)), decreasing=T)[1:1])
  new_LG_subset <- subset(LG_subset, as.character(LG_subset$X.2) == as.character(abundant_chrom))
  title = paste(as.character(LG_name), "vs", as.character(abundant_chrom))
  #plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
  #abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
  lm_stats <- lm(new_LG_subset$VN.0.7.17.r1188 ~ new_LG_subset$X.1)
  if (lm_stats$coefficients[2] < 0){
    new_LG_subset <- new_LG_subset[order(new_LG_subset$VN.0.7.17.r1188), ]
    new_LG_subset$X.1 <- ((max(new_LG_subset$X.1)+1)-new_LG_subset$X.1)
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1, spar = mean(sparvector))
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="l", lty = 2, col="blue",ylab="recombination rate [cM/Mb]", xlab = "bp", main = title)
    points(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec, col = "black", pch = 16)
  }
  if (lm_stats$coefficients[2] >= 0 ){
    new_LG_subset <- new_LG_subset[order(new_LG_subset$VN.0.7.17.r1188), ]
    d <- smooth.spline(new_LG_subset$VN.0.7.17.r1188 ,new_LG_subset$X.1, spar = mean(sparvector))
    plot(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188, main = title)
    abline(lm(new_LG_subset$X.1 ~ new_LG_subset$VN.0.7.17.r1188), col = "green")
    rec<-stats:::predict.smooth.spline(d,new_LG_subset$VN.0.7.17.r1188,deriv=1)$y*1000000
    rec[rec<0]<-0
    new_LG_subset$rec <- rec
    plot(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec,type="l", lty = 2, col="blue", ylab="recombination rate [cM/Mb]", xlab = "bp", main = title)
    points(new_LG_subset$VN.0.7.17.r1188, new_LG_subset$rec, col = "black", pch = 16)
  }
}

#--------------------------------
#now check coverage 
#--------------------------------

#load coverage data
cov_df <- read.csv("FalconR2_Phase30kb_windows_scaffolds.csv", header = FALSE)

mean(cov_df$V3)

max(cov_df$V3)
xmax <- max(cov_df$V2)

#and plot for each scaffold
# Create new column filled with default colour
cov_df$Colour="black"

#calulate mean and sd values
cov_mean <- mean(cov_df$V3)
#cov_test <- subset(cov_df, cov_df$V3 < 50)
cov_sd <- sd(cov_df$V3)
#cov_test_sd <- sd(cov_test$V3)
#cov_test_mean <- mean(cov_test$V3)

# Set new column values to appropriate colours
cov_df$Colour[cov_df$V3>=20]="red"
cov_df$Colour[cov_df$V3<=10]="blue"

cov_summary_low <- sum(cov_df$Colour == "blue")
cov_summary_mid <- sum(cov_df$Colour == "black")
cov_summary_high <- sum(cov_df$Colour == "red")
cov_summary_total <- cov_summary_low + cov_summary_mid + cov_summary_high
cov_total <- nrow(cov_df)

#plot coverage with windows coloured by coverage and a regression line
for (i in 0:39){
  scaffname <- paste("scaffold", as.character(i), "_", sep = "")
  cov_df_subset <- cov_df[grep(scaffname, cov_df$V1),]
  over40 <- sum(cov_df_subset$V3 > 40)
  under40 <- sum(cov_df_subset$V3 < 41)
  cov_title = paste(scaffname, " 0-40 = ", as.character(under40), "| 40+ = ", as.character(over40), "| total windows = ", as.character(nrow(cov_df_subset)),sep = "")
  plot(cov_df_subset$V3 ~ cov_df_subset$V2, ylim = c(0,40), main = cov_title, col=cov_df_subset$Colour, xlab = "bp", ylab = "coverage")
  abline(lm(cov_df_subset$V3 ~ cov_df_subset$V2), col = "green")
}

cov_df_noex <- subset(cov_df, cov_df$V3 < 60)
title_rows <- nrow(cov_df_noex)
hist(cov_df_noex$V3, main = paste("Falcon - Coverage < 60 ", " total = ", as.character(title_rows)))

cov_df_noex30 <- subset(cov_df, cov_df$V3 < 30)
title_rows30 <- nrow(cov_df_noex30)
hist(cov_df_noex30$V3, main = paste("Falcon - Coverage < 30 ", " total = ", as.character(title_rows30)))

cov_mean
cov_summary_low
cov_summary_mid
cov_summary_high
cov_total

#> cov_mean
#[1] 14.99976
#> cov_summary_low
#[1] 18141
#> cov_summary_mid
#[1] 54987
#> cov_summary_high
#[1] 6238
#> cov_total
#[1] 79366