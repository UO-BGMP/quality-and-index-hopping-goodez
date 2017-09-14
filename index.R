setwd("/Users/Zach/Desktop/qual_index_swap/Output_Real/")

### per-bp scores ###
#####################

R1_bp <- read.delim("R1_per_bp.tsv", sep = "\t", header = T)
R2_bp <- read.delim("R2_per_bp.tsv", sep = "\t", header = T)
R3_bp <- read.delim("R3_per_bp.tsv", sep = "\t", header = T)
R4_bp <- read.delim("R4_per_bp.tsv", sep = "\t", header = T)

par(mfrow=c(2,1))
plot(x=R1_bp$bp.position, y=R1_bp$mean.PHRED.score, type='l', pch=19,
     xlab="Position in read (bp)", ylab="Mean quality score", col="gold4", 
     main="Sequence Reads - Average quality score by position", ylim=c(30,40),
     cex.axis=0.8)
lines(R4_bp$bp.position, y=R4_bp$mean.PHRED.score, col="aquamarine4")
legend(75,32.5, c("R1 file","R4 file"), lty=c(1,1), col=c("gold4","aquamarine4"), cex=0.7,
       bty='n', yjust=0.5, y.intersp = 0.4)


plot(x=R2_bp$bp.position, y=R2_bp$mean.PHRED.score, type='l', pch=19,
     xlab="Position in read (bp)", ylab="Mean quality score", 
     main="Barcode Reads - Average quality score by position", ylim=c(30,40),
     col="deeppink3", cex.axis=0.8)
lines(R3_bp$bp.position, y=R3_bp$mean.PHRED.score, col="steelblue")
legend(6.5,32.5, c("R2 file","R3 file"), lty=c(1,1), col=c("deeppink3","steelblue"), cex=0.7,
       bty='n', yjust=0.5, y.intersp = 0.4)

### per-read scores ###
#######################
R1_read <- read.delim("R1_hist.tsv", sep = "\t", header=F)
R2_read <- read.delim("R2_hist.tsv", sep = "\t", header=F)
R3_read <- read.delim("R3_hist.tsv", sep = "\t", header=F)
R4_read <- read.delim("R4_hist.tsv", sep = "\t", header=F)


### Sequence reads
par(mfrow=c(1,2))
barplot(R1_read$V1, names.arg = R1_read$V2, cex.names = 0.8, xlab="Mean quality score of read",
        ylab="Frequency", ylim=c(0,2.5e+08), col="cyan", main="R1 file")
barplot(R4_read$V1, names.arg = R4_read$V2, cex.names = 0.8, xlab="Mean quality score of read",
        ylab="Frequency", ylim=c(0,2.5e+08), col="cyan", main="R4 file")

### Index reads
barplot(R2_read$V1, names.arg = R3_read$V2, cex.names = 0.8, xlab="Mean quality score of read",
        ylab="Frequency", ylim=c(0,2e+08), col="gold3", main="R2 file")
barplot(R3_read$V1, names.arg = R4_read$V2, cex.names = 0.8, xlab="Mean quality score of read",
        ylab="Frequency", ylim=c(0,2e+08), col="gold3", main="R3 file")

### demultiplex results ###
###########################

demulti <- read.delim("demultiplex.out.tsv", sep = "\t", header=T)

barplot(sort(demulti$PERCENT[1:24]), names.arg = demulti$ID[1:24], cex.names=0.75,
        ylab="Correct index pairs (% of all reads)", xlab="Index Pair", col="steelblue",
        main="Retained index pairs for each index")

############

demulti_summary <- data.frame(demulti[25:27,1:3])
summ_names= c("Error", "Index Hopped", "Matched")

summ_plot = barplot(demulti_summary$PERCENT, names.arg = summ_names, cex.names=0.8,
                    ylab="% of reads", xlab="", col="cyan",
                    main="Retained index pairs for each index")

summ_plot

midpoints <- summ_plot
text(midpoints, c(15,15,40), labels=c("6.5%","2.2%","91.3%"), 
     cex=0.7, col="darkred")


