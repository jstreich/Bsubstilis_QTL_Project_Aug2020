###############################################################################################################
# Bacterial Recombination Project; Collaboration between Jacobson Group and Michener Group
# Version: 1.0.5
# Created: 2020/06/02
# email: streich.jared@gmail.com, ju0@ornl.gov if still at ORNL
###############################################################################################################

###############################################################################################################
######################################### Load Packages #######################################################
###############################################################################################################

library(ecodist)
library(CircStats)
library(circlize)
library(mclust)
library(pracma)
library(TSA)
library(scatterplot3d)
library(rgl)
library(car)
library(plot3D)
library(ecodist)
library(vegan)
library(WaveletComp)
library(wmtsa)
library(scales)


###############################################################################################################
######################################## Set Color Palette ####################################################
###############################################################################################################

##### Colors
cols <- c("olivedrab3", "cadetblue2", "darkorange2", "orangered2", "aquamarine3", 
			"azure4", "gold3", "slateblue2", "pink3", "plum2")

##### Cols2
cols2 <- colorRampPalette(c("dodgerblue4", "white", "red2"))( 200 )

plot(c(1:length(cols)), col = cols, pch = 15, cex = 8)

###############################################################################################################
################################## Propreitary functions for script ###########################################
###############################################################################################################

##### Create function to convert normal ped file to haploid with only once cell per variant
pfix <- function(x){
		x <- x[,7:ncol(x)]
		y <- x[c(T,F)]
		return(y)
}




###############################################################################################################
###################################### Set Working Directory ##################################################
###############################################################################################################

setwd("/Users/ju0/Desktop/b.subtilis_recomb/2019-09_203_trmd_Samples/Bac_sub_ManuscriptPeds")

###############################################################################################################
######################################### Read Plot Data ######################################################
###############################################################################################################

##### 168E1_x_3A27K3 - Prototrophs 18x
map.data.1 <- read.table("168E1_x_3A27K3_Prototrophs_18x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.1 <- read.table("168E1_x_3A27K3_Prototrophs_18x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.1 <- pfix(ped.data.1)
dim(ped.data.1)


##### Run Wavelet on population 1 for recombination hotspot
ped.pop.sum <- colSums(ped.data.1)
ped.pop.sum <- rbind(colnames(ped.pop.sum), ped.pop.sum)
colnames(ped.data.1) <- map.data.1[,4]
ped.pop.sum[ped.pop.sum == 40] <- 0
par(mfrow = c(1,1))

##### Create empty vector of zeros to impute genomic features on
gnm <- rep(0, times = max(map.data.1[,4]))

## Impute GC marks into blank genome
for(i in 1:length(gnm)){
	pos <- map.data.1[i,4]
	gnm[pos] <- ped.pop.sum[i]
}
gnm <- gnm[complete.cases(gnm)]
dim(as.matrix(gnm))


##### Reduce the size of the vector by a power of 1000x, the function wvlt_bin_data2 is lower in the script
ped.sums <- wvlt_bin_data2(gnm, red  = 1000)
dim(as.matrix(ped.sums))
plot(ped.sums, type = "l", col = "olivedrab3")


# pdf(file = "pop1_snp_density_along_genome.pdf", width = 10, height = 4)
# note pdf plots do not seem to work with this plot type
plot(ped.sums, type = "l", col = "olivedrab3", xlab = "B. subtilis Genome, Binned Variants", 
	ylab = "Sum of Population Binned Markers", main = "Sum of Binned Population Markers of 168E1 x 3A27K3")
# dev.off()


wave.let <- wavCWT(ped.sums, wavelet = "gaussian2")
str(wave.let)
deltat(ped.sums) * c(1, length(ped.sums))
summary(wave.let)
png(file = "pop_1_wavelet_Variant_Density_binned_4kPos_July_2020.png", width = 3000, height = 1200)
plot(wave.let, col = cols2, xlab = "Genome Positions Scale 1:2", main = "Pop 1 Wavelet, Gaussian2 Wave")
dev.off()

plot(wave.let$)


set.seed(5362)
ped.sums.nv <- ped.sums
ped.sums.nv[3200:3400] <- sample(ped.sums, size = 201)
test.9 <- wavCWT(ped.sums.nv, wavelet = "gaussian2")


png(file = "pop_1_wavelet_Variant_Density_binned_4kPos_SelectionRemovedLoci3.2-3.4.png", width = 3000, height = 1200)
plot(test.9, col = cols2, xlab = "Genome Positions Scale 1:2", main = "Pop 1 Wavelet, Gaussian2 Wave")
dev.off()


pdf(file = "pop1_snp_density_along_genome_selectionRemovedloci3.2-3.4.pdf", width = 10, height = 4)
plot(ped.sums.nv, type = "l", col = "olivedrab3", xlab = "B. subtilis Genome, Binned Variants", 
	ylab = "Sum of Population Binned Markers", main = "Sum of Binned Population Markers of 168E1 x 3A27K3")
dev.off()


nrow(map.data.0)
test.7 <- pfix(ped.pop.sum)
plot(as.numeric(test.7))
test.7 <- pfix(test.7)
test.y <- wavCWT(test.7, wavelet = "gaussian2")

png(file = "pop_1_wavelet_Variant_Density.png", width = 3000, height = 1200)
plot(test.y, col = cols2, xlab = "Genome Positions Scale 1:2", main = "Pop 1 Wavelet, Gaussian2 Wave")
dev.off()



##### 3A27E1K3_x_2A11_EM_14x
map.data.2 <- read.table("3A27E1K3_x_2A11_EM_14x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.2 <- read.table("3A27E1K3_x_2A11_EM_14x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.2 <- pfix(ped.data.2)
ped.pop.sum <- colSums(ped.data.2)
ped.pop.sum <- rbind(colnames(ped.pop.sum), ped.pop.sum)
colnames(ped.data.2) <- map.data.2[,4]
ped.pop.sum[ped.pop.sum == 40] <- 0
colnames(ped.data.2) <- map.data.2[,4]

##### 3A27E1K3_x_28A5_HK_17x
map.data.3 <- read.table("3A27E1K3_x_28A5_HK_17x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.3 <- read.table("3A27E1K3_x_28A5_HK_17x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.3 <- pfix(ped.data.3)
colnames(ped.data.3) <- map.data.3[,4]

##### 168K3_x_3A27E1_18x_Double Resistant
map.data.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.4 <- pfix(ped.data.4)
colnames(ped.data.4) <- map.data.4[,4]

##### 3A27E1K3_x_3A1_EM_16x
map.data.5 <- read.table("3A27E1K3_x_3A1_EM_16x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.5 <- read.table("3A27E1K3_x_3A1_EM_16x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.5 <- pfix(ped.data.5)
colnames(ped.data.5) <- map.data.5[,4]

##### 168K3_x_3A27E1_Prototrophs_18x
map.data.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.6 <- pfix(ped.data.6)
colnames(ped.data.6) <- map.data.6[,4]


##### 3A27E1K3_x_2A11_HK_15x
map.data.7 <- read.table("3A27E1K3_x_2A11_HK_15x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.7 <- read.table("3A27E1K3_x_2A11_HK_15x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.7 <- pfix(ped.data.7)
colnames(ped.data.7) <- map.data.7[,4]
dim(map.data.7)

##### 3A27E1K3_x_28A5_EM_18x
map.data.8 <- read.table("3A27E1K3_x_28A5_EM_18x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.8 <- read.table("3A27E1K3_x_28A5_EM_18x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.8 <- pfix(ped.data.8)
colnames(ped.data.8) <- map.data.8[,4]

##### 168E1_x_3A27K3_Double_resistant_18x
map.data.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.9 <- pfix(ped.data.9)
colnames(ped.data.9) <- map.data.9[,4]

##### 3A27E1K3_x_3A1_HK_16x
map.data.0 <- read.table("3A27E1K3_x_3A1_HK_16x_bcf_DP20-GQ6.vcf_rcd012.map", header = F)
ped.data.0 <- read.table("3A27E1K3_x_3A1_HK_16x_bcf_DP20-GQ6.vcf_rcd012.ped", header = F)
ped.data.0 <- pfix(ped.data.0)
colnames(ped.data.0) <- map.data.0[,4]
dim(ped.data.0)



###############################################################################################################
#################################### Read Genome Information/Statistics #######################################
###############################################################################################################

##### Mean Depth Stats
depth.1 <- read.table("168E1_x_3A27K3_Prototrophs_18x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.2 <- read.table("3A27E1K3_x_2A11_EM_14x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.3 <- read.table("3A27E1K3_x_28A5_HK_17x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.5 <- read.table("3A27E1K3_x_3A1_EM_16x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.7 <- read.table("3A27E1K3_x_2A11_HK_15x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.8 <- read.table("3A27E1K3_x_28A5_EM_18x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
depth.0 <- read.table("3A27E1K3_x_3A1_HK_16x_bcf_DP20-GQ6.vcf_meanDepth.ldepth.mean", header = T)
head(depth.1)


##### Site Depth Stats
s.depth.1 <- read.table("168E1_x_3A27K3_Prototrophs_18x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.2 <- read.table("3A27E1K3_x_2A11_EM_14x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.3 <- read.table("3A27E1K3_x_28A5_HK_17x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.5 <- read.table("3A27E1K3_x_3A1_EM_16x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.7 <- read.table("3A27E1K3_x_2A11_HK_15x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.8 <- read.table("3A27E1K3_x_28A5_EM_18x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
s.depth.0 <- read.table("3A27E1K3_x_3A1_HK_16x_bcf_DP20-GQ6.vcf_siteDepth.ldepth", header = T)
head(s.depth.1)


##### Get seq depth for each individual in a shuffled population
seq.depth.1 <- read.table("168E1_x_3A27K3_Prototrophs_18x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.2 <- read.table("3A27E1K3_x_2A11_EM_14x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.3 <- read.table("3A27E1K3_x_28A5_HK_17x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.5 <- read.table("3A27E1K3_x_3A1_EM_16x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.7 <- read.table("3A27E1K3_x_2A11_HK_15x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.8 <- read.table("3A27E1K3_x_28A5_EM_18x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)
seq.depth.0 <- read.table("3A27E1K3_x_3A1_HK_16x_bcf_DP20-GQ6.vcf_depth.idepth", header = T)




###############################################################################################################
######################################## Get Recombination Length #############################################
###############################################################################################################

##### Process and plot ironed out markers between parents, binary to each parent A or B.
## This recodes variants to be parent A as '0', or as parent B as '2'
# par(mfrow = c(5,1))

recode_by_parent <- function(x, clm){
	x <- rbind(c(1:ncol(x)), x)
	y <- x[, x[(clm)+1,] == 0]
	cls2fx <- as.numeric(y[1,])
	for(i in 1:length(cls2fx)){
		j <- as.numeric(cls2fx[i])
		rs.1 <- x[,j]
		rs.1[rs.1 == 2] <- 1
		rs.1[rs.1 == 0] <- 2
		rs.1[clm + 1] <- 2
		rs.1[rs.1 == 1] <- 0
		x[,j] <- rs.1
	}
	return(x[2:nrow(x), ])
	}


##### Pop1
rowSums(ped.data.1)/2
p.1 <- recode_by_parent(ped.data.1, 1)
write.table(p.1, "168E1_x_3A27K3_Prototrophs_18x.txt", quote = F, col.names = T, row.names = T)
p.1 <- read.delim("168E1_x_3A27K3_Prototrophs_18x.txt", header = T, sep = " ")
clmn <- gsub("X","",colnames(p.1))
colnames(p.1) <- clmn
p.1 <- p.1[!is.na(colnames(p.1)), ]

##### Pop2
rowSums(ped.data.2)/2
p.2 <- recode_by_parent(ped.data.2, 1)
dim(p.2)
write.table(p.2, "3A27E1K3_x_2A11_HK_15x.txt", quote = F, col.names = T, row.names = T)
p.2 <- read.table("3A27E1K3_x_2A11_HK_15x.txt")
colnames(p.2) <- map.data.2[,4]

##### Pop3
rowSums(ped.data.3)
p.3 <- recode_by_parent(ped.data.3, 19)
write.table(p.3, "3A27E1K3_x_28A5_HK_17x.txt", col.names = T, row.names = T)
p.3 <- read.table("3A27E1K3_x_28A5_HK_17x.txt")
colnames(p.3) <- map.data.3[,4]

##### Pop4
rowSums(ped.data.4)
p.4 <- recode_by_parent(ped.data.4, 1)
write.table(p.4, "168K3_x_3A27E1_Double_resistant_18x.txt", quote = F, col.names = T, row.names = T)
p.4 <- read.table("168K3_x_3A27E1_Double_resistant_18x.txt")
colnames(p.4) <- map.data.4[,4]

##### Pop5
p.5 <- recode_by_parent(ped.data.5, 1)
write.table(p.5, "3A27E1K3_x_3A1_EM_16x.txt", quote = F, col.names = T, row.names = T)
p.5 <- read.table("3A27E1K3_x_3A1_EM_16x.txt")
colnames(p.5) <- map.data.5[,4]

##### Pop6 *Note, this population failed in the wetlab
# p.6 <- recode_by_parent(ped.data.6, 1)
# write.table(p.6, "168K3_x_3A27E1_Prototrophs_18x.txt", quote = F, col.names = T, row.names = T)
# p.6 <- read.table("168K3_x_3A27E1_Prototrophs_18x.txt")
# colnames(p.6) <- map.data.6[,4]


##### Pop7
p.7 <- recode_by_parent(ped.data.7, 1)
write.table(p.7, "3A27E1K3_x_2A11_HK_15x.txt", quote = F, col.names = T, row.names = T)
p.7 <- read.table("3A27E1K3_x_2A11_HK_15x.txt")
colnames(p.7) <- map.data.7[,4]

##### Pop8
rowSums(ped.data.8)
p.8 <- recode_by_parent(ped.data.8, 20)
write.table(p.8, "3A27E1K3_x_28A5_EM_18x.txt", col.names = T, row.names = T, quote = F)
p.8 <- read.table("3A27E1K3_x_28A5_EM_18x.txt")
colnames(p.8) <- map.data.8[,4]

##### Pop9
p.9 <- recode_by_parent(ped.data.9, 1)
write.table(p.9, "168E1_x_3A27K3_Double_resistant_18x.txt", quote = F, col.names = T, row.names = T)
p.9 <- read.table("168E1_x_3A27K3_Double_resistant_18x.txt")
colnames(p.9) <- map.data.9[,4]

##### Pop10
rowSums(ped.data.0)
p.0 <- recode_by_parent(ped.data.0, 1)
write.table(p.0, "3A27E1K3_x_3A1_HK_16x.txt", quote = F, col.names = T, row.names = T)
p.0 <- read.table("3A27E1K3_x_3A1_HK_16x.txt")
colnames(p.0) <- map.data.0[,4]

###############################################################################################################
############################## Get marker distances for each shuffling population #############################
###############################################################################################################
p.var.1 <- (map.data.1[2:nrow(map.data.1),4]) - (map.data.1[(1:nrow(map.data.1)-1),4])
p.var.2 <- (map.data.2[2:nrow(map.data.2),4]) - (map.data.2[(1:nrow(map.data.2)-1),4])
p.var.3 <- (map.data.3[2:nrow(map.data.3),4]) - (map.data.3[(1:nrow(map.data.3)-1),4])
p.var.4 <- (map.data.4[2:nrow(map.data.4),4]) - (map.data.4[(1:nrow(map.data.4)-1),4])
p.var.5 <- (map.data.5[2:nrow(map.data.5),4]) - (map.data.5[(1:nrow(map.data.5)-1),4])
p.var.6 <- (map.data.6[2:nrow(map.data.6),4]) - (map.data.6[(1:nrow(map.data.6)-1),4])
p.var.7 <- (map.data.7[2:nrow(map.data.7),4]) - (map.data.7[(1:nrow(map.data.7)-1),4])
p.var.8 <- (map.data.8[2:nrow(map.data.8),4]) - (map.data.8[(1:nrow(map.data.8)-1),4])
p.var.9 <- (map.data.9[2:nrow(map.data.9),4]) - (map.data.9[(1:nrow(map.data.9)-1),4])
p.var.0 <- (map.data.0[2:nrow(map.data.0),4]) - (map.data.0[(1:nrow(map.data.0)-1),4])


#######################################################################################################
################################ Get length of recombination block per indv. ##########################
#######################################################################################################

##### Create Insert length and start stop Function
blocklength_v3 <- function(z){
	rnms <- rownames(z)
	for(l in 1:nrow(z)){
		bp <- 1
		fp <- 1
		x <- as.numeric(z[l,])
		while(fp < length(x)){
			snps <- 0
			while(x[fp] != 2){
				fp <- fp + 1
				if(fp >= length(x)){
					break
				}
			}
			strt <- as.numeric(colnames(z)[fp])
			bp <- fp
			while(x[bp] != 0){
				bp <- bp + 1
				snps <- snps + 1
				if(bp >= length(x)){
					break
				}	
			}
			fp <- bp
			end <- as.numeric(colnames(z)[bp])
			rw <- c(rnms[l], strt, end, (end - strt), snps)
			if(!exists("rw.p")){
				rw.p <- rw
			}
			else{
				rw.p <- rbind(rw.p, rw)
			}
		}
	}
	rw.p <- rw.p[rw.p[,5] > 0, ]
	return(rw.p)
	}


##### Run each population through insertion length function
##### Pop1
testout.1.1 <- blocklength_v3(p.1[2:20,])
head(testout.1.1)
write.table(testout.1.1, "recomb_block_pop1_v3.txt", row.names = T, col.names = F, quote = F)
testout.1 <- read.table("recomb_block_pop1_v3.txt")
head(testout.1)

##### Pop2
p.2[1:20,1:20]
testout.2 <- blocklength_v3(p.2[2:16,])
write.table(testout.2, "recomb_block_pop2_v3.txt", row.names = T, col.names = F, quote = F)
testout.2 <- read.table("recomb_block_pop2_v3.txt")
head(testout.2)

##### Pop3
testout.3 <- blocklength_v3(p.3[2:18,])
write.table(testout.3, "recomb_block_pop3_v3.txt", row.names = T, col.names = F, quote = F)
testout.3 <- read.table("recomb_block_pop3_v3.txt")
testout.3 <- cbind(rep("rw", times = nrow(testout.3)), testout.3)

##### Pop4
p.4[1:20,1:20]
testout.4 <- blocklength_v3(p.4[2:20,])
write.table(testout.4, "recomb_block_pop4_v3.txt", row.names = T, col.names = F, quote = F)
testout.4 <- read.table("recomb_block_pop4_v3.txt")

##### Pop5
p.5[1:20,1:20]
testout.5 <- blocklength_v3(p.5[1:17,])
write.table(testout.5, "recomb_block_pop5_v3.txt", row.names = T, col.names = F, quote = F)
testout.5 <- read.table("recomb_block_pop5_v3.txt")

##### Pop6
# rowSums(p.6)
# p.6[1:21,1:20]
# testout.6 <- blocklength_v3(p.6[2:20,])
# write.table(testout.6, "recomb_block_pop6_v3.txt", row.names = T, col.names = F)
# testout.6 <- read.table("recomb_block_pop6_v3.txt")

##### Pop7
p.7[1:20,1:20]
testout.7 <- blocklength_v3(p.7[2:nrow(p.7),])
write.table(testout.7, "recomb_block_pop7_v3.txt", row.names = F, col.names = F, quote = F)
testout.7 <- read.table("recomb_block_pop7_v3.txt")
testout.7 <- testout.7[complete.cases(testout.7[,5]), ]
testout.7 <- cbind(rep("rw", times = nrow(testout.7)), testout.7)

##### Remove flagged samples from pop7
testout.7 <- testout.7[testout.7[,2] != 18, ]
testout.7 <- testout.7[testout.7[,2] != 14, ]
testout.7 <- testout.7[testout.7[,2] != 10, ]
testout.7 <- testout.7[testout.7[,2] != 9, ]
testout.7 <- testout.7[testout.7[,2] != 7, ]
# testout.7 <- testout.7[testout.7[,2] != 11,] # 

##### Pop8
testout.8 <- blocklength_v3(p.8[1:19,])
write.table(testout.8, "recomb_block_pop8_v3.txt", row.names = T, col.names = F, quote = F)
testout.8 <- read.table("recomb_block_pop8_v3.txt")
testout.8 <- cbind(rep("rw", times = nrow(testout.8)), testout.8)

##### Pop9
rowSums(p.9)/2
testout.9 <- blocklength_v3(p.9[2:20,])
write.table(testout.9, "recomb_block_pop9_v3.txt", row.names = F, col.names = F)
testout.9 <- read.table("recomb_block_pop9_v3.txt")
testout.9 <- testout.9[complete.cases(testout.9[,5]), ]
testout.9 <- cbind(rep("rw", times = nrow(testout.9)), testout.9)

##### Pop10
testout.0 <- blocklength_v3(p.0[2:18,])
write.table(testout.0, "recomb_block_pop10_v3.txt", row.names = T, col.names = F, quote = F)
testout.0 <- read.table("recomb_block_pop10_v3.txt")
testout.0 <- cbind(rep("rw", times = nrow(testout.0)), testout.0)

##### Add column names to insert length files
colnames(testout.1) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.2) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.3) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.4) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.5) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.6) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.7) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.8) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.9) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")
colnames(testout.0) <- c("row", "indv", "start", "end", "insertSize", "nSNPs")


###############################################################################################################
################################## Remove mutations between reference and parent ##############################
###############################################################################################################
##### Filter out dominant markers from pop1

table.5.3.1 <- table(testout.1[,2], paste(testout.1[,3], "_", testout.1[,4], sep = ""))
table.5.3.1[table.5.3.1 == 2] <- 1
crop.marks.1 <- colSums(table.5.3.1)
crop.marks.1 <- crop.marks.1[crop.marks.1 >= (max(crop.marks.1)-3)]
crop.marks.1 <- t(matrix(unlist(strsplit(names(crop.marks.1), split = "_")), nrow = 2))
test.marks.1.nmatch <- !match(testout.1[,3], crop.marks.1[,1])
test.marks.1.nmatch[is.na(test.marks.1.nmatch)] <- 1
test.marks.1.nmatch <- cbind(test.marks.1.nmatch, testout.1)
test.marks.1.nmatch <- test.marks.1.nmatch[test.marks.1.nmatch[,1] > 0, ]
test.marks.1.nmatch <- test.marks.1.nmatch[,2:ncol(test.marks.1.nmatch)]
dim(test.marks.1.nmatch)
mean(test.marks.1.nmatch[,5])
sd(test.marks.1.nmatch[,5])

write.table(test.marks.1.nmatch, "recomb_block_pop1_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop2
table.5.3.2 <- table(testout.2[,2], paste(testout.2[,3], "_", testout.2[,4], sep = ""))
table.5.3.2[table.5.3.2 == 2] <- 1
crop.marks.2 <- colSums(table.5.3.2)
crop.marks.2 <- crop.marks.2[crop.marks.2 >= (max(crop.marks.2)-3)]
crop.marks.2 <- t(matrix(unlist(strsplit(names(crop.marks.2), split = "_")), nrow = 2))
test.marks.2.nmatch <- !match(testout.2[,3], crop.marks.2[,1])
test.marks.2.nmatch[is.na(test.marks.2.nmatch)] <- 1
test.marks.2.nmatch <- cbind(test.marks.2.nmatch, testout.2)
test.marks.2.nmatch <- test.marks.2.nmatch[test.marks.2.nmatch[,1] > 0, ]
test.marks.2.nmatch <- test.marks.2.nmatch[,2:ncol(test.marks.2.nmatch)]
dim(test.marks.2.nmatch)
mean(test.marks.2.nmatch[,5])
sd(test.marks.2.nmatch[,5])

write.table(test.marks.2.nmatch, "recomb_block_pop2_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop3
table.5.3.3 <- table(testout.3[,2], paste(testout.3[,3], "_", testout.3[,4], sep = ""))
table.5.3.3[table.5.3.3 == 2] <- 1
crop.marks.3 <- colSums(table.5.3.3)
crop.marks.3 <- crop.marks.3[crop.marks.3 >= (max(crop.marks.3)-3)]
crop.marks.3 <- t(matrix(unlist(strsplit(names(crop.marks.3), split = "_")), nrow = 2))
test.marks.3.nmatch <- !match(testout.3[,3], crop.marks.3[,1])
test.marks.3.nmatch[is.na(test.marks.3.nmatch)] <- 1
test.marks.3.nmatch <- cbind(test.marks.3.nmatch, testout.3)
test.marks.3.nmatch <- test.marks.3.nmatch[test.marks.3.nmatch[,1] > 0, ]
test.marks.3.nmatch <- test.marks.3.nmatch[,2:ncol(test.marks.3.nmatch)]
dim(test.marks.3.nmatch)
mean(test.marks.3.nmatch[,5])
sd(test.marks.3.nmatch[,5])

write.table(test.marks.3.nmatch, "recomb_block_pop3_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop4
table.5.3.4 <- table(testout.4[,2], paste(testout.4[,3], "_", testout.4[,4], sep = ""))
table.5.3.4[table.5.3.4 == 2] <- 1
crop.marks.4 <- colSums(table.5.3.4)
crop.marks.4 <- crop.marks.4[crop.marks.4 >= (max(crop.marks.4)-3)]
crop.marks.4 <- t(matrix(unlist(strsplit(names(crop.marks.4), split = "_")), nrow = 2))
test.marks.4.nmatch <- !match(testout.4[,3], crop.marks.4[,1])
test.marks.4.nmatch[is.na(test.marks.4.nmatch)] <- 1
test.marks.4.nmatch <- cbind(test.marks.4.nmatch, testout.4)
test.marks.4.nmatch <- test.marks.4.nmatch[test.marks.4.nmatch[,1] > 0, ]
test.marks.4.nmatch <- test.marks.4.nmatch[,2:ncol(test.marks.4.nmatch)]
dim(test.marks.4.nmatch)
mean(test.marks.4.nmatch[,5])
sd(test.marks.4.nmatch[,5])

write.table(test.marks.4.nmatch, "recomb_block_pop4_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop5
table.5.3.5 <- table(testout.5[,2], paste(testout.5[,3], "_", testout.5[,4], sep = ""))
table.5.3.5[table.5.3.5 == 2] <- 1
crop.marks.5 <- colSums(table.5.3.5)
crop.marks.5 <- crop.marks.5[crop.marks.5 >= (max(crop.marks.5)-3)]
crop.marks.5 <- t(matrix(unlist(strsplit(names(crop.marks.5), split = "_")), nrow = 2))
test.marks.5.nmatch <- !match(testout.5[,3], crop.marks.5[,1])
test.marks.5.nmatch[is.na(test.marks.5.nmatch)] <- 1
test.marks.5.nmatch <- cbind(test.marks.5.nmatch, testout.5)
test.marks.5.nmatch <- test.marks.5.nmatch[test.marks.5.nmatch[,1] > 0, ]
test.marks.5.nmatch <- test.marks.5.nmatch[,2:ncol(test.marks.5.nmatch)]
dim(test.marks.5.nmatch)
mean(test.marks.5.nmatch[,5])
sd(test.marks.5.nmatch[,5])

write.table(test.marks.5.nmatch, "recomb_block_pop5_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop6
# table.5.3.6 <- table(testout.6[,2], paste(testout.6[,3], "_", testout.6[,4], sep = ""))
# table.5.3.6[table.5.3.6 == 2] <- 1
# crop.marks.6 <- colSums(table.5.3.6)
# crop.marks.6 <- crop.marks.6[crop.marks.6 >= (max(crop.marks.6)-3)]
# crop.marks.6 <- t(matrix(unlist(strsplit(names(crop.marks.6), split = "_")), nrow = 2))
# test.marks.6.nmatch <- !match(testout.6[,3], crop.marks.6[,1])
# test.marks.6.nmatch[is.na(test.marks.6.nmatch)] <- 1
# test.marks.6.nmatch <- cbind(test.marks.6.nmatch, testout.6)
# test.marks.6.nmatch <- test.marks.6.nmatch[test.marks.6.nmatch[,1] > 0, ]
# test.marks.6.nmatch <- test.marks.6.nmatch[,2:ncol(test.marks.6.nmatch)]
# dim(test.marks.6.nmatch)
# mean(test.marks.6.nmatch[,5])
# sd(test.marks.6.nmatch[,5])
# write.table(test.marks.6.nmatch, "recomb_block_pop6_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop7
table.5.3.7 <- table(testout.7[,2], paste(testout.7[,3], "_", testout.7[,4], sep = ""))
image(t(table.5.3.7))
table.5.3.7[table.5.3.7 == 2] <- 1
crop.marks.7 <- colSums(table.5.3.7)
max(crop.marks.7)
crop.marks.7 <- crop.marks.7[crop.marks.7 >= (max(crop.marks.7)-3)]
crop.marks.7 <- t(matrix(unlist(strsplit(names(crop.marks.7), split = "_")), nrow = 2))
test.marks.7.nmatch <- !match(as.numeric(testout.7[,4]), as.numeric(crop.marks.7[,2]))
test.marks.7.nmatch[is.na(test.marks.7.nmatch)] <- 1
test.marks.7.nmatch <- cbind(test.marks.7.nmatch, testout.7)
test.marks.7.nmatch <- test.marks.7.nmatch[test.marks.7.nmatch[,1] > 0, ]
test.marks.7.nmatch <- test.marks.7.nmatch[,2:ncol(test.marks.7.nmatch)]
dim(test.marks.7.nmatch)
mean(test.marks.7.nmatch[,5])
sd(test.marks.7.nmatch[,5])

write.table(test.marks.7.nmatch, "recomb_block_pop7_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop8
table.5.8 <- table(testout.8[,2], paste(testout.8[,3], "_", testout.8[,4], sep = ""))
table.5.8[table.5.8 == 2] <- 1
crop.marks.8 <- colSums(table.5.8)
crop.marks.8 <- crop.marks.8[crop.marks.8 >= (max(crop.marks.8)-3)]
crop.marks.8 <- t(matrix(unlist(strsplit(names(crop.marks.8), split = "_")), nrow = 2))
test.marks.8.nmatch <- !match(testout.8[,3], crop.marks.8[,1])
test.marks.8.nmatch[is.na(test.marks.8.nmatch)] <- 1
test.marks.8.nmatch <- cbind(test.marks.8.nmatch, testout.8)
test.marks.8.nmatch <- test.marks.8.nmatch[test.marks.8.nmatch[,1] > 0, ]
test.marks.8.nmatch <- test.marks.8.nmatch[,2:ncol(test.marks.8.nmatch)]
dim(test.marks.8.nmatch)
mean(test.marks.8.nmatch[,5])
sd(test.marks.8.nmatch[,5])

write.table(test.marks.8.nmatch, "recomb_block_pop8_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop9
table.5.3.9 <- table(testout.9[,2], paste(testout.9[,3], "_", testout.9[,4], sep = ""))
table.5.3.9[table.5.3.9 == 2] <- 1
crop.marks.9 <- colSums(table.5.3.9)
crop.marks.9 <- crop.marks.9[crop.marks.9 >= (max(crop.marks.9)-3)]
crop.marks.9 <- t(matrix(unlist(strsplit(names(crop.marks.9), split = "_")), nrow = 2))
test.marks.9.nmatch <- !match(testout.9[,3], crop.marks.9[,1])
test.marks.9.nmatch[is.na(test.marks.9.nmatch)] <- 1
test.marks.9.nmatch <- cbind(test.marks.9.nmatch, testout.9)
test.marks.9.nmatch <- test.marks.9.nmatch[test.marks.9.nmatch[,1] > 0, ]
test.marks.9.nmatch <- test.marks.9.nmatch[,2:ncol(test.marks.9.nmatch)]
dim(test.marks.9.nmatch)
mean(test.marks.9.nmatch[,5])
sd(test.marks.9.nmatch[,5])

write.table(test.marks.9.nmatch, "recomb_block_pop9_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Filter out dominant markers from pop10 aka 0
table.5.3.0 <- table(testout.0[,2], paste(testout.0[,3], "_", testout.0[,4], sep = ""))
table.5.3.0[table.5.3.0 == 2] <- 1
crop.marks.0 <- colSums(table.5.3.0)
crop.marks.0 <- crop.marks.0[crop.marks.0 >= (max(crop.marks.0)-3)]
crop.marks.0 <- t(matrix(unlist(strsplit(names(crop.marks.0), split = "_")), nrow = 2))
test.marks.0.nmatch <- !match(testout.0[,3], crop.marks.0[,1])
test.marks.0.nmatch[is.na(test.marks.0.nmatch)] <- 1
test.marks.0.nmatch <- cbind(test.marks.0.nmatch, testout.0)
test.marks.0.nmatch <- test.marks.0.nmatch[test.marks.0.nmatch[,1] > 0, ]
test.marks.0.nmatch <- test.marks.0.nmatch[,2:ncol(test.marks.0.nmatch)]
test.marks.0.nmatch <- test.marks.0.nmatch[complete.cases(test.marks.0.nmatch[,5]),]
dim(test.marks.0.nmatch)
mean(test.marks.0.nmatch[,5])
sd(test.marks.0.nmatch[,5])

write.table(test.marks.0.nmatch, "recomb_block_pop0_v3_copy_Check.txt", sep = "	", col.names = T, row.names = T, quote = F)

##### Create Mixed Population of 1, 4, 9
max(test.marks.1.nmatch[,2])
max(test.marks.4.nmatch[,2])
pop.4.t <- test.marks.4.nmatch
pop.9.t <- test.marks.9.nmatch
pop.4.t[,2] <- pop.4.t[,2] + max(test.marks.1.nmatch[,2])
pop.9.t[,2] <- pop.9.t[,2] + max(pop.4.t[,2])


pop149 <- rbind(test.marks.1.nmatch, pop.4.t, pop.9.t)



###############################################################################################################
###################### Create historgrams of insert size by frequency across populations ######################
###############################################################################################################
hist.1 <- test.marks.1.nmatch[!duplicated(paste(test.marks.1.nmatch[,2], "_", test.marks.1.nmatch[,3], sep = "")), ]
hist.2 <- test.marks.2.nmatch[!duplicated(paste(test.marks.2.nmatch[,2], "_", test.marks.2.nmatch[,3], sep = "")), ]
hist.3 <- test.marks.3.nmatch[!duplicated(paste(test.marks.3.nmatch[,2], "_", test.marks.3.nmatch[,3], sep = "")), ]
hist.4 <- test.marks.4.nmatch[!duplicated(paste(test.marks.4.nmatch[,2], "_", test.marks.4.nmatch[,3], sep = "")), ]
hist.5 <- test.marks.5.nmatch[!duplicated(paste(test.marks.5.nmatch[,2], "_", test.marks.5.nmatch[,3], sep = "")), ]
hist.6 <- test.marks.6.nmatch[!duplicated(paste(test.marks.6.nmatch[,2], "_", test.marks.6.nmatch[,3], sep = "")), ]
hist.7 <- test.marks.7.nmatch[!duplicated(paste(test.marks.7.nmatch[,2], "_", test.marks.7.nmatch[,3], sep = "")), ]
hist.8 <- test.marks.8.nmatch[!duplicated(paste(test.marks.8.nmatch[,2], "_", test.marks.8.nmatch[,3], sep = "")), ]
hist.9 <- test.marks.9.nmatch[!duplicated(paste(test.marks.9.nmatch[,2], "_", test.marks.9.nmatch[,3], sep = "")), ]
hist.0 <- test.marks.0.nmatch[!duplicated(paste(test.marks.0.nmatch[,2], "_", test.marks.0.nmatch[,3], sep = "")), ]
hist.149 <- pop.9.t[!duplicated(paste(pop.9.t[,2], "_", pop.9.t[,3], sep = "")), ]


####################################################################################################################################
################################################# Plot jitter and density plots ####################################################
####################################################################################################################################

##### Create jitter plot of insertion sizes
par(mfrow = c(1,1))
max.b.t <- max(c(max(hist.1[,5]), max(hist.2[,5]), max(hist.3[,5]), max(hist.4[,5]), max(hist.5[,5]),
			max(hist.6[,5]), max(hist.7[,5]), max(hist.8[,5]), max(hist.9[,5]), max(hist.0[,5])))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))
# Start plot
pdf("10_popsBsub_jitterplots.png", height = 800, width = 1300)
plot(jitter(rep(0,length(hist.1[,5])),amount=0.2), hist.1[,5],
     xlim=range(-0.5,9.5), ylim=range(-3,max(max.b.t)), main = "Insertion Size of 10 Populations",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[1], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
abline(h = as.integer(tpl), lwd = 0.8, col = "grey40", lty = 2)
abline(h = (as.integer(tpl)/2), lwd = 0.8, col = "grey40", lty = 3)
points(jitter(rep(1,length(hist.2[,5])), amount=0.2), hist.2[,5], cex = 0.75, col = cols[2], pch = 16)
points(jitter(rep(2,length(hist.3[,5])), amount=0.2), hist.3[,5], cex = 0.75, col = cols[3], pch = 16)
points(jitter(rep(3,length(hist.4[,5])), amount=0.2), hist.4[,5], cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(4,length(hist.5[,5])), amount=0.2), hist.5[,5], cex = 0.75, col = cols[5], pch = 16)
points(jitter(rep(5,length(hist.6[,5])), amount=0.2), hist.6[,5], cex = 0.75, col = cols[6], pch = 16)
points(jitter(rep(6,length(hist.7[,5])), amount=0.2), hist.7[,5], cex = 0.75, col = cols[7], pch = 16)
points(jitter(rep(7,length(hist.8[,5])), amount=0.2), hist.8[,5], cex = 0.75, col = cols[8], pch = 16)
points(jitter(rep(8,length(hist.9[,5])), amount=0.2), hist.9[,5], cex = 0.75, col = cols[9], pch = 16)
points(jitter(rep(9,length(hist.0[,5])), amount=0.2), hist.0[,5], cex = 0.75, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()


##### Only 3a27 x 168 populations
par(mfrow = c(1,1))
max.b.t <- max(c(log10(max(hist.1[,5])), log10(max(hist.4[,5])), log10(max(hist.9[,5]))))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))

pdf("ThreePops_3a27x168_InsertSizeBsub_jitterplots_log.pdf", height = 8, width = 13)
plot(jitter(rep(0,length(hist.1[,5])),amount=0.2), log10(hist.1[,5]),
     xlim=range(-0.5,2.5), ylim=range(-3,max(max.b.t)), main = "Insertion Size of Populations 1, 4, and 9",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[1], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
points(jitter(rep(1,length(hist.4[,5])), amount=0.2), log10(hist.4[,5]), cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(2,length(hist.9[,5])), amount=0.2), log10(hist.9[,5]), cex = 0.75, col = cols[9], pch = 16)
axis(2, labels = T, tick = T)
dev.off()

hist.0 <- hist.0[complete.cases(hist.0[,4]), ]

##### Pops 2-7, 3-8, 5-10
max.nna <- function(x){
	x <- x[!is.na(x)]
	return(x)
}


max.b.t <- max.nna(c(log10(max(hist.2[,5])), log10(max(hist.7[,5])), log10(max(hist.3[,5])), 
	log10(max(hist.8[,5])), log10(max(hist.5[,5])), log10(max(hist.0[,5]))))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))
boxplot(cbind(log10(hist.2[,5]), log10(hist.7[,5]), log10(hist.3[,5]), log10(hist.8[,5]), 
	log10(hist.5[,5]), log10(hist.0[,5])), col = cols)

boxplot(cbind(hist.2[,5], hist.7[,5], hist.3[,5], hist.8[,5], hist.5[,5], hist.0[,5]), pch = 16, cex = 0.8,
	col = c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10]), border = "grey40")

##### Jitter plots of non-3a27 and 168 parents
pdf("ThreePopPairss_2-7_3-8_5-10_InsertSizeBsub_jitterplots_Aug2020.pdf", height = 8, width = 13)
plot(jitter(rep(0,length(hist.2[,5])),amount=0.2), log10(hist.2[,5]),
     xlim=range(-0.5,4), ylim=range(-3,max(max.b.t)), main = "Insertion Size of Populations 2-7, 3-8, and 5-10",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[2], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
points(jitter(rep(0.5,length(hist.7[,5])), amount=0.2), log10(hist.7[,5]), cex = 0.75, col = cols[7], pch = 16)
points(jitter(rep(1.5,length(hist.3[,5])), amount=0.2), log10(hist.3[,5]), cex = 0.75, col = cols[3], pch = 16)
points(jitter(rep(2,length(hist.8[,5])), amount=0.2), log10(hist.8[,5]), cex = 0.75, col = cols[8], pch = 16)
points(jitter(rep(3,length(hist.5[,5])), amount=0.2), log10(hist.5[,5]), cex = 0.75, col = cols[5], pch = 16)
points(jitter(rep(3.5,length(hist.0[,5])), amount=0.2), log10(hist.0[,5]), cex = 0.75, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()


##### Compare two examples
max.b.t <- max(c(max(hist.2[,5]), max(hist.4[,5]), max(hist.5[,5])))
plot(jitter(rep(0,length(hist.2[,5])),amount=0.2), hist.2[,5],
     xlim=range(-0.5, 3.5), ylim=range(-3, max(max.b.t)),
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[2], pch = 16)
points(jitter(rep(2,length(hist.4[,5])), amount=0.2), hist.4[,5], cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(3,length(hist.5[,5])), amount=0.2), hist.5[,5], cex = 0.75, col = cols[5], pch = 16)


##### Make GGPLOTs of distributions
# install.packages("ggjoy")
library(ggridges)
library(ggplot2)

test.data <- rbind(
			cbind(rep(1, times = length(hist.2[,5])), log10(hist.2[,5])), 
			cbind(rep(2, times = length(hist.7[,5])), log10(hist.7[,5])),
			cbind(rep(4, times = length(hist.3[,5])), log10(hist.3[,5])),
			cbind(rep(5, times = length(hist.8[,5])), log10(hist.8[,5])),
			cbind(rep(7, times = length(hist.5[,5])), log10(hist.5[,5])),
			cbind(rep(8, times = length(hist.0[,5])), log10(hist.0[,5])))
cols4 <- c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10])
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Log_Values")
head(test.data)

test.data <- test.data[test.data[,2] > 0,]

pdf(file = "Fig_4B_Pops_2-7_3-8_5-10_LogBPDistances_Aug_2020.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Log_Values, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "º", point_size = 2, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()


##### Distribution and variance tests

##### Pop2 to pop7
wilcox.test(log10(hist.2[,5]),  log10(hist.7[,5]))
ks.test(log10(hist.2[,5]), log10(hist.7[,5]))
t.test(log10(hist.2[,5]), log10(hist.7[,5]))
var.test(log10(hist.2[,5]), log10(hist.7[,5]), alternative = "greater")


##### Pop3 to pop8
wilcox.test(hist.3[,5], hist.8[,5])
ks.test(hist.3[,5], hist.8[,5])
t.test(hist.8[,5], hist.3[,5])
var.test(log10(hist.3[,5]), log10(hist.8[,5]), alternative = "greater")


##### pop5 t0 pop10
wilcox.test(log10(hist.5[,5]),  log10(hist.0[,5]))
ks.test(log10(hist.5[,5]), log10(hist.0[,5]))
t.test(log10(hist.5[,5]), log10(hist.0[,5]))
var.test(log10(hist.5[,5]), log10(hist.0[,5]), alternative = "greater")

##### pop1 to pop4
wilcox.test(log10(hist.1[,5]),  log10(hist.4[,5]))
ks.test(log10(hist.1[,5]), log10(hist.4[,5]))
t.test(log10(hist.1[,5]), log10(hist.4[,5]))
var.test(log10(hist.1[,5]), log10(hist.4[,5]), alternative = "greater")

##### pop1 to pop9
wilcox.test(log10(hist.1[,5]),  log10(hist.9[,5]))
ks.test(log10(hist.1[,5]), log10(hist.9[,5]))
t.test(log10(hist.1[,5]), log10(hist.9[,5]))
var.test(log10(hist.1[,5]), log10(hist.9[,5]), alternative = "greater")

##### pop4 to pop9
wilcox.test(log10(hist.4[,5]),  log10(hist.9[,5]))
ks.test(log10(hist.4[,5]), log10(hist.9[,5]))
t.test(log10(hist.4[,5]), log10(hist.9[,5]))
var.test(log10(hist.4[,5]), log10(hist.9[,5]), alternative = "greater")

###############################################################################################################
####################### Table out the number of insertions per strain in populations ##########################
###############################################################################################################

max.ins.num <- max(c(rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5])), 
	rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5])),
  rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5])), rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5])),
  rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5])), rowSums(table(test.marks.6.nmatch[,2], test.marks.6.nmatch[,5])),
  rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5])), rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5])),
  rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5])), rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))))

t.1 <- rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5]))
t.2 <- rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5]))
t.3 <- rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5]))
t.4 <- rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5]))
t.5 <- rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5]))
t.6 <- rowSums(table(test.marks.6.nmatch[,2], test.marks.6.nmatch[,5]))
t.7 <- rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5]))
t.8 <- rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5]))
t.9 <- rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5]))
t.0 <- rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))
t.149 <- rowSums(table(pop149[,2], pop149[,5]))
hist(149)
##### Pops 149
pdf("Three_3a27x168_NumbOfRecombJitterplots.pdf", height = 8, width = 8)
plot(jitter(rep(0, times = length(t.149)), amount=0.2), t.149,
     	main = "Number of Cross-Over Events pops1,4 and 9 [replace text]",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = "Orangered3", pch = 16, xlab = "Individuals", 	
     	ylab = "Number of Cross-Over Events")
axis(2, labels = T, tick = T)
text(-0.13,70, labels = "Cross-Over Events n = 769", cex = 1.3)
dev.off()


pdf("10_popsBsub_NumbOfRecombJitterplots.png", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.1)), amount=0.2), t.1,
     xlim=range(-0.5,9.5), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[1], pch = 16, xlab = "10 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
abline(h = 190, lwd = 0.8, col = "grey40", lty = 2)
abline(h = (190/2), lwd = 0.8, col = "grey40", lty = 3)
points(jitter(rep(1,length(t.2)), amount=0.2), t.2, cex = 1.5, col = cols[2], pch = 16)
points(jitter(rep(2,length(t.3)), amount=0.2), t.3, cex = 1.5, col = cols[3], pch = 16)
points(jitter(rep(3,length(t.4)), amount=0.2), t.4, cex = 1.5, col = cols[4], pch = 16)
points(jitter(rep(4,length(t.5)), amount=0.2), t.5, cex = 1.5, col = cols[5], pch = 16)
points(jitter(rep(5,length(t.6)), amount=0.2), t.6, cex = 1.5, col = cols[6], pch = 16)
points(jitter(rep(6,length(t.7)), amount=0.2), t.7, cex = 1.5, col = cols[7], pch = 16)
points(jitter(rep(7,length(t.8)), amount=0.2), t.8, cex = 1.5, col = cols[8], pch = 16)
points(jitter(rep(8,length(t.9)), amount=0.2), t.9, cex = 1.5, col = cols[9], pch = 16)
points(jitter(rep(9,length(t.0)), amount=0.2), t.0, cex = 1.5, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()




max.ins.num <- max(c(rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5])), 
 rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5])),
  rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5]))))

pdf("Fig_2A_ThreePops_168x3a27_Bsub_NumbOfRecombJitterplots_2.pdf", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.1)), amount=0.2), t.1,
     xlim=range(-0.5,2.5), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[1], pch = 16, xlab = "168 x 3A27 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
points(jitter(rep(1,length(t.4)), amount=0.2), t.4, cex = 1.5, col = cols[4], pch = 16)
points(jitter(rep(2,length(t.9)), amount=0.2), t.9, cex = 1.5, col = cols[9], pch = 16)
axis(2, labels = T, tick = T)
dev.off()



##### Make dotplot for 2-7, 3-8,  and 5-10 shuffles
max.ins.num <- max(c(rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5])), 
  rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5])),
  rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5])),
  rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5])),
  rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5])),
  rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))))


pdf("Fig4A_ThreePops_2-7_3-8_5-10_Bsub_NumbOfRecombJitterplots_2_Aug_2020.pdf", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.2)), amount=0.1), t.2,
     xlim=range(-0.5,4), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[2], pch = 16, xlab = "168 x 3A27 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
points(jitter(rep(0.6,length(t.7)), amount=0.1), t.7, cex = 1.5, col = cols[7], pch = 16)
points(jitter(rep(1.6,length(t.3)), amount=0.1), t.3, cex = 1.5, col = cols[3], pch = 16)
points(jitter(rep(2.2,length(t.8)), amount=0.1), t.8, cex = 1.5, col = cols[8], pch = 16)
points(jitter(rep(3.2,length(t.5)), amount=0.1), t.5, cex = 1.5, col = cols[5], pch = 16)
points(jitter(rep(3.8,length(t.0)), amount=0.1), t.0, cex = 1.5, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()




wilcox.test(t.1 ~ t.4)


test.data <- rbind(
			cbind(rep(1, times = length(t.2)), log10(t.2)), 
			cbind(rep(2, times = length(t.7)), log10(t.7)),
			cbind(rep(4, times = length(t.3)), log10(t.3)),
			cbind(rep(5, times = length(t.8)), log10(t.8)),
			cbind(rep(7, times = length(t.5)), log10(t.5)),
			cbind(rep(8, times = length(t.0)), log10(t.0)))
cols4 <- c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10])
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Log_Insertions")
head(test.data)

library(ggplot2)
pdf(file = "Pops_2-7_3-8_5-10_LogBPDistances_Aug_2020.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Log_Insertions, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "º", point_size = 2, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()


test.data <- rbind(
			cbind(rep(1, times = length(t.2)), (t.2)), 
			cbind(rep(2, times = length(t.7)), (t.7)),
			cbind(rep(4, times = length(t.3)), (t.3)),
			cbind(rep(5, times = length(t.8)), (t.8)),
			cbind(rep(7, times = length(t.5)), (t.5)),
			cbind(rep(8, times = length(t.0)), (t.0)))
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Insertions")

pdf(file = "Pops_2-7_3-8_5-10_Insertions_Aug_2020.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Insertions, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "•", point_size = 5, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()



##############################################################################################################
################################################## Bio Circos ################################################
##############################################################################################################

# install.packages('BioCircos')
# if (!require('devtools')){install.packages('devtools')}
# devtools::install_github('lvulliard/BioCircos.R', build_vignettes = TRUE)


library(BioCircos)



##### Set Directory
setwd("/Users/ju0/Desktop/b.subtilis_recomb/2019-09_203_trmd_Samples/Bac_sub_ManuscriptPeds")


########################################################################################
################################## Build and Parse Data sets ###########################
########################################################################################

BioCircos()

chr <- rep("chr1", times = nrow(testout.7))

head(testout.5)
plot(testout.5[,2])
max(testout.5[,4])


myGenome = list("B. subtilis" = 4010000)
hist(testout.3[,2], breaks = max(testout.3[,2]))
myGenome1 <- testout.3[testout.3[,2] == 1,]
# myGenome = list(myGenome1[,3], myGenome1[,5])

BioCircos(genome = myGenome, genomeFillColor = c("pink4", "pink4"), genomeBorderColor = "pink3", 
	genomeTextColor = "pink4", genomeTicksScale = 1e+6, chrPad = 0.002, genomeTicksTextSize = "0.6em", 
	genomeBorderSize = 0.2)

hist(testout.3[,2])
myGenome1 <- testout.3[testout.3[,2] == 6,]
arcs_chromosomes = "B. subtilis"
arcs_begin = myGenome1[,3]
arcs_end = myGenome1[,4]
tracklist = BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, minRadius = 1, 
	maxRadius = 1.14, colors = c("darkorange3", "darkorange3"))

BioCircos(tracklist, genomeFillColor = c("grey20", "grey20"), chrPad = 0.0, 
		displayGenomeBorder = F, genomeTicksDisplay = F, genomeLabelTextSize = 0, 
		genome = myGenome, opacities = rep(1, times = 1000))


#################################################################################
########################## Loop to build many tracks ############################
#################################################################################

##### modify the names in this section to get each population's circos plot
# To get a PDF run plot and opens in browswer then print, but save as pdf

# Set pop color
pop.numb <- 7

# Set pop filtered markers
pop <- test.marks.7.nmatch

# Set pop unfiltered markers
pop <- testout.7

map.data <- map.data.7[,4]
posVert = p.var.9
dpth <- s.depth.9
m.depth <- depth.9


##### Population seven individuals to remove
testout.7 <- testout.7[testout.7[,2] != 18, ]
testout.7 <- testout.7[testout.7[,2] != 14, ]
testout.7 <- testout.7[testout.7[,2] != 10, ]
testout.7 <- testout.7[testout.7[,2] != 9, ]
testout.7 <- testout.7[testout.7[,2] != 7, ]


########## create circos plot
max.i <- pop[!duplicated(pop[,2]), 2]
hist(pop[,2])
j <- 1
##### Start BioCircos
pop.col <- "dodgerblue1"
opvar <- 1/(max(as.numeric(pop[,2])))
for(l in 1:length(max.i)){
	i <- max.i[l]
	myGenome1 <- pop[pop[,2] == i,]
	arcs_chromosomes = "B. subtilis"
	arcs_begin = myGenome1[,3]
	arcs_end = myGenome1[,4]
	minr <- 1.00
	maxr <- 1.14
	if(i == min(pop[,2])){
		tracklist = BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, opacities = rep(0.1, times = 100), 
			color = c(pop.col, pop.col))
	}
	else{
		tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, opacities = rep(opvar, times = 100), 
			color = c(pop.col, pop.col))
	}
}

j <- 1.02
for(l in 1:length(max.i)){
	i <- max.i[l]
	myGenome1 <- pop[pop[,2] == i,]
	arcs_chromosomes = "B. subtilis"
	arcs_begin = myGenome1[,3]
	arcs_end = myGenome1[,4]
	maxr <- j - 0.05
	minr <- maxr - 0.04
	j <- j - 0.05
	if(i == 2000){
		tracklist = BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, color = c(pop.col, pop.col))
	}
	else{
		tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, color = c(pop.col, pop.col))
	}
}


BioCircos(tracklist, genomeFillColor = c("grey30", "grey30"), chrPad = 0.0, 
		displayGenomeBorder = T, genomeTicksDisplay = F, genomeLabelTextSize = 0, 
		genome = myGenome, opacities = rep(1, times = 1000))



#########################################################################
######################## Create Extra Data Tracks #######################
#########################################################################

##### Read in files for methylation and gene postion
genepositions <- read.table("3a27_genome_anno_reduced.txt", header = T)
methyl_sites <- read.table("3a27_gaygnnnnnnctt_frwrd.txt", header = F)

##### Set Genome Size
myGenome = list("B. subtilis" = 4010000)

##### Set sample tracks to population
map.data <- map.data.1[,4]
m.depth.p <- m.depth[,3]
dpth.p <- dpth[,3]
posVert = p.var.1


##### Get a randomly sampled and ordered set of data points
# Reduction power
pwr <- 10
samp.track <- cbind(posVert, map.data[1:length(posVert)], dpth.p[1:length(posVert)], m.depth.p[1:length(posVert)])
samp.track <- samp.track[sample(1:nrow(samp.track), size = round((length(posVert)/pwr)), replace = F), ]
samp.track <- samp.track[order(samp.track[,2]), ]
map.data.r <- samp.track[,2]
posVert.r <- samp.track[,1]
m.depth.p <- samp.track[,4]
dpth.p <- samp.track[,3]

##### Space Between Markers of Parent A and B
minr = 1.14
maxr = 1.22
tracklist = tracklist + BioCircosLineTrack("B.subtilis", arcs_chromosomes, position = map.data.r, values = posVert.r, 
		minRadius = minr, maxRadius = maxr, color = "#a1887a", width = 0.8)
tracklist = tracklist + BioCircosBackgroundTrack("B.subtilis", arcs_chromosomes, fillColors = '#f0f0f2', borderSize = 0, 
		minRadius = minr, maxRadius = maxr)

##### Site Read Density
minr = 1.22
maxr = 1.30
tracklist = tracklist + BioCircosLineTrack("B.subtilis", arcs_chromosomes, position = map.data.r, values = m.depth.p, 
		minRadius = minr, maxRadius = maxr, color = "#498df2", width = 0.8)
tracklist = tracklist + BioCircosBackgroundTrack("B.subtilis", arcs_chromosomes, fillColors = '#fafaff', borderSize = 0, 
		minRadius = minr, maxRadius = maxr)

##### Mean Read Depth
minr = 1.30
maxr = 1.33
tracklist = tracklist + BioCircosLineTrack("B.subtilis", arcs_chromosomes, position = map.data.r, values = dpth.p, 
		minRadius = minr, maxRadius = maxr, color = "#a0de6a", width = 0.8)
tracklist = tracklist + BioCircosBackgroundTrack("B.subtilis", arcs_chromosomes, fillColors = '#f0f0f2', borderSize = 0, 
		minRadius = minr, maxRadius = maxr)

##### Positive Strand Genes
genepositions.p <- genepositions[genepositions[,7] == "+", ]
arcs_begin <- genepositions.p[,4]
arcs_end   <- genepositions.p[,4] + genepositions.p[,6]
minr = 1.33
maxr = 1.36
tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
	minRadius = minr, maxRadius = maxr, color = "#aaad39")

##### Positive Strand Genes
minr = 1.36
maxr = 1.39
genepositions.n <- genepositions[genepositions[,7] == "-", ]
arcs_begin <- genepositions.n[,4] - genepositions.n[,6]
arcs_end   <- genepositions.n[,4]
tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
	minRadius = minr, maxRadius = maxr, color = "#e3b520")

##### Positive Strand Genes
minr = 1.39
maxr = 1.42
methyl_sites <- methyl_sites[methyl_sites[,6] == "+", ]
arcs_begin <- methyl_sites[,2]
arcs_end   <- methyl_sites[,3] + 200
tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
	minRadius = minr, maxRadius = maxr, color = "#050e3b")

BioCircos(tracklist, displayGenomeBorder = T, LINEMouseOutDisplay = F, genomeLabelTextSize = 0, genome = myGenome)


##############################################################################################################
################################# Sample Tracks with Density on Genome Track #################################
##############################################################################################################

pop.numb <- 9
pop <- test.marks.9.nmatch
map.data <- map.data.1[,4]
posVert = p.var.9
dpth <- s.depth.9
m.depth <- depth.9

##### Start BioCircos
pop.col <- cols[pop.numb]
opvar <- 1/(max(as.numeric(pop[,2])))
for(i in (min(pop[,2])):max(pop[,2])){
	myGenome1 <- pop[pop[,2] == i,]
	arcs_chromosomes = "B. subtilis"
	arcs_begin = myGenome1[,3]
	arcs_end = myGenome1[,4]
	minr <- 1.00
	maxr <- 1.14
	if(i == min(pop[,2])){
		tracklist = BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, opacities = rep(0.1, times = 100), 
			color = c(pop.col, pop.col))
	}
	else{
		tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, opacities = rep(opvar, times = 100), 
			color = c(pop.col, pop.col))
	}
}

j <- 1.02
for(i in (min(pop[,2])):max(pop[,2])){
	myGenome1 <- pop[pop[,2] == i,]
	arcs_chromosomes = "B. subtilis"
	arcs_begin = myGenome1[,3]
	arcs_end = myGenome1[,4]
	maxr <- j - 0.05
	minr <- maxr - 0.04
	j <- j - 0.05
	if(i == 2000){
		tracklist = BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, color = c(pop.col, pop.col))
	}
	else{
		tracklist = tracklist + BioCircosArcTrack("B. subtilis", arcs_chromosomes, arcs_begin, arcs_end, 
			minRadius = minr, maxRadius = maxr, color = c(pop.col, pop.col))
	}
}


BioCircos(tracklist, genomeFillColor = c("grey30", "grey30"), chrPad = 0.0, 
		displayGenomeBorder = T, genomeTicksDisplay = F, genomeLabelTextSize = 0, 
		genome = myGenome, opacities = rep(1, times = 1000))




###########################################################################################
################################ Genome Track with Line Chart #############################
###########################################################################################


head(p.2[1:5, 1:5])


chrVert =  rep(c(1, 3, 5), c(20,10,5))
posVert = c(249250621*log(c(20:1, 10:1, 5:1), base = 20))


plot(sin(posVert[1:80000]))


redux <- function(x){
	x <- x[c(T, F, F, F, F, F, F)]
	return(x)
}

p.3.r <- redux(p.3)
dim(p.3.r)
p.3.r[1:19, 1:30]
chrVert = "B. subtilis"
vals = (colSums(p.3.r)/(max(p.3.r)))
range(vals)
pos <- as.numeric(colnames(p.3.r))
range()
pos[1:200]

tracks = BioCircosLineTrack('LineTrack1', as.character(chrVert), vals, values = cos(pos))

tracks = tracks + BioCircosLineTrack('LineTrack2', as.character(chrVert), 0.95*vals, 
  values = sin(posVert), color = "#40D4B9")

tracks = tracks + BioCircosBackgroundTrack('Bg', fillColors = '#FFEEBB', borderSize = 0)

BioCircos(tracks, genome = myGenome, chrPad = 0.00, genomeFillColor = c("pink4", "pink4"))





###############################################################################################################
############################################ Permutation testing ##############################################
###############################################################################################################
par(mfrow = c(1,1))
setwd("~/Desktop/b.subtilis_recomb/2019-09_203_trmd_Samples/Bac_sub_ManuscriptPeds/")

##### Read in files for methylation and gene postion
# Gene postions
genepositions <- read.table("3a27_genome_anno_reduced.txt", header = T)
gene.sizes <- genepositions[,6]
genepositions <- c(genepositions[,4], genepositions[,5])

# Methylation sites
motif.f <- read.delim("3a27_gaygnnnnnnctt_frwrd.txt", header = F)
motif.n <- read.delim("bsub3a27_aagnnnnnncrtc.txt", header = F)
colnames(motif.f) <- c("chr", "start", "end", "dot", "1k", "strand")
colnames(motif.n) <- c("chr", "start", "end", "dot", "1k", "strand")
motif.f <- round((motif.f[,2] + motif.f[,3])/2)
motif.n <- round((motif.n[,2] + motif.n[,3])/2)
length(motif.f)
length(motif.n)
motif <- c(motif.f, motif.n)
motif <- motif[complete.cases(motif)]
motif <- motif[order(motif)]

length(motif)
# Load population recombination positions

pop <- 1
pop.c <- 1
pop.rcmb.sts <- test.marks.1.nmatch
pop.rcmb.frt <- pop.rcmb.sts[,3]
pop.rcmb.end <- pop.rcmb.sts[,4]
pop.rcmb.rnd <- round((pop.rcmb.frt + pop.rcmb.end)/2)

pop.rcmb.frt <- pop.rcmb.frt[!duplicated(pop.rcmb.frt)]
pop.rcmb.end <- pop.rcmb.end[!duplicated(pop.rcmb.end)]
pop.rcmb.rnd <- pop.rcmb.rnd[!duplicated(pop.rcmb.rnd)]

population <- pop.rcmb.frt
dim(as.matrix(population))
population <- population[complete.cases(population)]
##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?


##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 46

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Start Loop for permutation test, Get list of Random distances to population insertion sites
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	random.set <- sample(r.num, size = sample.size, replace = F)
	random.set <- random.set[order(random.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:sample.size){
		rand.data.point <- which(abs(random.set-population[j])==min(abs(random.set-population[j])))
		closest.random.dist <- ((random.set[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.random.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.random.dist)
		}
	}
	rand.mean <- mean(clst.rand.dist.p)
	rand.stdv <- sd(clst.rand.dist.p)
	if(i == 1){
		rand.mean.p <- rand.mean
		rand.stdv.p <- rand.stdv
	}
	else{
		rand.mean.p <- c(rand.mean.p, rand.mean)
		rand.stdv.p <- c(rand.stdv.p, rand.stdv)		
	}
}
mean(rand.mean.p)
mean(rand.stdv.p)
range(rand.mean.p)
hist(rand.mean.p, breaks = 100, col = cols[pop], border = "grey40")

##### Methylation variance test
##### Start Loop for permutation test
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	motif.set <- sample(motif, size = sample.size, replace = F)
	motif.set <- motif.set[order(motif.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:sample.size){
		motif.data.point <- which(abs(motif.set-population[j])==min(abs(motif.set-population[j])))
		closest.motif.dist <- ((motif.set[motif.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.motif.dist.p <- closest.motif.dist
		}
		else{
			clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
		}
	}
	motif.mean <- mean(clst.motif.dist.p)
	motif.stdv <- sd(clst.motif.dist.p)
	if(i == 1){
		motif.mean.p <- motif.mean
		motif.stdv.p <- motif.stdv
	}
	else{
		motif.mean.p <- c(motif.mean.p, motif.mean)
		motif.stdv.p <- c(motif.stdv.p, motif.stdv)		
	}
}
mean(motif.mean.p)
mean(motif.stdv.p)
range(motif.mean.p)
hist(motif.mean.p, breaks = 100, col = cols[pop], border = "grey40")

plot(density(motif.mean.p), type = "l")
lines(density(rand.mean.p), type = "l")


mean(rand.mean.p)
mean(rand.stdv.p)
range(rand.mean.p)

for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}


##### Motif Data
mean(clst.motif.dist.p)
sd(clst.motif.dist.p)
range(clst.motif.dist.p)

##### Random Data
mean(motif.mean.p)
mean(motif.stdv.p)
range(motif.mean.p)

##### Create Stats for pop
popstats <- rbind(c(mean(motif.mean.p), mean(motif.stdv.p), range(motif.mean.p)),
	mean(clst.motif.dist.p), sd(clst.motif.dist.p), range(clst.motif.dist.p))
write.table(popstats, file = paste(pop,"pop_stats.txt", sep = "_"))

##### Random Data
plot(rand.mean.p, cex = 0.01)
lines(rand.mean.p, col = "aquamarine3", lwd = 0.7)
lines(rand.stdv.p, col = "pink3", lwd = 0.7)
mean(rand.mean.p)
mean(rand.stdv.p)

##### histograms of random and methylated distances to recombination sites
pdf(file = paste(pop, "PopOverRand_hist.pdf",sep = "_"), width = 10, height = 8)
hist(rand.mean.p, breaks = 100, col = "azure3")
rug(rand.mean.p, col = "grey40")
dev.off()

pdf(file = paste(pop, "PopOverRand_hist.pdf",sep = "_"), width = 10, height = 8)
hist(clst.motif.dist.p, breaks = 70, col = "goldenrod3")
rug(clst.motif.dist.p, col = "grey40")
dev.off()



par(mfrow = c(1,1))
plot(density(motif.mean.p), type = "l", col = cols[pop])
lines(density(rand.mean.p), type = "l", col = "grey40")




wilcox.test(motif.mean.p,  rand.mean.p)
ks.test(motif.mean.p, rand.mean.p)
t.test(motif.mean.p, rand.mean.p)

var.test(motif.mean.p, rand.mean.p, alternative = "greater")


################################################################################################
#################### Check for population recombination vs. gene position ######################
################################################################################################

pop <- 1
pop.rcmb.sts <- test.marks.1.nmatch
pop.rcmb.frt <- pop.rcmb.sts[,3]
pop.rcmb.end <- pop.rcmb.sts[,4]
pop.rcmb.rnd <- round((pop.rcmb.frt + pop.rcmb.end)/2)
pop.rcmb.frt <- pop.rcmb.frt[!duplicated(pop.rcmb.frt)]
pop.rcmb.end <- pop.rcmb.end[!duplicated(pop.rcmb.end)]
pop.rcmb.rnd <- pop.rcmb.rnd[!duplicated(pop.rcmb.rnd)]
population <- pop.rcmb.frt

##### Number of tests
data.times <- 1000

##### Set genome_size to genome size
genome_size <- (max(population))

##### Set an initial random seed for further randomize future seeds
s.seed <- 3

##### Get number of tests per population
data.times <- 1000

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(0:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of tests
ntst <- 1000

##### Test if gene position is correlated with recombination location
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	genepositions.set <- sample(genepositions, size = sample.size, replace = T)
	genepositions.set <- genepositions.set[order(genepositions.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:sample.size){
		genepositions.data.point <- which(abs(genepositions.set-population[j])==min(abs(genepositions.set-population[j])))
		closest.genepositions.dist <- ((genepositions.set[genepositions.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.genepositions.dist.p <- closest.genepositions.dist
		}
		else{
			clst.genepositions.dist.p <- c(clst.genepositions.dist.p, closest.genepositions.dist)
		}
	}
	genepositions.mean <- mean(clst.genepositions.dist.p)
	genepositions.stdv <- sd(clst.genepositions.dist.p)
	if(i == 1){
		genepositions.mean.p <- genepositions.mean
		genepositions.stdv.p <- genepositions.stdv
	}
	else{
		genepositions.mean.p <- c(genepositions.mean.p, genepositions.mean)
		genepositions.stdv.p <- c(genepositions.stdv.p, genepositions.stdv)		
	}
}

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Test if gene position is correlated with recombination location
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	r.num.set <- sample(r.num, size = sample.size, replace = T)
	r.num.set <- r.num.set[order(r.num.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:sample.size){
		r.num.set.data.point <- which(abs(r.num.set-population[j])==min(abs(r.num.set-population[j])))
		closest.r.num.dist <- ((r.num.set[r.num.set.data.point] - population[j])^2)^0.5
		if(j == 1){
			closest.r.num.dist.p <- closest.r.num.dist
		}
		else{
			closest.r.num.dist.p <- c(closest.r.num.dist.p, closest.r.num.dist)
		}
	}
	closest.r.num.dist.mean <- mean(closest.r.num.dist.p)
	closest.r.num.dist.stdv <- sd(closest.r.num.dist.p)
	if(i == 1){
		closest.r.num.dist.mean.p <- closest.r.num.dist.mean
		closest.r.num.dist.stdv.p <- closest.r.num.dist.stdv
	}
	else{
		closest.r.num.dist.mean.p <- c(closest.r.num.dist.mean.p, closest.r.num.dist.mean)
		closest.r.num.dist.stdv.p <- c(closest.r.num.dist.stdv.p, closest.r.num.dist.stdv)		
	}
}


##### Create Histogram of distance to gene start/stop
pdf(file = paste(pop, "dist_2_geneSrtStp.pdf"), width = 10, height = 8)
hist(genepositions.mean.p, breaks = 100, col = cols[pop])
rug(genepositions.mean.p, col = "grey40")
dev.off()

hist(genepositions.stdv.p, breaks = 100, col = "pink3")
rug(genepositions.stdv.p, col = "grey40")

##### Plot Random and Gene Position distributions
pdf(file = paste("Pop", pop, "GenePosition_Vs_Insertion_July2020.pdf",sep = "_"), width = 10, height = 8)
plot(density(genepositions.mean.p), col = cols[pop], lwd = 1.5, main = paste("Pop_", pop, "_GenePos_Vs_Insertion", sep = ""),
	ylim = c(0,max(c(density(genepositions.mean.p)$y, density(closest.r.num.dist.mean.p)$y))))
lines(density(closest.r.num.dist.mean.p), col = "azure3", lwd = 1.5)
dev.off()



################################################################################################
######################## Check for population recombination vs. random #########################
################################################################################################

pop <- 1
pop.rcmb.sts <- test.marks.1.nmatch
pop.rcmb.frt <- pop.rcmb.sts[,3]
pop.rcmb.end <- pop.rcmb.sts[,4]
pop.rcmb.rnd <- round((pop.rcmb.frt + pop.rcmb.end)/2)
pop.rcmb.frt <- pop.rcmb.frt[!duplicated(pop.rcmb.frt)]
pop.rcmb.end <- pop.rcmb.end[!duplicated(pop.rcmb.end)]
pop.rcmb.rnd <- pop.rcmb.rnd[!duplicated(pop.rcmb.rnd)]
population <- pop.rcmb.frt

##### Number of tests
data.times <- 1000

##### Set genome_size to genome size
genome_size <- 4010000

##### Set an initial random seed for further randomize future seeds
s.seed <- 3

##### Get number of tests per population
data.times <- 1000

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(0:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get random numbers to length of genome
r.num <- c(1:genome_size)

##### Get randome numbers to length of genome
r.num.t <- c(1:genome_size)
set.seed(0570)
r.num.t <- sample(r.num, size = sample.size*4, replace = F)
r.num.t <- r.num.t[order(r.num.t)]
##### Number of tests
ntst <- 1000
i <- 1
j <- 1
##### Test if gene position is correlated with recombination location
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	r.num.set <- sample(r.num, size = sample.size, replace = T)
	r.num.set <- r.num.set[order(r.num.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:length(r.num.set)){
		r.num.data.point <- which(abs(r.num.set-r.num.t[j])==min(abs(r.num.set-r.num.t[j])))
		r.num.data.point.dist <- ((r.num.set[r.num.data.point] - r.num.t[j])^2)^0.5
		if(j == 1){
			r.num.data.point.dist.p <- r.num.data.point.dist
		}
		else{
			r.num.data.point.dist.p <- c(r.num.data.point.dist.p, r.num.data.point.dist)
		}
	}
	r.num.data.point.dist.p.mean <- mean(r.num.data.point.dist.p)
	r.num.data.point.dist.p.stdv <- sd(r.num.data.point.dist.p)
	if(i == 1){
		r.num.data.point.dist.p.mean.p <- r.num.data.point.dist.p.mean
		r.num.data.point.dist.p.stdv.p <- r.num.data.point.dist.p.stdv
	}
	else{
		r.num.data.point.dist.p.mean.p <- c(r.num.data.point.dist.p.mean.p, r.num.data.point.dist.p.mean)
		r.num.data.point.dist.p.stdv.p <- c(r.num.data.point.dist.p.stdv.p, r.num.data.point.dist.p.stdv)		
	}
}
hist(r.num.data.point.dist.p.mean.p, breaks = 100, col = cols[pop], border = "grey40")
dev.off()



##### Get random numbers to length of genome
r.num <- c(1:genome_size)
population <- population[order(population)]
##### Test if gene position is correlated with recombination location
for(i in 1:as.numeric(ntst)){
	set.seed(random.pop.nums[i])
	r.num.set <- sample(r.num, size = sample.size, replace = T)
	r.num.set <- r.num.set[order(r.num.set)]
	# which(abs([vector] - [singl value])==min(abs([vector] - [single value])))
	for(j in 1:sample.size){
		r.num.set.data.point <- which(abs(r.num.set-population[j])==min(abs(r.num.set-population[j])))
		closest.r.num.dist <- ((r.num.set[r.num.set.data.point] - population[j])^2)^0.5
		if(j == 1){
			closest.r.num.dist.p <- closest.r.num.dist
		}
		else{
			closest.r.num.dist.p <- c(closest.r.num.dist.p, closest.r.num.dist)
		}
	}
	closest.r.num.dist.mean <- mean(closest.r.num.dist.p)
	closest.r.num.dist.stdv <- sd(closest.r.num.dist.p)
	if(i == 1){
		closest.r.num.dist.mean.p <- closest.r.num.dist.mean
		closest.r.num.dist.stdv.p <- closest.r.num.dist.stdv
	}
	else{
		closest.r.num.dist.mean.p <- c(closest.r.num.dist.mean.p, closest.r.num.dist.mean)
		closest.r.num.dist.stdv.p <- c(closest.r.num.dist.stdv.p, closest.r.num.dist.stdv)		
	}
}

pdf(file = paste("Pop",pop,"_random_Vs_population_july2020.pdf"), width = 10, height = 8)
plot(density(closest.r.num.dist.mean.p), col = cols[pop], type = "l")
lines(density(r.num.data.point.dist.p.mean.p))
dev.off()

##### Create Histogram of distance to gene start/stop
pdf(file = paste(pop, "dist_2_geneSrtStp.pdf"), width = 10, height = 8)
hist(genepositions.mean.p, breaks = 100, col = cols[pop])
rug(genepositions.mean.p, col = "grey40")
dev.off()

hist(genepositions.stdv.p, breaks = 100, col = "pink3")
rug(genepositions.stdv.p, col = "grey40")

##### Plot Random and Gene Position distributions
pdf(file = paste("Pop", pop, "GenePosition_Vs_Insertion_july_2020.pdf",sep = "_"), width = 10, height = 8)
plot(density(genepositions.mean.p), col = cols[pop], lwd = 1.5, main = paste("Pop_", pop, "_GenePos_Vs_Insertion", sep = ""),
	ylim = c(0,max(c(density(genepositions.mean.p)$y, density(closest.r.num.dist.mean.p)$y))))
lines(density(closest.r.num.dist.mean.p), col = "azure3", lwd = 1.5)
dev.off()







##########################################################################################
################################## Create GC Content Data ################################
##########################################################################################

##### Read in GC content file
gc.pos <- read.table("bsub_3a27_GC_positions.txt", header = F)

##### Create blank gc genome
gcnm <- rep(0, times = max(map.data.1[,4]))

##### Impute GC marks into blank genome
for(i in 1:nrow(gc.pos)){
	pos <- gc.pos[i,]
	gcnm[pos] <- 1
}
gcnm <- gcnm[complete.cases(gcnm)]
plot(gcnm[1:1000])
pop <- 1

i <- 1
j <- 12
for(j in 2:12){
	scl <- 2^j
	for(i in 1:length(population)){
		p <- population[i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}
hist(j.m.p[,4], col = "grey60", breaks = 30)
pdf(file = paste("Pop_", pop, "_NumberOfInsertions_ToGC.pdf"), width = 11, height = 6)
plot(j.m.p[,4], type = "l", lwd = 1.4, col = cols[pop], ylab = "Percent GC Content", 
	ylim = c(0.40,0.47), xlab = "Exponential Increase in Bases Considered 2^n")
##### Set first seed
set.seed(8375)
ntst <- 1000
for(i in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	for(j in 2:12){
		scl <- round(2^j)
		for(k in 1:length(rand.smp)){
			p.r <- rand.smp[k]
			p.pre.r <- p.r - scl
			p.pst.r <- p.r + scl
			gc.scn.r <- mean(gcnm[p.pre.r:p.pst.r])
			if(k == 1){
				gc.scn.p.r <- gc.scn.r
			}
			else{
				gc.scn.p.r <- c(gc.scn.p.r, gc.scn.r)
			}
		}
		j.m <- gc.scn.p.r
		if(j == 2){
			j.m.p.r <- j.m
		}
		else{
			j.m.p.r <- c(j.m.p.r, j.m)
		}
	}
	if(i == 1){
		j.m.p.r.p <- j.m.p.r
		lines((j.m.p.r.p[i]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
	else{
		j.m.p.r.p <- rbind(j.m.p.r.p, j.m.p.r)
		lines((j.m.p.r.p[i,]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
}
lines(j.m.p, col = "darkorange3", lwd = 4)
dev.off()

pdf(file = paste("Pop_", pop, "_Percent_Of_Window_GC_Content.pdf", sep = ""))
hist(j.m.p.r.p[,4], breaks = 75, col = cols[pop], border = "grey30", xlab = "Distribution of GC Content",
	ylab = "Density", main = "Percent GC Content at Insertion Site to Random", freq = F)
abline(v = j.m.p[,4], lty = 2)
dev.off()

i <- 1
p.t.p <- rep(0, size = ncol(j.m.p.r.p))
for(i in 1:ncol(j.m.p.r.p)){
	p.test <- c(j.m.p.r.p[,i], j.m.p[i])
	p.t.p[i] <- p.adjust(pnorm(abs(scale(p.test)))[length(p.test)], method = "fdr")
}

pdf(file = paste("Pop_", pop, "_boxplots_Of_Window_GC_Content.pdf", sep = ""), height = 3, width = 8)
boxplot(j.m.p.r.p, 	ylab = "% GC Content", col = cols[pop], border = "grey40",
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:12; 4:4096bp)")
dev.off()


pdf(file = paste("Pop_", pop, "_p-values_Of_Window_GC_Content.pdf", sep = ""), height = 3, width = 8)
plot(1-p.t.p, type = "l", col = cols[pop], lwd = 1.9, ylab = "p-values", 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()




###################################
###### Permutation test of outward 
###################################

pdf(file = paste("Pop_", pop, "_NumberOfInsertions_ToGC_July_2020.pdf"), width = 11, height = 6)
plot(j.m.p[,4], type = "l", lwd = 1.4, col = cols[pop], ylab = "Percent GC Content", 
	ylim = c(0.40,0.47), xlab = "Exponential Increase in Bases Considered 2^n")
##### Set first seed
set.seed(8375)
ntst <- 1000
for(i in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	for(j in 2:12){
		scl <- round(2^j)
		for(k in 1:length(rand.smp)){
			p.r <- rand.smp[k]
			p.pre.r <- p.r - scl
			p.pst.r <- p.r + scl
			gc.scn.r <- mean(gcnm[p.pre.r:p.pst.r])
			if(k == 1){
				gc.scn.p.r <- gc.scn.r
			}
			else{
				gc.scn.p.r <- c(gc.scn.p.r, gc.scn.r)
			}
		}
		j.m <- mean(gc.scn.p.r)
		if(j == 2){
			j.m.p.r <- j.m
		}
		else{
			j.m.p.r <- c(j.m.p.r, j.m)
		}
	}
	if(i == 1){
		j.m.p.r.p <- j.m.p.r
		lines((j.m.p.r.p[i]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
	else{
		j.m.p.r.p <- rbind(j.m.p.r.p, j.m.p.r)
		lines((j.m.p.r.p[i,]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
}
lines(j.m.p, col = "darkorange3", lwd = 4)
dev.off()




##############################
###### 
##############################





set.seed(722)
gene.lengths <- sample(genepositions, size = length(genepositions), replace = F)
plot(as.numeric(gene.lengths), cex = 0.08)
dev.off()

gene.length.wv <- wavCWT(genepositions, wavelet = "gaussian2")
pop <- 1

# note **** get cols2 color palette from below **** note #

pdf(file = "3A27_GeneSize_acrossGnome.pdf", width = 8, height = 3)
plot(gene.length.wv, col = cols2, phase = T, xlab = "Genes Along Genome", main = "Wavelet Analysis of Gene Size along the 3A27 Genome")
dev.off()


##### Set first seed
set.seed(0736)
i <- 1
j <- 1
k <- 1
ntst <- 1000
for(i in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	gene.lengths <- sample(as.numeric(gene.lengths), size = length(rand.smp), replace = F)
	rand.smp <- gene.lengths + rand.smp
	for(j in 2:12){
		scl <- round(2^j)
		for(k in 1:length(rand.smp)){
			p.r <- rand.smp[k]
			p.pre.r <- p.r - scl
			p.pst.r <- p.r + scl
			gc.scn.r <- mean(gcnm[p.pre.r:p.pst.r])
			if(k == 1){
				gc.scn.p.r <- gc.scn.r
			}
			else{
				gc.scn.p.r <- c(gc.scn.p.r, gc.scn.r)
			}
		}
		j.m <- mean(gc.scn.p.r)
		if(j == 2){
			j.m.p.r <- j.m
		}
		else{
			j.m.p.r <- c(j.m.p.r, j.m)
		}
	}
	if(i == 1){
		j.m.p.r.p <- j.m.p.r
		lines((j.m.p.r.p[i]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
	else{
		j.m.p.r.p <- rbind(j.m.p.r.p, j.m.p.r)
		lines((j.m.p.r.p[i,]), col = alpha(rgb(0,0,0),0.05), lwd = 1.1)
	}
}
pdf(file = paste("Pop_", pop, "_Percent_Of_Window_GC_Content_withGeneSizes_window=2048bases.pdf", sep = ""))
hist(j.m.p.r.p[,11], breaks = 75, col = cols[pop], border = "grey30", xlab = "Distribution of GC Content",
	ylab = "Density", main = "Percent GC Content at Insertion Site to Random", freq = F)
abline(v = j.m.p[11], lty = 2)
dev.off()

i <- 1
p.t.p <- rep(0, size = ncol(j.m.p.r.p))
for(i in 1:ncol(j.m.p.r.p)){
	p.test <- c(j.m.p.r.p[,i], j.m.p[i])
	pnorm(abs(scale(p.test)))
	p.t.p[i] <- pnorm(abs(scale(p.test)))[length(p.test)]
}

pdf(file = paste("Pop_", pop, "_boxplots_Of_Window_GC_Content_withGeneSizes.pdf", sep = ""), height = 3, width = 8)
boxplot(j.m.p.r.p, 	ylab = "% GC Content", col = cols[pop], 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()


pdf(file = paste("Pop_", pop, "_p-values_Of_Window_GC_Content_withGeneSizes.pdf", sep = ""), height = 3, width = 8)
plot(1-p.t.p, type = "l", col = cols[pop], lwd = 1.9, ylab = "p-values", 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()








########## Check for number of SNPs to insertion positions ##########
gnm <- rep(0, times = max(map.data.1[,4]))
length(gnm)
nrow(map.data.1)
i <- 1
for(i in 1:nrow(map.data.1)){
	pos <- map.data.1[i,4]
	gnm[pos] <- 1
}

gnm <- gnm[complete.cases(gnm)]
plot(gnm[1:10000])


for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(population)){
		p <- population[i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- mean(gc.scn.p)
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- c(j.m.p, j.m)
	}
}

pdf(file = paste("Pop_", pop, "_SNP-Density_test.pdf", sep = ""), height = 8, width = 11)
plot(j.m.p, type = "l", lwd = 1.4, col = cols[pop], ylab = "Percent Of Window Can Be Polymorphic Variants", 
	ylim = c(0.00,0.15), xlab = "Exponential Increase in Bases of Window 2^i")
dev.off()


##### Set first seed
set.seed(8375)
ntst <- 1000
for(i in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	for(j in 2:12){
		scl <- round(2^j)
		for(k in 1:length(rand.smp)){
			p.r <- rand.smp[k]
			p.pre.r <- p.r - scl
			p.pst.r <- p.r + scl
			gc.scn.r <- mean(gnm[p.pre.r:p.pst.r])
			if(k == 1){
				gc.scn.p.r <- gc.scn.r
			}
			else{
				gc.scn.p.r <- c(gc.scn.p.r, gc.scn.r)
			}
		}
		j.m <- mean(gc.scn.p.r)
		if(j == 2){
			j.m.p.r <- j.m
		}
		else{
			j.m.p.r <- c(j.m.p.r, j.m)
		}
	}
	if(i == 1){
		j.m.p.r.p <- j.m.p.r
		lines((j.m.p.r.p[i]), col = alpha(rgb(0,0,0),0.03), lwd = 1.1)
	}
	else{
		j.m.p.r.p <- rbind(j.m.p.r.p, j.m.p.r)
		lines((j.m.p.r.p[i,]), col = alpha(rgb(0,0,0),0.03), lwd = 1.1)
	}
}

dev.off()

dim(j.m.p.r.p)

n <- 4
pdf(file = paste("Pop_", pop, "_Percent_Of_Window_Variant_Content.pdf", sep = ""))
hist(j.m.p.r.p[,n], breaks = 75, col = "honeydew2", border = "grey30", xlab = "Distribution of GC Content",
	ylab = "Density", main = "Percent Variants Content at Insertion Site to Random", freq = F)
lines(density(j.m.p.r.p[,n]), col = cols[pop])
abline(v = j.m.p[3], lty = 2)
dev.off()

i <- 1
p.t.p <- rep(0, size = ncol(j.m.p.r.p){)
for(i in 1:ncol(j.m.p.r.p)){
	p.test <- c(j.m.p.r.p[,i], j.m.p[i])
	pnorm(abs(scale(p.test)))
	p.t.p[i] <- pnorm(abs(scale(p.test)))[length(p.test)]
}

pdf(file = paste("Pop_", pop, "_boxplots_Of_Window_Polymorphic.pdf", sep = ""), height = 3, width = 8)
boxplot(j.m.p.r.p, 	ylab = "% GC Content", col = cols[pop], 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()


pdf(file = paste("Pop_", pop, "_p-values_Of_Window_Polymorphic.pdf", sep = ""), height = 3, width = 8)
plot(1-p.t.p, type = "l", col = cols[pop], lwd = 1.9, ylab = "p-values", 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()



########## Check for number of SNPs to insertion positions ##########
gnm <- rep(0, times = max(map.data.1[,4]))
length(gnm)
nrow(map.data.1)
i <- 1
for(i in 1:nrow(map.data.1)){
	pos <- map.data.1[i,4]
	gnm[pos] <- 1
}

gnm <- gnm[complete.cases(gnm)]
plot(gnm[1:10000])
length(gnm)

for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(population)){
		p <- population[i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}

pdf(file = paste("Pop_", pop, "_SNP-Density_test.pdf", sep = ""), height = 8, width = 11)
plot(j.m.p[,4], type = "l", lwd = 1.4, col = cols[pop], ylab = "Percent Of Window Can Be Polymorphic Variants", 
	ylim = c(0.00,0.15), xlab = "Exponential Increase in Bases of Window 2^i")
dev.off()


##### Set first seed
set.seed(8375)
ntst <- 1000
for(i in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	for(j in 2:12){
		scl <- round(2^j)
		for(k in 1:length(rand.smp)){
			p.r <- rand.smp[k]
			p.pre.r <- p.r - scl
			p.pst.r <- p.r + scl
			gc.scn.r <- mean(gnm[p.pre.r:p.pst.r])
			if(k == 1){
				gc.scn.p.r <- gc.scn.r
			}
			else{
				gc.scn.p.r <- c(gc.scn.p.r, gc.scn.r)
			}
		}
		j.m <- mean(gc.scn.p.r)
		if(j == 2){
			j.m.p.r <- j.m
		}
		else{
			j.m.p.r <- c(j.m.p.r, j.m)
		}
	}
	if(i == 1){
		j.m.p.r.p <- j.m.p.r
		# lines((j.m.p.r.p[i]), col = alpha(rgb(0,0,0),0.03), lwd = 1.1)
	}
	else{
		j.m.p.r.p <- rbind(j.m.p.r.p, j.m.p.r)
		# lines((j.m.p.r.p[i,]), col = alpha(rgb(0,0,0),0.03), lwd = 1.1)
	}
}
2^8
n <- 8
pdf(file = paste("Pop_", pop, "_Percent_Variants", 2^n, "bpWindow_Var_Content_July2020.pdf", sep = ""))
hist(j.m.p[,n], breaks = 9, col = "honeydew2", border = "grey30", xlab = "Local SNP Desnity",
	ylab = "Density", main = "Percent Variants at 256bp Window at Insertion Site to Random", freq = F, ylim = c(0,50))
lines(density(j.m.p[,n]), col = cols[pop])
lines(density(j.m.p.p[,n], adjust = 5), col = "grey40")
# abline(v = j.m.p[3], lty = 2)
dev.off()

for(l in 1:ntst){
	ran.seed <- c(0:9999)
	ran.seed <- sample(ran.seed, size = ntst, replace = F)
	set.seed(ran.seed[i])
	rand.nums <- min(population):(length(gcnm)-min(population))
	rand.smp <- sample(rand.nums, size = length(population), replace = F)
	for(j in 2:12){
		scl <- 2^j
		for(i in 1:length(population)){
			p <- rand.smp[j]
			p.pre <- p - scl
			p.pst <- p + scl
			gc.scn <- mean(gnm[p.pre:p.pst])
			if(i == 1){
				gc.scn.p <- gc.scn
			}
			else{
				gc.scn.p <- c(gc.scn.p, gc.scn)
			}
		}
		j.m <- gc.scn.p
		if(j == 2){
			j.m.p <- j.m	
		}
		else{
			j.m.p <- cbind(j.m.p, j.m)
		}
	}
	if(l == 1){
		j.m.p.p <- j.m.p
	}
	else{
		j.m.p.p <- rbind(j.m.p.p, j.m.p)
	}
}








test5 <- sample((c(1:1000)/100), replace = T)
test4 <- sin(2*(c(1:1000)))
test3 <- sin(c(1:1000))
test6 <- rep(0, times = 1000)
test6[400:600] <- 1
test6[800:900] <- -1
test.2 <- wavCWT(c(test6 + test3, test6 + test3, test6 + test3), wavelet = "gaussian2")
plot(test.2, col = cols2, phase = T)
print(test.2)
plot(c(test6 + test3, test6 + test3, test6 + test3), type = "l")


i <- 1
p.t.p <- rep(0, size = ncol(j.m.p.r.p){)
for(i in 1:ncol(j.m.p.r.p)){
	p.test <- c(j.m.p.r.p[,i], j.m.p[i])
	pnorm(abs(scale(p.test)))
	p.t.p[i] <- pnorm(abs(scale(p.test)))[length(p.test)]
}

pdf(file = paste("Pop_", pop, "_boxplots_Of_Window_Polymorphic.pdf", sep = ""), height = 3, width = 8)
boxplot(j.m.p.r.p, 	ylab = "% GC Content", col = cols[pop], border = "grey40",
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()


pdf(file = paste("Pop_", pop, "_p-values_Of_Window_Polymorphic.pdf", sep = ""), height = 3, width = 8)
plot(1-p.t.p, type = "l", col = cols[pop], lwd = 1.9, ylab = "p-values", 
	xlab = "Number of Bases Examined for GC Content in Exponential Scale, 2^n, (n = 2:11)")
dev.off()


known.pos <- map.data.1[,4]

for i in iterations{
	Get a list of known SNP positions and select a subset based on pop dynamics, 
	for each random position add a randomly chosen actual insertion lengh, 
	get nearest and closest snp
	expand out to 
}


############################################################################################
######################## Permutation Bidirection to 5' 3' Insertion ########################
############################################################################################

##### Set initial Seed
set.seed(4332)
ntst <- 1000
i <- 1
j <- 1
l <- 2
m <- 1
pop.dyn <- nrow(test.marks.1.nmatch)
ran.nms <- sample(c(0:9999), size = ntst, replace = F)
insertion <- test.marks.1.nmatch[,5]
insertion <- insertion[insertion > 512]
for(i in 1:ntst){
	##### Set random seed
	set.seed(ran.nms[i])
	##### Get Random Set of SNPs
	pos.nms <- sample(map.data.1[,4], size = pop.dyn, replace = F)
	##### Get random set of Insertion Lengths
	ins.lnth <- sample(insertion, size = pop.dyn, replace = T)
	pop.3p <- pos.nms + ins.lnth
	pop.3p <- pop.3p[pop.3p < max(map.data.1[,4])]
	ins.lnth <- ins.lnth[1:length(pop.3p)]
	pos.nms <- pos.nms[1:length(pop.3p)]
	for(j in 1:length(pop.3p)){
		pop.3.j <- which(abs(map.data.1[,4] - pop.3p[j])==min(abs(map.data.1[,4] - pop.3p[j])))
		pop.3.p <- map.data.1[pop.3.j,4]
		if(j == 1){
			pop.3.p.p <- pop.3.p
		}
		else{
			pop.3.p.p <- c(pop.3.p.p, pop.3.p)
		}
	}
	pop.3.p.p <- pop.3.p.p[1:length(pop.3p)]
	for(l in 2:12){
		scl <- 2^l
		pos.3.pst <- pop.3.p.p + scl
		pos.3.pre <- pop.3.p.p - scl
		pos.5.pst <- pos.nms + scl
		pos.5.pre <- pos.nms - scl
		pos.53.p <- cbind(pos.5.pre, pos.5.pst, pos.3.pre, pos.3.pst)
		length(pos.5.pst)
		pos.53.p <- pos.53.p[pos.53.p[,1] > 0, ]
		pos.53.p <- pos.53.p[pos.53.p[,3] > 0, ]
		for(m in 1:nrow(pos.53.p)){
			mn.5 <- mean(gcnm[pos.53.p[m,1]:pos.53.p[m,2]])
			mn.3 <- mean(gcnm[pos.53.p[m,3]:pos.53.p[m,4]])
			if(m == 1){
				mn.5.p <- mn.5
				mn.3.p <- mn.3
			}
			else{
				mn.5.p <- c(mn.5.p, mn.5)
				mn.3.p <- c(mn.3.p, mn.3)
			}	
		}
		if(l == 2){
			mn.5.p.p <- mean(mn.5.p)
			mn.3.p.p <- mean(mn.3.p)
		}
		else{
			mn.5.p.p <- cbind(mn.5.p.p, mean(mn.5.p))
			mn.3.p.p <- cbind(mn.3.p.p, mean(mn.3.p))
		}
	}
	if(i == 1){
		mn.5.p.p.p <- mn.5.p.p
		mn.3.p.p.p <- mn.3.p.p
	}
	else{
		mn.5.p.p.p <- rbind(mn.5.p.p.p, mn.5.p.p)
		mn.3.p.p.p <- rbind(mn.3.p.p.p, mn.3.p.p)	
	}
}


n = 3
# pdf(file = paste("Pop_", pop, "In-Silico.pdf", sep = ""), width = 8, height = 6)
plot.5 <- mn.5.p.p.p[complete.cases(mn.5.p.p.p[,n]), ]
plot.3 <- mn.3.p.p.p[complete.cases(mn.3.p.p.p[,n]), ]
plot(density(plot.5[,n]), type = "l", col = "grey40")	
lines(density(plot.3[,n]), col = "olivedrab3")
# dev.off()


n = 5
pdf(file = paste("Pop_", pop, "_WindowSize=",4^n , "In-Silico.pdf", sep = ""), width = 8, height = 6)
par(mfrow = c(2,2))
plot.5 <- mn.5.p.p.p[complete.cases(mn.5.p.p.p[,n]), ]
plot.3 <- mn.3.p.p.p[complete.cases(mn.3.p.p.p[,n]), ]
hist(j.m.p.5.1ko[,n], breaks = 50, col = "olivedrab3", border = "grey80", 
	main = paste("GC% at Actual 5' Positions with ", 4^n, "bp Windows", sep=""),
	xlab = "Range of %GC Content")

hist(plot.5[,n], col = "green4", breaks = 50, border = "grey80", 
	main = paste("GC% at Random 5' Positions with ", 4^n, "bp Windows", sep=""),
	xlab = "Range of %GC Content")

hist(j.m.p.3.1ko[,n], col = "grey60", breaks = 50, border = "grey80", 
	main = paste("GC% at Actual 3' Positions with ", 4^n, "bp Windows", sep=""),
	xlab = "Range of %GC Content")

hist(plot.3[,n], breaks = 50, col = "grey40", border = "grey80", 
	main = paste("GC% at Random 3' Positions with ", 4^n, "bp Windows", sep=""),
	xlab = "Range of %GC Content")
dev.off()



hist(mn.m.p.p.p[,6], breaks = 500, col = "olivedrab3", border = "grey40")
hist(j.m.p.5.1ko[,6], breaks = 500, col = "olivedrab3", border = "grey40")


lines(plot)
abline(v = j.m.p.5[n], lty = 2, col = "grey40")
abline(v = j.m.p.3[n], lty = 2, col = "olivedrab3")

n=5
wilcox.test(mn.5.p.p.p[,n],  j.m.p.5.1ko[,n])
ks.test(mn.5.p.p.p[,n],  j.m.p.5.1ko[,n])
t.test(mn.5.p.p.p[,n],  j.m.p.5.1ko[,n])

wilcox.test(mn.3.p.p.p[,n], j.m.p.3.1ko[,n])
ks.test(mn.3.p.p.p[,n], j.m.p.3.1ko[,n])
t.test(mn.3.p.p.p[,n], j.m.p.3.1ko[,n])





var.test(log10(hist.1[,5]), log10(hist.4[,5]), alternative = "greater")



################## 5' get mean across tests ##################
j <- 1
for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,3][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- mean(gc.scn.p)
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- c(j.m.p, j.m)
	}
}
j.m.p.5 <- j.m.p

################## 3' get mean across tests ##################
for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,4][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- mean(gc.scn.p)
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- c(j.m.p, j.m)
	}
}
j.m.p.3 <- j.m.p

################## 5' get lists across tests ##################
i <- 12
j <- 1
for(j in 2:12){
	scl <- 2^j
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,3][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}
j.m.p.5.1ko <- j.m.p
hist(j.m.p.5.1ko[,10], breaks = 50)

################## 3' get lists across tests ##################
for(j in 2:12){
	scl <- 2^j
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,4][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}
j.m.p.3.1ko <- j.m.p

################## 5' get lists across tests ##################
for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,3][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}
j.m.p.5.1ko <- j.m.p


################## 3' get lists across tests ##################
for(j in 2:12){
	scl <- round(2^j)
	for(i in 1:length(test.marks.1.nmatch[,5])){
		p <- test.marks.1.nmatch[,4][i]
		p.pre <- p - scl
		p.pst <- p + scl
		gc.scn <- mean(gcnm[p.pre:p.pst])
		if(i == 1){
			gc.scn.p <- gc.scn
		}
		else{
			gc.scn.p <- c(gc.scn.p, gc.scn)
		}
	}
	j.m <- gc.scn.p
	if(j == 2){
		j.m.p <- j.m	
	}
	else{
		j.m.p <- cbind(j.m.p, j.m)
	}
}
j.m.p.3.1ko <- j.m.p



n = 3
par(mfrow = c(2,1))


n = 5
mx.5 <- max(density(mn.5.p.p.p[,n])$y)
mx.3 <- max(density(mn.3.p.p.p[,n])$y)
hist(j.m.p.5.1k[,n], breaks = 50, col = "olivedrab3", border = "grey30", freq = F, xlim = c(0,1), ylim = c(0,mx.5))
lines(density(j.m.p.5.1k[,n]),  col = "grey50", xlim = c(0,1))
lines(density(mn.5.p.p.p[,n]),  col = "grey20", xlim = c(0,1))
hist(j.m.p.3.1k[,n], breaks = 50, col = "grey80", border = "grey30", freq = F, xlim = c(0,1), ylim = c(0,mx.5))
lines(density(j.m.p.3.1k[,n]), col = "grey20", xlim = c(0,1))
lines(density(mn.3.p.p.p[,n]),  col = "grey20", xlim = c(0,1))


n = 3
hist(mn.5.p.p.p[,n], breaks = 100, col = "olivedrab3", border = "grey40",  
	xlab = paste("bp window = ", 4^n, ", n = 4^",n,  sep = ""),
	main = paste("Five Prime of Random Insertion, n = ", n, sep = ""), freq = F)
	abline(v = j.m.p.5[n], lwd = 0.8, lty = 2, col = "grey40")


n = 3
hist(mn.3.p.p.p[,n], breaks = 100, col = "olivedrab3", border = "grey40",  
	xlab = paste("bp window = ", 4^n, ", n = 4^",n, sep = ""), 
	main = paste("Three Prime of Random Insertion, n = ", n, sep = ""))
	abline(v = j.m.p.3[n], lwd = 0.8, lty = 2, col = "grey40")


mean(test.marks.1.nmatch[,5])


##### Set initial Seed
set.seed(4332)
ntst <- 10000
i <- 1
j <- 1
l <- 2
m <- 1
pop.dyn <- nrow(test.marks.1.nmatch)
ran.nms <- sample(c(0:9999), size = ntst, replace = F)
for(i in 1:ntst){
	##### Set random seed
	set.seed(ran.nms[i])
	##### Get Random Set of SNPs
	pos.nms <- sample(map.data.1[,4], size = pop.dyn, replace = F)
	for(l in 2:12){
		scl <- 2^l
		pos.pst <- pos.nms + scl
		pos.pre <- pos.nms - scl
		pos.mean <- cbind(pos.pre, pos.pst)
		pos.mean <- pos.mean[pos.mean[,1] < max(pos.mean[,1]), ]
		pos.mean <- pos.mean[pos.mean[,2] < max(pos.mean[,2]), ]
		pos.mean <- pos.mean[pos.mean[,1] > 0, ]
		pos.mean <- pos.mean[pos.mean[,2] > 0, ]		
		for(m in 1:nrow(pos.mean)){
			mn.m <- mean(gcnm[pos.mean[m,1]:pos.mean[m,2]])
			if(m == 1){
				mn.m.p <- mn.m
			}
			else{
				mn.m.p <- c(mn.m.p, mn.m)
			}
		}
		if(l == 2){
			mn.m.p.p <- mean(mn.m.p)
		}
		else{
			mn.m.p.p <- cbind(mn.m.p.p, mean(mn.m.p))
		}
	}
	if(i == 1){
		mn.m.p.p.p <- mn.m.p.p
	}
	else{
		mn.m.p.p.p <- rbind(mn.m.p.p.p, mn.m.p.p)
	}
}
hist(mn.m.p.p.p[,6], breaks = 500, col = "olivedrab3", border = "grey40")



###############################################################################################################
################################################ Poisson Analysis #############################################
###############################################################################################################

ppois(length(population),population)
plot(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5], cex = 0.75, pch = 16, col = "grey30")
points(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5], cex = 0.5, pch = 16, col = "olivedrab3")

##### read name file
names.file <- read.table("168E1_x_3A27K3_Prototrophs_18x_names.txt")
paste(as.character(names.file[i,]), "_trmd.bam_unique_mapped.bam_sorted_reads.bam_read_groups.bam", sep = "")

##### Dir to bam files on ORNL BRUT
"/home/ju0/BactAssembly/2019-08-22_2ndRound_Analysis/allsamples_fastq/final_bams"

positions.raw <- paste(test.marks.1.nmatch[,3], test.marks.1.nmatch[,4], sep = "-")
positions.unq <- positions.raw[!duplicated(positions.raw)]
line.commands <- paste("samtools view -o 3a27_2_168_selectLoci.sam" ,"3A27E1_S1_L001_R1_001_trmd.bam_unique_mapped.bam_sorted_reads.bam_read_groups.bam" ,"chr1:", positions.unq, sep = "")

length(positions.unq)

line.commands[1]

pop <- 1
i <- 1
j <- 560
for(j in 1:nrow(test.marks.1.nmatch)){
	i <- as.matrix(test.marks.1.nmatch[j,2])
	samName <- as.matrix(names.file)[i-2]
	line.command <- paste("samtools view -o ", paste("Pop", pop, "_Ind_", samName,"_pedNumb", i-2,
	"_Pos", test.marks.1.nmatch[j,3], "-", test.marks.1.nmatch[j,4], ".sam ", sep = ""),
	paste(samName, "_R1_001_trmd.bam_unique_mapped.bam_sorted_reads.bam_read_groups.bam ", sep = ""),
	paste("chr1:", test.marks.1.nmatch[j,3], "-", test.marks.1.nmatch[j,4], sep = ""), sep = " ")
	if(j == 1){
		line.command.p <- as.character(line.command)
	}
	else{
		line.command.p <- rbind(line.command.p, as.character(line.command))
	}
}
line.command.p[1:4]
dim(line.command.p)
write.table(line.command, "pop1_read_extractionFor_Ref168.sh", quote = F, col.names = F, row.names = F, sep = "	")


###############################################################################################################
######################################## Wavelet Analysis #####################################################
###############################################################################################################

##### Bining Function
wvlt_bin_data <- function(x, bins){
	mx <- max(gnm)
	mxb <- round((mx)/bins)
	for(i in 1:(mxb-1)){
		if(i == 1){
			i.s <- 1
			i.e <- bins
			bn.p <- sum(x[i.s:i.e])
		}
		else{
			i.s <- i*bins
			i.e	<- i.s+bins
			bn <- sum(x[i.s:i.e])
			bn.p <- c(bn.p, bn)
		}
	}
	return(bn.p)
}

wvlt_bin_data2 <- function(x, red){
	mx <- length(x)
	mxb <- round((mx)/red)
	for(i in 1:(mxb-1)){
		if(i == 1){
			i.s <- 1
			i.e <- red
			bn.p <- sum(x[i.s:i.e])
		}
		else{
			i.s <- i*red
			i.e	<- i.s+red
			bn <- sum(x[i.s:i.e])
			bn.p <- c(bn.p, bn)
		}
	}
	return(bn.p)
}

head(depth.1)
plot(s.depth.1[1:10000,4], type = "l", lwd = 1.2)
dim(s.depth.1)

cols2 <- colorRampPalette(c("dodgerblue4", "white", "red2"))( 200 )
var.depth <- depth.1[,4]
var.depth[is.na(var.depth)] <- 0
plot(var.depth, type = "l")
var.depth.t <- var.depth[c(T,F)]
length(var.depth.t)
test.1 <- wavCWT(var.depth.t, wavelet = "gaussian2")
plot(test.1, col = cols2, series = F, phase = T, xlab = "chr",main = paste("Pop", pop, "Variant Coverage Wavelet"))

########## Create GC Content Data ##########
##### Read in GC content file
gc.pos <- read.table("bsub_3a27_GC_positions.txt", header = F)

##### Create blank gc genome
gnm <- rep(0, times = 4010000)

##### Impute GC marks into blank genome
for(i in 1:nrow(gc.pos)){
	pos <- gc.pos[i,]
	gnm[pos] <- 1
}
gnm <- gnm[complete.cases(gnm)]

#### bin gc content 
test.1 <- wvlt_bin_data2(gnm, red = 100)
test.1 <- test.1[complete.cases(test.1)]
test.1 <- test.1/max(test.1)

########## Try Wavelet Comp R Package ##########
my.data <- as.data.frame(test.1)
my.w <- analyze.wavelet(my.data, "test.1", loess.span = 0, dt = 2, dj = 1/250, lowerPeriod = 1, 
		upperPeriod = 128, make.pval = TRUE, n.sim = 10)

##### Create a color Palette
cols3 <- colorRampPalette(c("steelblue3", "lightblue3", "white", "thistle3", "indianred3"))( 200 )

##### Plot row sums
test.1.sum <- rowSums(as.matrix(my.w$Power.pval))
plot(test.1.sum, cex = 0.3, pch = 16, col = "grey30")

##### Create png file of output
png(file = "3a27_GC_Content_3988Bins_501_levels16-64.png", width = 2500, height = 800)
image(t(as.matrix(my.w$Power.pval)), col = cols3, main = "3a27 GC Content 3988 Bins, 501 Levels")
axis(1, tick = T)
dev.off()


##### Population specific Wavelets Plus GC Content
pop.wvt <- colSums(p.1[2:20,])
length(pop.wvt)

##### try wavelet package wmtsa
cols2 <- colorRampPalette(c("dodgerblue4", "white", "red2"))( 200 )
wcwt <- wavCWT(c(test.1, test.1, test.1))
wcwt <- wavCWT(test.1, wavelet = "gaussian2")
plot(wcwt, col = cols2, series = F, phase = T, xlab = "chr",main = paste("Pop", pop, "Recombination_Wavelet"))
par(mfrow = c(1,1))
##### get object structure
str(wcwt)


########## Population specific Wavelets Plus GC Content ##########
pop.wvt <- colSums(p.1[2:20,])
pop.wvt <- pop.wvt[!is.na(names(pop.wvt))]
pop.wvt <- pop.wvt[pop.wvt != 0]
pop.wvt <- pop.wvt/2
plot(pop.wvt)
mx <- names(pop.wvt[length(pop.wvt)])
gnm <- rep(0, times = mx)
names(pop.wvt)

##### Impute markers to genome of zeros
for(i in 1:length(pop.wvt)){
	pos <- as.numeric(as.character(names(pop.wvt[i])))
	gnm[pos] <- pop.wvt[i]
}

##### Examine characteristics of imputed genome
length(gnm)
max(gnm)
which.max(gnm)
dim(as.matrix(gnm))
plot(gnm[1:10000])





########## Wavelet for variants between genomes
parents.delta <- wvlt_bin_data2(gnm, bins = 100)
length(parents.delta)
parents.delta <- parents.delta[complete.cases(parents.delta)]

par.delta.cwt <- wavCWT(parents.delta, wavelet = "gaussian2")
png(file = paste("Parents_Variants_wavelet_Gauss2.png", sep = ""), width = 2900, height = 1500)
plot(par.delta.cwt, col = cols2, phase = T, xlab = paste("Pop ", pop, sep = ""))
dev.off()





##### Bin genome down
markers.1 <- wvlt_bin_data2(gnm, red = 10)
markers.1 <- markers.1[markers.1 <= 19]
markers.1 <- markers.1/max(markers.1)
plot(markers.1)

##### Run Gaussian 2 Wavelet Transform, Continuous Wavelet
markers.1.cwt <- wavCWT(markers.1, wavelet = "gaussian2")
# Plot
pdf(file = paste("Pop_", pop, "_wavelet_Gauss2.pdf", sep = ""), width = 15, height = 5)
plot(markers.1.cwt, col = cols2, phase = T, xlab = paste("Pop ", pop, sep = ""))
dev.off()

##### mak GC and population markers the same length for wavelet analysis
length(test.1)
length(markers.1)
test.1 <- test.1[1:length(markers.1)]

##### Run wavelet on GC plus population markers
gc.pop <- wavCWT((test.1 + markers.1), wavelet = "gaussian2")

##### Create Figure of Population plus GC wavelet
png(file = paste("Pop_", pop, "_And_GC_wavelet_Gauss2.png", sep = ""), width = 900, height = 300)
plot(gc.pop, col = cols2, phase = T, xlab = paste("Pop ", pop, sep = ""))
dev.off()

##### Create Population only wavelet
png(file = paste("Pop_", pop, "_wavelet_Gauss2.png", sep = ""), width = 900, height = 300)
plot(markers.1.cwt, col = cols2, phase = T, xlab = paste("Pop ", pop, sep = ""), type = "image", 
	grid.size = 300)
dev.off()
print(markers.1.cwt)

##### Examine structure of marker file
str(markers.1.cwt)
image(as.matrix(markers.1.cwt))

##### Create Plot of Wavelet peaks
cwt.peak <- rowSums(as.matrix(markers.1.cwt))
plot(cwt.peak, type = "l", lwd = 1.2, col = cols[pop])

pdf(file = paste("Pop_", pop, "_wavelet_topography.pdf", sep = ""), width = 15, height = 5)
barplot(cwt.peak, col =cols[pop], border = F, main = paste("Pop_", pop, "_wavelet_topography", sep = ""),
	xlab = "Genome Position", ylab = "Sum of Coefficients")
lines(cbind(c(0:48000), c(0,0)))
text(cbind(c(0,48000), c(-8, -8)), label = c("0","40100000"))
dev.off()

cwt.peak[cwt.peak < 35] <- 0
cwt.peak <- cwt.peak[20000:length(cwt.peak)]
(which.max(cwt.peak)+20000)*100

cwt.peak <- rowSums(as.matrix(markers.1.cwt))
cwt.peak[cwt.peak < 90] <- 0

which.max(cwt.peak)*100

image(as.matrix(markers.1.cwt[nrow(markers.1.cwt):1,]), col = cols2)

markers.1[markers.1 == 1] <- 0

png(file = paste("Pop_", pop, "_wavelet_peaks_Gauss2.png", sep = ""), width = 1500, height = 1200)
par(mfrow = c(2,1))
plot(markers.1[length(markers.1):1], lwd = 1.2, type ="l", col = cols[pop], xlab = "Positions in Genome (/100bp)")
plot(markers.1.cwt, col = cols2, phase = T, xlab = paste("Pop ", pop, sep = ""))
dev.off()

###### Create Motif Wavelet 
# Create blank gc genome
gnm <- rep(0, times = max(motif))

for(i in 1:length(motif)){
	pos <- motif[i]
	gnm[pos] <- 1
}
gnm <- c(gnm)
plot(gnm[1:120000])
test.1 <- wvlt_bin_data2(gnm, bins = 1000)
test.1 <- test.1[complete.cases(test.1)]
dim(as.matrix(test.1))

my.data <- as.data.frame(test.1)
my.w <- analyze.wavelet(my.data, "test.1", loess.span = 0, dt = 1, dj = 1/250, lowerPeriod = 1, 
		upperPeriod = 128, make.pval = TRUE, n.sim = 10)

cols3 <- colorRampPalette(c("steelblue2", "lightblue", "white", "thistle2", "indianred3"))( 200 )

png(file = "3a27_MotifSites_4002Bins_501_levels1-128.png", width = 2500, height = 800)
image(t(as.matrix(my.w$Power.pval)), col = cols3, main = "3a27 Motif Sites 3988 Bins, 501 Levels")
axis(1, tick = T)
dev.off()



###### Create Population Marker Wavelet 
# Create blank gc genome
dim(p.9)
p.9[,1:20]
p1.wvl <- cbind(as.numeric(colnams(p.9)), colSums(p.9[2:20,]))
p1.wvl <- p1.wvl[complete.cases(p1.wvl[,1]), ]
rw.nms.pls <- max(as.numeric(p1.wvl[,1]))
# buff <- p1.wvl[1:100,]
# buff[1,] <- buff[,1] + rw.nms.pls
# p1.wvl <- rbind(p1.wvl, buff)


gnm <- rep(0, times = max(map.data.1[,4]))
length(gnm)

# Impute GC marks into blank genome
i <- 1
for(i in 1:nrow(p1.wvl)){
	pos <- as.numeric(p1.wvl[i,1])
	pos.val <- as.numeric(p1.wvl[i,2])
	gnm[pos] <- pos.val
}

test.1 <- wvlt_bin_data2(gnm, bins = 1000)
test.1 <- test.1[complete.cases(test.1)]
my.data <- as.data.frame(test.1)
my.w <- analyze.wavelet(my.data, "test.1", loess.span = 0, dt = 1, dj = 1/250, lowerPeriod = 1, 
		upperPeriod = 128, make.pval = TRUE, n.sim = 10)

# wt.image(my.w, color.key = "quantile", n.levels = 250, legend.params = list(lab = "wavelet power levels", mar = 4.7))
cols3 <- colorRampPalette(c("indianred2", "white", "steelblue1"))( 300 )

pop <- 9
par(mfrow = c(1,1))
# # png(file = paste("pop", pop, "_wavelet.png", sep = ""), width = 4000, height = 600)
pdf(file = paste("pop_",pop , "period1-128_.pdf", sep = ""), width = 16, height = 4)
wt.image(my.w, color.key = "quantile", n.levels = 250, legend.params = list(lab = "wavelet power levels", mar = 4.7))
dev.off()


image(t(as.matrix(my.w$Power.pval)), col = cols3)
# dev.off()

image(matrix(c(1:100), nrow = 10, ncol = 10), col = cols3)
pwr.lev <- my.w$Power.pval[,1:100]

max(my.w$Power.pval)
min(my.w$Power.pval)

gradient.1 <- seq(0, 1, 0.01)

image(t(cbind(gradient.1, gradient.1, gradient.1, gradient.1, gradient.1, gradient.1, gradient.1)), col = cols)





#########################################################################################################
####################################### GAGGAC Methylation sites ########################################
#########################################################################################################


motif.f <- read.delim("GAGGAC_3A27_05-26-2020.txt", header = F)
motif.n <- read.delim("GAGGAC_NegStrnd_3a27_05-26-2020.txt", header = F)
colnames(motif.f) <- c("chr", "start", "end", "dot", "1k", "strand")
colnames(motif.n) <- c("chr", "start", "end", "dot", "1k", "strand")
motif.f <- round((motif.f[,2] + motif.f[,3])/2)
motif.n <- round((motif.n[,2] + motif.n[,3])/2)
length(motif.f)
length(motif.n)
motif <- c(motif.f, motif.n)
motif <- motif[complete.cases(motif)]
motif <- motif[order(motif)]

length(motif)
motif.gaggac <- motif
##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 46

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}




##### Get actual distance to methylation sites from population markers
j <- 1
for(j in 1:length(population)){
	mtf.tmp.5 <- motif[motif >= (population[j]+20)]
	mtf.tmp.3 <- motif[motif <= (population[j]-20)]
	motif.data.point.5 <- which(abs(mtf.tmp.5-population[j])==min(abs(mtf.tmp.5-population[j])))
	closest.motif.dist.5 <- ((mtf.tmp.5[motif.data.point] - population[j])^2)^0.5
	motif.data.point.3 <- which(abs(mtf.tmp.3-population[j])==min(abs(mtf.tmp.3-population[j])))
	closest.motif.dist.3 <- ((mtf.tmp.3[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.5.p <- closest.motif.dist.5
		clst.motif.dist.3.p <- closest.motif.dist.3
	}
	else{
		clst.motif.dist.5.p <- c(clst.motif.dist.5.p, closest.motif.dist.5)
		clst.motif.dist.3.p <- c(clst.motif.dist.3.p, closest.motif.dist.3)
	}
}












##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_GGAGAC-Methyl2Insrt_vs_Rand2Inst_pop_July2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)


motif.outlier.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (sd(clst.motif.dist.p)*4)]
random.outlier.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (sd(clst.rand.dist.p.p)*4)]

motif.iqr.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (IQR(clst.motif.dist.p)*1.5)]
random.iqr.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (IQR(clst.rand.dist.p.p)*1.5)]



wilcox.test(motif.outlier.rmvd, random.outlier.rmvd)
t.test(motif.outlier.rmvd, random.outlier.rmvd)
ks.test(motif.outlier.rmvd, random.outlier.rmvd)


hist(motif.outlier.rmvd, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(motif.outlier.rmvd), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(random.outlier.rmvd), col = "bisque4", lwd = 3)




####################################################################################################################
################################## AAGNNNNNNCRTC Methylation site to random ########################################
####################################################################################################################


motif.f <- read.delim("AAGNNNNNNCRTC_3a27_05-26-2020.txt", header = F)
motif.n <- read.delim("AAGNNNNNNCRTC_Negstrand_3a27_05-26-2020.txt", header = F)
colnames(motif.f) <- c("chr", "start", "end", "dot", "1k", "strand")
colnames(motif.n) <- c("chr", "start", "end", "dot", "1k", "strand")
motif.f <- round((motif.f[,2] + motif.f[,3])/2)
motif.n <- round((motif.n[,2] + motif.n[,3])/2)
length(motif.f)
length(motif.n)
motif <- c(motif.f, motif.n)
motif <- motif[complete.cases(motif)]
motif <- motif[order(motif)]

length(motif)
motif.aagnnnnnncrtc <- motif
##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 0909

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}



##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_AAGNNNNNNCRTC-Methyl2Insrt_vs_Rand2Inst_pop_July_2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)


motif.outlier.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (sd(clst.motif.dist.p)*4)]
random.outlier.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (sd(clst.rand.dist.p.p)*4)]

motif.iqr.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (IQR(clst.motif.dist.p)*1.5)]
random.iqr.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (IQR(clst.rand.dist.p.p)*1.5)]



wilcox.test(motif.outlier.rmvd, random.outlier.rmvd)
t.test(motif.outlier.rmvd, random.outlier.rmvd)
ks.test(motif.outlier.rmvd, random.outlier.rmvd)


hist(motif.outlier.rmvd, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(motif.outlier.rmvd), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(random.outlier.rmvd), col = "bisque4", lwd = 3)


#########################################################################################################
####################################### Both Methylation Sites ##########################################
#########################################################################################################

motif <- c(motif.gaggac, motif.aagnnnnnncrtc)

##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 0909

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}



##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_GAGGAC-AAGNNNNNNCRTC-Methyl2Insrt_vs_Rand2Inst_pop_July_2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to GAGGAC and AAGNNNNNNCRTC Methylation and Insertion to Random", breaks = 60)
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)




#########################################################################################################################################
################ Permutation testing for 3a27 to 168 reference markers for recombination, unedited and not complete #####################
#########################################################################################################################################










#########################################################################################################
####################################### GAGGAC Methylation sites ########################################
#########################################################################################################


motif.f <- read.delim("ref168_gaggac_FwdStrnd.txt", header = F)
motif.n <- read.delim("ref168_gaggac_NegStrnd.txt", header = F)
colnames(motif.f) <- c("chr", "start", "end", "dot", "1k", "strand")
colnames(motif.n) <- c("chr", "start", "end", "dot", "1k", "strand")
motif.f <- round((motif.f[,2] + motif.f[,3])/2)
motif.n <- round((motif.n[,2] + motif.n[,3])/2)
length(motif.f)
length(motif.n)
motif <- c(motif.f, motif.n)
motif <- motif[complete.cases(motif)]
motif <- motif[order(motif)]

length(motif)
motif.gaggac <- motif

##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(ref3a27.ref168)

i <- 1
gnm <- rep(0, times = genome_size)
for(i in 1:length(motif.gaggac)){
	pos <- motif.gaggac[i]
	gnm[pos] <- 1
}
plot(gnm[1:10000])


est.n.gaggac <- genome_size/length(motif)
plot(c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)/length(motif.gaggac), type = "l")
##### Set an initial random seed for further randomize future seeds
s.seed <- 46

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- nrow(ref3a27.ref168)

##### Get random numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000

head(ref3a27.ref168)
population <- c(ref3a27.ref168[,3], ref3a27.ref168[,4])


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}
hist(clst.motif.dist.p, breaks = 100, col = "olivedrab3", border = "grey40", freq = F)
lines(density(clst.motif.dist.p), col = "grey30", lwd = 1.8)


##### Get actual distance to methylation sites from population markers
j <- 1
for(j in 1:length(population)){
	mtf.tmp.5 <- motif[motif >= (population[j]+20)]
	mtf.tmp.3 <- motif[motif <= (population[j]-20)]
	motif.data.point.5 <- which(abs(mtf.tmp.5-population[j])==min(abs(mtf.tmp.5-population[j])))
	closest.motif.dist.5 <- ((mtf.tmp.5[motif.data.point] - population[j])^2)^0.5
	motif.data.point.3 <- which(abs(mtf.tmp.3-population[j])==min(abs(mtf.tmp.3-population[j])))
	closest.motif.dist.3 <- ((mtf.tmp.3[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.5.p <- closest.motif.dist.5
		clst.motif.dist.3.p <- closest.motif.dist.3
	}
	else{
		clst.motif.dist.5.p <- c(clst.motif.dist.5.p, closest.motif.dist.5)
		clst.motif.dist.3.p <- c(clst.motif.dist.3.p, closest.motif.dist.3)
	}
}





pop <- 1



##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_GGAGAC-Methyl2Insrt_vs_Rand2Inst_pop_July2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60, )
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)


motif.outlier.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (sd(clst.motif.dist.p)*4)]
random.outlier.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (sd(clst.rand.dist.p.p)*4)]

motif.iqr.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (IQR(clst.motif.dist.p)*1.5)]
random.iqr.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (IQR(clst.rand.dist.p.p)*1.5)]



wilcox.test(motif.outlier.rmvd, random.outlier.rmvd)
t.test(motif.outlier.rmvd, random.outlier.rmvd)
ks.test(motif.outlier.rmvd, random.outlier.rmvd)


hist(motif.outlier.rmvd, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(motif.outlier.rmvd), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(random.outlier.rmvd), col = "bisque4", lwd = 3)



i <- 1
j <- 1
l <- 1
for(j in 1:10){
	j.wnd <- j*1000
	for(i in 1:ntst){
		set.seed(random.pop.nums[i])
		rnd.nms <- c(1:4010000)
		rnd.nms <- sample(rnd.nms, size = length(motif.gaggac), replace = F)
		rnd.nms <- rnd.nms[order(rnd.nms)]
		rnd.nms <- rnd.nms[rnd.nms > j.wnd]
		rnd.nms <- rnd.nms[rnd.nms < length(gnm)-j.wnd]
		for(l in 1:rnd.nms)
			l.range <- c((rnd.nms[l]):(rnd.nms[l]+j.wnd), (rnd.nms[l]):(rnd.nms[l]-j.wnd))

	
	}
}















####################################################################################################################
################################## AAGNNNNNNCRTC Methylation site to random ########################################
####################################################################################################################


motif.f <- read.delim("AAGNNNNNNCRTC_3a27_05-26-2020.txt", header = F)
motif.n <- read.delim("AAGNNNNNNCRTC_Negstrand_3a27_05-26-2020.txt", header = F)
colnames(motif.f) <- c("chr", "start", "end", "dot", "1k", "strand")
colnames(motif.n) <- c("chr", "start", "end", "dot", "1k", "strand")
motif.f <- round((motif.f[,2] + motif.f[,3])/2)
motif.n <- round((motif.n[,2] + motif.n[,3])/2)
length(motif.f)
length(motif.n)
motif <- c(motif.f, motif.n)
motif <- motif[complete.cases(motif)]
motif <- motif[order(motif)]

length(motif)
motif.aagnnnnnncrtc <- motif
##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 0909

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}



##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_AAGNNNNNNCRTC-Methyl2Insrt_vs_Rand2Inst_pop_July_2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)


motif.outlier.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (sd(clst.motif.dist.p)*4)]
random.outlier.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (sd(clst.rand.dist.p.p)*4)]

motif.iqr.rmvd <- clst.motif.dist.p[clst.motif.dist.p < (IQR(clst.motif.dist.p)*1.5)]
random.iqr.rmvd <- clst.rand.dist.p.p[clst.rand.dist.p.p < (IQR(clst.rand.dist.p.p)*1.5)]



wilcox.test(motif.outlier.rmvd, random.outlier.rmvd)
t.test(motif.outlier.rmvd, random.outlier.rmvd)
ks.test(motif.outlier.rmvd, random.outlier.rmvd)


hist(motif.outlier.rmvd, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to Methylation and Insertion to Random", breaks = 60)
lines(density(motif.outlier.rmvd), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(random.outlier.rmvd), col = "bisque4", lwd = 3)


#########################################################################################################
####################################### Both Methylation Sites ##########################################
#########################################################################################################

motif <- c(motif.gaggac, motif.aagnnnnnncrtc)

##### Are methylation sites closer or farther to recombination spots than random?
##### Are gene postions closer or farther to recombination sites than random?

##### Set genome_size to genome size
genome_size <- max(population)

##### Set an initial random seed for further randomize future seeds
s.seed <- 0909

##### Get number of tests per population
data.times <- 9999

##### List of random numbers to choose one for seed set
set.seed(s.seed)
s.num <- c(0:9999)
s.num <- sample(s.num, size = 9999, replace = F)
s.seed <- s.num[1]
set.seed(s.seed)
random.pop.nums <- c(1:9999)
random.pop.nums <- sample(random.pop.nums, size = data.times, replace = F)

##### Get dimensions of data set
sample.size <- length(population)

##### Get randome numbers to length of genome
r.num <- c(1:genome_size)

##### Number of Tests
ntst <- 1000


##### Get actual distance to methylation sites from population markers
for(j in 1:length(motif)){
	motif.data.point <- which(abs(motif-population[j])==min(abs(motif-population[j])))
	closest.motif.dist <- ((motif[motif.data.point] - population[j])^2)^0.5
	if(j == 1){
		clst.motif.dist.p <- closest.motif.dist
	}
	else{
		clst.motif.dist.p <- c(clst.motif.dist.p, closest.motif.dist)
	}
}



##### Generate list of 1066 random positions
plot(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 1.5, ylim = c(0,0.00065))
for(i in 1:ntst){
	set.seed(random.pop.nums[i])
	rnd.nms <- c(1:4010000)
	rnd.nms <- sample(rnd.nms, size = length(motif), replace = F)
	rnd.nms <- rnd.nms[order(rnd.nms)]
	for(j in 1:length(rnd.nms)){
		rand.data.point <- which(abs(rnd.nms-population[j])==min(abs(rnd.nms-population[j])))
		closest.rand.dist <- ((rnd.nms[rand.data.point] - population[j])^2)^0.5
		if(j == 1){
			clst.rand.dist.p <- closest.rand.dist
		}
		else{
			clst.rand.dist.p <- c(clst.rand.dist.p, closest.rand.dist)
		}
	}
	if(i == 1){
		clst.rand.dist.p.p <- clst.rand.dist.p
	}
	else{
	clst.rand.dist.p.p <- c(clst.rand.dist.p.p, clst.rand.dist.p)		
	}
	lines(density(clst.rand.dist.p.p), col = alpha("grey50",0.05), lwd = 1.1)
}
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 5)

##### Write out pdf of methylation to insertion overlayed with random positions to insertions
pdf(file = paste("Actual_GAGGAC-AAGNNNNNNCRTC-Methyl2Insrt_vs_Rand2Inst_pop_July_2020", pop, ".pdf", sep = ""), height = 6, width = 11)
hist(clst.motif.dist.p, col = "honeydew2", border = "grey40", freq = F, xlab = "bases", 
	main = "Histogram of Insertion to GAGGAC and AAGNNNNNNCRTC Methylation and Insertion to Random", breaks = 60)
lines(density(clst.motif.dist.p), type = "l", col = cols[pop], lwd = 4, ylim = c(0,0.0005))
lines(density(clst.rand.dist.p.p), col = "bisque4", lwd = 3)
dev.off()

wilcox.test(clst.motif.dist.p, clst.rand.dist.p.p)
t.test(clst.motif.dist.p, clst.rand.dist.p.p)
ks.test(clst.motif.dist.p, clst.rand.dist.p.p)



################################################################################################################################################
######################################################### Alternate Set of Figure Plots ########################################################
################################################################################################################################################


###############################################################################################################
###################### Create historgrams of insert size by frequency across populations ######################
###############################################################################################################
hist.1 <- test.marks.1.nmatch[!duplicated(paste(test.marks.1.nmatch[,2], "_", test.marks.1.nmatch[,3], sep = "")), ]
hist.2 <- test.marks.2.nmatch[!duplicated(paste(test.marks.2.nmatch[,2], "_", test.marks.2.nmatch[,3], sep = "")), ]
hist.3 <- test.marks.3.nmatch[!duplicated(paste(test.marks.3.nmatch[,2], "_", test.marks.3.nmatch[,3], sep = "")), ]
hist.4 <- test.marks.4.nmatch[!duplicated(paste(test.marks.4.nmatch[,2], "_", test.marks.4.nmatch[,3], sep = "")), ]
hist.5 <- test.marks.5.nmatch[!duplicated(paste(test.marks.5.nmatch[,2], "_", test.marks.5.nmatch[,3], sep = "")), ]
hist.6 <- test.marks.6.nmatch[!duplicated(paste(test.marks.6.nmatch[,2], "_", test.marks.6.nmatch[,3], sep = "")), ]
hist.7 <- test.marks.7.nmatch[!duplicated(paste(test.marks.7.nmatch[,2], "_", test.marks.7.nmatch[,3], sep = "")), ]
hist.8 <- test.marks.8.nmatch[!duplicated(paste(test.marks.8.nmatch[,2], "_", test.marks.8.nmatch[,3], sep = "")), ]
hist.9 <- test.marks.9.nmatch[!duplicated(paste(test.marks.9.nmatch[,2], "_", test.marks.9.nmatch[,3], sep = "")), ]
hist.0 <- test.marks.0.nmatch[!duplicated(paste(test.marks.0.nmatch[,2], "_", test.marks.0.nmatch[,3], sep = "")), ]
hist.149 <- pop.9.t[!duplicated(paste(pop.9.t[,2], "_", pop.9.t[,3], sep = "")), ]

##### Histogram of insertion size
pdf("pop1_bact_insrtsize.png", height = 8, width = 13)
hist(hist.1[,5], breaks = 100, col = cols[1], main = "Pop 1", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop2_bact_insrtsize.png", height = 8, width = 13)
hist(hist.2[,5], breaks = 100, col = cols[2], main = "Pop 2", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop3_bact_insrtsize.png", height = 8, width = 13)
hist(hist.3[,5], breaks = 100, col = cols[3], main = "Pop 3", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop4_bact_insrtsize.png", height = 8, width = 13)
hist(hist.4[,5], breaks = 100, col = cols[4], main = "Pop 4", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop5_bact_insrtsize.png", height = 8, width = 13)
hist(hist.5[,5], breaks = 100, col = cols[5], main = "Pop 5", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop6_bact_insrtsize.png", height = 8, width = 13)
hist(hist.6[,5], breaks = 100, col = cols[6], main = "Pop 6", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop7_bact_insrtsize.png", height = 8, width = 13)
hist(hist.7[,5], breaks = 100, col = cols[7], main = "Pop 7", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop8_bact_insrtsize.png", height = 8, width = 13)
hist(hist.8[,5], breaks = 100, col = cols[8], main = "Pop 8", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop9_bact_insrtsize.png", height = 8, width = 13)
hist(hist.9[,5], breaks = 100, col = cols[9], main = "Pop 9", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop10_bact_insrtsize.png", height = 8, width = 13)
hist(hist.0[,5], breaks = 100, col = cols[10], main = "Pop 10", xlab = "Insert Size (Bases)")
dev.off()

pdf("pop149_bact_insrtsize_InDels.pdf", height = 8, width = 13)
par(mfrow = c(2,1))
hist(log10(hist.149[,5]), breaks =40, col = "Orangered2", main = "Frequency - Pop: 1,4,9", xlab = "Insert Size (Bases)")
hist(log10(hist.149[,5]), breaks = 40, col = "Orangered2", main = "Frequency - Pop :1,4,9", 
	xlab = "Insert Size (Bases)", prob = T, border = "grey40")
lines(density(log10(hist.149[,5]), adjust = 0.5), col = "grey40")
rug(log10(hist.149[,5]), col = "grey40")
dev.off()

par(mfrow = c(1,1))
plot(density(log10(hist.149[,5])))
lines(density(log10(hist.149[,5])), lwd = 1.5, col = "orangered3")

##### Plot density plots
# lines(density(test.marks.1.nmatch[,5]), col = cols[1])
plot(density(test.marks.5.nmatch[,5]), col = cols[5], xlim=range(-4000,150000), ylim=range(0,0.00025), lwd = 2)
lines(density(test.marks.2.nmatch[,5]), col = cols[2], lwd = 2)
lines(density(test.marks.4.nmatch[,5]), col = cols[4], lwd = 2)
lines(density(test.marks.3.nmatch[,5]), col = cols[3], lwd = 2)
lines(density(test.marks.1.nmatch[,5]), col = cols[1], lwd = 2)
# lines(density(test.marks.6.nmatch[,5]), col = cols[6], lwd = 2)
lines(density(test.marks.7.nmatch[,5]), col = cols[7], lwd = 2)
lines(density(test.marks.8.nmatch[,5]), col = cols[8], lwd = 2)
# lines(density(test.marks.9.nmatch[,5]), col = cols[9], lwd = 2)
lines(density(test.marks.0.nmatch[,5]), col = cols[0], lwd = 2)

##### Create jitter plot of insertion sizes
par(mfrow = c(1,1))
max.b.t <- max(c(max(hist.1[,5]), max(hist.2[,5]), max(hist.3[,5]), max(hist.4[,5]), max(hist.5[,5]),
			max(hist.6[,5]), max(hist.7[,5]), max(hist.8[,5]), max(hist.9[,5]), max(hist.0[,5])))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))
# Start plot
pdf("10_popsBsub_jitterplots.png", height = 800, width = 1300)
plot(jitter(rep(0,length(hist.1[,5])),amount=0.2), hist.1[,5],
     xlim=range(-0.5,9.5), ylim=range(-3,max(max.b.t)), main = "Insertion Size of 10 Populations",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[1], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
abline(h = as.integer(tpl), lwd = 0.8, col = "grey40", lty = 2)
abline(h = (as.integer(tpl)/2), lwd = 0.8, col = "grey40", lty = 3)
points(jitter(rep(1,length(hist.2[,5])), amount=0.2), hist.2[,5], cex = 0.75, col = cols[2], pch = 16)
points(jitter(rep(2,length(hist.3[,5])), amount=0.2), hist.3[,5], cex = 0.75, col = cols[3], pch = 16)
points(jitter(rep(3,length(hist.4[,5])), amount=0.2), hist.4[,5], cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(4,length(hist.5[,5])), amount=0.2), hist.5[,5], cex = 0.75, col = cols[5], pch = 16)
points(jitter(rep(5,length(hist.6[,5])), amount=0.2), hist.6[,5], cex = 0.75, col = cols[6], pch = 16)
points(jitter(rep(6,length(hist.7[,5])), amount=0.2), hist.7[,5], cex = 0.75, col = cols[7], pch = 16)
points(jitter(rep(7,length(hist.8[,5])), amount=0.2), hist.8[,5], cex = 0.75, col = cols[8], pch = 16)
points(jitter(rep(8,length(hist.9[,5])), amount=0.2), hist.9[,5], cex = 0.75, col = cols[9], pch = 16)
points(jitter(rep(9,length(hist.0[,5])), amount=0.2), hist.0[,5], cex = 0.75, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()


##### Only 3a27 x 168 populations
par(mfrow = c(1,1))
max.b.t <- max(c(log10(max(hist.1[,5])), log10(max(hist.4[,5])), log10(max(hist.9[,5]))))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))

pdf("ThreePops_3a27x168_InsertSizeBsub_jitterplots_log.pdf", height = 8, width = 13)
plot(jitter(rep(0,length(hist.1[,5])),amount=0.2), log10(hist.1[,5]),
     xlim=range(-0.5,2.5), ylim=range(-3,max(max.b.t)), main = "Insertion Size of Populations 1, 4, and 9",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[1], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
points(jitter(rep(1,length(hist.4[,5])), amount=0.2), log10(hist.4[,5]), cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(2,length(hist.9[,5])), amount=0.2), log10(hist.9[,5]), cex = 0.75, col = cols[9], pch = 16)
axis(2, labels = T, tick = T)
dev.off()


##### Pops 2-7, 3-8, 5-10
max.b.t <- max(c(log10(max(hist.2[,5])), log10(max(hist.7[,5])), log10(max(hist.3[,5])), log10(max(hist.8[,5])), log10(max(hist.5[,5])), log10(max(hist.0[,5]))))
# Define dash line positions
tpl <- nchar(max.b.t) - 1
tpl <- as.numeric(substr(max.b.t, 1, tpl-(tpl-1)))
tpl <- tpl*(10^(as.numeric(nchar(max.b.t)-1)))
boxplot(cbind(log10(hist.2[,5]), log10(hist.7[,5]), log10(hist.3[,5]), log10(hist.8[,5]), log10(hist.5[,5]), log10(hist.0[,5])), col = cols)

boxplot(cbind(hist.2[,5], hist.7[,5], hist.3[,5], hist.8[,5], hist.5[,5], hist.0[,5]), pch = 16, cex = 0.8,
	col = c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10]), border = "grey40")


pdf("ThreePopPairss_2-7_3-8_5-10_InsertSizeBsub_jitterplots.pdf", height = 8, width = 13)
plot(jitter(rep(0,length(hist.2[,5])),amount=0.2), log10(hist.2[,5]),
     xlim=range(-0.5,4), ylim=range(-3,max(max.b.t)), main = "Insertion Size of Populations 2-7, 3-8, and 5-10",
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[2], pch = 16, xlab = "Shuffled Populations", ylab = "Number of Insertions")
points(jitter(rep(0.5,length(hist.7[,5])), amount=0.2), log10(hist.7[,5]), cex = 0.75, col = cols[7], pch = 16)
points(jitter(rep(1.5,length(hist.3[,5])), amount=0.2), log10(hist.3[,5]), cex = 0.75, col = cols[3], pch = 16)
points(jitter(rep(2,length(hist.8[,5])), amount=0.2), log10(hist.8[,5]), cex = 0.75, col = cols[8], pch = 16)
points(jitter(rep(3,length(hist.5[,5])), amount=0.2), log10(hist.5[,5]), cex = 0.75, col = cols[5], pch = 16)
points(jitter(rep(3.5,length(hist.0[,5])), amount=0.2), log10(hist.0[,5]), cex = 0.75, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()


##### Compare two examples
max.b.t <- max(c(max(hist.2[,5]), max(hist.4[,5]), max(hist.5[,5])))
plot(jitter(rep(0,length(hist.2[,5])),amount=0.2), hist.2[,5],
     xlim=range(-0.5, 3.5), ylim=range(-3, max(max.b.t)),
     axes=FALSE,frame.plot=TRUE, cex = 0.75, col = cols[2], pch = 16)
points(jitter(rep(2,length(hist.4[,5])), amount=0.2), hist.4[,5], cex = 0.75, col = cols[4], pch = 16)
points(jitter(rep(3,length(hist.5[,5])), amount=0.2), hist.5[,5], cex = 0.75, col = cols[5], pch = 16)


##### Make GGPLOTs of distributions
install.packages("ggjoy")
library(ggridges)
library(ggplot2)

test.data <- rbind(
			cbind(rep(1, times = length(hist.2[,5])), log10(hist.2[,5])), 
			cbind(rep(2, times = length(hist.7[,5])), log10(hist.7[,5])),
			cbind(rep(4, times = length(hist.3[,5])), log10(hist.3[,5])),
			cbind(rep(5, times = length(hist.8[,5])), log10(hist.8[,5])),
			cbind(rep(7, times = length(hist.5[,5])), log10(hist.5[,5])),
			cbind(rep(8, times = length(hist.0[,5])), log10(hist.0[,5])))
cols4 <- c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10])
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Log_Values")
head(test.data)

pdf(file = "Pops_2-7_3-8_5-10_LogBPDistances.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Log_Values, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "º", point_size = 2, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()


test.data <- rbind(
			cbind(rep(1, times = length(hist.2[,5])), (hist.2[,5])), 
			cbind(rep(2, times = length(hist.7[,5])), (hist.7[,5])),
			cbind(rep(4, times = length(hist.3[,5])), (hist.3[,5])),
			cbind(rep(5, times = length(hist.8[,5])), (hist.8[,5])),
			cbind(rep(7, times = length(hist.5[,5])), (hist.5[,5])),
			cbind(rep(8, times = length(hist.0[,5])), (hist.0[,5])))
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","BasePairs")
ggplot(test.data, aes(x = BasePairs, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 6, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "•", point_size = 5, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()

hist(log10(hist.2[,5]), breaks = 100)
hist(log10(hist.7[,5]), breaks = 100)
wilcox.test(log10(hist.2[,5]),  log10(hist.7[,5]))
ks.test(log10(hist.2[,5]), log10(hist.7[,5]))
t.test(log10(hist.2[,5]), log10(hist.7[,5]))

hist(log10(hist.3[,5]), breaks = 100)
hist(log10(hist.8[,5]), breaks = 100)
plot(density(log10(hist.3[,5])))
plot(density(log10(hist.8[,5])))
wilcox.test(hist.3[,5], hist.8[,5])
ks.test(hist.3[,5], hist.8[,5])
t.test(hist.8[,5], hist.3[,5])

plot(density(log10(as.numeric(hist.0[,5]))))
plot(density(log10(as.numeric(hist.5[,5]))))
wilcox.test(log10(hist.5[,5]),  log10(hist.0[,5]))
ks.test(log10(hist.5[,5]), log10(hist.0[,5]))
t.test(log10(hist.5[,5]), log10(hist.0[,5]))


wilcox.test(log10(hist.1[,5]),  log10(hist.4[,5]))
ks.test(log10(hist.1[,5]), log10(hist.4[,5]))
t.test(log10(hist.1[,5]), log10(hist.4[,5]))


wilcox.test(log10(hist.1[,5]),  log10(hist.9[,5]))
ks.test(log10(hist.1[,5]), log10(hist.9[,5]))
t.test(log10(hist.1[,5]), log10(hist.9[,5]))


wilcox.test(log10(hist.4[,5]),  log10(hist.9[,5]))
ks.test(log10(hist.4[,5]), log10(hist.9[,5]))
t.test(log10(hist.4[,5]), log10(hist.9[,5]))

var.test(log10(hist.1[,5]), log10(hist.4[,5]), alternative = "greater")

###############################################################################################################
####################### Table out the number of insertions per strain in populations ##########################
###############################################################################################################

max.ins.num <- max(c(rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5])), 
	rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5])),
  rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5])), rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5])),
  rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5])), rowSums(table(test.marks.6.nmatch[,2], test.marks.6.nmatch[,5])),
  rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5])), rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5])),
  rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5])), rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))))

t.1 <- rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5]))
t.2 <- rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5]))
t.3 <- rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5]))
t.4 <- rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5]))
t.5 <- rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5]))
t.6 <- rowSums(table(test.marks.6.nmatch[,2], test.marks.6.nmatch[,5]))
t.7 <- rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5]))
t.8 <- rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5]))
t.9 <- rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5]))
t.0 <- rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))
t.149 <- rowSums(table(pop149[,2], pop149[,5]))
hist(149)
##### Pops 149
pdf("Three_3a27x168_NumbOfRecombJitterplots.pdf", height = 8, width = 8)
plot(jitter(rep(0, times = length(t.149)), amount=0.2), t.149,
     	main = "Number of Cross-Over Events pops1,4 and 9 [replace text]",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = "Orangered3", pch = 16, xlab = "Individuals", 	
     	ylab = "Number of Cross-Over Events")
axis(2, labels = T, tick = T)
text(-0.13,70, labels = "Cross-Over Events n = 769", cex = 1.3)
text(-0.13,70, labels = "Cross-Over Events = ", cex = 1.3)
dev.off()


pdf("10_popsBsub_NumbOfRecombJitterplots.png", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.1)), amount=0.2), t.1,
     xlim=range(-0.5,9.5), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[1], pch = 16, xlab = "10 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
abline(h = 190, lwd = 0.8, col = "grey40", lty = 2)
abline(h = (190/2), lwd = 0.8, col = "grey40", lty = 3)
points(jitter(rep(1,length(t.2)), amount=0.2), t.2, cex = 1.5, col = cols[2], pch = 16)
points(jitter(rep(2,length(t.3)), amount=0.2), t.3, cex = 1.5, col = cols[3], pch = 16)
points(jitter(rep(3,length(t.4)), amount=0.2), t.4, cex = 1.5, col = cols[4], pch = 16)
points(jitter(rep(4,length(t.5)), amount=0.2), t.5, cex = 1.5, col = cols[5], pch = 16)
points(jitter(rep(5,length(t.6)), amount=0.2), t.6, cex = 1.5, col = cols[6], pch = 16)
points(jitter(rep(6,length(t.7)), amount=0.2), t.7, cex = 1.5, col = cols[7], pch = 16)
points(jitter(rep(7,length(t.8)), amount=0.2), t.8, cex = 1.5, col = cols[8], pch = 16)
points(jitter(rep(8,length(t.9)), amount=0.2), t.9, cex = 1.5, col = cols[9], pch = 16)
points(jitter(rep(9,length(t.0)), amount=0.2), t.0, cex = 1.5, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()




max.ins.num <- max(c(rowSums(table(test.marks.1.nmatch[,2], test.marks.1.nmatch[,5])), 
 rowSums(table(test.marks.4.nmatch[,2], test.marks.4.nmatch[,5])),
  rowSums(table(test.marks.9.nmatch[,2], test.marks.9.nmatch[,5]))))

pdf("ThreePops_168x3a27_InDels_Bsub_NumbOfRecombJitterplots_2.pdf", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.1)), amount=0.2), t.1,
     xlim=range(-0.5,2.5), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[1], pch = 16, xlab = "168 x 3A27 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
points(jitter(rep(1,length(t.4)), amount=0.2), t.4, cex = 1.5, col = cols[4], pch = 16)
points(jitter(rep(2,length(t.9)), amount=0.2), t.9, cex = 1.5, col = cols[9], pch = 16)
axis(2, labels = T, tick = T)
dev.off()

t.1, t.4, t.9

##### Make dotplot for 2-7, 3-8,  and 5-10 shuffles
max.ins.num <- max(c(rowSums(table(test.marks.2.nmatch[,2], test.marks.2.nmatch[,5])), 
  rowSums(table(test.marks.7.nmatch[,2], test.marks.7.nmatch[,5])),
  rowSums(table(test.marks.3.nmatch[,2], test.marks.3.nmatch[,5])),
  rowSums(table(test.marks.8.nmatch[,2], test.marks.8.nmatch[,5])),
  rowSums(table(test.marks.5.nmatch[,2], test.marks.5.nmatch[,5])),
  rowSums(table(test.marks.0.nmatch[,2], test.marks.0.nmatch[,5]))))


pdf("ThreePops_2-7_3-8_5-10_Bsub_NumbOfRecombJitterplots_2.pdf", height = 8, width = 13)
plot(jitter(rep(0, times = length(t.2)), amount=0.1), t.2,
     xlim=range(-0.5,4), ylim=range(-3, max.ins.num), 
     	main = "Number of Cross-Over Events per Shuffled Population",
     axes=FALSE,frame.plot=TRUE, cex = 1.5, col = cols[2], pch = 16, xlab = "168 x 3A27 Shuffled Populations", 	
     	ylab = "Number of Cross-Over Events")
points(jitter(rep(0.6,length(t.7)), amount=0.1), t.7, cex = 1.5, col = cols[7], pch = 16)
points(jitter(rep(1.6,length(t.3)), amount=0.1), t.3, cex = 1.5, col = cols[3], pch = 16)
points(jitter(rep(2.2,length(t.8)), amount=0.1), t.8, cex = 1.5, col = cols[8], pch = 16)
points(jitter(rep(3.2,length(t.5)), amount=0.1), t.5, cex = 1.5, col = cols[5], pch = 16)
points(jitter(rep(3.8,length(t.0)), amount=0.1), t.0, cex = 1.5, col = cols[10], pch = 16)
axis(2, labels = T, tick = T)
dev.off()




wilcox.test(t.1 ~ t.4)


test.data <- rbind(
			cbind(rep(1, times = length(t.2)), log10(t.2)), 
			cbind(rep(2, times = length(t.7)), log10(t.7)),
			cbind(rep(4, times = length(t.3)), log10(t.3)),
			cbind(rep(5, times = length(t.8)), log10(t.8)),
			cbind(rep(7, times = length(t.5)), log10(t.5)),
			cbind(rep(8, times = length(t.0)), log10(t.0)))
cols4 <- c(cols[2], cols[7], cols[3], cols[8], cols[5], cols[10])
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Log_Insertions")
head(test.data)

pdf(file = "Pops_2-7_3-8_5-10_LogBPDistances.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Log_Insertions, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "º", point_size = 2, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()


test.data <- rbind(
			cbind(rep(1, times = length(t.2)), (t.2)), 
			cbind(rep(2, times = length(t.7)), (t.7)),
			cbind(rep(4, times = length(t.3)), (t.3)),
			cbind(rep(5, times = length(t.8)), (t.8)),
			cbind(rep(7, times = length(t.5)), (t.5)),
			cbind(rep(8, times = length(t.0)), (t.0)))
test.data <- data.frame(test.data)
colnames(test.data) <- c("Shuffle_Population","Insertions")

pdf(file = "Pops_2-7_3-8_5-10_Insertions.pdf", width = 8, height = 5)
ggplot(test.data, aes(x = Insertions, y = Shuffle_Population, group = Shuffle_Population, fill = factor(Shuffle_Population))) + 
	stat_density_ridges(quantile_lines = TRUE, scale = 0.6, quantiles = 4, color = "grey90",
	jittered_points = TRUE,
	position = position_points_jitter(width = 0.05, height = 0.05),
	point_shape = "•", point_size = 5, point_alpha = 2, alpha = 1.7, point_color = "grey40") + 
	theme_ridges(grid = FALSE, center_axis_labels = TRUE) + 
	scale_fill_cyclical(values = c(cols4)) +
	coord_flip()
dev.off()

