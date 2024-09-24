# samrtpca
## smartpca和flashpca结果完全一致，norm pca 的pve 略小

library(smartsnp) #smartpca r实
library(algatr)
library(vcfR)
library(SNPRelate)
library(SeqArray)
library(here)
library(flashpcaR)





vcf_ld <- read.vcfR("../ld/snp_removeIndividual_ldpruned.vcf")
vcf_ld_dose <- vcf_to_dosage(vcf_ld)
vcf_ld_dose <- simple_impute(vcf_ld_dose, FUN = median)
pca <- prcomp(vcf_ld_dose,scale=T)
plot(pca$x,main ="norm_pca")
std_dev <- pca$sdev
# 计算每个主成分的方差
variance <- std_dev^2
# 计算总方差
total_variance <- sum(variance)
# 计算每个主成分的方差解释比例（PVE）
pve <- variance / total_variance
barplot(pve, main = "Proportion of Variance Explained (PVE)", xlab = "Principal Components", ylab = "PVE")



system("~/biosoft/plink2 --vcf ../ld/snp_removeIndividual_ldpruned.vcf --recode A-transpose --out ../ld/snp_removeIndividual_ldpruned")
numSamples = nrow(read.table("../ld/snp_removeIndividual.fam"))
# There is just a single group in this data
group_id <- rep(c("pop"), length.out = numSamples)
# Running smart_pca
sm.pca <- smart_pca(snp_data = "../ld/snp_removeIndividual_ldpruned.traw", 
                    sample_group = group_id,
                    missing_value = NA,program_svd = "bootSVD") # 计算所有方差
sm.pve <- as.data.frame(sm.pca[["pca.eigenvalues"]])
sm.pve <- as.data.frame(t(sm.pve))
colnames(sm.pve) <- c("eig","pve","cum_pve")
sm.pca_bar <- barplot(sm.pve$pve[1:10],main = "smartpca pve(%)")
text(x= sm.pca_bar,y = sm.pve$pve[1:10], label = round(sm.pve$pve[1:10],2), pos = 1, cex = 0.8, col = "red")
# Here is a plot of the first two components:
plot(sm.pca$pca.sample_coordinates[, c(3,4)],main="sm.pca")


fh.pca <-flashpcaR::flashpca(X = "../ld/snp_removeIndividual_ldpruned",stand = "binom2")
plot(fh.pca$vectors,main = "flashpca")

flashpca_bar <- barplot(fh.pca$pve,main="flashpca pve(%)")
text(x= flashpca_bar,y = fh.pca$pve, label = round(100*fh.pca$pve,2), pos = 1, cex = 0.8, col = "red")





