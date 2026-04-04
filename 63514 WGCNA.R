library(WGCNA)
library(flashClust)
library(tidyverse)
library(dplyr)

setwd("I:/ข้าวหอม (CAFs)/WGCNA 2026")
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#1.Preparing data and removing outlier
data <- read.csv("I:/ข้าวหอม (CAFs)/WGCNA 2026/63514 raw.csv")
data  =as.data.frame(data )
rownames(data ) <- data $X
data $X=NULL
datExpr= as.data.frame(t(data[,]))
names(datExpr)= row.names(data)
rownames(datExpr)=names(data)
dim(datExpr)

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 160, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
keepSamples[c(1,2)]
clust
datExpr = datExpr[keepSamples,]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


datTraits <- read.csv("I:/ข้าวหอม (CAFs)/WGCNA 2026/traits.csv")
rownames(datTraits) <- datTraits$X
datTraits$X <- NULL
datTraits <- t(datTraits)
datTraits <- as.data.frame(datTraits)
datTraits <- datTraits[keepSamples, ]
datTraits


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


save(datExpr, datTraits, file = "Samplepreparing.RData")


#2.Creating Network
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##############################
png("SoftThreshold_WGCNA.png",
    width = 9 * 300,    # 9 inches × 300 dpi
    height = 5 * 300,   # 5 inches × 300 dpi
    res = 300)

par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")

text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     cex = cex1,
     col = "red")

abline(h = 0.90, col = "red")

# Mean connectivity
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")

text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     cex = cex1,
     col = "red")

mtext("A",
      side = 3,
      outer = TRUE,
      line = -2,
      adj = 0.05,
      cex = 2,
      font = 2)

dev.off()
############################################################################

softPower = 14
adjacency = adjacency(datExpr, power = softPower, type = "unsigned") 

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 20
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


MEDissThres = 0.0
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts

save(dynamicMods, moduleLabels,moduleColors,geneTree, MEList, MEs, MEDiss, METree, file= "Clusteringgene2.RData")

#3.Association module and phenotype
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot

###########################################################
png("ModuleTraitHeatmap63514.png",
    width = 2000 ,    
    height = 5000 ,   
    res = 300)
par(mar = c(5, 10, 4, 2))   
par(mfrow = c(1,1))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


mtext("B",
      side = 3,
      outer = TRUE,
      line = -2,
      adj = 0.05,
      cex = 2,
      font = 2)

dev.off()

#############################################################################################
# Calculate module eigengenes
MEs = moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# Choose a module
module = "greenyellow"

# Get the module eigengene for the chosen module
# used to extract the module eigengene for a specific module from the MEs data frame in the context of WGCNA (Weighted Gene Co-expression Network Analysis).
ME = MEs[, paste("ME", module, sep = "")]

# Calculate Module Membership (kME) for the chosen module
kME = as.data.frame(cor(datExpr, ME, use = "p"))

# Calculate Gene Significance (GS) for the trait of interest
GS = as.data.frame(cor(datExpr, datTraits, use = "p"))

# Extract MM and GS for the chosen module
genesInModule = which(moduleColors == module)
MM = kME[genesInModule, 1]
GS = GS[genesInModule, 1]

# Calculate correlation and p-value
correlation = cor(MM, GS, use = "p")
p_value = cor.test(MM, GS)$p.value

# Format p-value for display
formatted_p_value = format(p_value, scientific = TRUE)

# Add the correlation and p-value to the plot
Text <- text(x = max(MM) * 0.8, y = max(GS) * 0.8, 
             labels = paste("correlation =", round(correlation, 2), "\np-value =", formatted_p_value), 
             pos = 4)

# Create a scatter plot

png("Greenyellow.png",
    width = 6 * 300,
    height = 6 * 300,
    res = 300)

par(mfrow = c(1, 1))
par(mar = c(6, 8.5, 3, 3))

plot(MM, GS, 
     xlab = paste("Module Membership in", module, "module"), 
     ylab = "Gene Significance", 
     main = paste("Module Membership vs. Gene Significance \n",
                  "cor =", round(correlation, 2), ", p =", formatted_p_value),
     pch = 20, col = module)

# Add a regression line
abline(lm(GS ~ MM), col = "black")

# Add vertical and horizontal lines
abline(v =0.6, col = "red", lty = 2) # Vertical line at mean of MM
abline(h = 0.5, col = "red", lty = 2) # Horizontal line at mean of GS


mtext("C",
      side = 3,
      outer = TRUE,
      line = -2,
      adj = 0.05,
      cex = 2,
      font = 2)

dev.off()


# Calculate Module Membership (kME) for the chosen module (ต้องรัน2โค้ดนี้เสมอก่อนรันแยกดูแต่ละmodule)
kME1 = as.data.frame(cor(datExpr, ME, use = "p"))
# Calculate Gene Significance (GS) for the trait of interest (if trait is available)
# Replace 'trait' with your actual trait data if applicable
GS1 = as.data.frame(cor(datExpr, datTraits, use = "p"))

# Subset the data to include only genes in the chosen module
greenyellow_indices <- which(moduleColors == "greenyellow")
greenyellow <- names(datExpr)[greenyellow_indices]
# Subset kME and GS to include only genes in the "greenyellow" module
kME_greenyellow <- kME1[greenyellow_indices, , drop = FALSE]
GS_greenyellow <- GS1[greenyellow_indices, , drop = FALSE]

# Filter genes based on MM > 0.5 and GS > 0.03
selected_genes_greenyellow<- greenyellow[kME_greenyellow[, 1] > 0.6 & GS_greenyellow[, 1] > 0.5]

# Print the selected genes
selected_genes_greenyellow <- as.data.frame(selected_genes_greenyellow)
write.table(selected_genes_greenyellow, file = "greenyellow.txt", row.names = FALSE, sep = "\t")


#####################################################################################

# Calculate module eigengenes
MEs = moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# Choose a module
module = "lightcyan"

# Get the module eigengene for the chosen module
# used to extract the module eigengene for a specific module from the MEs data frame in the context of WGCNA (Weighted Gene Co-expression Network Analysis).
ME = MEs[, paste("ME", module, sep = "")]

# Calculate Module Membership (kME) for the chosen module
kME = as.data.frame(cor(datExpr, ME, use = "p"))

# Calculate Gene Significance (GS) for the trait of interest
GS = as.data.frame(cor(datExpr, datTraits, use = "p"))

# Extract MM and GS for the chosen module
genesInModule = which(moduleColors == module)
MM = kME[genesInModule, 1]
GS = GS[genesInModule, 1]

# Calculate correlation and p-value
correlation = cor(MM, GS, use = "p")
p_value = cor.test(MM, GS)$p.value

# Format p-value for display
formatted_p_value = format(p_value, scientific = TRUE)

# Add the correlation and p-value to the plot
Text <- text(x = max(MM) * 0.8, y = max(GS) * 0.8, 
             labels = paste("correlation =", round(correlation, 2), "\np-value =", formatted_p_value), 
             pos = 4)

# Create a scatter plot
png("Lightcyan.png",
    width = 6 * 300,
    height = 6 * 300,
    res = 300)

par(mfrow = c(1, 1))
par(mar = c(6, 8.5, 3, 3))

plot(MM, GS, 
     xlab = paste("Module Membership in", module, "module"), 
     ylab = "Gene Significance", 
     main = paste("Module Membership vs. Gene Significance \n",
                  "cor =", round(correlation, 2), ", p =", formatted_p_value),
     pch = 20, col = "blue")

# Add a regression line
abline(lm(GS ~ MM), col = "black")

# Add vertical and horizontal lines
abline(v =0.6, col = "red", lty = 2) # Vertical line at mean of MM
abline(h = 0.5, col = "red", lty = 2) # Horizontal line at mean of GS

dev.off()


# Calculate Module Membership (kME) for the chosen module (ต้องรัน2โค้ดนี้เสมอก่อนรันแยกดูแต่ละmodule)
kME1 = as.data.frame(cor(datExpr, ME, use = "p"))
# Calculate Gene Significance (GS) for the trait of interest (if trait is available)
# Replace 'trait' with your actual trait data if applicable
GS1 = as.data.frame(cor(datExpr, datTraits, use = "p"))
#1.Black
# Subset the data to include only genes in the chosen module
lightcyan_indices <- which(moduleColors == "lightcyan")
lightcyan <- names(datExpr)[lightcyan_indices]
# Subset kME and GS to include only genes in the "greenyellow" module
kME_lightcyan <- kME1[lightcyan_indices, , drop = FALSE]
GS_lightcyan <- GS1[lightcyan_indices, , drop = FALSE]

# Filter genes based on MM > 0.5 and GS > 0.03
selected_genes_lightcyan<- lightcyan[kME_lightcyan[, 1] > 0.6 & GS_lightcyan[, 1] > 0.5]

# Print the selected genes
selected_genes_lightcyan<- as.data.frame(selected_genes_lightcyan)
write.table(selected_genes_lightcyan, file = "lightcyan.txt", row.names = FALSE, sep = "\t")

#######################################################################################################
# Calculate module eigengenes
MEs = moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# Choose a module
module = "salmon"

# Get the module eigengene for the chosen module
# used to extract the module eigengene for a specific module from the MEs data frame in the context of WGCNA (Weighted Gene Co-expression Network Analysis).
ME = MEs[, paste("ME", module, sep = "")]

# Calculate Module Membership (kME) for the chosen module
kME = as.data.frame(cor(datExpr, ME, use = "p"))

# Calculate Gene Significance (GS) for the trait of interest
GS = as.data.frame(cor(datExpr, datTraits, use = "p"))

# Extract MM and GS for the chosen module
genesInModule = which(moduleColors == module)
MM = kME[genesInModule, 1]
GS = GS[genesInModule, 1]

# Calculate correlation and p-value
correlation = cor(MM, GS, use = "p")
p_value = cor.test(MM, GS)$p.value

# Format p-value for display
formatted_p_value = format(p_value, scientific = TRUE)

# Add the correlation and p-value to the plot
Text <- text(x = max(MM) * 0.8, y = max(GS) * 0.8, 
             labels = paste("correlation =", round(correlation, 2), "\np-value =", formatted_p_value), 
             pos = 4)

# Create a scatter plot
png("Salmon.png",
    width = 6 * 300,
    height = 6 * 300,
    res = 300)

par(mfrow = c(1, 1))
par(mar = c(6, 8.5, 3, 3))

plot(MM, GS, 
     xlab = paste("Module Membership in", module, "module"), 
     ylab = "Gene Significance", 
     main = paste("Module Membership vs. Gene Significance \n",
                  "cor =", round(correlation, 2), ", p =", formatted_p_value),
     pch = 20, col = module)

# Add a regression line
abline(lm(GS ~ MM), col = "black")

# Add vertical and horizontal lines
abline(v =0.6, col = "red", lty = 2) # Vertical line at mean of MM
abline(h = 0.5, col = "red", lty = 2) # Horizontal line at mean of GS

dev.off()

# Calculate Module Membership (kME) for the chosen module (ต้องรัน2โค้ดนี้เสมอก่อนรันแยกดูแต่ละmodule)
kME1 = as.data.frame(cor(datExpr, ME, use = "p"))
# Calculate Gene Significance (GS) for the trait of interest (if trait is available)
# Replace 'trait' with your actual trait data if applicable
GS1 = as.data.frame(cor(datExpr, datTraits, use = "p"))
#3.Pink
# Subset the data to include only genes in the chosen module
salmon_indices <- which(moduleColors == "salmon")
salmon <- names(datExpr)[salmon_indices]
# Subset kME and GS to include only genes in the "greenyellow" module
kME_salmon<- kME1[salmon_indices, , drop = FALSE]
GS_salmon <- GS1[salmon_indices, , drop = FALSE]

# Filter genes based on MM > 0.5 and GS > 0.03
selected_genes_salmon<-salmon[kME_salmon[, 1] > 0.6 & GS_salmon[, 1] > 0.6]

# Print the selected genes
selected_genes_salmon<- as.data.frame(selected_genes_salmon)
write.table(selected_genes_salmon, file = "salmon.txt", row.names = FALSE, sep = "\t")
