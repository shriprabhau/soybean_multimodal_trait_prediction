##############
###rr-BLUP###
#############

#for calculating genomic heritability without GxE
library(rrBLUP)

#library(snpReady)
#matrix <- G.matrix(genotype, method = "VanRaden", format = "wide", plot = FALSE)
#K <- matrix$Ga
data = read.csv("seed_oil.csv", header = TRUE)

grm <- as.matrix(read.table("seed_oil.rel", header = FALSE))
ids <- read.table("seed_oil.rel.id")
clean_ids <- sub("_(.*)", "", ids$V2)
rownames(grm) <- clean_ids
colnames(grm) <- clean_ids

K <- grm
#resultFT <- kin.blup(data = data, geno = Data.storage.ID, pheno = data$value, K = K)
resultFT <- kin.blup(
  data = data,
  geno = "Data.storage.ID",
  pheno = "value",
  K = K)
resultFT <- kin.blup(
  data = data,
  geno = "Data.storage.ID",
  pheno = "value",
  K = K,
  fixed = "method_name") #with env as fixed effect.
resultFT$Vg
resultFT$Ve
H2_FT <- (resultFT$Vg)/((resultFT$Vg)+(resultFT$Ve))
H2_FT

##################
####Using lme4####
##################
#Custom path for RTools
rtools_path <- "C:/Users/22888361/OneDrive - The University of Western Australia/Documents/rtools44"
# Add RTools bin directories to PATH
Sys.setenv(PATH = paste(
  file.path(rtools_path, "x86_64-w64-mingw32.static.posix", "bin"),
  file.path(rtools_path, "usr", "bin"),
  Sys.getenv("PATH"),
  sep = ";"
))
# Calculate braod sense and narrow sense heritability

library(lme4)
library(tidyverse)

data = read.table("seed_oil.csv", header = T, sep = ",")

#Remove missing values
data = data %>% filter(!is.na(value))

#make a new variable location_year
#data$ENV = paste0(data$Location,"_",data$Year)

#Run model


model <- lmer(value ~ method_name + (1 | Data.storage.ID), data = data) #without GxE effect with env as fixed effect


# Print the model summary to view the variance components
summary(model)

# Extract the variance components
varcomp <- as.data.frame(VarCorr(model))

## Broad sense heritability ##

# Calculate genetic variance (VG) which includes variance due to Line and Location
VG <- varcomp[1, "vcov"] + varcomp[2, "vcov"]  # Genotypic variance (Line + Location)

# Extract the residual (error) variance (Ve)
Ve <- attr(VarCorr(model), "sc")^2  # Residual variance

# Calculate the total phenotypic variance (VP)
VP <- VG + Ve

# Calculate broad-sense heritability (H2)
H2 <- VG / VP

# Output the broad-sense heritability
H2


## narrow sense heritability ##


# Additive genetic variance (from Line)
additive_variance <- varcomp[varcomp$grp == "Data.storage.ID", "vcov"]

# Residual variance (from the model's scale)
residual_variance <- attr(VarCorr(model), "sc")^2  # Residual variance

# Phenotypic variance is the sum of genetic variance and residual variance
phenotypic_variance <- additive_variance + residual_variance

# Calculate narrow-sense heritability (h^2)
h2 <- additive_variance / phenotypic_variance

# Print narrow-sense heritability
print(paste("Narrow-sense heritability: ", round(h2, 3)))