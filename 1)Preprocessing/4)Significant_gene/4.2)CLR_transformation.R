library(plyr)
library(dplyr)
library(scales)
library(phyloseq)

# Set working directory to this file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Species level
Metatable <- read.csv("../3)Merge_OTUtables/Metadata_Merged_rm_duplicates.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../3)Merge_OTUtables/TAXtable_merged_family.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../3)Merge_OTUtables/OTUtable_merged_family.csv",sep=",")
names(OTUtable)[1] = "OTU"

# Transform into matrices
OTUtable <- OTUtable %>%
  tibble::column_to_rownames("OTU")
TAXtable <- TAXtable %>%
  tibble::column_to_rownames("OTU")
Metatable <- Metatable %>%
  tibble::column_to_rownames("Name")
tax_mat <- as.matrix(TAXtable)
otu_mat <- as.matrix(OTUtable)

# Make phyloseq object
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Metatable)
Phylo <- phyloseq(OTU, TAX, samples)

OTUtable <- otu_table(Phylo)

# CLR Transformation
# -------------------
logx = log(OTUtable+1)
xclr = logx - outer(rowMeans(logx),rep(1,length=ncol(logx)))
# centring xclr
cxclr = xclr - outer(rep(1,nrow(xclr)), colMeans(xclr))  

# Save file
# ----------
write.csv(cxclr, paste("../3)Merge_OTUtables/OTUtable_merged_family_clr.csv", sep=""))
