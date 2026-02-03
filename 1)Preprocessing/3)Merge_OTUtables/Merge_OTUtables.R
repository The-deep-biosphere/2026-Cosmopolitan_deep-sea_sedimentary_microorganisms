library(plyr)
library(dplyr)
library(scales)
library(phyloseq)

# Set working directory to this file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## --------------------
## PHYLOSEQ OBJECT 1
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/AMOR/Metadata_AMOR.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/AMOR/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/AMOR/OTUtable_decontaminated.csv",sep=",")
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
Phylo1 <- phyloseq(OTU, TAX, samples)
Phylo1 <- subset_samples(Phylo1, Depth!="B")
Phylo1 <- subset_samples(Phylo1, Depth!="R")

## --------------------
## PHYLOSEQ OBJECT 2
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/NWPO/Metadata_NWPO.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/NWPO/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/NWPO/OTUtable_NWPO_decontaminated.csv",sep=",",check.names = FALSE)
names(OTUtable)[1] = "OTU"


# Make sure the OTU numbers are not the same as Phylo object 1
OTU_number = c()
for (i in 1:nrow(OTUtable)){
  OTU = strsplit(OTUtable[i,1], "_")[[1]][2]
  OTU = strtoi(OTU)+100000
  OTU_number = c(OTU_number,paste(c("OTU_", OTU), collapse = ""))
}

OTUtable[,1] = OTU_number
TAXtable[,1] = OTU_number

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
Phylo2 <- phyloseq(OTU, TAX, samples)

## --------------------
## PHYLOSEQ OBJECT 3
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/WNAG/Metadata_WNAG.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/WNAG/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/WNAG/OTUtable_WNAG_decontaminated.csv",sep=",",check.names = FALSE)
names(OTUtable)[1] = "OTU"


# Make sure the OTU numbers are not the same as Phylo object 1
OTU_number = c()
for (i in 1:nrow(OTUtable)){
  OTU = strsplit(OTUtable[i,1], "_")[[1]][2]
  OTU = strtoi(OTU)+200000
  OTU_number = c(OTU_number,paste(c("OTU_", OTU), collapse = ""))
}

OTUtable[,1] = OTU_number
TAXtable[,1] = OTU_number

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
Phylo3 <- phyloseq(OTU, TAX, samples)

## --------------------
## PHYLOSEQ OBJECT 4
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/MAR/Metadata_MAR.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/MAR/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/MAR/OTUtable_MAR_decontaminated.csv",sep=",",check.names = FALSE)
names(OTUtable)[1] = "OTU"


# Make sure the OTU numbers are not the same as Phylo object 1
OTU_number = c()
for (i in 1:nrow(OTUtable)){
  OTU = strsplit(OTUtable[i,1], "_")[[1]][2]
  OTU = strtoi(OTU)+300000
  OTU_number = c(OTU_number,paste(c("OTU_", OTU), collapse = ""))
}

OTUtable[,1] = OTU_number
TAXtable[,1] = OTU_number

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
Phylo4 <- phyloseq(OTU, TAX, samples)


## --------------------
## PHYLOSEQ OBJECT 5
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/SCS/Metadata_SCS.csv",sep=",")
names(Metatable)[1] <- "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/SCS/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/SCS/OTUtable_SCS_decontaminated.csv",sep=",",check.names = FALSE)
names(OTUtable)[1] = "OTU"


# Make sure the OTU numbers are not the same as Phylo object 1
OTU_number = c()
for (i in 1:nrow(OTUtable)){
  OTU = strsplit(OTUtable[i,1], "_")[[1]][2]
  OTU = strtoi(OTU)+400000
  OTU_number = c(OTU_number,paste(c("OTU_", OTU), collapse = ""))
}

OTUtable[,1] = OTU_number
TAXtable[,1] = OTU_number

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
Phylo5 <- phyloseq(OTU, TAX, samples)


## --------------------
## PHYLOSEQ OBJECT 6
##--------------------
# Species level
Metatable <- read.csv("../2)OTUtable_finalization/SPO/Metadata_SPO.csv",sep=",")
names(Metatable)[1] = "Name"
TAXtable <- read.csv("../2)OTUtable_finalization/SPO/taxtable_decontaminated.csv",sep=",")
names(TAXtable)[1] = 'OTU'
OTUtable <- read.csv("../2)OTUtable_finalization/SPO/OTUtable_SPO_decontaminated.csv",sep=",",check.names = FALSE)
names(OTUtable)[1] = "OTU"


# Make sure the OTU numbers are not the same as Phylo object 1
OTU_number = c()
for (i in 1:nrow(OTUtable)){
  OTU = strsplit(OTUtable[i,1], "_")[[1]][2]
  OTU = strtoi(OTU)+500000
  OTU_number = c(OTU_number,paste(c("OTU_", OTU), collapse = ""))
}

OTUtable[,1] = OTU_number
TAXtable[,1] = OTU_number


# Get percentage data
OTUtable <- OTUtable %>%
  tibble::column_to_rownames("OTU")
# OTU_perc <- prop.table(as.matrix(OTUtable), margin = 2)

# Transform into matrices
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
Phylo6 <- phyloseq(OTU, TAX, samples)

## ---------------------------------
## MERGE OBJECT 2,3,4,5,6 with OBJECT 1
##----------------------------------

Phylo <- merge_phyloseq(Phylo1,Phylo2,Phylo3,Phylo4,Phylo5,Phylo6)

nsamples(Phylo1)+nsamples(Phylo2)+nsamples(Phylo3)+nsamples(Phylo4)+nsamples(Phylo5)+nsamples(Phylo6) == nsamples(Phylo)
ntaxa(Phylo1)+ntaxa(Phylo2)+ntaxa(Phylo3)+ntaxa(Phylo4)+ntaxa(Phylo5)+ntaxa(Phylo6) == ntaxa(Phylo)

write.csv(otu_table(Phylo), "../3)Merge_OTUtables/OTUtable_merged.csv")
write.csv(tax_table(Phylo), "../3)Merge_OTUtables/TAXtable_merged.csv")
