## Libraries
library(decontam)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(Biostrings)
library(compositions)
library(zCompositions)
library(ggforce)
library(tibble)       # Needed for converting column to row names

# Set working directory to this file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --------------------------------
# Choose the directory in pre-decontamination folders for one of the sites.
# --------------------------------

# AMOR
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/AMOR/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/AMOR/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/AMOR/Metadata_AMOR.csv",sep=",")
names(Metatable)[1] = "Name"

####
# FIX PROBLEM WITH UNPUBLISHED DATA OUT OF HERE
###

# MAR
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/MAR/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/MAR/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/MAR/Metadata_MAR.csv",sep=",")
names(Metatable)[1] = "Name"

# NWPO
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/NWPO/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/NWPO/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/NWPO/Metadata_NWPO.csv",sep=",")
names(Metatable)[1] = "Name"

# SPO
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/SPO/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/SPO/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/SPO/Metadata_SPO.csv",sep=",")
names(Metatable)[1] = "Name"

# WNAG
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/WNAG/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/WNAG/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/WNAG/Metadata_WNAG.csv",sep=",")
names(Metatable)[1] = "Name"

# SCS
# ------
OTUtable <- read.csv("../2)OTUtable_finalization/SCS/Pre-decontamination/OTUtable.csv",sep="\t")
names(OTUtable)[1] = "OTU"
TAXtable <- read.csv("../2)OTUtable_finalization/SCS/Pre-decontamination/TAXtable.csv",sep=",")
names(TAXtable)[1] = "OTU"
Metatable <- read.csv("../2)OTUtable_finalization/SCS/Metadata_SCS.csv",sep=",")
names(Metatable)[1] = "Name"

# convert the samples names
for (i in 2:length(names(OTUtable))){
  if (names(OTUtable)[i] %in% Metatable$SeqRun){
    row_name = which(Metatable$SeqRun == names(OTUtable)[i])
    names(OTUtable)[i] <- Metatable[row_name,1]
  }
}

# Transform into matrices
OTUtable <- OTUtable %>%
  tibble::column_to_rownames("OTU") 
TAXtable <- TAXtable %>% 
  tibble::column_to_rownames("OTU")
Metatable <- Metatable %>% 
  tibble::column_to_rownames("Name") 
otu_mat <- as.matrix(OTUtable)
tax_mat <- as.matrix(TAXtable)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(Metatable)
AR <- phyloseq(OTU, TAX, samples)

AR_rel <- transform_sample_counts(AR, function(x) x/sum(x)) 
AR
sum(sample_sums(AR)) 

# check which samples are not present in your metadata
for (i in 1:length(sample_names(OTU))){
  if (sample_names(OTU)[i] %in% sample_names(AR)){
  }else {
    print('error')
    print(sample_names(OTU)[i])
  }
}

counts <- data.frame(sample_sums(AR))
counts$Sample <- rownames(counts)
counts$Location <- sample_data(AR)$Core
counts$depth <- sample_data(AR)$Depth
counts$Nreads <- (counts[,1] <= 50)
counts$PercTotal <- 100

## Remove OTUs with no domain assignments
AR_BA <- subset_taxa(AR, DOMAIN %in% c(" Archaea", " Bacteria"))
AR_rel_BA  <- subset_taxa(AR_rel, DOMAIN %in% c(" Archaea", " Bacteria"))
AR_BA
sum(sample_sums(AR_BA)) # Total reads  118619252

#Visualize the amount of reads lost (relatively).
counts$percBA <- sample_sums(AR_BA)*100/counts$sample_sums.AR.
Counts_long <- melt(counts[,2:6])
Counts_long <- Counts_long[-which(Counts_long['variable'] == 'depth'),]
ggplot(Counts_long, aes(x = variable, y = value, group = Sample, color = Location)) + geom_path() + labs(x = 'Step', y = 'Relative amount of reads remaining')

# What was not assigned to the Archaea and Bacteria domain?
AR_rel_nonBA <- subset_taxa(AR_rel, !DOMAIN %in% c(" Archaea", " Bacteria"))
table(as.data.frame(tax_table(AR_rel_nonBA))$DOMAIN)

# Plot the most abundant ones that are not assigned.
plot(sort(taxa_sums(AR_rel_nonBA)), ylab = 'Cumulative relative abundance', xlab = 'OTUs')
# most of the are normally Eukaryota (Chloroplast or Mitochondria)

rm(AR_rel_nonBA)

# Remove list of known contaminant
# This is based on the paper from Eisenhofer in 2019, where there is a list of known lab contaminants. 
# The main idea is to remove all OTUs which genus taxonomic assignment is present in this list.

Eisenhofer <- as.list(read.csv("Pre-decontamination/Contaminants_Eisenhofer.csv", sep = ",", header = FALSE)[,1])
AR_Eisen <- subset_taxa(AR_BA, !GENUS %in% Eisenhofer)
AR_Eisen
sum(sample_sums(AR_Eisen))

counts$percEisen <- sample_sums(AR_Eisen)*100/counts$sample_sums.AR.
Counts_long <- melt(counts[,2:7])
ggplot(Counts_long, aes(x = variable, y = value, group = Sample, color = Location)) + geom_path()

# This often has a major impact, we can have a look which samples are most impacted.

rownames(counts)[which(counts$percEisen <= 85)]

# Which OTUs are the cause?
AR_Eisen_rel_out <- subset_taxa(AR_BA, GENUS %in% Eisenhofer)
AR_Eisen_rel_out_melt <- psmelt(AR_Eisen_rel_out) # Get a long version of the phyloseq object
ggplot(AR_Eisen_rel_out_melt, aes(Sample, Abundance, color = OTU)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~OTU) + theme(legend.position="none")

# Remove all from dataset
AR_Eisen <- subset_taxa(AR_BA, !GENUS %in% Eisenhofer)
AR_rel_Eisen <- subset_taxa(AR_BA, !GENUS %in% Eisenhofer)
AR_Eisen
sum(sample_sums(AR_Eisen))


counts$percEisenNew <- sample_sums(AR_Eisen)*100/counts$sample_sums.AR.
Counts_long <- melt(counts[,2:8])
ggplot(Counts_long, aes(x = variable, y = value, group = Sample, color = Location)) + geom_path()


rm(Counts_long, Eisenhofer, AR_Eisen_rel_out)

write.csv(otu_table(AR_Eisen), "OTUtable_decontaminated.csv")
write.csv(tax_table(AR_Eisen), "taxtable_decontaminated.csv")

