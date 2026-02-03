
# I would like to replace the empty spaces in your tax table with "Unclassified"
library(stringr)

tax <- read.csv("assignments.csv",sep = "\t")

for (i in 1:nrow(tax)){
  if (tax[i,2] == "No hits"){
    tax[i,2:12] <- "Unclassified"
  } else if (tax[i,3] == ""){
    tax[i,3:12] <- "Unclassified"
  } else if (tax[i,4] == ""){
    tax[i,4:12] <- "Unclassified"
  } else if (tax[i,5] == ""){
    kingdom <- paste("Unclassified_", str_squish(tax[i,4]), sep = "")
    tax[i,5:12] <- kingdom
  } else if (tax[i,6] == ""){
    superphylum <- paste("Unclassified_", str_squish(tax[i,5]), sep = "")
    tax[i,6:12] <- superphylum
  } else if (tax[i,7] == ""){
    phylum <- paste("Unclassified_", str_squish(tax[i,6]), sep = "")
    tax[i,7:12] <- phylum
  } else if (tax[i,8] == ""){
    class <- paste("Unclassified_", str_squish(tax[i,7]), sep = "")
    tax[i,8:12] <- class
  } else if (tax[i,9] == ""){
    order <- paste("Unclassified_", str_squish(tax[i,8]), sep = "")
    tax[i,9:12] <- order
  } else if (tax[i,10] == ""){
    family <- paste("Unclassified_", str_squish(tax[i,9]), sep = "")
    tax[i,10:12] <- family
  } else if (tax[i,11] == ""){
    genus <- paste("Unclassified_", str_squish(tax[i,10]), sep = "")
    tax[i,11:12] <- genus
  } else if (tax[i,12] == ""){
    species <- paste("Unclassified_", str_squish(tax[i,11]), sep = "")
    tax[i,12] <- species
  }
}

write.csv(tax,"assignments_with_unclassified.csv", row.names = FALSE)

# for (i in 1:nrow(tax)){
#   if (tax[i,2] == "No hits"){
#     tax[i,2:12] <- "Unclassified"
#   } 
#   if (tax[i,6] == ""){
#     tax[i,6] <- "Unclassified"
#   } 
#   if (tax[i,7] == ""){
#     tax[i,7] <- "Unclassified"
#   } 
#   if (tax[i,8] == ""){
#     tax[i,8] <- "Unclassified"
#   } 
#   if (tax[i,9] == ""){
#     tax[i,9] <- "Unclassified"
#   } 
#   if (tax[i,10] == ""){
#     tax[i,10] <- "Unclassified"
#   } 
#   if (tax[i,11] == ""){
#     tax[i,11] <- "Unclassified"
#   }
#   if (tax[i,12] == ""){
#     tax[i,12] <- "Unclassified"
#   }
# }