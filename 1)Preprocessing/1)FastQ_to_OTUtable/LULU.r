

library(devtools)
# install_github("tobiasgf/lulu")
require(lulu) # Find here how to install: https://github.com/tobiasgf/lulu
require(methods)

matchlist = read.delim("LULU_match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus.all = read.delim("SCS_all.otutab.sorted.txt",row.names=1,header=T,sep="\t")
curated_result <- lulu(otus.all,matchlist, minimum_match = 97)
lulus = curated_result$curated_table
write.table(data.frame("OTU"=rownames(lulus),lulus),"SCS_table_curated.tsv", row.names=FALSE,
            quote=F, sep="\t")

curated_result$curated_table


