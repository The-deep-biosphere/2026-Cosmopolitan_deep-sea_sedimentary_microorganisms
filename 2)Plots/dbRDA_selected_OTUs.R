######## FUNCTION mBiplext ####################
# Compositional biplot, colored by population,
# with additional real variables incorporated.
# Programmed by J.J. Egozcue (2020) based on
# previous function mBPPOP (Dec. 2014)
##### draws a CoDa-biplot with data coming from
# populations coded by a number xpop indicating
# color from a color sequence in colist.
# carries out clr of data set.
# centres the clr and the additional real variables
# in extr. Then they are added to the centered clr.
# carries out svd of centred clr and added real
# extra variables in extr
# plots biplot with data colored by xpop (number of color)
##### input:
# x compositional data by rows (matrix)
# xpop factor indicating the population of each row in x
# extr a vector or a matrix of real variables to be added
#    to the biplot.
# namextr name of a variable (only for single variable added)
# biscale    = 1 covariance biplot
#            = 0 form biplot
# circ       = FALSE (default) ; TRUE plots a unit circle
#              in form biplots
# punch      = 0 do not plot data points
#            = 1 plot symbols for data points
#              following pchpoint
#            = 2 plot numbers for data points following
#              pchpoints
#            = 3 lines between consecutive data points
#              are drawn
# choice[1:2] = PC's to be plotted, eg c(1,2), c(1,3)...
# colist a color sequence for populations.
# pchpoint integer sequence determining the plot pch symbol
#     of each point; when punch=2 is the number to be plotted.
# optpdf = 1  prints on pdf
# filename (optional) defines the name of the output pdf file.
#         by default formBPPOP.pdf or covBPPOP.pdf are used
#### output: a list containing
# the svd matrices: U, V, and singular values in D
# explained variance in the biplot explvar
# additionally
# a pdf file containing the biplot is printed in working dir.
###################################################
if (T){
  mBiplext <- function(x,xpop=NULL,extr=NULL,choice=c(1,2),
                       biscale=1,punch=1,colist=1:10,
                       circ=FALSE,colcirc="grey70",
                       optpdf=0,filename=NULL,
                       namextr=c("tot"),
                       colray="red",colextr="darkgreen",
                       cextext=1,lwdray=1,pchpoint=1){
    # point colors
    colpoint = rep(1,length(colist))
    if(!is.null(xpop)){
      colpoint = colist
    }
    # clr of x
    logx= log(x)
    xclr = logx - outer(rowMeans(logx),rep(1,length=ncol(logx)))
    # centring xclr
    cxclr = xclr - outer(rep(1,nrow(xclr)), colMeans(xclr))  
    # centering real variables extr, if any
    extrue=FALSE
    nextr=NULL
    if(!is.null(extr)){
      if(is.vector(extr)){
        cextr = extr - mean(extr)
        mextr = matrix(cextr,nrow=length(extr),ncol=1)
        colnames(mextr)=namextr
        nextr=1
        extrue = TRUE
      }
      if(is.matrix(extr)){
        namextr = colnames(extr)
        mextr = extr-outer(rep(1,nrow(extr)),colMeans(extr))
        nextr=ncol(mextr) }
      # append real variables in extr  
      cxclr1= cbind(cxclr,mextr)
      colnames(cxclr1)=c(colnames(xclr),namextr)
      cxclr=cxclr1
      extrue = TRUE
    }
    # svd (cxclr)
    SVDxclr = svd(cxclr)
    U = SVDxclr$u
    V = SVDxclr$v
    rownames(V)=colnames(cxclr)
    D = SVDxclr$d
    # scores and loadings
    ## covariance biplot
    if(biscale==1){
      ld=t(diag(D)%*%t(V))/sqrt(nrow(cxclr))
      mainT="covariance biplot"
      fileT="covBiplext"
      mld=max(abs(ld[,choice]))
      msc=max(abs(U[,choice]))
      sc=U *(mld/msc)  # scaling scores
      xylimit=c(-mld,mld)
    }
    
    ## form biplot
    ## scaling: unit norm of V-vectors
    if(biscale==0){
      sc=U%*%diag(D)
      ld=V
      mainT="form biplot"
      fileT="formBiplext"
      mld = max(abs(ld[,choice]))     # scaling basis vectors
      msc = max(abs(sc[,choice]))
      sc = sc*(mld/msc)
      xylimit = c(-mld,mld)
    }
    
    # numeric output
    variances = D^2
    totvar=sum(variances)/(nrow(x)-1)
    extrvar=0
    if(extrue==TRUE){
      extrvar = var(extr) }
    explvar = (variances[choice[1]]+variances[choice[2]])/sum(variances)
    # names
    #  clrnames=paste("clr.",colnames(x),sep="")
    clrnames=colnames(cxclr)
    if(choice[1] == 1){
      expl=100*variances[choice[1]]/sum(variances)
      xlabel=paste("first axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[1] == 2){
      expl=100*variances[choice[1]]/sum(variances)
      xlabel=paste("second axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[1] == 3){
      expl=100*variances[choice[1]]/sum(variances)
      xlabel=paste("third axis,", " var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[1] >= 4){
      expl=100*variances[choice[1]]/sum(variances)
      xlabel=paste(paste(choice[1],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[2] == 1){
      expl=100*variances[choice[2]]/sum(variances)
      ylabel=paste("first axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[2] == 2){
      expl=100*variances[choice[2]]/sum(variances)
      ylabel=paste("second axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[2] == 3){
      expl=100*variances[choice[2]]/sum(variances)
      ylabel=paste("third axis,"," var% ",format(expl,nsmall=1,digits=3),sep="")}
    if(choice[2] >= 4){
      expl=100*variances[choice[2]]/sum(variances)
      ylabel=paste(paste(choice[2],"th axis",sep=""), " var% ",format(expl,nsmall=1,digits=3),sep="")}
    
    if(punch==0){pun="n"}
    if(punch==1){pun="p"}
    if(punch==2){pun="n"}
    if(punch==3){pun="b"}
    
    # pdf output
    filenam=paste(fileT,".pdf",sep="")
    if(optpdf==1){
      if(is.null(filename)==FALSE){filenam=filename}
      pdf(filenam, width=5, height=5, fam="Times")}
    
    plot(sc[,choice],col=colist,type=pun,cex=0.8,asp=1,
         xlim=xylimit,ylim=xylimit,main=mainT,
         xlab=xlabel,ylab=ylabel,pch=pchpoint)
    # only form biplot: unit circle on variables
    if(circ==TRUE & biscale==0){
      theta = seq(from=0, to=(2*pi),length=150)
      xc=cos(theta)
      yc=sin(theta)
      lines(xc,yc,col="grey70")
    }
    
    #  this is for changing punch in place of color
    #  plot(sc[,choice],col="black",type=pun,cex=0.8,asp=1,
    #       pch=colpoint,
    #       xlim=xylimit,ylim=xylimit,main=mainT,
    #       xlab=xlabel,ylab=ylabel)
    
    if(punch==2){
      #    text(sc[,choice],labels=(1:nrow(sc)),col=colpoint,cex=0.8)       
      text(sc[,choice],labels=pchpoint,col=colpoint,cex=0.8)
    }
    for(i in 1:ncol(x)){
      xx=rbind(c(0,0),ld[i,choice])
      lines(xx, col=colray,lwd=lwdray)
      xtext = ld[i,choice[1]]
      ytext = ld[i,choice[2]]
      text(xtext,ytext,labels=clrnames[i],col=colray,
           pos=2,offset=0.3,cex=cextext)
    }
    if(!is.null(nextr)){
      for(iex in 1:nextr){
        nnrow = ncol(x)+iex
        xxetr = rbind(c(0,0),ld[nnrow,choice])
        lines(xxetr,col=colextr,lwd=lwdray)
        xtextr = ld[nnrow,choice[1]]
        ytextr = ld[nnrow,choice[2]]
        text(xtextr,ytextr,labels=clrnames[nnrow],col=colextr,
             pos=2,offset=0.3,cex=cextext)
      } }
    
    if(optpdf==1){
      dev.off()
    }
    lout = list("U"=U,"V"=V,"D"=D,"explvar"=explvar,"totvar"=totvar,
                "extrvar"=extrvar)
    return(lout)
  }
  
  library(tibble)       # Needed for converting column to row names
  library(phyloseq)
  library(ggplot2)
  library(plotly)
  library(vegan)
  library(matrixStats)
  library(ggpubr)
  library(cowplot)
  library(ggforce)
  library(ggrepel)
  library(stringr)

  # Set working directory to this file directory
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  # Species level
  Metatable <- read.csv("../1)Preprocessing/4)Significant_gene/MetaGeoSeqData_with_interpolation_MergedData.csv",sep=",")
  names(Metatable)[1] <- "Name"
  TAXtable <- read.csv("../1)Preprocessing/3)Merge_OTUtables/TAXtable_merged_family.csv",sep=",")
  names(TAXtable)[1] = 'OTU'
  OTUtable <- read.csv("../1)Preprocessing/3)Merge_OTUtables/OTUtable_merged_family.csv",sep=",")
  names(OTUtable)[1] = "OTU"
  
  if (T){
    OTUtable <- OTUtable %>%
      tibble::column_to_rownames("OTU") 
    TAXtable <- TAXtable %>% 
      tibble::column_to_rownames("OTU")
    Metatable <- Metatable %>% 
      tibble::column_to_rownames("Name") 
    tax_mat <- as.matrix(TAXtable)
    otu_mat <- as.matrix(OTUtable)
    
    OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
    TAX = tax_table(tax_mat)
    samples = sample_data(Metatable)
    Phylo <- phyloseq(OTU, TAX, samples)
    
    Phylo <- subset_samples(Phylo, Depth!="B")
    Phylo <- subset_samples(Phylo, Depth!="R")

    # remove fjord and continental shelf cores
    Phylo <- subset_samples(Phylo, Core!=22)
    Phylo <- subset_samples(Phylo, Core!=21)
    Phylo <- subset_samples(Phylo, Core!=19)
    Phylo <- subset_samples(Phylo, Core!=18)
    
    # OTUs based on binning family level using OWN CODE (py)
    # ------------------------------------------------------
    OTU_selection <- c('OTU_1010', 'OTU_10037', 'OTU_10', 'OTU_11724', 'OTU_10033', 'OTU_10146',
                       'OTU_10100', 'OTU_1', 'OTU_1000', 'OTU_10511', 'OTU_10002', 'OTU_10260',
                       'OTU_10512', 'OTU_10054', 'OTU_10022', 'OTU_10246','OTU_11074','OTU_10014','OTU_1033','OTU_10379')   #'OTU_10391','OTU_130'
    
    OTU_selection <- unique(OTU_selection)
    Phylo.subset <- subset_taxa(Phylo, rownames(tax_table(Phylo)) %in% OTU_selection)
    
    # remove unclassified from name
    for (i in 1:length(tax_table(Phylo.subset)[,9])) {
      str <- tax_table(Phylo.subset)[i,9]
      new_str <- gsub('Unclassified_','',str)
      tax_table(Phylo.subset)[i,9] <- new_str
    }
    # replace names with Family taxa
    taxa_names(Phylo.subset) = tax_table(Phylo.subset)[,9]
    
    Phylo_OTU = as(otu_table(Phylo.subset), "matrix")
    OTUtable_no0s <- Phylo_OTU + 1

    # So let's make a PCA on the clr transformed data (ie. aitchison distance).
    PCjuanjo <- mBiplext(t(OTUtable_no0s), extr = NULL, biscale = 0)
    
    # Variance explained by each PC
    expl_var <- (PCjuanjo$D)^2/sum((PCjuanjo$D)^2)
    head((PCjuanjo$D)^2/sum((PCjuanjo$D)^2))
    sum(((PCjuanjo$D)^2/sum((PCjuanjo$D)^2))[1:2])
    # Get the scores for the 3 first PCs
    PCA_scores <- as.data.frame(PCjuanjo$U%*%diag(PCjuanjo$D))[,1:3]
    
    PCA_scores$sample <- row.names(sample_data(Phylo))
    PCA_scores$Loc <- sample_data(Phylo)$Location
    PCA_scores$Depth <- strtoi(sample_data(Phylo)$Depth)
    PCA_scores$core <- sample_data(Phylo)$Core
    PCA_scores$CN <- sample_data(Phylo)$CNratio
    PCA_scores$Corg <- sample_data(Phylo)$Corg
    PCA_scores$Norg <- sample_data(Phylo)$Norg
    PCA_scores$NH4 <- sample_data(Phylo)$Ammonia_interpolated
    PCA_scores$NO2 <- sample_data(Phylo)$NO2_interpolated
    PCA_scores$NO3 <- sample_data(Phylo)$NO3_interpolated
    PCA_scores$O2 <- sample_data(Phylo)$Oxygen_interpolated
    PCA_scores$Mn <- sample_data(Phylo)$Mn_interpolated
    PCA_scores$SO4 <- sample_data(Phylo)$Sulfate_interpolated
    
    # Remember that ggplot shuffles things around, so you may want to factorise some variables in order depending on your needs. See example for the barplot.
    # If you want the loadings, and plot only the top ones, you could make a new data frame with them, and then select only a subselection (eg. the longest arrows). Here is an example for the loadings of the form biplot (which as we saw is not the good plot for looking at the loadings).
    
    # Get the loadings for the 3 firsy PCs
    Loadings <- as.data.frame(PCjuanjo$V[,1:3])
    # Compute the length over the 3 PCs using Pythagorus
    Loadings$length <- sqrt(Loadings$V1^2+Loadings$V2^2+Loadings$V3^2)
    # And select the 20 longest ones (for example).
    Loadings <- Loadings[order(Loadings$length, decreasing = TRUE), ]
    
  }
  
  PCA_scores$site_colour <- PCA_scores$Loc
  
  PCA_scores$site_colour[which(PCA_scores$site_colour == 'SCS')] = 'red'
  PCA_scores$site_colour[which(PCA_scores$site_colour == 'AMOR')] = 'green' 
  PCA_scores$site_colour[which(PCA_scores$site_colour == 'Pacific')] = 'blue'
  PCA_scores$site_colour[which(PCA_scores$site_colour == 'North Pond')] = 'lightblue'
  PCA_scores$site_colour[which(PCA_scores$site_colour == 'Atlantic')] = 'orange'  
}
# ----------
# PCA plot
# ----------
# 
# plt <-ggplot(PCA_scores, aes(x = V1, y = V2, text=sample))+
#   geom_point(aes(col = Loc)) 
# 
# ggplotly(plt)
# 
# plot_ly(x=PCA_scores$V1, y=PCA_scores$V2, z=PCA_scores$V3,
#         type="scatter3d", mode="markers", color=PCA_scores$Loc, size =4)

f <- 30
text_f <- 1.1

aerobic_nitrification <- c("Nitrosopumilaceae","Nitrospinaceae")
aerobic_heterotrophy <- c("Kiloniellaceae","SAR202 clade","Woeseiaceae","JTB23",
                          'Magnetospiraceae','Defluviicoccales')
anearobic_amm_ox <- c("Scalinduaceae","GWA2-50-13")
anearobic_heterotrophy <- c("JS1","Dehalococcoidia","AKAU3564 sediment group",
                            "SG8-4","Aerophobales","Sh765B-AG-111","Desulfosarcinaceae","Desulfatiglandaceae",
                            "GIF3","Lokiarchaeia")

Loadings$col <- 'red'
rownames(Loadings) <- str_squish(rownames(Loadings))
names <- rownames(Loadings)

Loadings$col[which(names %in% aerobic_nitrification)] = '#c9a409'#'deepskyblue'
Loadings$col[which(names %in% aerobic_heterotrophy)] = '#218423'#'green'
Loadings$col[which(names %in% anearobic_amm_ox)] = '#ed2e3b'#'purple'
Loadings$col[which(names %in% anearobic_heterotrophy)] = '#162699'#'brown'

plt1 <-ggplot(PCA_scores, aes(x = V1, y = V2))+
  geom_point(aes(shape = Loc),alpha=0.2) +
  labs(x=paste('PC1 (',round(expl_var[1]*100, digits=1),'%)', sep = ''), y=paste('PC2 (',round(expl_var[2]*100, digits=1),'%)', sep = ''))+
  geom_segment(data = Loadings, aes(x = 0, y = 0, xend = V1*f, yend = V2*f,colour = col), alpha=0.4, size = 0.5, arrow = arrow(length = unit(0.2, "cm")))+
  geom_text_repel(data=Loadings, aes(x= V1*f*text_f, y= V2*f*text_f,  label = rownames(Loadings),colour = col))+
  # xlim(-18,20)+
  theme_bw(base_size = 15) +
  theme(legend.background = element_rect(linetype="solid", 
                                         colour ="black"))+
  labs(colour = "Groups", shape='Sites') +
  scale_colour_manual(values = c("#c9a409", "#218423", "#ed2e3b","#162699",'red'),
                      name="Groups",
                      breaks=c("#c9a409", "#218423", "#ed2e3b","#162699",'red'),
                      labels=c("Aerobic nitrification", "Facultative anaerobic heterotrophy     ", 
                               "Anearobic ammox","Anearobic heterotrophy     ",'random'))
plt1


# ====================
# FAMILY vs Chemistry
# ====================

NH4 <- PCA_scores$NH4
NO3 <- PCA_scores$NO3
NO2 <- PCA_scores$NO2
O2 <- PCA_scores$O2
Mn <- PCA_scores$Mn
SO4 <- PCA_scores$SO4

loc_na <- which(is.na(O2))
loc_na <- c(loc_na,which(is.na(NO3)))
loc_na <- c(loc_na,which(is.na(NO2)))
loc_na <- c(loc_na,which(is.na(NH4)))
loc_na <- c(loc_na,which(is.na(Mn)))
loc_na <- c(loc_na,which(is.na(SO4)))

logx= log(otu_table(Phylo.subset)+1)
xclr = logx - outer(rowMeans(logx),rep(1,length=ncol(logx)))
OTU_matrix <- t(xclr + outer(rep(1,nrow(xclr)), -colMins(xclr)))
OTU_matrix <- OTU_matrix[-loc_na,]

O2 <- O2[-loc_na]
NO3 <- NO3[-loc_na]
NO2 <- NO2[-loc_na]
NH4 <- NH4[-loc_na]
Mn <- Mn[-loc_na]
NOx <- NO3+NO2
SO4 <-SO4[-loc_na]

if (F){
  chem <- as.data.frame(NO3)
  chem <- cbind(chem,NH4,Mn,O2)
  colnames(chem) <- c('NO3', 'NH4','Mn','O2')#
  
  logx= log(chem+1e-130)
  xclr = logx - outer(rowMeans(logx),rep(1,length=ncol(logx)))
  chem <- t(xclr + outer(rep(1,nrow(xclr)), -colMins(as.matrix(xclr))))
  chem <- t(chem)
  
  O2 <- chem[,4]
  NO3 <- chem[,1]
  NH4 <- chem[,2]
  Mn <- chem[,3]
}

colnames(OTU_matrix) <- str_squish(colnames(OTU_matrix))

dbRDA2 <- rda(OTU_matrix ~ O2 + NO3 + NH4 + Mn, dist="euc")
dbRDA2_var <- summary(dbRDA2)$cont$importance[2,]*100
dbRDA2_var

scores(dbRDA2)$sites

plot(dbRDA2, type='n', scaling=3,xlim=c(-4,4),ylim=c(-4,4))
points(dbRDA2, display="sites", pch=20, cex=1, col=PCA_scores$site_colour[-loc_na], scaling=3) # the samples
points(dbRDA2, display="species", pch=4, cex=0.7, col="grey", scaling=3)  # the OTUs
text(dbRDA2, scaling=3, display="bp", col="#5c0b6b", cex=1.1)                           # the predictors
legend("bottomright", legend=c('SCS','AMOR'), bty="n", col=c('red','green'), pch=20, cex=1)

smry <- scores(dbRDA2)
df1  <- data.frame(smry$sites[,1:2])       # PC1 and PC2
taxa  <- data.frame(smry$species[,1:2])     # loadings for PC1 and PC2
biplot <- data.frame(smry$biplot)

aerobic_nitrification <- c("Nitrosopumilaceae","Nitrospinaceae")
aerobic_heterotrophy <- c("Kiloniellaceae","SAR202 clade","Woeseiaceae","JTB23",
                          'Magnetospiraceae','Defluviicoccales')
anearobic_amm_ox <- c("Scalinduaceae","GWA2-50-13")
anearobic_heterotrophy <- c("JS1","Dehalococcoidia","AKAU3564 sediment group",
                            "SG8-4","Aerophobales","Sh765B-AG-111","Desulfosarcinaceae","Desulfatiglandaceae",
                            "GIF3","Lokiarchaeia")
families <- c(aerobic_nitrification,aerobic_heterotrophy,anearobic_amm_ox,anearobic_heterotrophy)

taxa$number <-  1:nrow(taxa)
taxa$col <- 'red'
taxa$cat <- 0

taxa$col[which(rownames(taxa) %in% aerobic_nitrification)] = '#c9a409'#'deepskyblue'
taxa$cat[which(rownames(taxa) %in% aerobic_nitrification)] = 1
taxa$col[which(rownames(taxa) %in% aerobic_heterotrophy)] = '#218423'#'green'
taxa$cat[which(rownames(taxa) %in% aerobic_heterotrophy)] = 2
taxa$col[which(rownames(taxa) %in% anearobic_amm_ox)] = '#ed2e3b'#'purple'
taxa$cat[which(rownames(taxa) %in% anearobic_amm_ox)] = 3
taxa$col[which(rownames(taxa) %in% anearobic_heterotrophy)] = '#162699'#'brown'
taxa$cat[which(rownames(taxa) %in% anearobic_heterotrophy)] = 4

taxa <- taxa[match(families, rownames(taxa)),]
taxa$number <-  1:nrow(taxa)

fac_RDA1 <- 0.5*biplot$RDA1/abs(biplot$RDA1)
fac_RDA2 <- 0.1*biplot$RDA2/abs(biplot$RDA2)

scaling <- 2.
scaling1 <- 3.5
scaling2 <- 1.5
plt3 <-ggplot(taxa, aes(x = RDA1, y = RDA2))+
labs(x=paste0('RDA1 (',round(dbRDA2_var[1],1),'%)'), y=paste0('RDA2 (',round(dbRDA2_var[2],1),'%)'))+
  geom_text(data=biplot, aes(x= RDA1*scaling1*scaling+fac_RDA1, y= RDA2*scaling2*scaling+fac_RDA2,  label = rownames(biplot)))+ # ,  label = rownames(taxa) ,  label = taxa$number , label.size = 0.15
  geom_segment(data = biplot, aes(x = 0, y = 0, xend = RDA1*scaling1*scaling, yend = RDA2*scaling2*scaling), alpha=0.4, size = 0.5, arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(data=taxa, aes(x= RDA1*scaling1*.4, y= RDA2*scaling2*1.5, label = number,colour = col))+ # ,  label = rownames(taxa) ,  label = taxa$number
  theme_bw(base_size = 15)+
  scale_colour_manual(values = c("#c9a409", "#218423", "#ed2e3b","#162699"),
                      name="Groups",
                      breaks=c("#c9a409", "#218423", "#ed2e3b","#162699"),
                      labels=c("Aerobic nitrification", "Facultative Anaerobic heterotrophy    ",
                               "Anearobic ammox","Anearobic heterotrophy    "))+
  ylim(-6,8) +
  scale_x_continuous(breaks = seq(-6, 6.5, by = 2))

plt3 <- plt3 + theme(legend.position="none")
plt3

final_plt <- ggarrange(plt1,plt3, 
                       ncol=2,
                       widths = c(1, .75),
                       labels=c("A","B"),
                       font.label = list(size = 23, color = "black"))


ggsave(filename = "PCA_RDA_notrepelled_withMn.png",
       plot = final_plt,
       width = 40,
       height = 18,
       units = "cm",
       dpi = 300)

