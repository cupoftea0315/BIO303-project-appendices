# GEND000181: Transcriptional response to SARS-CoV-2 infection.
# GEND000481: A safe inhalational treatment prevents SARS–CoV-2 viral replication in human airway epithelial cells.

# Read data
rm(list=ls())
pkgs <- c('sva', 'ggplot2', 'RUVSeq','edgeR',"dplyr","pheatmap","ggsci",
          "ggplotify","reshape2","EDASeq","RColorBrewer","cluster","ggvenn")
ins <- lapply(pkgs, library, character.only = TRUE)
load("GEN.Rdata")

# RUVg adjustment---------------------------------------------------------------
# ref: negative control genes
# k: number of surrogate variables
RUVg_res <- function(Data,group,ref = NULL,k = NULL){
  x <- as.factor(group)
  # least significantly DE genes based 
  # on a first-pass DE analysis performed prior to RUVg normalization
  # DE analysis ,top 500 genes
  if(is.null(ref)){
    design <- model.matrix(~x,data=Data)
    y <- DGEList(counts=Data,group = x)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)
    top <- topTags(lrt, n=nrow(Data))$table
    ref <- rownames(Data)[which(!(rownames(Data) %in% rownames(Data)[1:500]))]
  }
  if(is.null(k)){
    # calculating the number of surrogate variables
    k = sva::num.sv(dat = Data, mod = design,method = "be") 
  }
  
  adjusted <- RUVg(as.matrix(Data),ref,k=k,isLog = F)
  return(adjusted)
}
rawdata_combined_filt <- rawdata_combined[which(apply(rawdata_combined,1,
                                                              function(x)length(x[x>5])>=2)),] # Filtering
rawdata_combined_RUVg <- RUVg_res(rawdata_combined_filt,group,ref = NULL,k = 2)
# rawdata_combined_RUVg$normalizedCounts

# RUVs adjustment---------------------------------------------------------------
RUVs_res <- function(Data, group, ref = NULL, k = NULL){
  x <- as.factor(group)
  if(is.null(ref)){
    design <- model.matrix(~x, data=Data)
    y <- DGEList(counts=Data, group = x)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)
    top <- topTags(lrt, n=nrow(Data))$table
    ref <- rownames(Data)[which(!(rownames(Data) %in% rownames(Data)[1:500]))]
  }
  if(is.null(k)){
    k = sva::num.sv(dat = Data, mod = design, method = "be") 
  }
  differences <- makeGroups(as.factor(group))  # a matrix to specifying the replicates
  adjusted <- RUVs(as.matrix(Data), ref, k, differences, isLog = F)
  return(adjusted)
}
rawdata_combined_filt <- rawdata_combined[which(apply(rawdata_combined,1,
                                                                function(x)length(x[x>5])>=2)),] # Filtering
rawdata_combined_RUVs <- RUVs_res(rawdata_combined_filt, group, ref = NULL, k = 2)
# rawdata_combined_RUVs$normalizedCounts

# ComBat adjustment-------------------------------------------------------------
ComBat_res <- function(Data, batch, group){
  Data <- as.matrix(Data)
  Data <- Data[apply(Data, 1, var) > 1, ]   # Remove genes with too low variance
  Data <- Data[complete.cases(Data), ]   # Remove missing value
  batch <- as.factor(batch)
  mod <- model.matrix(~group)
  adjusted <- ComBat(dat = Data, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
  return(adjusted)
}
rawdata_combined_ComBat <- ComBat_res(rawdata_combined, batch, group)

# ComBat-seq adjustment---------------------------------------------------------
ComBatseq_res <- function(Data, batch, group){
  Data <- as.matrix(Data)
  batch <- as.factor(batch)
  adjusted <- ComBat_seq(counts = Data, batch = batch, group = group)
  return(adjusted)
}
rawdata_combined_ComBatseq <- ComBatseq_res(rawdata_combined, batch, group)

# RLE plot----------------------------------------------------------------------
RLEplotfunction <- function(Data,group,batch,method,version){
  Data <- as.matrix(Data)
  plotname <-paste0(getwd(),"/",version,"_",method,"_","RLEplot",".pdf",sep = "")
  pdf(plotname,width = 12.0, height = 6.0)
  n.batch <- nlevels(as.factor(batch))
  nbatch <- as.factor(batch)
  colors <- brewer.pal(n.batch, "Accent")
  plotRLE(Data,col = colors[nbatch],ann=T,xaxt="n",main=method)
  dev.off()
}
version <- "GEN_lung_cov"
RLEplotfunction(rawdata_combined,group,batch,"rawcount",version)
RLEplotfunction(rawdata_combined_RUVg$normalizedCounts,group,batch,"RUVg",version)
RLEplotfunction(rawdata_combined_RUVs$normalizedCounts,group,batch,"RUVs",version)
RLEplotfunction(FPKMdata_combined,group,batch,"FPKM",version)
RLEplotfunction(TPMdata_combined,group,batch,"TPM",version)
RLEplotfunction(rawdata_combined_ComBat, group, batch, "ComBat", version)
RLEplotfunction(rawdata_combined_ComBatseq, group, batch, "ComBat_seq", version)

# PCA analysis------------------------------------------------------------------
pca_correctedfunction <- function(Data,batch,group,method,version){
  pca_corrected_obj = prcomp(Data)
  pca.var = pca_corrected_obj$sdev^2 %>% as.data.frame()
  pca.var$var = round(pca.var$. / sum(pca.var) * 100, 2)
  pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)
  #assign labels to the data frame
  pca_corrected[,"group"] = group
  pca_corrected[,"batch"] = batch
  plotname <-paste0(getwd(),"/",version,"_",method,"_","pca",".png",sep = "")
  png(file=plotname,width=2000,height=2000,res=72*5)
  batch <- as.factor(batch)

  p2 = ggplot(data=pca_corrected, aes(x=PC1,
                                      y=PC2, 
                                      color=batch, shape=group))
  p2 = p2 + geom_point(size=3)
  titlename <- paste0(method)
  p2 = p2 + labs(title=titlename, color="batch", shape="group",
                 x= paste('PC1(', pca.var$var[1],'%)', sep = ''),
                 y= paste('PC2(', pca.var$var[2],'%)', sep = ''))
  plot(p2)
  dev.off()
  return(p2)
}

# conduct pca analysis on samples
version <- "GEN_lung_cov"
rawcount_pca <- pca_correctedfunction(rawdata_combined,batch,group,"rawcount",version)
rawdata_combined_RUVg <- rawdata_combined_RUVg$normalizedCounts
RUVg_pca <- pca_correctedfunction(rawdata_combined_RUVg,batch,group,"RUVg",version)
rawdata_combined_RUVs <- rawdata_combined_RUVs$normalizedCounts
RUVs_pca <- pca_correctedfunction(rawdata_combined_RUVs,batch,group,"RUVs",version)
FPKM_pca <- pca_correctedfunction(FPKMdata_combined,batch,group,"FPKM",version)
TPM_pca <- pca_correctedfunction(TPMdata_combined,batch,group,"TPM",version)
rawcount_ComBat_pca <- pca_correctedfunction(rawdata_combined_ComBat, batch, group, "ComBat", version)
rawcount_ComBatseq_pca <- pca_correctedfunction(rawdata_combined_ComBatseq, batch, group, "ComBat_seq", version)

# Save data as Rdata
save(rawdata_combined, rawdata_combined_RUVg, rawdata_combined_RUVs, rawdata_combined_ComBat, 
     rawdata_combined_ComBatseq, group, batch, file = "GEN_after.Rdata")