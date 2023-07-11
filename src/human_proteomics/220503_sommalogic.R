library(SomaDataIO)
library(openxlsx)
library(heatmaply)
library(proBatch)
library(limma)
library(ggrepel)
#library(structToolbox)

selector <- function(x){
  to.solve <- all.soma[which(all.soma$PhenoID == x),]
  if(nrow(to.solve) > 1){
    rownames(to.solve) <- to.solve$Barcode
    reduced.to.solve <- to.solve[,seq(30,ncol(to.solve))]
    cvs <- apply(reduced.to.solve,1,function(y){sd(y)/mean(y)})
    return(all.soma[which(all.soma$Barcode == names(which(cvs == min(cvs)))),])
  }else if(nrow(to.solve) == 1){
    return(to.solve)
  }else{
    stop("Inconsequent multisample tracking")
  }
  
}
pareto_normalize <- function(data){
  pareto.result <- apply(data, 2, function(x){x / sqrt(sd(x,na.rm = T))})
  return(pareto.result)
}
autoscaling <- function(data){
  autoscaling.result <- apply(data, 2, function(x){x / sd(x,na.rm = T)})
  return(autoscaling.result)
}

amines <- read.xlsx("~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Amines")
acyl <- read.xlsx("~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Acylcarnitines")
tgs <- read.xlsx("~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Positive lipids TG")
non.tgs <- read.xlsx("~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Positive lipids non TG")
soma <- read.xlsx("~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Raw SomaScan")

t <- read_adat("Documents/BPD_info/BPD_patient_clustering/rawData/raw_somascan/CPC-17-191_20171113/CPC-17-191.hybNorm.medNorm.20171113.adat")
tt <- read.table("Documents/BPD_info/BPD_patient_clustering/rawData/raw_somascan/CPC-17-191_20171113/CPC-17-191.hybNorm.plateScale.medNorm.calibrate.bridge.20171113.adat",
                 sep = "\t",skip=10,header = T)

id.table <- read.xlsx("Documents/BPD_info/Protein_IDs_20201012.xlsx")
raw.clinical <- read.xlsx(xlsxFile = "~/Documents/BPD_info/20191001/OMICs/All_omics_BPD.xlsx",sheet = "Clinical data")
clinical.data <- read.xlsx(xlsxFile = "~/Documents/BPD_info/clinical_from_theresa.xlsx")

met.samples <- rbind(
  amines[,1:2],
  acyl[,1:2],
  tgs[,1:2],
  non.tgs[,1:2]
)

met.samples <- met.samples[!duplicated(met.samples),]

reduced.soma <- do.call(rbind,apply(met.samples,1, function(x){
  soma[which(soma$PhenoID %in% x[1] & soma$Zeitpunkt %in% x[2]),]
}))

remain.soma <- soma[which(!soma$PhenoID %in% reduced.soma$PhenoID),]

raw.soma <- rbind(reduced.soma,remain.soma)

reduced.tt <- tt[which(tt$Barcode %in% raw.soma$Barcode),]

reduced.tt <- reduced.tt[order(reduced.tt$Barcode),]
raw.soma <- raw.soma[order(raw.soma$Barcode),]

reduced.tt$Barcode
raw.soma$Barcode

all.soma <- cbind(raw.soma,reduced.tt[,seq(26,ncol(reduced.tt))])
all.soma <- all.soma[order(all.soma$PhenoID),]
colnames(all.soma) <- c(colnames(raw.soma),attributes(t)$Col.Meta$Target)

id.table[7,grep("^Human$",x = id.table[7,],invert = T)]
all.soma <- all.soma[,grep(pattern = "HPV|HIV",x = colnames(all.soma),invert = T)]

list.soma <- do.call(rbind,lapply(unique(all.soma$PhenoID),selector))
soma <- list.soma[,seq(30,ncol(list.soma))]
rownames(soma) <- list.soma$PhenoID

######

any(duplicated(rownames(soma)))

soma.clinical <- do.call(rbind,apply(list.soma[,1:2],1, function(x){
  raw.clinical[which(raw.clinical$PhenoID %in% x[1] & raw.clinical$Sample.time %in% x[2]),]
}))
soma.clinical <- unique(soma.clinical)
soma.clinical <- soma.clinical[match(rownames(soma),soma.clinical$PhenoID),]

soma <- apply(soma, 2, function(x) as.numeric(as.character(x)))

######

soma <- as.data.frame(t(soma))
colnames(soma) <- list.soma$PhenoID

# normalization and plots

norm.soma <- pareto_normalize(soma)
norm.soma <- log2(norm.soma)
boxplot(norm.soma, las=2)
heatmaply(norm.soma,showticklabels = c(TRUE,FALSE))

# batch effect detection
## annotation table creation

clinical.data <- clinical.data[!duplicated(clinical.data$PhenoID),]
clinical.data <- clinical.data[which(clinical.data$PhenoID %in% colnames(soma)),]

list.soma <- list.soma[!duplicated(list.soma),]

clinical.data <- clinical.data[order(clinical.data$PhenoID),]
list.soma <- list.soma[order(list.soma$PhenoID),]

clinical.data$sample.times <- as.numeric(gsub(pattern = "Tag | - Geburt",replacement = "",x = list.soma$Zeitpunkt))
clinical.data$gender <- ifelse(soma.clinical$Gender == "weiblich",0,1)
clinical.data$batch <- as.numeric(gsub(pattern = "Somalogic_Neo_",replacement = "",x = list.soma$Entnahmeantrag))

annot.soma <- clinical.data
colnames(annot.soma)[1] <- "FullRunName"

annot.soma$gender <- as.factor(annot.soma$gender)
annot.soma$BPD <- as.factor(ifelse(annot.soma$BPD_grad > 0,1,0))
annot.soma$BPD_Grad <- as.factor(annot.soma$BPD_grad)
annot.soma$batch <- as.factor(annot.soma$batch)


## ploting batch effec detection - hierarchichal clustering and PCA

array.comp <- c("gender","BPD","BPD_Grad")
plot_hierarchical_clustering(norm.soma,
                             sample_annotation = annot.soma,
                             factors_to_plot = array.comp,
                             distance = 'euclidean', agglomeration = 'complete',
                             label_samples = FALSE)

plot_PCA(data_matrix = norm.soma,sample_annotation = annot.soma,color_by = "batch")

# Detection of outlier samples and remove them

pca.soma <- prcomp(t(norm.soma))$x
plot(pca.soma[,1],pca.soma[,2],cex=0.1)
text(pca.soma[,1],pca.soma[,2],labels=rownames(pca.soma),cex=0.6)

new.soma <- norm.soma[,which(!colnames(norm.soma) %in% c("Mvsm676","Mkuh511"))]

# Differential abundance analysis

annot.soma <- annot.soma[which(!annot.soma$FullRunName %in% c("Mvsm676","Mkuh511")),]
row.names(annot.soma) <- annot.soma$FullRunName

model.soma <- model.matrix(~0+BPD+batch+sample.times+gender+Gestational_age,annot.soma)
fit.soma <- lmFit(new.soma,model.soma)
esoma <- eBayes(fit.soma) 
result.soma <- topTable(esoma, coef=2, adjust="BH", sort.by= "P", number="Inf")
feature.list.soma <- rownames(result.soma[which(result.soma$adj.P.Val < 0.05),])
length(feature.list.soma)

new.soma <- new.soma[order(rownames(new.soma)),]
result.soma <- result.soma[order(rownames(result.soma)),]

m.soma <- cbind(result.soma,new.soma)
write.table(m.soma,"~/Documents/BPD_info/BPD_patient_clustering/omics_tables/2205_omic_tables/220503_sommalogic_info.tsv",
            quote = F,sep = "\t")

# volvanoplot

new.result.soma <- result.soma
new.result.soma$threshold <- new.result.soma$adj.P.Val < 0.05
ggplot(new.result.soma) +
  geom_point(aes(x=logFC, y=-log10(adj.P.Val),colour=threshold)) +
  scale_color_manual(values=c("grey","#E7298A")) +
  geom_text_repel(data =subset(new.result.soma[1:10,], threshold == TRUE), 
                  mapping =aes(x= logFC, y = -log10(adj.P.Val), label = rownames(new.result.soma)[1:10])) +
  xlab("log2 fold change")  +
  ylab("-log10 adjusted p-value") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black") +
  theme_classic(base_size = 22) + theme(legend.position = "none")

######

batch.corrected.soma <- as.data.frame(removeBatchEffect(new.soma,batch = annot.soma$batch))
plot_PCA(data_matrix = batch.corrected.soma,sample_annotation = annot.soma,color_by = "batch")
proteomics <- batch.corrected.soma[feature.list.soma,which(colnames(batch.corrected.soma) %in% annot.soma$FullRunName[which(annot.soma$BPD == 1)])]
