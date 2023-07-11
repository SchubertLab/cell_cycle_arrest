library(ggplot2)
library(gprofiler2)
library(UpSetR)
library(igraph)
library(graphlayouts)
library(ggraph)
library(ComplexHeatmap)
library(gridExtra)
library(grid)
library(lattice)

set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")

#### FUNCTIONS

binox.creator <- function(pathways,DEG){
  to.binox <- data.frame()
  for(a in 1:nrow(pathways)){
    path.genes <- gost(query = pathways[a,"term_id"], organism = "mmusculus")$meta$genes_metadata$query$query_1$ensgs
    to.binox <- rbind(to.binox,data.frame(Node=DEG$ensembl_gene_id[DEG$ensembl_gene_id %in% path.genes],
                                          GroupNames=gsub(pattern = "[[:blank:]]",replacement = "_",
                                                          x =  pathways[a,"term_name"])))
    
  }
  return(to.binox)
}

GOandKEGG <- function(GO,KEGG){
  count.matrix <- read.table("~/Documents/PhD/in_vitro/Count_matrix/Count.matrix.xls",sep = "\t",header = T,stringsAsFactors = F,
                             row.names = 1)
  fisher.matrix <- matrix(0,length(levels(as.factor(KEGG$GroupNames))),length(levels(as.factor(GO$GroupNames))))
  for(i in 1:nrow(fisher.matrix)){
    for(j in 1:ncol(fisher.matrix)){
      kegg.gene.set <- KEGG$Node[KEGG$GroupNames == levels(as.factor(KEGG$GroupNames))[j]]
      go.gene.set <- GO$Node[GO$GroupNames == levels(as.factor(GO$GroupNames))[i]]
      fisher.matrix[i,j] <- fisher.test(matrix(c(length(which(kegg.gene.set %in% go.gene.set)),
                                                 length(go.gene.set) - length(which(kegg.gene.set %in% go.gene.set)),
                                                 length(kegg.gene.set) - length(which(kegg.gene.set %in% go.gene.set)),
                                                 nrow(count.matrix) - length(go.gene.set) - length(kegg.gene.set) + length(which(kegg.gene.set %in% go.gene.set))
      ),2,2),alternative = "greater")$p.value
    }
  }
  
  pre.final <- p.adjust(fisher.matrix,method = "fdr")
  final <- matrix(pre.final,length(levels(as.factor(KEGG$GroupNames))),length(levels(as.factor(GO$GroupNames))))
  rownames(final) <- levels(as.factor(KEGG$GroupNames))
  colnames(final) <- levels(as.factor(GO$GroupNames))
  
  return(final)
}

enrichment_profile = function(pathways, DEG){
  # data is the dataframe containing only significant genes for a specific cell type
  # pathway is a df of pathway ids/go ids that are statistically enriched from geneprofile
  to.plot <- data.frame()
  for(i in 1:nrow(pathways)){
    res = gost(query = pathways[i,"term_id"], organism = "mmusculus")
    gene_set = res$meta$genes_metadata$query$query_1$ensgs
    
    to.plot <- rbind(to.plot,
                     data.frame(
                       Name=pathways[i,"term_name"],
                       mean_FC=c(mean(DEG[DEG$ensembl_gene_id %in% gene_set,]$log2FoldChange[DEG[DEG$ensembl_gene_id %in% gene_set,]$log2FoldChange > 0]),
                                 mean(DEG[DEG$ensembl_gene_id %in% gene_set,]$log2FoldChange[DEG[DEG$ensembl_gene_id %in% gene_set,]$log2FoldChange < 0])),
                       Type=c("up-regulated","down-regulated")
                     ))
  }
  
  to.plot[is.nan(to.plot$mean_FC),"mean_FC"] <- 0
  
  ggplot(to.plot, aes(x = mean_FC, y = reorder(Name, abs(mean_FC)), fill = Type)) +
    geom_bar(stat = "identity")
}

graph.KEGG.GO <- function(KEGG.GO){
  net.kegg.go <- data.frame()
  for(i in seq(nrow(KEGG.GO))){
    for(j in seq(ncol(KEGG.GO))){
      if(KEGG.GO[i,j] < 0.05){
        net.kegg.go <- rbind(net.kegg.go,c(rownames(KEGG.GO)[i],colnames(KEGG.GO)[j],KEGG.GO[i,j]))
      }
    }
  }
  colnames(net.kegg.go) <- c("source","target","p-value")
  
  net.kegg.go$source <- toupper(net.kegg.go$source)
  net.kegg.go$target <- toupper(net.kegg.go$target)
  
  graph <- graph_from_data_frame(net.kegg.go,directed = F,)
  
  node.atr <- rbind(data.frame(name=levels(as.factor(net.kegg.go$source)), source="KEGG"),
                    data.frame(name=levels(as.factor(net.kegg.go$target)), source="GO"))
  rep.node <- node.atr$name[duplicated(node.atr$name)]
  node.atr <- node.atr[!duplicated(node.atr$name),]
  node.atr$source[node.atr$name %in% rep.node] <- "KEGG/GO"
  
  
  graph <- set_vertex_attr(graph = graph,name = "source",value = as.factor(node.atr$source))
  
  V(graph)$source[V(graph)$source == 1] <- "GO"
  V(graph)$source[V(graph)$source == 2] <- "KEGG"
  V(graph)$source[V(graph)$source == 3] <- "KEGG/GO"
  
  ggraph(graph,layout = "focus",
         focus=which(V(graph)$name == names(sort(degree(graph),decreasing = T)[1]))) + 
    geom_edge_link0(edge_colour = "grey66") + 
    geom_node_point(aes(colour=source,size=centrality_degree())) +
    geom_node_text(aes(label = name),size = 3,
                   family = "serif",repel = TRUE) +
    theme_graph()
}

#### in vitro

raw.AT <- read.table("~/Documents/PhD/in_vitro/diffExp/Hyperoxia_AT.xls",sep = "\t",header = T,stringsAsFactors = F)
DEG.AT <- raw.AT[which(raw.AT$padj < 0.05),]
AT.gost <- gost(DEG.AT$ensembl_gene_id,organism = "mmusculus",sources=c("GO:MF", "GO:BP","GO:CC", "KEGG"))
AT.KEGG <- AT.gost$result[AT.gost$result$source == "KEGG",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
AT.BP <- AT.gost$result[AT.gost$result$source == "GO:BP",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
AT.MF <- AT.gost$result[AT.gost$result$source == "GO:MF",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
AT.CC <- AT.gost$result[AT.gost$result$source == "GO:CC",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]

#write.table(x = AT.KEGG,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/AT_KEGG.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = AT.BP,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/AT_BP.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = AT.MF,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/AT_MF.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = AT.CC,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/AT_CC.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)

raw.EC <- read.table("~/Documents/PhD/in_vitro/diffExp/Hyperoxia_EC.xls",sep = "\t",header = T,stringsAsFactors = F)
DEG.EC <- raw.EC[which(raw.EC$padj < 0.05),]
EC.gost <- gost(DEG.EC$ensembl_gene_id,organism = "mmusculus",sources=c("GO:MF", "GO:BP","GO:CC", "KEGG"))
EC.KEGG <- EC.gost$result[EC.gost$result$source == "KEGG",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
EC.BP <- EC.gost$result[EC.gost$result$source == "GO:BP",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
EC.MF <- EC.gost$result[EC.gost$result$source == "GO:MF",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
EC.CC <- EC.gost$result[EC.gost$result$source == "GO:CC",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]

#write.table(x = EC.KEGG,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/EC_KEGG.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = EC.BP,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/EC_BP.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = EC.MF,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/EC_MF.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = EC.CC,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/EC_CC.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)

raw.MFB <- read.table("~/Documents/PhD/in_vitro/diffExp/Hyperoxia_MFB.xls",sep = "\t",header = T,stringsAsFactors = F)
DEG.MFB <- raw.MFB[which(raw.MFB$padj < 0.05),]
MFB.gost <- gost(DEG.MFB$ensembl_gene_id,organism = "mmusculus",sources=c("GO:MF", "GO:BP","GO:CC", "KEGG"))
MFB.KEGG <- MFB.gost$result[MFB.gost$result$source == "KEGG",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
MFB.BP <- MFB.gost$result[MFB.gost$result$source == "GO:BP",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
MFB.MF <- MFB.gost$result[MFB.gost$result$source == "GO:MF",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
MFB.CC <- MFB.gost$result[MFB.gost$result$source == "GO:CC",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]

#write.table(x = MFB.KEGG,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/MFB_KEGG.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = MFB.BP,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/MFB_BP.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = MFB.MF,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/MFB_MF.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = MFB.CC,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/MFB_CC.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)

###### in vivo

raw.IVC <- read.table("~/Documents/in_vitro/diffExp/inVivo_hyperoxia.txt",sep = "\t",header = T,stringsAsFactors = F)
DEG.IVC <- raw.IVC[which(raw.IVC$padj < 0.05),]
IVC.gost <- gost(DEG.IVC$ensembl_gene_id,organism = "mmusculus",sources=c("GO:MF", "GO:BP","GO:CC", "KEGG"))
IVC.KEGG <- IVC.gost$result[IVC.gost$result$source == "KEGG",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVC.BP <- IVC.gost$result[IVC.gost$result$source == "GO:BP",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVC.MF <- IVC.gost$result[IVC.gost$result$source == "GO:MF",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVC.CC <- IVC.gost$result[IVC.gost$result$source == "GO:CC",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]

#write.table(x = IVC.KEGG,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Hyperoxia_KEGG.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVC.BP,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Hyperoxia_BP.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVC.MF,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Hyperoxia_MF.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVC.CC,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Hyperoxia_CC.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)

raw.IVM <- read.table("~/Documents/in_vitro/diffExp/inVivo_mixed.txt",sep = "\t",header = T,stringsAsFactors = F)
DEG.IVM <- raw.IVM[which(raw.IVM$padj < 0.05),]
IVM.gost <- gost(DEG.IVM$ensembl_gene_id,organism = "mmusculus",sources=c("GO:MF", "GO:BP","GO:CC", "KEGG"))
IVM.KEGG <- IVM.gost$result[IVM.gost$result$source == "KEGG",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVM.BP <- IVM.gost$result[IVM.gost$result$source == "GO:BP",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVM.MF <- IVM.gost$result[IVM.gost$result$source == "GO:MF",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]
IVM.CC <- IVM.gost$result[IVM.gost$result$source == "GO:CC",c("term_id", "source", "term_name", "p_value", "query_size", "term_size", "intersection_size")]

#write.table(x = IVM.KEGG,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Mixed_KEGG.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVM.BP,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Mixed_BP.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVM.MF,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Mixed_MF.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)
#write.table(x = IVM.CC,file = "~/Documents/PhD/My_Results/microBulk/report/RNAseq2021/KEGG_and_GO_tables/Mixed_CC.txt",quote = F,
 #           sep = "\t",row.names = F,col.names = T)

###### UpSet plots

# Common genes - Total

frame <- fromList(list(DEG.EC$ensembl_gene_id,
                       DEG.MFB$ensembl_gene_id,
                       DEG.AT$ensembl_gene_id))
                       #DEG.IVM$ensembl_gene_id))

colnames(frame) <- c(paste0("EC"," (N = ",nrow(DEG.EC),")"),
                     paste0("MFB"," (N = ",nrow(DEG.MFB),")"),
                     paste0("ATll"," (N = ",nrow(DEG.AT),")"))
                     #paste0("in vivo"," (N = ",nrow(DEG.IVM),")"))

upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#A79AD4","#DF758C","#ACD49F"),
            sets = colnames(frame),keep.order = T)

# Common genes - UP and DOWN

UP.AT <- DEG.AT[DEG.AT$log2FoldChange > 0,]
DOWN.AT <- DEG.AT[DEG.AT$log2FoldChange < 0,]

UP.EC <- DEG.EC[DEG.EC$log2FoldChange > 0,]
DOWN.EC <- DEG.EC[DEG.EC$log2FoldChange < 0,]

UP.MFB <- DEG.MFB[DEG.MFB$log2FoldChange > 0,]
DOWN.MFB <- DEG.MFB[DEG.MFB$log2FoldChange < 0,]

UP.IV <- DEG.IVM[DEG.IVM$log2FoldChange > 0,]
DOWN.IV <- DEG.IVM[DEG.IVM$log2FoldChange < 0,]

frame <- fromList(list(UP.EC$ensembl_gene_id,
                       DOWN.EC$ensembl_gene_id,
                       UP.MFB$ensembl_gene_id,
                       DOWN.MFB$ensembl_gene_id,UP.AT$ensembl_gene_id,
                       DOWN.AT$ensembl_gene_id))
                       #UP.IV$ensembl_gene_id,
                       #DOWN.IV$ensembl_gene_id))

colnames(frame) <- c(paste0("EC UP"," (N = ",nrow(UP.EC),")"),
                     paste0("EC DOWN"," (N = ",nrow(DOWN.EC),")"),
                     paste0("MFB UP"," (N = ",nrow(UP.MFB),")"),
                     paste0("MFB DOWN"," (N = ",nrow(DOWN.MFB),")"),
                     paste0("ATll UP"," (N = ",nrow(UP.AT),")"),
                     paste0("ATll DOWN"," (N = ",nrow(DOWN.AT),")"))
                     #paste0("in vivo UP"," (N = ",nrow(UP.IV),")"),
                     #paste0("in vivo DOWN"," (N = ",nrow(DOWN.IV),")"))

print(upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#A79AD4","#A79AD4",
                                                                                              "#DF758C","#DF758C",
                                                                                              "#ACD49F","#ACD49F"),
            sets = colnames(frame),keep.order = T,mb.ratio = c(0.35,0.65)))

# KEGG

frame <- fromList(list(AT.KEGG$term_id,
                       EC.KEGG$term_id,
                       MFB.KEGG$term_id,
                       IVM.KEGG$term_id))

colnames(frame) <- c(paste0("ATll"," (N = ",nrow(AT.KEGG),")"),
                     paste0("EC"," (N = ",nrow(EC.KEGG),")"),
                     paste0("MFB"," (N = ",nrow(MFB.KEGG),")"),
                     paste0("in vivo"," (N = ",nrow(IVM.KEGG),")"))

print(upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#ACD49F","#A79AD4","#DF758C","orange"),
            sets = colnames(frame),keep.order = T))


# GO:BP

frame <- fromList(list(AT.BP$term_id,
                       EC.BP$term_id,
                       MFB.BP$term_id,
                       IVM.BP$term_id))

colnames(frame) <- c(paste0("ATll"," (N = ",nrow(AT.BP),")"),
                     paste0("EC"," (N = ",nrow(EC.BP),")"),
                     paste0("MFB"," (N = ",nrow(MFB.BP),")"),
                     paste0("in vivo"," (N = ",nrow(IVM.BP),")"))

print(upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#ACD49F","#A79AD4","#DF758C","orange"),
            sets = colnames(frame),keep.order = T))

# GO:MF

frame <- fromList(list(AT.MF$term_id,
                       EC.MF$term_id,
                       MFB.MF$term_id,
                       IVM.MF$term_id))

colnames(frame) <- c(paste0("ATll"," (N = ",nrow(AT.MF),")"),
                     paste0("EC"," (N = ",nrow(EC.MF),")"),
                     paste0("MFB"," (N = ",nrow(MFB.MF),")"),
                     paste0("in vivo"," (N = ",nrow(IVM.MF),")"))

print(upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#ACD49F","#A79AD4","#DF758C","orange"),
            sets = colnames(frame),keep.order = T))

# CC

frame <- fromList(list(AT.CC$term_id,
                       EC.CC$term_id,
                       MFB.CC$term_id,
                       IVM.CC$term_id))

colnames(frame) <- c(paste0("ATll"," (N = ",nrow(AT.CC),")"),
                     paste0("EC"," (N = ",nrow(EC.CC),")"),
                     paste0("MFB"," (N = ",nrow(MFB.CC),")"),
                     paste0("in vivo"," (N = ",nrow(IVM.CC),")"))

print(upset(frame,nsets = 10, order.by = "freq",empty.intersections = NULL,sets.bar.color = c("#ACD49F","#A79AD4","#DF758C","orange"),
            sets = colnames(frame),keep.order = T))

# Enrichment levels

enrichment_profile(pathways = AT.KEGG,DEG = DEG.AT)
enrichment_profile(pathways = EC.KEGG,DEG = DEG.EC)
enrichment_profile(pathways = MFB.KEGG,DEG = DEG.MFB)
enrichment_profile(pathways = IVC.KEGG,DEG = DEG.IVC)
enrichment_profile(pathways = IVM.KEGG,DEG = DEG.IVM)

enrichment_profile(pathways = head(AT.BP,10),DEG = DEG.AT)
enrichment_profile(pathways = head(EC.BP,10),DEG = DEG.EC)
enrichment_profile(pathways = head(MFB.BP,10),DEG = DEG.MFB)
enrichment_profile(pathways = head(IVC.BP,10),DEG = DEG.IVC)
enrichment_profile(pathways = head(IVM.BP,10),DEG = DEG.IVM)

enrichment_profile(pathways = head(AT.MF,10),DEG = DEG.AT)
enrichment_profile(pathways = head(EC.MF,10),DEG = DEG.EC)
enrichment_profile(pathways = head(MFB.MF,10),DEG = DEG.MFB)
enrichment_profile(pathways = head(IVC.MF,10),DEG = DEG.IVC)
enrichment_profile(pathways = head(IVM.MF,10),DEG = DEG.IVM)

##################### TABLES - BINOX

# Binox files

AT.binox <- binox.creator(pathways = AT.KEGG,DEG = DEG.AT)
EC.binox <- binox.creator(pathways = EC.KEGG,DEG = DEG.EC)
MFB.binox <- binox.creator(pathways = MFB.KEGG,DEG = DEG.MFB)
IVC.binox <- binox.creator(pathways = IVC.KEGG,DEG = DEG.IVC)
IVM.binox <- binox.creator(pathways = IVM.KEGG,DEG = DEG.IVM)

#write.table(AT.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/AT_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(EC.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/EC_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(MFB.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/MFB_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(IVC.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/IVH_Hyperoxia.tsv",quote = F,sep = "\t",row.names = F)
#write.table(IVM.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/IVM_Hyperoxia.tsv",quote = F,sep = "\t",row.names = F)

AT.GO.binox <- binox.creator(pathways = AT.BP[1:20,],DEG = DEG.AT)
EC.GO.binox <- binox.creator(pathways = EC.BP[1:20,],DEG = DEG.EC)
MFB.GO.binox <- binox.creator(pathways = MFB.BP[1:20,],DEG = DEG.MFB)
IVC.GO.binox <- binox.creator(pathways = IVC.BP[1:20,],DEG = DEG.IVC)
IVM.GO.binox <- binox.creator(pathways = IVM.BP,DEG = DEG.IVM)

#write.table(AT.GO.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO/AT_GO_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(EC.GO.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO/EC_GO_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(MFB.GO.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO/MFB_GO_Stretch.tsv",quote = F,sep = "\t",row.names = F)
#write.table(IVC.GO.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO/IVH_GO_Hyperoxia.tsv",quote = F,sep = "\t",row.names = F)
#write.table(IVM.GO.binox,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO/IVM_GO_Hyperoxia.tsv",quote = F,sep = "\t",row.names = F)

# GO and KEGG

AT.KEGG.GO <- GOandKEGG(GO = AT.GO.binox,KEGG = AT.binox)
EC.KEGG.GO <- GOandKEGG(GO = EC.GO.binox,KEGG = EC.binox)
MFB.KEGG.GO <- GOandKEGG(GO = MFB.GO.binox,KEGG = MFB.binox)
IVC.KEGG.GO <- GOandKEGG(GO = IVC.GO.binox,KEGG = IVC.binox)
IVM.KEGG.GO <- GOandKEGG(GO = IVM.GO.binox,KEGG = IVM.binox)

#write.table(AT.KEGG.GO,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO2KEGG/AT_K2G_Stretch.tsv",quote = F,sep = "\t",row.names = T)
#write.table(EC.KEGG.GO,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO2KEGG/EC_K2G_Stretch.tsv",quote = F,sep = "\t",row.names = T)
#write.table(MFB.KEGG.GO,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO2KEGG/MFB_K2G_Stretch.tsv",quote = F,sep = "\t",row.names = T)
#write.table(IVC.KEGG.GO,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO2KEGG/IVH_K2G_Hyperoxia.tsv",quote = F,sep = "\t",row.names = T)
#write.table(IVM.KEGG.GO,"~/Documents/PhD/My_Results/microBulk/binox_files/newGO2KEGG/IVM_K2G_Hyperoxia.tsv",quote = F,sep = "\t",row.names = T)

graph.KEGG.GO(AT.KEGG.GO)
graph.KEGG.GO(EC.KEGG.GO)
graph.KEGG.GO(MFB.KEGG.GO)
graph.KEGG.GO(IVC.KEGG.GO)
graph.KEGG.GO(IVM.KEGG.GO)

# Binox results - graphs

AT.KEGG.BINOX <- read.table("~/Documents/PhD/My_Results/microBulk/binox_files/newKEGG/AT_vs_IVH_Hyperoxia.tsv",
                            header = T, row.names = 1,comment.char = "")

AT.KEGG.BINOX[which(AT.KEGG.BINOX$X.5.FDR < 0.05 & AT.KEGG.BINOX$X.6.relationType == "+"),]











