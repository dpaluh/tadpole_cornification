setwd("")
Path_toCount=""

library(GO.db)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(Mfuzz)
library(data.table)
library(ggrepel)
library(pheatmap)
library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)

GTF_file <- "GCF_027579735.1_aBomBom1.pri_genomic.gtf"
genome_annotation.df <- read.delim(GTF_file, header=FALSE, comment.char="#", stringsAsFactors=FALSE)
gtf_headers <- c("sequence", "source", "feature", "start", "end", "score","strand", "phase", "attributes")
colnames(genome_annotation.df) <- gtf_headers
myfeatures<-c("gene")
genes_only<-subset(genome_annotation.df, genome_annotation.df$feature %in% myfeatures )
genenames = unlist(lapply(genes_only$attributes,function(x) strsplit(x,';')[[1]][1]))
genenames<-gsub("gene_id ", "", genenames)
genes_only$gene<-genenames
genes_only$gene_pos<-paste(genes_only$gene, genes_only$sequence, sep="_")


# Run ABC model -----------------------------------------------------------
count_name = list.files(path = Path_toCount, full.names=F, pattern="ReadsPerGene.out.tab")
sample_name = gsub("ReadsPerGene.out.tab","",count_name)
stage = unlist(lapply(sample_name,function(x) strsplit(x,'_')[[1]][3]))
batch = unlist(lapply(sample_name,function(x) strsplit(x,'_')[[1]][4]))
sampleTable = data.frame(
  sampleName = sample_name,
  fileName = count_name,
  condition_stage = stage,
  condition_batch = batch)

files<-list.files(Path_toCount,full.names=T)
files<-files[1:9]
countData <- data.frame(data.table::fread(files[1]))[, c(1, 4)]
colnames(countData)[colnames(countData)=="V1"]<-"Gene"
colnames(countData)[colnames(countData)=="V4"]<-gsub("ReadsPerGene.out.tab","",gsub(paste0(Path_toCount,"/"),"",files[1]))

for(i in 2:length(files)){
  countData = cbind(countData,data.frame(data.table::fread(files[i]))[4])
  colnames(countData)[colnames(countData)=="V4"]<-gsub("ReadsPerGene.out.tab","",gsub(paste0(Path_toCount,"/"),"",files[i]))
}

#Skip first 4 rows that contain output stats
countData = countData[c(5:nrow(countData)),]
rownames(countData) = countData$Gene
countData$Gene<-NULL
nrow(countData)

sampleTable<-sampleTable[sampleTable$sampleName!="BO_MO_C_2_26",] ## This sample was chosen to be removed from doing all reps PCA first
colnames(countData) <- basename(colnames(countData))
countData <- countData[, colnames(countData) != "BO_MO_C_2_26"]

ddsHTSeq = DESeqDataSetFromMatrix(countData = countData, colData = sampleTable, design = ~  condition_batch + condition_stage )
ddsHTSeq_keep <- ddsHTSeq[ rowSums(counts(ddsHTSeq)>=5) >= 3, ]

ddsHTSeq<-DESeq(ddsHTSeq_keep)
Bombina_ABC_ddsHTSeq_results<-results(ddsHTSeq)
Bombina_ABC_ddsHTSeq<-DESeq(ddsHTSeq)


# Run PCA -----------------------------------------------------------------
ddsHTSeq_trans <- vst(Bombina_ABC_ddsHTSeq , blind=FALSE) #blind=FALSE takes design into account, default is blind=T
data <- plotPCA(ddsHTSeq_trans, intgroup=c("condition_stage"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
mycolors<-c("greenyellow","green4","blue3")
names(mycolors)<-c("A","B","C")
plot<-ggplot(data, aes(PC1, PC2, color=condition_stage)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+
  theme(aspect.ratio=1, axis.text=element_text(size=18), axis.title=element_text(size=18),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        plot.title = element_text(size=20), 
        legend.key.size = unit(1, 'cm'),legend.text = element_text(size = 18),legend.title = element_text(size = 20),
        legend.key=element_blank(),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = 1))+
  labs(colour = "Stage")+ggtitle("Bombina bombina PCA")+scale_color_manual(values=mycolors)
ggsave(plot,file="Bombina_PCA_C26rmd.pdf",device=pdf,width = 210, height=297, units="mm")
ggsave(plot,file="Bombina_PCA_C26rmd.png",device=png,width = 210, height=297, units="mm")


# Correct modeling for multiple comparisons -------------------------------

Bombina_AB_ddsHTSeq_results<-results(Bombina_ABC_ddsHTSeq,
                 contrast = c("condition_stage","B","A"))

Bombina_BC_ddsHTSeq_results<-results(Bombina_ABC_ddsHTSeq,
                 contrast = c("condition_stage","C","B"))

Bombina_AC_ddsHTSeq_results<-results(Bombina_ABC_ddsHTSeq,
                                     contrast = c("condition_stage","C","A"))


# Filtering DEGs ----------------------------------------------------------

##Here looking getting deferentially expressed genes, making two lists, one for B/A and one for C/B
Bombina_AB<-as.data.frame(Bombina_AB_ddsHTSeq_results)
Bombina_BC<-as.data.frame(Bombina_BC_ddsHTSeq_results)
Bombina_AC<-as.data.frame(Bombina_AC_ddsHTSeq_results)

# Volcano Plots -----------------------------------------------------------
mycolors<-c("greenyellow","green4","blue3","grey50")
names(mycolors)<-c("A","B","C","unbiased")
plot_df<-Bombina_AB
plot_df$bias<-"unbiased"
plot_df$bias[plot_df$padj<0.05&is.na(plot_df$padj)!=T&plot_df$log2FoldChange<(-1)]<-"A"
plot_df$bias[plot_df$padj<0.05&is.na(plot_df$padj)!=T&plot_df$log2FoldChange>1]<-"B"
table(plot_df$bias)
volcano_plot<-ggplot()+geom_point(data=plot_df,aes(x=log2FoldChange, y=-log10(padj),color=bias))+
  theme_light() +
  theme(plot.title=element_text(size=10),legend.position = "none",
        axis.text.x=element_text(size=23),
        axis.title.x=element_text(size=23),
        axis.text.y=element_text(size=23),
        axis.title.y=element_text(size=23),
        aspect.ratio = 9/8, axis.title=element_text(size=14),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  geom_hline(yintercept=-log10(0.05), col="darkred",linetype = 2) +  
  geom_vline(xintercept=-1, col="darkred",linetype = 2) +
  geom_vline(xintercept=1, col="darkred",linetype = 2) +
  ggtitle("B vs A")+
  ylab("-log10(p-adj)")+xlab("Log2 Fold-change")+scale_color_manual(values=mycolors)
ggsave(volcano_plot,file="Bombina_AB_Volcano.pdf",device=pdf, width = 105, height=148, units="mm")
ggsave(volcano_plot,file="Bombina_AB_Volcano.png",device=png, width = 105, height=148, units="mm")

plot_df<-Bombina_BC
plot_df$bias<-"unbiased"
plot_df$bias[plot_df$padj<0.05&is.na(plot_df$padj)!=T&plot_df$log2FoldChange<(-1)]<-"B"
plot_df$bias[plot_df$padj<0.05&is.na(plot_df$padj)!=T&plot_df$log2FoldChange>1]<-"C"
table(plot_df$bias)


# significantly DE genes -----------------------------------------------------------
AB_DEGs<-rownames(Bombina_AB[Bombina_AB$padj<0.05 & is.na(Bombina_AB$padj)!=T & abs(Bombina_AB$log2FoldChange)>1,])
BC_DEGs<-rownames(Bombina_BC[Bombina_BC$padj<0.05 & is.na(Bombina_BC$padj)!=T & abs(Bombina_BC$log2FoldChange)>1,])
AC_DEGs<-rownames(Bombina_AC[Bombina_AC$padj<0.05 & is.na(Bombina_AC$padj)!=T & abs(Bombina_AC$log2FoldChange)>1,])
AB_DEGs_df<-Bombina_AB[Bombina_AB$padj<0.05 & is.na(Bombina_AB$padj)!=T & abs(Bombina_AB$log2FoldChange)>1,]
BC_DEGs_df<-Bombina_BC[Bombina_BC$padj<0.05 & is.na(Bombina_BC$padj)!=T & abs(Bombina_BC$log2FoldChange)>1,]
AC_DEGs_df<-Bombina_AC[Bombina_AC$padj<0.05 & is.na(Bombina_AC$padj)!=T & abs(Bombina_AC$log2FoldChange)>1,]

write.csv(AB_DEGs_df,file = "AB_DEGs.csv")
write.csv(BC_DEGs_df,file = "BC_DEGs.csv")


##The top biased genes in the AB model are all B biased

# mFuzz -------------------------------------------------------------------

AB_DEGs<-rownames(Bombina_AB[Bombina_AB$padj<0.05 & is.na(Bombina_AB$padj)!=T & abs(Bombina_AB$log2FoldChange)>1,])
BC_DEGs<-rownames(Bombina_BC[Bombina_BC$padj<0.05 & is.na(Bombina_BC$padj)!=T & abs(Bombina_BC$log2FoldChange)>1,])
AC_DEGs<-rownames(Bombina_AC[Bombina_AC$padj<0.05 & is.na(Bombina_AC$padj)!=T & abs(Bombina_AC$log2FoldChange)>1,])
AB_DEGs_df<-Bombina_AB[Bombina_AB$padj<0.05 & is.na(Bombina_AB$padj)!=T & abs(Bombina_AB$log2FoldChange)>1,]
BC_DEGs_df<-Bombina_BC[Bombina_BC$padj<0.05 & is.na(Bombina_BC$padj)!=T & abs(Bombina_BC$log2FoldChange)>1,]
AC_DEGs_df<-Bombina_AC[Bombina_AC$padj<0.05 & is.na(Bombina_AC$padj)!=T & abs(Bombina_AC$log2FoldChange)>1,]

##DEGs are obtained from pairwise comparisons of the stages, 
DEGs<-unique(c(AB_DEGs,BC_DEGs,AC_DEGs))


ABC_DEGs_normcounts_df<-counts(Bombina_ABC_ddsHTSeq,normalized=T)
ABC_DEGs_normcounts_df_DEGs<-ABC_DEGs_normcounts_df[rownames(ABC_DEGs_normcounts_df)%in%DEGs,]
nrow(ABC_DEGs_normcounts_df_DEGs)

BO_A<-as.data.frame(ABC_DEGs_normcounts_df_DEGs[,colnames(ABC_DEGs_normcounts_df)%like%"_A_"])
BO_B<-as.data.frame(ABC_DEGs_normcounts_df_DEGs[,colnames(ABC_DEGs_normcounts_df)%like%"_B_"])
BO_C<-as.data.frame(ABC_DEGs_normcounts_df_DEGs[,colnames(ABC_DEGs_normcounts_df)%like%"_C_"])

BO_A$A_Median<-apply(BO_A,1,function(x){median(x[x>0])})
BO_A$A_Median[is.na(BO_A$A_Median)]<-0
BO_B$B_Median<-apply(BO_B,1,function(x){median(x[x>0])})
BO_B$B_Median[is.na(BO_B$B_Median)]<-0
BO_C$C_Median<-apply(BO_C,1,function(x){median(x[x>0])})
BO_C$C_Median[is.na(BO_C$C_Median)]<-0

BO_Medians<-cbind(BO_A,BO_B,BO_C)
BO_Medians<-BO_Medians[,colnames(BO_Medians)%like%"_Median"]
ex.m <- as.matrix(BO_Medians) 
eset <- new('ExpressionSet', exprs=ex.m) 
AB.r <- filter.NA(eset, thres=0.25) 
AB.f <- fill.NA(AB.r,mode="mean")
tmp <- filter.std(AB.f,min.std=0)
AB.s <- standardise(AB.f)
m1 <- mestimate(AB.s)  ##m=3.6803
cl <- mfuzz(AB.s,c=6,m=m1)
cl$cluster
write.csv(cl$cluster,file="cl$cluster.csv")
citation("clusterProfiler")

for (i in 1:6){
  pdf(paste0("~/Bombina/Cluster6_",i,"_mFuzz.pdf"))
  mfuzz.plot2(AB.s,cl=cl,x11=F,Xwidth=12,Xheight=10,xlab="Stage",single=i,time.labels=c("Before","During","After"))
  mtext(side=3, text=paste0(length(names(cl$cluster[cl$cluster==i]))," genes"),cex=1.3)
  dev.off()
}



# GO Term Analysis --------------------------------------------------------

DEGs<-unique(c(AB_DEGs,BC_DEGs,AC_DEGs))
DEGs_for_GO<-DEGs
universe<-rownames(Bombina_ABC_ddsHTSeq_results)

Go_annotations_df<-read.csv("GCF_027579735.1-RS_2023_03_gene_ontology.gaf",sep="\t",skip=8)
Go_annotations_df$GO_ID
Go_annotations_df$Term<-Term(Go_annotations_df$GO_ID)
Go_annotations_df$Ontology<-Ontology(Go_annotations_df$GO_ID)
##First column gene ID second column GO ID
GeneID_Symbol_Translater<-dplyr::select(genes_only,gene,attributes)
GeneID_Symbol_Translater$GeneID<-str_split_i(str_split_i(GeneID_Symbol_Translater$attributes,"db_xref GeneID:",2),";",1)
GeneID_Symbol_Translater<-dplyr::select(GeneID_Symbol_Translater,gene,GeneID)

##Now translate your gene symbols to gene IDs
DEGs_GeneIDs_GO<-GeneID_Symbol_Translater[GeneID_Symbol_Translater$gene%in%DEGs_for_GO,]
DEGs_GeneIDs_GO<-DEGs_GeneIDs_GO$GeneID
Universe_GeneIDs_GO<-GeneID_Symbol_Translater[GeneID_Symbol_Translater$gene%in%universe,]
Universe_GeneIDs_GO<-Universe_GeneIDs_GO$GeneID

cl_i<-names(cl$cluster[cl$cluster==2])
for (i in 1:6) {
  
  test_cat<-"BP"
  message("Begin cluster ", i)
  cl_i<-names(cl$cluster[cl$cluster==i])
  cl_i<-gsub("gene-","",cl_i) 
  DEGs_for_GO<-gsub("gene-","",cl_i)
  DEGs_GeneIDs_GO<-GeneID_Symbol_Translater[GeneID_Symbol_Translater$gene%in%DEGs_for_GO,]$GeneID
  message("Running Test for Cluster ", i)
  try(results1<-enricher(
    gene=DEGs_GeneIDs_GO,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = Universe_GeneIDs_GO,
    minGSSize = 10,
    maxGSSize = length(DEGs_GeneIDs_GO),
    gson = NULL,
    TERM2GENE = Go_annotations_df[Go_annotations_df$Ontology==test_cat,c("GO_ID","GeneID")],
    TERM2NAME = Go_annotations_df[Go_annotations_df$Ontology==test_cat,c("GO_ID","Term")]
  ))
  
  assign(paste0("Cluster_",i,"_",test_cat,"_results"),results1)
  Nr_significant<-nrow(as.data.frame(results1))
  message("Plotting Cluster ", i)
  
  try(GO_dotplot<-dotplot(results1,showCategory=Nr_significant,font.size = 10)+
        ggtitle(paste0("Cluster ",i,"  (",str_split_i(as.data.frame(results1)$GeneRatio[1],"/",2),"/",length(cl_i)," genes analyzed)"))+
        theme(aspect.ratio = 1,
              axis.text.x=element_text(size=15),
              axis.text.y=element_text(size=13),
              axis.title.x = element_text(size=15),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12))
              )
  
  message("Saving Cluster ", i)
  
  try(ggsave(plot=GO_dotplot,file=paste0("~/Bombina",test_cat,"_6_Clusters_",i,".pdf"),device=pdf,width=210,height=297,units = "mm"))
  try(ggsave(plot=GO_dotplot,file=paste0("~/Bombina",test_cat,"_6_Clusters_",i,".png"),device=png,width=210,height=297,units = "mm"))
  rm(results1)
  rm(GO_dotplot)
}


for (i in c(1,2,4,5,6)){
  GOTerm_results<-as.data.frame(get(paste0("Cluster_",i,"_",test_cat,"_results")))
  GOTerm_results$Ontology_category<-NA
  GOTerm_results$Ontology_category<-Ontology(as.data.frame(GOTerm_results)$ID)
  GOTerm_results<-dplyr::select(GOTerm_results,Description,GeneRatio,Count,Ontology_category,geneID,p.adjust)
  GOTerm_results$Gene_Symbols<-NA
  for (k in 1:nrow(GOTerm_results)){
    GOTerm_results$Gene_Symbols[k]<-paste(GeneID_Symbol_Translater$gene[GeneID_Symbol_Translater$GeneID%in%unlist(strsplit(GOTerm_results$geneID[k],"/"))],collapse = ",")
  }
  BP<-GOTerm_results[GOTerm_results$Ontology_category=="BP",]
  BP<-BP[order(BP$p.adjust,decreasing = F),]
  CC<-GOTerm_results[GOTerm_results$Ontology_category=="CC",]
  CC<-CC[order(CC$p.adjust,decreasing = F),]
  MF<-GOTerm_results[GOTerm_results$Ontology_category=="MF",]
  MF<-MF[order(MF$p.adjust,decreasing = F),]
  GOTerm_results_ordered<-rbind(BP,CC,MF)
  GOTerm_results_ordered$Cluster<-i
  GOTerm_results_ordered$GO_ID<-rownames(GOTerm_results_ordered)
  GOTerm_results_ordered<-dplyr::select(GOTerm_results_ordered,Ontology_category,Cluster,GO_ID,Description,GeneRatio,p.adjust,Gene_Symbols)
  write.csv(GOTerm_results_ordered,file=paste0("~/Bombina/GOTerm_results_cluster_",test_cat,"_",i,".csv"))
}

for (i in c(1,2,4,5,6)){
  obj <- paste0("Cluster_",i,"_",test_cat,"_results")
  x <- as.data.frame(get(obj))
  cat(obj, "rows =", nrow(x), "\n")
}
# End ---------------------------------------------------------------------



##########
##########
#Line Plots of keratin and AEDC genes
##########
##########
setwd("")

library(dplyr)
library(DESeq2)
library(tidyr)
library(ggrepel)



# All STAR count files
files <- list.files(pattern = "ReadsPerGene.out.tab$")

# Read the first file
counts <- read.delim(files[1], header = FALSE)
head(counts)
# Keep only gene rows

counts <- counts[-(1:4), c(1,4)]
colnames(counts) <- c("gene", tools::file_path_sans_ext(basename(files[1])))
# Add remaining samples
for(i in 2:length(files)) {
  
  tmp <- read.delim(files[i], header = FALSE)
  tmp <- tmp[-(1:4), c(1,4)]
  
  # Make sure genes are in the same order
  if (!identical(tmp$V1, counts$gene)) {
    stop(paste("Gene order differs in", files[i]))
  }
  
  counts[[tools::file_path_sans_ext(basename(files[i]))]] <- tmp$V4
}

counts

counts <- counts %>%
  dplyr::select(-BO_MO_C_2_26ReadsPerGene.out)

rownames(counts) <- counts$gene
counts <- counts[,-1]

coldata <- data.frame(
  stage = factor(c("A","A","A","B","B","B","C","C"),
                 levels = c("A","B","C"))
)
rownames(coldata) <- colnames(counts)

str(counts)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ stage
)

dds <- estimateSizeFactors(dds)

#matrix of variance-stabilized expression values
vsd <- vst(dds, blind = TRUE)

vsd_matrix <- assay(vsd)
rownames(vsd_matrix) <- sub("^gene-", "", rownames(vsd_matrix))
colnames(vsd_matrix) <- sub("ReadsPerGene\\.out$", "", colnames(vsd_matrix))

krt_typeI <- c(
  "LOC128637710",
  "LOC128641512",
  "LOC128641528",
  "LOC128641544",
  "LOC128641560",
  "LOC128641577",
  "LOC128641588",
  "LOC128641597",
  "LOC128643914",
  "KRT222",
  "LOC128645469",
  "LOC128645489",
  "LOC128645496",
  "LOC128645504",
  "LOC128645512",
  "LOC128645523",
  "LOC128645530",
  "LOC128645537",
  "LOC128645544",
  "LOC128645554",
  "LOC128645561",
  "LOC128645566",
  "LOC128645577",
  "LOC128645584",
  "LOC128645590",
  "LOC128645597",
  "LOC128645604",
  "LOC128645610",
  "LOC128645617",
  "LOC128645623",
  "LOC128645636",
  "LOC128645644",
  "KRT18",
  "LOC128663343",
  "LOC128663371",
  "LOC128663387",
  "LOC128663395",
  "LOC128663404",
  "LOC128663413",
  "LOC128663422",
  "LOC128663439",
  "LOC128663448",
  "LOC128663456"
)

krt_typeII <- c(
  "LOC128645268",
  "LOC128652065",
  "LOC128652066",
  "LOC128652614",
  "LOC128652616",
  "LOC128652617",
  "LOC128652618",
  "LOC128652619",
  "KRT8",
  "LOC128654423",
  "LOC128654424",
  "LOC128654425",
  "LOC128654427",
  "LOC128654428",
  "LOC128654429",
  "LOC128654430",
  "LOC128654431",
  "LOC128654432",
  "LOC128654433",
  "LOC128654434",
  "LOC128654435",
  "LOC128654436",
  "LOC128654438",
  "LOC128654439"
)

aedc <- c(
  "LOC128639475",
  "LOC128639476",
  "LOC128639561",
  "LOC128639477",
  "LOC128639632",
  "LOC128639479",
  "LOC128639480",
  "LOC128639481",
  "LOC128640370",
  "LOC128639482",
  "LOC128640372",
  "LOC128640373",
  "LOC128639483",
  "LOC128639485",
  "LOC128639486",
  "LOC128639487",
  "LOC128639488",
  "LOC128639489",
  "LOC128639490",
  "LOC128639491",
  "LOC128639492",
  "LOC128639493",
  "LOC128639494",
  "LOC128639845",
  "LOC128639496"
)

krt_typeI_matrix <- vsd_matrix[rownames(vsd_matrix) %in% krt_typeI, ]
krt_typeII_matrix <- vsd_matrix[rownames(vsd_matrix) %in% krt_typeII, ]
aedc_matrix <- vsd_matrix[rownames(vsd_matrix) %in% aedc, ]


plot_matrix <- rbind(
  krt_typeI_matrix,
  krt_typeII_matrix,
  aedc_matrix
)

gene_class <- data.frame(
  gene = rownames(plot_matrix),
  class = c(
    rep("Type I keratin", nrow(krt_typeI_matrix)),
    rep("Type II keratin", nrow(krt_typeII_matrix)),
    rep("AEDC", nrow(aedc_matrix))
  )
)

plot_data <- as.data.frame(plot_matrix) %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  ) %>%
  left_join(gene_class, by = "gene")

aes(color = class)

plot_data <- plot_data %>%
  mutate(
    stage = case_when(
      grepl("_A_", sample) ~ "A",
      grepl("_B_", sample) ~ "B",
      grepl("_C_", sample) ~ "C"
    )
  )

plot_data$stage <- factor(plot_data$stage,
                          levels = c("A", "B", "C"))

# Krt34 paralogs to highlight
krt34_genes <- c(
  "LOC128641588",
  "LOC128643914",
  "LOC128645584",
  "LOC128645590",
  "LOC128663371",
  "LOC128663387",
  "LOC128663395",
  "LOC128663404",
  "LOC128663413",
  "LOC128663422",
  "LOC128663439",
  "LOC128663448",
  "LOC128663456"
)

# Subset type I keratins
typeI_plot <- plot_data %>%
  filter(class == "Type I keratin") %>%
  mutate(
    highlight = ifelse(gene %in% krt34_genes,
                       gene,
                       "Other Type I")
  )

ggplot(typeI_plot,
       aes(x = stage,
           y = expression,
           group = gene,
           color = highlight)) +
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2,
    alpha = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "Other Type I" = "grey80",
      setNames(
        rainbow(length(krt34_genes)),
        krt34_genes
      )
    )
  ) +
  theme_classic() +
  labs(
    x = "Developmental stage",
    y = "Variance-stabilized expression",
    color = NULL
  )


krt59_genes <- c(
  "LOC128652065",
  "LOC128652066",
  "LOC128652616",
  "LOC128652617",
  "LOC128652618",
  "LOC128652619",
  "LOC128654431",
  "LOC128654432",
  "LOC128654433",
  "LOC128654434",
  "LOC128654435",
  "LOC128654436"
)

# Subset type II keratins
typeII_plot <- plot_data %>%
  filter(class == "Type II keratin") %>%
  mutate(
    highlight = ifelse(gene %in% krt59_genes,
                       gene,
                       "Other Type II")
  )

ggplot(typeII_plot,
       aes(x = stage,
           y = expression,
           group = gene,
           color = highlight)) +
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2,
    alpha = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  scale_color_manual(
    values = c(
      "Other Type II" = "grey80",
      setNames(
        rainbow(length(krt59_genes)),
        krt59_genes
      )
    )
  ) +
  theme_classic() +
  labs(
    x = "Developmental stage",
    y = "Variance-stabilized expression",
    color = NULL
  )




aedc_names <- c(
  "LOC128639475" = "aedc01",
  "LOC128639476" = "aedc02",
  "LOC128639561" = "aedc03",
  "LOC128639477" = "aedc04",
  "LOC128639632" = "aedc05",
  "LOC128639479" = "aedc06",
  "LOC128639480" = "aedc07",
  "LOC128639481" = "aedc08",
  "LOC128640370" = "aedc09",
  "LOC128639482" = "aedc10",
  "LOC128640372" = "aedc11",
  "LOC128640373" = "aedc12",
  "LOC128639483" = "aedc13",
  "LOC128639485" = "aedc14",
  "LOC128639486" = "aedc15",
  "LOC128639487" = "aedc16",
  "LOC128639488" = "aedc17",
  "LOC128639489" = "aedc18",
  "LOC128639490" = "aedc19",
  "LOC128639491" = "aedc20",
  "LOC128639492" = "aedc21",
  "LOC128639493" = "aedc22",
  "LOC128639494" = "aedc23",
  "LOC128639845" = "aedc24",
  "LOC128639496" = "aedc25"
)

aedc_plot <- plot_data %>%
  filter(class == "AEDC") %>%
  mutate(
    gene = recode(gene, !!!aedc_names)
  )

# Plot AEDC expression trajectories
ggplot(aedc_plot,
       aes(x = stage,
           y = expression,
           group = gene,
           color = gene)) +
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2,
    alpha = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  theme_classic() +
  labs(
    x = "Developmental stage",
    y = "Variance-stabilized expression",
    color = "AEDC gene"
  )


label_data <- aedc_plot %>%
  group_by(gene, stage) %>%
  summarize(expression = mean(expression), .groups = "drop") %>%
  group_by(gene) %>%
  slice_tail(n = 1) %>%
  ungroup()

ggplot(aedc_plot,
       aes(x = stage,
           y = expression,
           group = gene,
           color = gene)) +
  geom_point(
    position = position_jitter(width = 0.08),
    size = 2,
    alpha = 0.6
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    direction = "y",
    hjust = 0,
    nudge_x = 0.2,
    segment.color = "grey70",
    show.legend = FALSE
  ) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 80, 5.5, 5.5)
  ) +
  labs(
    x = "Developmental stage",
    y = "Variance-stabilized expression"
  )