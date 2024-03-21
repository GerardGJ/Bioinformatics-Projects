library(dplyr)
library(tidyr)
library(naniar)
library(PhosR)
library(ggplot2)
library(ggpubr)
library(ggfortify)
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
#BiocManager::install("fgsea")
#BiocManager::install("piano")
library(enrichplot)
library(fgsea)
library(sgof)
library(stats)
library(piano)
setwd("/Users/Gerard/Desktop/Perseus_R_Analysis/")


table_in = read.table("proteinGroups_obob_liver_exHEPA.txt", sep = "\t", header = T)
matrix1 <- table_in %>% select(LFQ.intensity.lean_1,LFQ.intensity.lean_2,LFQ.intensity.lean_3,LFQ.intensity.lean_4,
                               LFQ.intensity.ob_1, LFQ.intensity.ob_2, LFQ.intensity.ob_3, LFQ.intensity.ob_4,
                               Peptides, Razor...unique.peptides, Unique.peptides, Sequence.coverage....,
                               Unique...razor.sequence.coverage...., Unique.sequence.coverage...., Mol..weight..kDa.,
                               Q.value, Score, Only.identified.by.site, Reverse, Potential.contaminant, Protein.IDs,
                               Majority.protein.IDs, Protein.names, Gene.names, id)
#Remove the identified by site, reverse ,and Potential contaminant

  #Using which
matrix2 <- matrix1[-which(matrix1$Only.identified.by.site == "+"),]
matrix3 <- matrix2[-which(matrix2$Reverse == "+"),]
matrix4 <- matrix3[-which(matrix3$Potential.contaminant == "+"),]
  
  #Using filter
matrix5 <- matrix1 %>% filter(Only.identified.by.site !="+",Reverse !="+",Potential.contaminant !="+")

#Log2 transform the data on the LFQ cols
matrix5[,1:8] <- log2(matrix5[1:8])
matrix6 <- matrix5 %>% replace_with_na(replace = list(LFQ.intensity.lean_1 = -Inf,LFQ.intensity.lean_2 = -Inf,LFQ.intensity.lean_3 = -Inf,LFQ.intensity.lean_4 = -Inf,
                                                      LFQ.intensity.ob_1 = -Inf, LFQ.intensity.ob_2 = -Inf, LFQ.intensity.ob_3 = -Inf, LFQ.intensity.ob_4 = -Inf))

#Filtering by valid values
  #matrix7 <- matrix6[rowSums(is.na(matrix6[ , 1:8])) <= 4,] Not very useful, doesn't take into account the groups

grps <- as.factor(c("Lean","Lean","Lean","Lean","Ob","Ob","Ob","Ob"))#Annotate the groups

matrix_to_keep <- selectGrps(matrix6[1:8],grps,0.75,n=1)

matrix7 <- matrix6 %>% right_join(matrix_to_keep, by = c("LFQ.intensity.lean_1","LFQ.intensity.lean_2","LFQ.intensity.lean_3","LFQ.intensity.lean_4",
                                                         "LFQ.intensity.ob_1", "LFQ.intensity.ob_2", "LFQ.intensity.ob_3", "LFQ.intensity.ob_4"))
#Imputation
set.seed(123)

matrix_imputed <- tImpute(matrix7[1:8], 1.8, 0.3)
matrix_imputed <- cbind(matrix_imputed,matrix7[,-(1:8)])

#multiple Histograms
ggarrange(ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.lean_1)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.lean_2)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.lean_3)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.lean_4)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.ob_1)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.ob_2)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.ob_3)),
          ggplot(matrix_imputed[,1:8]) + geom_histogram(aes(LFQ.intensity.ob_4)),
          ncol = 2, nrow =4)

#Looking at the correlation
pairs(matrix_imputed[,1:8])

#Perform and plot the PCA
pca <- prcomp(matrix_imputed[,1:8], scale. = T)

autoplot(pca)

ggplot(mapping = aes(pca[["x"]][,1], pca[["x"]][,2])) + 
  geom_point()

ggplot(mapping = aes(pca[["rotation"]][,1], pca[["rotation"]][,2], color = grps)) + 
  geom_point()

#T-test

p_vals <- data.frame()
foldchange <- data.frame()
matrix_to_Ttest_GOBP <- unique(matrix_imputed[,-26])
for(i in 1:nrow(matrix_imputed)){
  
  row = matrix_imputed[i,1:8]
  
  t_test <- t.test(row[5:8], row[1:4])
  
  p_vals <- rbind(p_vals,t_test[["p.value"]])
  foldchange <- rbind(foldchange,t_test[["estimate"]][["mean of x"]] - t_test[["estimate"]][["mean of y"]])
}

p_vals_adjusted <- p.adjust(as.matrix(p_vals),method = "fdr")
t_test_sig <- data.frame()

for(p_val in p_vals_adjusted){
  if(p_val <= 0.05){
    t_test_sig <- rbind(t_test_sig, 0) # 0 significant
  } else {
    t_test_sig <- rbind(t_test_sig,1) # 1 not significant
  }
}
pval_sig <- cbind(as.data.frame(p_vals_adjusted),t_test_sig,foldchange)
colnames(pval_sig) <- c("p_vals_adjusted", "t_test_sig","foldchange")
matrix_with_pvals <- cbind(matrix_imputed,pval_sig)

length(which(matrix_with_pvals$p_vals_adjusted <= 0.05))

#Annotaion #2
annotation_table <- read.table(gzfile("mainAnnot.mus_musculusR.txt.gz"), header = T, sep = "\t", fill = T)
annotations <- c("GOBP.name", "GOMF.name", "GOCC.name", "KEGG.name", "Keywords")
matrix_annotations_GOBP = data.frame()
matrix_annotations_GOMF = data.frame()
matrix_annotations_GOCC = data.frame()
matrix_annotations_KEGG = data.frame()


for(Protein.IDs in matrix_with_pvals$Protein.IDs){
  id = strsplit(x = Protein.IDs, split = ";")[[1]][1]
  row_uni = which(annotation_table$UniProt == id)
    if(length(row_uni) != 0){
      row_to_annotate = annotation_table[row_uni,]
      
      for(listPathways in strsplit(x = row_to_annotate$GOBP.name, split = ";")){
        for(pathway in listPathways){
          matrix_annotations_GOBP = rbind(matrix_annotations_GOBP,c(pathway,id,Protein.IDs))
        }
      }
      for(listPathways in strsplit(x = row_to_annotate$GOMF.name, split = ";")){
        for(pathway in listPathways){
          matrix_annotations_GOMF = rbind(matrix_annotations_GOMF,c(pathway,id,Protein.IDs))
        }
      }
      for(listPathways in strsplit(x = row_to_annotate$GOCC.name, split = ";")){
        for(pathway in listPathways){
          matrix_annotations_GOCC = rbind(matrix_annotations_GOCC,c(pathway,id,Protein.IDs))
        }
      }
      for(listPathways in strsplit(x = row_to_annotate$KEGG.name, split = ";")){
        for(pathway in listPathways){
          matrix_annotations_KEGG = rbind(matrix_annotations_KEGG,c(pathway,id,Protein.IDs))
        }
      }
      
    }
}
matrix_annotated_GOBP <- matrix_with_pvals %>% inner_join(matrix_annotations_GOBP, by = "Protein.IDs")
matrix_annotated_GOMF <- matrix_with_pvals %>% inner_join(matrix_annotations_GOMF, by = "Protein.IDs")
matrix_annotated_GOCC <- matrix_with_pvals %>% inner_join(matrix_annotations_GOCC, by = "Protein.IDs")
matrix_annotated_KEGG <- matrix_with_pvals %>% inner_join(matrix_annotations_KEGG, by = "Protein.IDs")

matrix_annotated_KEGG <- unique(matrix_annotated_KEGG)


N_GOBP <- nrow(unique(matrix_annotated_GOBP[,-29]))
n_GOBP <- length(which(unique(matrix_annotated_GOBP[,-29])$t_test_sig == 0))
N_GOMF <- nrow(unique(matrix_annotated_GOMF[,-29]))
n_GOMF <- length(which(unique(matrix_annotated_GOMF[,-29])$t_test_sig == 0))
N_GOCC <- nrow(unique(matrix_annotated_GOCC[,-29]))
n_GOCC <- length(which(unique(matrix_annotated_GOCC[,-29])$t_test_sig == 0))
N_KEGG <- nrow(unique(matrix_annotated_KEGG[,-29]))
n_KEGG <- length(which(unique(matrix_annotated_KEGG[,-29])$t_test_sig == 0))

#Volcano Plot
ggplot(matrix_with_pvals, aes(foldchange, -log10(p_vals_adjusted))) + 
  geom_point(aes(color = as.factor(t_test_sig)))

#Fisher exact test

  #GOBP
matrix_of_annotations_GOBP <- matrix_annotated_GOBP %>% select("pathway","t_test_sig")
table_GOBP <- as.data.frame(table(matrix_of_annotations_GOBP))
table_GOBP <- reshape(table_GOBP, idvar = "pathway", timevar = "t_test_sig", direction = "wide")

pVals_Fisher <- data.frame()
for(i in 1:nrow(table_GOBP)){
  row = table_GOBP[i,]
  
  K = sum(row[2:3]) #Total number of genes in that pathway
  k = row[2][1,1] #Total number of significant genes in the pathway
  
  c11 = k
  c21 = K-k
  c12 = n_GOBP-k
  c22 = N_GOBP-K-n_GOBP+k
  
  m <- matrix(c(c11,c21,c12,c22),2, 2)
  fisher_res = fisher.test(m,alternative ="greater")
  pVals_Fisher = rbind(pVals_Fisher, fisher_res[["p.value"]])
}
table_GOBP <- cbind(table_GOBP, pVals_Fisher)

  #GOMF

matrix_of_annotations_GOMF <- matrix_annotated_GOMF %>% select("pathway","t_test_sig")
table_GOMF <- as.data.frame(table(matrix_of_annotations_GOMF))
table_GOMF <- reshape(table_GOMF, idvar = "pathway", timevar = "t_test_sig", direction = "wide")

pVals_Fisher <- data.frame()
for(i in 1:nrow(table_GOMF)){
  row = table_GOMF[i,]
  
  K = sum(row[2:3]) #Total number of genes in that pathway
  k = row[2][1,1] #Total number of significant genes in the pathway
  
  c11 = k
  c21 = K-k
  c12 = n_GOMF-k
  c22 = N_GOMF-K-n_GOMF+k
  
  m <- matrix(c(c11,c21,c12,c22),2, 2)
  fisher_res = fisher.test(m)
  pVals_Fisher = rbind(pVals_Fisher, fisher_res[["p.value"]])
}
table_GOMF <- cbind(table_GOMF, pVals_Fisher)

  #GOCC

matrix_of_annotations_GOCC <- matrix_annotated_GOCC %>% select("pathway","t_test_sig")
table_GOCC <- as.data.frame(table(matrix_of_annotations_GOCC))
table_GOCC <- reshape(table_GOCC, idvar = "pathway", timevar = "t_test_sig", direction = "wide")

pVals_Fisher <- data.frame()
for(i in 1:nrow(table_GOCC)){
  row = table_GOCC[i,]
  
  K = sum(row[2:3]) #Total number of genes in that pathway
  k = row[2][1,1] #Total number of significant genes in the pathway
  
  c11 = k
  c21 = K-k
  c12 = n_GOCC-k
  c22 = N_GOCC-K-n_GOCC+k
  
  m <- matrix(c(c11,c21,c12,c22),2, 2)
  fisher_res = fisher.test(m)
  pVals_Fisher = rbind(pVals_Fisher, fisher_res[["p.value"]])
}
table_GOCC <- cbind(table_GOCC, pVals_Fisher)

  #KEGG

matrix_of_annotations_KEGG <- matrix_annotated_KEGG %>% select("pathway","t_test_sig")
table_KEGG <- as.data.frame(table(matrix_of_annotations_KEGG))
table_KEGG <- reshape(table_KEGG, idvar = "pathway", timevar = "t_test_sig", direction = "wide")

pVals_Fisher <- data.frame()
for(i in 1:nrow(table_KEGG)){
  row = table_KEGG[i,]
  
  K = sum(row[2:3]) #Total number of genes in that pathway
  k = row[2][1,1] #Total number of significant genes in the pathway
  
  c11 = k
  c21 = K-k
  c12 = n_KEGG-k
  c22 = N_KEGG-K-n_KEGG+k
  
  m <- matrix(c(c11,c21,c12,c22),2, 2)
  fisher_res = fisher.test(m)
  pVals_Fisher = rbind(pVals_Fisher, fisher_res[["p.value"]])
}
table_KEGG <- cbind(table_KEGG, pVals_Fisher)

#GSEA with fgsea (gene set enrichment analysis)

  #GOBP

GOBP_pathways <- matrix_annotated_GOBP %>% select(pathway,id.y)
GOBP_pathway_table <- aggregate(GOBP_pathways$id.y, list(GOBP_pathways$pathway), paste, collapse=" ")
for(i in 1:nrow(GOBP_pathway_table)){
  GOBP_pathway_table$x[[i]] = strsplit(x = GOBP_pathway_table$x[[i]], split = " ")
  GOBP_pathway_table[i,2] <- as.vector(GOBP_pathway_table$x[[i]])
}

GOBP_pathway_table <- t(GOBP_pathway_table)
colnames(GOBP_pathway_table) <- GOBP_pathway_table[1,]
GOBP_pathway_table <- GOBP_pathway_table[-1,]
GOBP_pathway_table <- as.list(GOBP_pathway_table)

#The ranks are the fold change in order
matrix_ranks <- matrix_annotated_GOBP %>% select(foldchange)
matrix_ranks <- t(as.vector(matrix_ranks))
names(matrix_ranks) <- matrix_annotated_GOBP$id.y
barplot(sort(matrix_ranks,decreasing = T))

fgseaRes_GOBP <- fgsea(pathways = GOBP_pathway_table, 
                  stats    = matrix_ranks,
                  eps      = 0.0,
                  minSize  = 1,
                  maxSize  = 500)
#As the names of the pathway are to long need to create a table where we connect the number to the pathway
head(fgseaRes_GOBP[order(pval), ])

  #GOMF

GOMF_pathways <- matrix_annotated_GOMF %>% select(pathway,id.y)
GOMF_pathway_table <- aggregate(GOMF_pathways$id.y, list(GOMF_pathways$pathway), paste, collapse=" ")
for(i in 1:nrow(GOMF_pathway_table)){
  GOMF_pathway_table$x[[i]] = strsplit(x = GOMF_pathway_table$x[[i]], split = " ")
  GOMF_pathway_table[i,2] <- as.vector(GOMF_pathway_table$x[[i]])
}

GOMF_pathway_table <- t(GOMF_pathway_table)
colnames(GOMF_pathway_table) <- GOMF_pathway_table[1,]
GOMF_pathway_table <- GOMF_pathway_table[-1,]
GOMF_pathway_table <- as.list(GOMF_pathway_table)

#The ranks are the fold change in order
matrix_ranks <- matrix_annotated_GOMF %>% select(foldchange)
matrix_ranks <- t(as.vector(matrix_ranks))
names(matrix_ranks) <- matrix_annotated_GOMF$id.y
barplot(sort(matrix_ranks,decreasing = T))

fgseaRes_GOMF <- fgsea(pathways = GOMF_pathway_table, 
                       stats    = matrix_ranks,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)
#As the names of the pathway are to long need to create a table where we connect the number to the pathway
head(fgseaRes_GOMF[order(pval), ])

  #GOCC
GOCC_pathways <- matrix_annotated_GOCC %>% select(pathway,id.y)
GOCC_pathway_table <- aggregate(GOCC_pathways$id.y, list(GOCC_pathways$pathway), paste, collapse=" ")
for(i in 1:nrow(GOCC_pathway_table)){
  GOCC_pathway_table$x[[i]] = strsplit(x = GOCC_pathway_table$x[[i]], split = " ")
  GOCC_pathway_table[i,2] <- as.vector(GOCC_pathway_table$x[[i]])
}

GOCC_pathway_table <- t(GOCC_pathway_table)
colnames(GOCC_pathway_table) <- GOCC_pathway_table[1,]
GOCC_pathway_table <- GOCC_pathway_table[-1,]
GOCC_pathway_table <- as.list(GOCC_pathway_table)

#The ranks are the fold change in order
matrix_ranks <- matrix_annotated_GOCC %>% select(foldchange)
matrix_ranks <- t(as.vector(matrix_ranks))
names(matrix_ranks) <- matrix_annotated_GOCC$id.y
barplot(sort(matrix_ranks,decreasing = T))

fgseaRes_GOCC <- fgsea(pathways = GOCC_pathway_table, 
                       stats    = matrix_ranks,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)
#As the names of the pathway are to long need to create a table where we connect the number to the pathway
head(fgseaRes_GOCC[order(pval), ])

  #KEGG

KEGG_pathways <- matrix_annotated_KEGG %>% select(pathway,id.y)
KEGG_pathway_table <- aggregate(KEGG_pathways$id.y, list(KEGG_pathways$pathway), paste, collapse=" ")
for(i in 1:nrow(KEGG_pathway_table)){
  KEGG_pathway_table$x[[i]] = strsplit(x = KEGG_pathway_table$x[[i]], split = " ")
  KEGG_pathway_table[i,2] <- as.vector(KEGG_pathway_table$x[[i]])
}

KEGG_pathway_table <- t(KEGG_pathway_table)
colnames(KEGG_pathway_table) <- KEGG_pathway_table[1,]
KEGG_pathway_table <- KEGG_pathway_table[-1,]
KEGG_pathway_table <- as.list(KEGG_pathway_table)

#The ranks are the fold change in order
matrix_ranks <- matrix_annotated_KEGG %>% select(foldchange)
matrix_ranks <- t(as.vector(matrix_ranks))
names(matrix_ranks) <- matrix_annotated_KEGG$id.y
barplot(sort(matrix_ranks,decreasing = T))

fgseaRes_KEGG <- fgsea(pathways = KEGG_pathway_table, 
                       stats    = matrix_ranks,
                       eps      = 0.0,
                       minSize  = 1,
                       maxSize  = 500)
head(fgseaRes_KEGG[order(pval), ])

#Extra Analysis

#Which genes are in all 3 GO dataframes?

sig_GOBP <- fgseaRes_GOBP[which(fgseaRes_GOBP$pval <= 0.05)]
sig_GOBP <- fgseaRes_GOMF[which(fgseaRes_GOMF$pval <= 0.05)]
sig_GOCC <- fgseaRes_GOCC[which(fgseaRes_GOCC$pval <= 0.05)]
sig_KEGG <- fgseaRes_KEGG[which(fgseaRes_KEGG$pval <= 0.05)]

intersect(intersect(sig_GOBP$leadingEdge,sig_GOMF$leadingEdge),sig_GOCC$leadingEdge)

index_GOBP <- sig_GOBP$pathway[which(sig_GOBP$leadingEdge == intersect(intersect(sig_GOBP$leadingEdge,sig_GOMF$leadingEdge),sig_GOCC$leadingEdge)[[1]])]
GOBP_Sig_pathway <- number_pathway_GOBP[index_GOBP,]

index_GOMF <- sig_GOMF$pathway[which(sig_GOMF$leadingEdge == intersect(intersect(sig_GOBP$leadingEdge,sig_GOMF$leadingEdge),sig_GOCC$leadingEdge)[[1]])]
GOMF_Sig_pathway <- number_pathway_GOMF[index_GOMF,]

index_GOCC <- sig_GOCC$pathway[which(sig_GOCC$leadingEdge == intersect(intersect(sig_GOBP$leadingEdge,sig_GOMF$leadingEdge),sig_GOCC$leadingEdge)[[1]])]
GOCC_Sig_pathway <- number_pathway_GOCC[index_GOCC,]

for(i in nrow(sig_GOBP_MF_CC)){
  sig_GOBP_MF_CC[[i,1]] <- number_pathway_GOBP[sig_GOBP_MF_CC[[i,1]],2]
  sig_GOBP_MF_CC[[i,2]] <- number_pathway_GOMF[sig_GOBP_MF_CC[[i,2]],2]
  sig_GOBP_MF_CC[[i,3]] <- number_pathway_GOCC[sig_GOBP_MF_CC[[i,3]],2]
}

#Which GO are in common in all 3 fisher
sig_GOBP_F <- table_GOBP[which(table_GOBP$X1 <= 0.05),]
sig_GOMF_F <- table_GOMF[which(table_GOMF$x1 <= 0.05),]
sig_GOCC_F <- table_GOCC[which(table_GOCC$x1 <= 0.05),]





a <- as.data.frame(sort(matrix_ranks, decreasing = T))
