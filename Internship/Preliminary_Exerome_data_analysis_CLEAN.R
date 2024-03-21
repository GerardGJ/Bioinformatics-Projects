# Loading libraries
library('limma')
library('statmod')
library('PhosR')
library('dplyr')
library('tidyverse')
library('naniar')
library('ggpubr')
library('PhosR')
library('ggrepel')
library('Rtsne')
library('umap')
library('readxl')
library('bayestestR')
library("fgsea")
setwd("/Users/Gerard/Desktop/Gerard_ExerOME/")
##### Functions #####
Calculate_AUC <- function(Intensity,subjects, grps){
  T24 <- grep("T24",colnames(Intensity))
  Intensity <- Intensity[,-T24]
  grps <- grps[-T24]
  subjects <- subjects[-T24]
  
  T180 <- grep("T180",colnames(Intensity))
  Intensity <- Intensity[,-T180]
  grps <- grps[-T180]
  subjects <- subjects[-T180]
  
  Time_points <- c("pre", "T0","T30","T60")
  x_values <- c(-45,0,30,60) + 45
  num_prot <- nrow(Intensity)
  auc_df <- matrix(nrow = num_prot, ncol = 25)
  index_y = 1
  for(sbj in levels(subjects)){
    
    cols_subject <- Intensity[,which(subjects == sbj)]
    order_cols <- order(match(grps[which(subjects == sbj)],Time_points))
    ordered_cols <- cols_subject[,order_cols]
    
    if(ncol(ordered_cols) == 4){
      for(i in 1:nrow(ordered_cols)){
        area <- area_under_curve(x_values,ordered_cols[i,])
        auc_df[i,index_y] <- area
      }
    }
    index_y <- index_y+1
    
  }
  auc_df <- as.data.frame(auc_df)
  colnames(auc_df) <- levels(subjects)
  rownames(auc_df) <- rownames(Intensity)
  return(t(auc_df))
}

Correlation_OneMethod <- function(input_clinical, input_expression, method_cor, filter){
  names1 <- colnames(input_clinical)
  names2 <- colnames(input_expression)
  Out_mat <- data.frame()
  for(i in 1:ncol(input_clinical)){
    for(j in 1:ncol(input_expression)){
      if(length(input_expression[,j][complete.cases(input_expression[,j])]) >= filter){
        correlation  <- cor.test(as.numeric(input_clinical[,i]),as.numeric(input_expression[,j]), method = method_cor)
        Out_mat <- rbind(Out_mat,c(names1[i], names2[j],correlation$estimate, correlation$p.value))
      }
    }
  }

  colnames(Out_mat) <- c("Clinical", "Protein", "correlation", "pVal")
  Out_mat$pVal <- as.numeric(Out_mat$pVal)
  Out_mat$p.adjust <- p.adjust(Out_mat$pVal, method = "BH")
  return(Out_mat)
}

Correlation_TwoMethod <- function(input_clinical, input_expression,filter){
  shapiro_res = data.frame()
  for(i in 1:ncol(input_clinical)){
    shap_res = shapiro.test(as.numeric(input_clinical[,i]))
    shapiro_res <- rbind(shapiro_res, shap_res[["p.value"]])
  }
  not_normal <- which(shapiro_res <= 0.05)
  
  kendall = Correlation_OneMethod(input_clinical[which(shapiro_res >= 0.05)], input_expression, "kendall",filter)

  pearson = Correlation_OneMethod(input_clinical[which(shapiro_res <= 0.05)], input_expression, "pearson",filter)
  
  Out_mat <- rbind(pearson,kendall)
  Out_mat$p.adjust <- p.adjust(Out_mat$pVal, method = "BH")
  
  return(Out_mat)
  
}

Delta_calculator <- function( Intensity_matrix, subjects,groups, timepoint_1, timepoint_2){
  idx <- which(groups == timepoint_1)
  idx <- append(idx, which(groups == timepoint_2))
  Intensity_matrix <- Intensity_matrix[,idx]
  subjects_small <- unique(subjects)
  delta_matrix <- data.frame()
  for(sbjct in subjects_small){
    tps <- which(subjects[idx] == sbjct)
    tp1 <- Intensity_matrix[,tps[1]]
    tp2 <- Intensity_matrix[,tps[2]]
    delta <- tp2-tp1
    delta_matrix <- rbind(delta_matrix, delta)
  }
  rownames(delta_matrix) <- subjects_small
  colnames(delta_matrix) <- rownames(Intensity_matrix)
  return(t(delta_matrix))
}

Znorm <- function(input_mat){
  Output_mat = data.frame()
  for(i in 1:ncol(input_mat)){
    col = as.numeric(input_mat[,i])
    mean_col = mean(col,na.rm = T)
    sd_col = sd(col, na.rm = T)
    Output_mat <- rbind(Output_mat,(col-mean_col)/sd_col)
  }
  return(t(Output_mat))
} 

##### Clinical Data ####
Clinical_data <- as.data.frame(t(na.omit(t(read_excel("clinicalData_ExerOME.xlsx")))))
#Correcting typos
Clinical_data[,55]
Clinical_data[22,55] <- 70.791
Clinical_data[,34]
Clinical_data[15,34] <- 0.486

##### Plasma Proteome #####
#### Data Preprocessing ####

#Importing data 
plasma_proteome <- read.delim("20210107_104622_20200106_ExerOME_plasma_DIA_updated_Report_protein_quant_pivot_Atul.txt", header = T)
#plasma_proteome <- read.csv("20211129_140619_20191230_Report_JKL_Q_value_filtering.csv", header = T)

# Extracting PG.Quantity columns
Intensity_plasma <- plasma_proteome[,grep("PG.Quantity", colnames(plasma_proteome))]

rownames(Intensity_plasma) <- plasma_proteome$PG.Genes

#Changing "Filtered" missing values to NAs
Intensity_plasma[] <- lapply(Intensity_plasma, function(x) as.numeric(as.character(x)))
dim(Intensity_plasma)

#Log2 transform
Intensity_plasma = as.matrix(log2(Intensity_plasma))

# Cleaning up sample labels
colnames(Intensity_plasma) <- sub("20191121_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Nov_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191122_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Nov_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191123_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Nov_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191126_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Nov_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191127_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Nov_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191202_QE8_nLC9_ASD_SA_ExerOME_plasma_", "Dec_", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub("20191230_QE8_nLC14_ASD_SA_ExerOME_plasma_", "Dec_", colnames(Intensity_plasma))

#Sample 2 had weird label. Changing it.
colnames(Intensity_plasma) <- sub("_20191121095600","", colnames(Intensity_plasma))
# One sample was mislabelled "Pre" instead of "pre"
colnames(Intensity_plasma) <- sub("Pre","pre", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- sub(".raw.PG.Quantity","", colnames(Intensity_plasma))
colnames(Intensity_plasma) <- substring(colnames(Intensity_plasma), 6)
colnames(Intensity_plasma)[10:99] <- substring(colnames(Intensity_plasma)[10:99], 2)
colnames(Intensity_plasma)[100:ncol(Intensity_plasma)] <- substring(colnames(Intensity_plasma)[100:ncol(Intensity_plasma)], 3)


Intensity_plasma[,13] <- Intensity_plasma[,150] #Replacing duplicate 
Intensity_plasma[,37] <- Intensity_plasma[,149] #Replacing duplicate
Intensity_plasma <- Intensity_plasma[,1:148] 

#Defining groups
sample_name = strsplit(gsub("^_", "", colnames(Intensity_plasma)), "_")
df = S4Vectors::DataFrame(
  time = sapply(sample_name, "[[", 2),
  subject = sapply(sample_name, "[[", 3),
  nLC = sapply(sample_name, "[[", 1))
rownames(df) = colnames(Intensity_plasma)
grps_plasma <- paste0(df$time)
subj_plasma <- paste0(df$subject)
nlc_plasma <- paste0(df$nLC)


plotQC(Intensity_plasma, labels = colnames(Intensity_plasma), panel = "quantify", grps = grps_plasma)


# Filtering criteria. n= in number of groups the percentage valid values should be found
Intensity_plasma <- selectGrps(Intensity_plasma, grps_plasma, 0.5, n=2)
Intensity_plasma_to_PCA <- selectGrps(Intensity_plasma, grps_plasma, 1, n=6)
dim(Intensity_plasma)

pca_plasma <- prcomp(t(Intensity_plasma_to_PCA))
ggplot(mapping = aes(pca_plasma$x[,1],pca_plasma$x[,2], color = df$subject)) + geom_point() #We can see there is no batch effect, but we can see some outliers

tsne_plasma <- Rtsne(t(Intensity_plasma_to_PCA), perplexity = 5)
ggplot(mapping = aes(tsne_plasma$Y[,1],tsne_plasma$Y[,2], color = df$subject)) + geom_point() #We can see there is no batch effect

#Remove batch effects
  #No batch effect

set.seed(123)
#Imp_Intensity_plasma <- scImpute(Intensity_plasma, 0.55, grps) # This one can be left out.
Intensity_plasma_NI <- Intensity_plasma #Save the not imputed version of the Intensity plasma
Intensity_plasma <- tImpute(Intensity_plasma) # Tailed-based imputation as in Perseus

#Median Scaling
Intensity_plasma_mediascale <- medianScaling(Intensity_plasma)

# Median normalize
data_median <- apply(Intensity_plasma, 2, median, na.rm=TRUE)
Intensity_plasma <- Intensity_plasma[] - data_median[col(Intensity_plasma)[]]
data_median <- apply(Intensity_plasma_NI, 2, median, na.rm=TRUE)
Intensity_plasma_NI <- Intensity_plasma_NI[] - data_median[col(Intensity_plasma)[]]



#### LIMMA ####
Group <- factor(grps_plasma, levels=c("pre","T0","T30","T60","T180","T24"))
design <- model.matrix(~0+Group)

# Calculating covariance for 
# On imputation data
corfit_plasma_II <- duplicateCorrelation(Intensity_plasma, design, block=df$subject)
corfit_plasma_II$consensus

# On non-imputed data
corfit_plasma_INI <- duplicateCorrelation(Intensity_plasma_NI, design, block=df$subject)
corfit_plasma_INI$consensus


# lmFit which in this case uses the gls.series function.
#On imputated data
fit_plasma_II <- lmFit(Intensity_plasma, design, block = df$subject, correlation=corfit_plasma_II$consensus)

# On non-imputed data
fit_plasma_INI <- lmFit(Intensity_plasma_NI, design, block = df$subject, correlation=corfit_plasma_INI$consensus)


# Making contrast for proteins differentially expressed at each time points compared to basal
cm <- makeContrasts("GroupT0-Grouppre", "GroupT30-Grouppre", "GroupT60-Grouppre", "GroupT180-Grouppre", "GroupT24-Grouppre", levels=design)

fit2_plasma <- contrasts.fit(fit_plasma_II, cm)
# Performing moderated t-test
fit2_plasma <- eBayes(fit2_plasma)

results <- decideTests(fit2_plasma,method="global")
summary(results)

# Extracting individuel comparisons
topTable(fit2_plasma, coef = 1, p.value = 0.05, number = Inf) # T0 - Pre
topTable(fit2_plasma, coef = 2, p.value = 0.05) # T30 - Pre
topTable(fit2_plasma, coef = 3, p.value = 0.05) # T60 - Pre
topTable(fit2_plasma, coef = 4, p.value = 0.05) # T180 - Pre
topTable(fit2_plasma, coef = 5, p.value = 0.05) # T24 - Pre

F_statistics_table_plasma <- topTable(fit2_plasma, coef = NULL, p.value = 0.05, number = Inf)
nrow(F_statistics_table_plasma)
write.table(F_statistics_table_plasma,file = "Log_Plasma",row.names = rownames(F_statistics_table_plasma), col.names = colnames(F_statistics_table_plasma))
# VOLCANO PLOTS
  #T0 vs Pre
T0_Pre_plasma <- topTable(fit2_plasma, coef = 1, number = Inf)
T0_Pre_plasma$diffexpressed <- "NO"
T0_Pre_plasma$diffexpressed[T0_Pre_plasma$logFC > 0 & T0_Pre_plasma$adj.P.Val < 0.05] <- "UP"
T0_Pre_plasma$diffexpressed[T0_Pre_plasma$logFC < 0 & T0_Pre_plasma$adj.P.Val < 0.05] <- "DOWN"
T0_Pre_plasma$delabel <- NA
T0_Pre_plasma$delabel[T0_Pre_plasma$diffexpressed != "NO"] <- rownames(T0_Pre_plasma)[T0_Pre_plasma$diffexpressed != "NO"]
plot_T0_Pre <- ggplot(data=T0_Pre_plasma, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T0 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T0_Pre

table(T0_Pre_plasma$logFC > 0 & T0_Pre_plasma$adj.P.Val < 0.05)["TRUE"]
table(T0_Pre_plasma$logFC < 0 & T0_Pre_plasma$adj.P.Val < 0.05)["TRUE"]

#T30 vs Pre
T30_Pre_plasma <- topTable(fit2_plasma, coef = 2, number = Inf)
T30_Pre_plasma$diffexpressed <- "NO"
T30_Pre_plasma$diffexpressed[T30_Pre_plasma$logFC > 0 & T30_Pre_plasma$adj.P.Val < 0.05] <- "UP"
T30_Pre_plasma$diffexpressed[T30_Pre_plasma$logFC < 0 & T30_Pre_plasma$adj.P.Val < 0.05] <- "DOWN"
T30_Pre_plasma$delabel <- NA
T30_Pre_plasma$delabel[T30_Pre_plasma$diffexpressed != "NO"] <- rownames(T30_Pre_plasma)[T30_Pre_plasma$diffexpressed != "NO"]
plot_T30_Pre_plasma <- ggplot(data=T30_Pre_plasma, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T30 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T30_Pre_plasma

table(T30_Pre_plasma$logFC > 0 & T30_Pre_plasma$adj.P.Val < 0.05)["TRUE"]
table(T30_Pre_plasma$logFC < 0 & T30_Pre_plasma$adj.P.Val < 0.05)["TRUE"]

#T60 vs Pre
T60_Pre_plasma <- topTable(fit2_plasma, coef = 3, number = Inf)
T60_Pre_plasma$diffexpressed <- "NO"
T60_Pre_plasma$diffexpressed[T60_Pre_plasma$logFC > 0 & T60_Pre_plasma$adj.P.Val < 0.05] <- "UP"
T60_Pre_plasma$diffexpressed[T60_Pre_plasma$logFC < 0 & T60_Pre_plasma$adj.P.Val < 0.05] <- "DOWN"
T60_Pre_plasma$delabel <- NA
T60_Pre_plasma$delabel[T60_Pre_plasma$diffexpressed != "NO"] <- rownames(T60_Pre_plasma)[T60_Pre_plasma$diffexpressed != "NO"]
plot_T60_Pre_plasma <- ggplot(data=T60_Pre_plasma, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T60 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T60_Pre_plasma

table(T60_Pre_plasma$logFC > 0 & T60_Pre_plasma$adj.P.Val < 0.05)["TRUE"]
table(T60_Pre_plasma$logFC < 0 & T60_Pre_plasma$adj.P.Val < 0.05)["TRUE"]

#T180 vs Pre
T180_Pre_plasma <- topTable(fit2_plasma, coef = 4, number = Inf)
T180_Pre_plasma$diffexpressed <- "NO"
T180_Pre_plasma$diffexpressed[T180_Pre_plasma$logFC > 0 & T180_Pre_plasma$adj.P.Val < 0.05] <- "UP"
T180_Pre_plasma$diffexpressed[T180_Pre_plasma$logFC < 0 & T180_Pre_plasma$adj.P.Val < 0.05] <- "DOWN"
T180_Pre_plasma$delabel <- NA
T180_Pre_plasma$delabel[T180_Pre_plasma$diffexpressed != "NO"] <- rownames(T180_Pre_plasma)[T180_Pre_plasma$diffexpressed != "NO"]
plot_T180_Pre_plasma <- ggplot(data=T180_Pre_plasma, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3.5, 3.5) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T180 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T180_Pre_plasma

table(T180_Pre_plasma$logFC > 0 & T180_Pre_plasma$adj.P.Val < 0.05)["TRUE"]
table(T180_Pre_plasma$logFC < 0 & T180_Pre_plasma$adj.P.Val < 0.05)["TRUE"]

#T24 vs Pre
T24_Pre_plasma <- topTable(fit2_plasma, coef = 5, number = Inf)
T24_Pre_plasma$diffexpressed <- "NO"
T24_Pre_plasma$diffexpressed[T24_Pre_plasma$logFC > 0 & T24_Pre_plasma$adj.P.Val < 0.05] <- "UP"
T24_Pre_plasma$diffexpressed[T24_Pre_plasma$logFC < 0 & T24_Pre_plasma$adj.P.Val < 0.05] <- "DOWN"
T24_Pre_plasma$delabel <- NA
T24_Pre_plasma$delabel[T24_Pre_plasma$diffexpressed != "NO"] <- rownames(T24_Pre_plasma)[T24_Pre_plasma$diffexpressed != "NO"]
plot_T24_Pre_plasma <- ggplot(data=T24_Pre_plasma, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3.5, 3.5) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T24 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T24_Pre_plasma

table(T24_Pre_plasma$logFC > 0 & T24_Pre_plasma$adj.P.Val < 0.05)["TRUE"]
table(T24_Pre_plasma$logFC < 0 & T24_Pre_plasma$adj.P.Val < 0.05)["TRUE"]

hej_LF <- topTable(fit2_plasma, coef = 2, number = 100) # T30 - Pre
hej_LF$AveExpr

F_statistics_table_plasma <- topTable(fit2_plasma, coef = NULL, p.value = 0.05, number = Inf)
nrow(F_statistics_table_plasma)

#### Clusters ####

Sigs1 <- T0_Pre_plasma[T0_Pre_plasma$diffexpressed != "NO",]
Sigs2 <- T30_Pre_plasma[T30_Pre_plasma$diffexpressed != "NO",]
Sigs3 <- T60_Pre_plasma[T60_Pre_plasma$diffexpressed != "NO",]
Sigs4 <- T180_Pre_plasma[T180_Pre_plasma$diffexpressed != "NO",]
Sigs5 <- T24_Pre_plasma[T24_Pre_plasma$diffexpressed != "NO",]
Sigs_list <- list(Sigs1 = Sigs1, Sigs2 = Sigs2, Sigs3 = Sigs3, Sigs4 = Sigs4, Sigs5 = Sigs5)

gene_names <- c()
for(i in seq_along(Sigs_list)){
  gene_name = rownames(Sigs_list[[i]])
  gene_names <- append(gene_names, gene_name)
}
gene_names <- unique(gene_names)

List_behavior <- list()
Time_points <- c("T0", "T30", "T60", "T180", "T24")
for(gene in gene_names){
  gene_behavior <- c()
  for(i in seq_along(Sigs_list)){
    index <- which(rownames(Sigs_list[[i]]) == gene)
    if(length(index) != 0){
      gene_behavior <- append(gene_behavior, Sigs_list[[i]]$diffexpressed[index])
    } else{
      gene_behavior <- append(gene_behavior, NA)
    }
  }
  names(gene_behavior) <- Time_points
  List_behavior[[gene]] <- gene_behavior
}


##### Urine #####
#### Data Preprocessing ####

#Import the data
urine_proteome <- read.delim("20210107_144259_20200107_ExerOME_urine_DIA_singleshot_Report_protein_quant_pivot_Atul.txt", header = T)

#Extract the PG.quantity cols
Intensity_urine <- urine_proteome[,grep("PG.Quantity", colnames(urine_proteome))]
rownames(Intensity_urine) <- urine_proteome$PG.ProteinAccessions

#Change filtered for NA
Intensity_urine[] <- lapply(Intensity_urine, function(x) as.numeric(as.character(x)))
#Log2 transform data
Intensity_urine <- as.matrix(log2(Intensity_urine))

#Cleaning colnames
colnames(Intensity_urine) <- sub("20200310_QE8_nLC14_ASD_SA_ExerOME_urine_", "", colnames(Intensity_urine))
colnames(Intensity_urine) <- sub(".raw.PG.Quantity", "", colnames(Intensity_urine))
colnames(Intensity_urine) <- substring(colnames(Intensity_urine),6)
colnames(Intensity_urine)[10:75] <- substring(colnames(Intensity_urine)[10:75],2)

#Defining Groups
sample_name = strsplit(gsub("^_", "", colnames(Intensity_urine)), "_")
df_urine = S4Vectors::DataFrame(
  time = sapply(sample_name, "[[", 1),
  subject = sapply(sample_name, "[[", 2))
rownames(df_urine) = colnames(Intensity_urine)
grps_urine <- paste0(df_urine$time)
subj_urine <- paste0(df_urine$subject)


plotQC(Intensity_urine, labels = colnames(Intensity_urine), panel = "quantify", grps = grps_urine)

# Filtering criteria. n= in number of groups the percentage valid values should be found
Intensity_urine <- selectGrps(Intensity_urine, grps_urine, 0.5, n=2)
Intensity_tsne_U <- selectGrps(Intensity_urine, grps_urine, 1, n=3)

tsne_urine <- Rtsne(t(Intensity_tsne_U), perplexity = 5)
ggplot(mapping = aes(tsne_urine$Y[,1],tsne_urine$Y[,2], color = subj_urine)) + geom_point()
#Remove Batch Effect
#Which are the different batches?

# Median normalize
data_median <- apply(Intensity_urine, 2, median, na.rm=TRUE)
Intensity_urine <- Intensity_urine[] - data_median[col(Intensity_urine)[]]

#Imputation
set.seed(123)
Intensity_urine_NI <- Intensity_urine #Save the not imputed version of the Intensity plasma
Intensity_urine <- tImpute(Intensity_urine) # Tailed-based imputation as in Perseus

#### LIMMA ####
Group <- factor(grps_urine, levels=c("pre","T0","T24"))
design <- model.matrix(~0+Group)

#Imputed
corfit_urine_I <- duplicateCorrelation(Intensity_urine,design,block = df_urine$subject)
corfit_urine_I$consensus
#Not Imputed
corfit_urine_NI <- duplicateCorrelation(Intensity_urine_NI,design,block = df_urine$subject)
corfit_urine_NI$consensus

fit_urine_I <- lmFit(Intensity_urine,design, block = df_urine$subject, correlation = corfit_urine_I$consensus)
fit_urine_NI <- lmFit(Intensity_urine_NI,design, block = df_urine$subject, correlation = corfit_urine_NI$consensus)

# Making contrast for proteins differentially expressed at each time points compared to basal
cm <- makeContrasts("GroupT0-Grouppre", "GroupT24-Grouppre", levels=design)

fit2_urine <- contrasts.fit(fit_urine_I, cm)
# Performing moderated t-test
fit2_urine <- eBayes(fit2_urine)

results <- decideTests(fit2_urine,method="global")
summary(results)

# Extracting individuel comparisons
topTable(fit2_urine, coef = 1, p.value = 0.05, number = Inf) # T0 - Pre
topTable(fit2_urine, coef = 2, p.value = 0.05) #T24 - Pre

F_statistics_table_urine <- topTable(fit2_urine, coef = NULL, p.value = 0.05, number = Inf)
nrow(F_statistics_table_urine)

#Volcano Plots
#T0 vs Pre
T0_Pre <- topTable(fit2_urine, coef = 1, number = Inf)
T0_Pre$diffexpressed <- "NO"
T0_Pre$diffexpressed[T0_Pre$logFC > 0 & T0_Pre$adj.P.Val < 0.05] <- "UP"
T0_Pre$diffexpressed[T0_Pre$logFC < 0 & T0_Pre$adj.P.Val < 0.05] <- "DOWN"
T0_Pre$delabel <- NA
T0_Pre$delabel[T0_Pre$diffexpressed != "NO"] <- rownames(T0_Pre)[T0_Pre$diffexpressed != "NO"]
plot_T0_Pre <- ggplot(data=T0_Pre, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T0 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T0_Pre


table(T0_Pre$logFC > 0 & T0_Pre$adj.P.Val < 0.05)["TRUE"]
table(T0_Pre$logFC < 0 & T0_Pre$adj.P.Val < 0.05)["TRUE"]

#T24 vs Pre
T24_Pre <- topTable(fit2_urine, coef = 2, number = Inf)
T24_Pre$diffexpressed <- "NO"
T24_Pre$diffexpressed[T24_Pre$logFC > 0 & T24_Pre$adj.P.Val < 0.05] <- "UP"
T24_Pre$diffexpressed[T24_Pre$logFC < 0 & T24_Pre$adj.P.Val < 0.05] <- "DOWN"
T24_Pre$delabel <- NA
T24_Pre$delabel[T24_Pre$diffexpressed != "NO"] <- rownames(T24_Pre)[T24_Pre$diffexpressed != "NO"]
plot_T24_Pre <- ggplot(data=T24_Pre, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3.5, 3.5) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T24 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T24_Pre

table(T24_Pre$logFC > 0 & T24_Pre$adj.P.Val < 0.05)["TRUE"]
table(T24_Pre$logFC < 0 & T24_Pre$adj.P.Val < 0.05)["TRUE"]


##### Saliva Proteome #####
#### Data Preprocessing ####

#Load the data
saliva_proteome <- read.delim("20210108_101357_20210107_ExerOME_saliva_DIA_singleshot_Report_protein_quant_pivot_Atul.txt", header =  T)

#Extracting the PG.Quantity colunms
Intensity_saliva <- saliva_proteome[, grep("PG.Quantity", colnames(saliva_proteome))]
rownames(Intensity_saliva) <- saliva_proteome$PG.ProteinAccessions

#Change the filtered to NA
Intensity_saliva[] <- lapply(Intensity_saliva, function(x) as.numeric(as.character(x)))
dim(Intensity_saliva)

#Log2 transform
Intensity_saliva = as.matrix(log2(Intensity_saliva))

#Correcting column names
colnames(Intensity_saliva) <- sub("20200203_QE8_", "Feb_", colnames(Intensity_saliva))
colnames(Intensity_saliva) <- sub("ASD_SA_ExerOME_saliva_", "", colnames(Intensity_saliva))
colnames(Intensity_saliva) <- sub("20200103_QE8_", "Jen_", colnames(Intensity_saliva))

  #Change Pre to pre
colnames(Intensity_saliva) <- sub("Pre","pre", colnames(Intensity_saliva))

colnames(Intensity_saliva) <- sub(".raw.PG.Quantity","",colnames(Intensity_saliva))
colnames(Intensity_saliva) <- sub("R","",colnames(Intensity_saliva))
colnames(Intensity_saliva) <- substring(colnames(Intensity_saliva), 6)
colnames(Intensity_saliva)[10:150] <- substring(colnames(Intensity_saliva)[10:150], 2)
colnames(Intensity_saliva)[100:150] <- substring(colnames(Intensity_saliva)[100:150], 2)

#Defining groups
sample_name = strsplit(gsub("^_", "", colnames(Intensity_saliva)), "_")
df = S4Vectors::DataFrame(
  time = sapply(sample_name, "[[", 3),
  subject = sapply(sample_name, "[[", 4),
  nLC = sapply(sample_name, "[[", 2),
  date = sapply(sample_name, "[[", 1))
rownames(df) = colnames(Intensity_saliva)
grps_saliva <- paste0(df$time)
nlc <- paste0(df$nLC)
dates <- paste0(df$date)
subj_saliva <- paste0(df$subject)


plotQC(Intensity_saliva, labels = colnames(Intensity_saliva), panel = "quantify", grps = grps_urine)


# Filtering criteria. n= in number of groups the percentage valid values should be found
Intensity_saliva <- selectGrps(Intensity_saliva, grps_saliva, 0.5, n=2)
Intensity_saliva_to_pca <- selectGrps(Intensity_saliva, grps_saliva, 1, n=6)


pca <- prcomp(t(Intensity_saliva_to_pca))
ggplot(mapping = aes(pca$x[,1], pca$x[,2], color = dates)) + geom_point() 

tsne <- Rtsne(t(Intensity_saliva_to_pca))
ggplot(mapping = aes(tsne$Y[,1], tsne$Y[,2], color = dates)) + geom_point()

umap <- umap(t(Intensity_saliva_to_pca))
ggplot(mapping = aes(umap$layout[,1], umap$layout[,2], color = dates)) + geom_point()

#Remove batch effect
df
df$batch <- "NA"

df$batch[which(nlc == "nLC2")] <- "Cluster 1"
df$batch[which(nlc == "nLC14")] <- "Cluster 2"

Group <- factor(grps_saliva, levels=c("pre","T0","T30","T60","T180","T24"))
design <- model.matrix(~0+Group+df$batch)
colnames(design) <- make.names(colnames(design))
hej_saliva <- removeBatchEffect(Intensity_saliva, batch=as.factor(df$batch), covariates=NULL,design=design)
hej_saliva_to_pca <- selectGrps(hej_saliva, grps_saliva, 1, n=6)


#pca <- prcomp(t(hej_saliva_to_pca))
#ggplot(mapping = aes(pca$x[,1], pca$x[,2], color = dates)) + geom_point() 

tsne <- Rtsne(t(hej_saliva_to_pca), perplexity = 6, theta = 0)
#ggplot(mapping = aes(tsne$Y[,1], tsne$Y[,2], color = dates)) + geom_point() 
ggplot(mapping = aes(tsne$Y[,1], tsne$Y[,2], color = subj_saliva)) + geom_point() #We can see clustering by subject

#umap <- umap(t(hej_saliva_to_pca),)
#ggplot(mapping = aes(umap$layout[,1], umap$layout[,2], color = dates)) + geom_point()
#ggplot(mapping = aes(umap$layout[,1], umap$layout[,2], color = subj_saliva)) + geom_point() #We can see clustering by subject

# Median scaling
Intensity_saliva_medianscal <- medianScaling(Intensity_saliva)

# Median normalize
data_median <- apply(Intensity_saliva, 2, median, na.rm=TRUE)
Intensity_saliva <- Intensity_saliva[] - data_median[col(Intensity_saliva)[]]
#data_median <- apply(hej_saliva, 2, median, na.rm=TRUE)
#hej_saliva <- hej_saliva[] - data_median[col(hej_saliva)[]]

#Imputation
Intensity_saliva_NI <- Intensity_saliva #Save the not imputed version of the Intensity plasma
Intensity_saliva <- tImpute(Intensity_saliva) # Tailed-based imputation as in Perseus
Intensity_saliva_medianscal <- tImpute(Intensity_saliva_medianscal)
#hej_saliva_NI <- hej_saliva#Save the not imputed version of the hej
#hej_saliva <- tImpute(hej_saliva) # Tailed-based imputation as in Perseus

#### LIMMA ####
#Imputed
corfit_saliva_II <- duplicateCorrelation(Intensity_saliva, design, block = df$subject)
corfit_saliva_II$consensus
#corfit_saliva_hI <- duplicateCorrelation(hej_saliva, design, block = df$subject)
#corfit_saliva_hI$consensus

#Not Imputed
corfit_saliva_INI <- duplicateCorrelation(Intensity_saliva_NI, design, block = df$subject)
corfit_saliva_INI$consensus
#corfit_saliva_hNI <- duplicateCorrelation(hej_saliva_NI, design, block = df$subject)
#corfit_saliva_hNI$consensus

fit_saliva_II <- lmFit(Intensity_saliva, design, block = df$subject, correlation = corfit_saliva_II$consensus)
#fit_saliva_hI <- lmFit(hej_saliva, design, block = df$subject, correlation = corfit_saliva_hI$consensus)
fit_saliva_INI <- lmFit(Intensity_saliva_NI, design, block = df$subject, correlation = corfit_saliva_INI$consensus)
#fit_saliva_hNI <- lmFit(hej_saliva_NI, design, block = df$subject, correlation = corfit_saliva_hNI$consensus)

# Making contrast for proteins differentially expressed at each time points compared to basal
cm <- makeContrasts("GroupT0-Grouppre", "GroupT30-Grouppre", "GroupT60-Grouppre", "GroupT180-Grouppre", "GroupT24-Grouppre", levels=design)

fit2_saliva <- contrasts.fit(fit_saliva_II, cm)
# Performing moderated t-test
fit2_saliva <- eBayes(fit2_saliva)

results <- decideTests(fit2_saliva,method="global")
summary(results)

# Extracting individual comparisons
topTable(fit2_saliva, coef = 1, p.value = 0.05, number = Inf) # T0 - Pre
topTable(fit2_saliva, coef = 2, p.value = 0.05) # T30 - Pre
topTable(fit2_saliva, coef = 3, p.value = 0.05) # T60 - Pre
topTable(fit2_saliva, coef = 4, p.value = 0.05) # T180 - Pre
topTable(fit2_saliva, coef = 5, p.value = 0.05) # T24 - Pre

F_statistics_table_saliva <- topTable(fit2_saliva, coef = NULL, p.value = 0.05, number = Inf)
nrow(F_statistics_table_saliva)
write.table(F_statistics_table_saliva, file = "Log_saliva.csv", sep = ",")
# VOLCANO PLOTS
#T0 vs Pre
T0_Pre_saliva <- topTable(fit2_saliva, coef = 1, number = Inf)
T0_Pre_saliva$diffexpressed <- "NO"
T0_Pre_saliva$diffexpressed[T0_Pre_saliva$logFC > 0 & T0_Pre_saliva$adj.P.Val < 0.05] <- "UP"
T0_Pre_saliva$diffexpressed[T0_Pre_saliva$logFC < 0 & T0_Pre_saliva$adj.P.Val < 0.05] <- "DOWN"
T0_Pre_saliva$delabel <- NA
T0_Pre_saliva$delabel[T0_Pre_saliva$diffexpressed != "NO"] <- rownames(T0_Pre_saliva)[T0_Pre_saliva$diffexpressed != "NO"]
plot_T0_Pre <- ggplot(data=T0_Pre_saliva, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T0 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T0_Pre

table(T0_Pre_saliva$logFC > 0 & T0_Pre_saliva$adj.P.Val < 0.05)["TRUE"]
table(T0_Pre_saliva$logFC < 0 & T0_Pre_saliva$adj.P.Val < 0.05)["TRUE"]

#T30 vs Pre
T30_Pre_saliva <- topTable(fit2_saliva, coef = 2, number = Inf)
T30_Pre_saliva$diffexpressed <- "NO"
T30_Pre_saliva$diffexpressed[T30_Pre_saliva$logFC > 0 & T30_Pre_saliva$adj.P.Val < 0.05] <- "UP"
T30_Pre_saliva$diffexpressed[T30_Pre_saliva$logFC < 0 & T30_Pre_saliva$adj.P.Val < 0.05] <- "DOWN"
T30_Pre_saliva$delabel <- NA
T30_Pre_saliva$delabel[T30_Pre_saliva$diffexpressed != "NO"] <- rownames(T30_Pre_saliva)[T30_Pre_saliva$diffexpressed != "NO"]
plot_T30_Pre <- ggplot(data=T30_Pre_saliva, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T30 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T30_Pre

table(T30_Pre_saliva$logFC > 0 & T30_Pre_saliva$adj.P.Val < 0.05)["TRUE"]
table(T30_Pre_saliva$logFC < 0 & T30_Pre_saliva$adj.P.Val < 0.05)["TRUE"]

#T60 vs Pre
T60_Pre_saliva <- topTable(fit2_saliva, coef = 3, number = Inf)
T60_Pre_saliva$diffexpressed <- "NO"
T60_Pre_saliva$diffexpressed[T60_Pre_saliva$logFC > 0 & T60_Pre_saliva$adj.P.Val < 0.05] <- "UP"
T60_Pre_saliva$diffexpressed[T60_Pre_saliva$logFC < 0 & T60_Pre_saliva$adj.P.Val < 0.05] <- "DOWN"
T60_Pre_saliva$delabel <- NA
T60_Pre_saliva$delabel[T60_Pre_saliva$diffexpressed != "NO"] <- rownames(T60_Pre_saliva)[T60_Pre_saliva$diffexpressed != "NO"]
plot_T60_Pre <- ggplot(data=T60_Pre_saliva, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3, 3) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T60 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T60_Pre

table(T60_Pre_saliva$logFC > 0 & T60_Pre_saliva$adj.P.Val < 0.05)["TRUE"]
table(T60_Pre_saliva$logFC < 0 & T60_Pre_saliva$adj.P.Val < 0.05)["TRUE"]

#T180 vs Pre
T180_Pre_saliva <- topTable(fit2_saliva, coef = 4, number = Inf)
T180_Pre_saliva$diffexpressed <- "NO"
T180_Pre_saliva$diffexpressed[T180_Pre_saliva$logFC > 0 & T180_Pre_saliva$adj.P.Val < 0.05] <- "UP"
T180_Pre_saliva$diffexpressed[T180_Pre_saliva$logFC < 0 & T180_Pre_saliva$adj.P.Val < 0.05] <- "DOWN"
T180_Pre_saliva$delabel <- NA
T180_Pre_saliva$delabel[T180_Pre_saliva$diffexpressed != "NO"] <- rownames(T180_Pre_saliva)[v$diffexpressed != "NO"]
plot_T180_Pre <- ggplot(data=T180_Pre_saliva, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3.5, 3.5) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T180 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T180_Pre

table(T180_Pre_saliva$logFC > 0 & T180_Pre_saliva$adj.P.Val < 0.05)["TRUE"]
table(T180_Pre_saliva$logFC < 0 & T180_Pre_saliva$adj.P.Val < 0.05)["TRUE"]

#T24 vs Pre
T24_Pre_saliva <- topTable(fit2_saliva, coef = 5, number = Inf)
T24_Pre_saliva$diffexpressed <- "NO"
T24_Pre_saliva$diffexpressed[T24_Pre_saliva$logFC > 0 & T24_Pre_saliva$adj.P.Val < 0.05] <- "UP"
T24_Pre_saliva$diffexpressed[T24_Pre_saliva$logFC < 0 & T24_Pre_saliva$adj.P.Val < 0.05] <- "DOWN"
T24_Pre_saliva$delabel <- NA
T24_Pre_saliva$delabel[T24_Pre_saliva$diffexpressed != "NO"] <- rownames(T24_Pre_saliva)[T24_Pre_saliva$diffexpressed != "NO"]
plot_T24_Pre <- ggplot(data=T24_Pre_saliva, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label = delabel)) + 
  geom_point(alpha = 0.4, size = 3)+ 
  xlim(-3.5, 3.5) + 
  scale_color_manual(values=c("blue", "grey43", "red")) + 
  ggtitle("T24 vs Pre") + 
  theme_minimal() + 
  geom_text_repel()
plot_T24_Pre

table(T24_Pre$logFC > 0 & T24_Pre$adj.P.Val < 0.05)["TRUE"]
table(T24_Pre$logFC < 0 & T24_Pre$adj.P.Val < 0.05)["TRUE"]

hej_saliva_LF <- topTable(fit2_saliva, coef = 2, number = 100) # T30 - Pre
hej_saliva_LF$AveExpr

#### Clusters ####
Sigs1 <- T0_Pre_saliva[T0_Pre_saliva$diffexpressed != "NO",]
Sigs2 <- T30_Pre_saliva[T30_Pre_saliva$diffexpressed != "NO",]
Sigs3 <- T60_Pre_saliva[T60_Pre_saliva$diffexpressed != "NO",]
Sigs4 <- T180_Pre_saliva[T180_Pre_saliva$diffexpressed != "NO",]
Sigs5 <- T24_Pre_saliva[T24_Pre_saliva$diffexpressed != "NO",]
Sigs_list_saliva <- list(Sigs1 = Sigs1, Sigs2 = Sigs2, Sigs3 = Sigs3, Sigs4 = Sigs4, Sigs5 = Sigs5)

gene_names_saliva <- c()
for(i in seq_along(Sigs_list_saliva)){
  gene_name = rownames(Sigs_list_saliva[[i]])
  gene_names_saliva <- append(gene_names_saliva, gene_name)
}
gene_names_saliva <- unique(gene_names_saliva)

List_behavior_saliva <- list()
Time_points <- c("T0", "T30", "T60", "T180", "T24")
for(gene in gene_names_saliva){
  gene_behavior <- c()
  for(i in seq_along(Sigs_list_saliva)){
    index <- which(rownames(Sigs_list_saliva[[i]]) == gene)
    if(length(index) != 0){
      gene_behavior <- append(gene_behavior, Sigs_list_saliva[[i]]$diffexpressed[index])
    } else{
      gene_behavior <- append(gene_behavior, NA)
    }
  }
  names(gene_behavior) <- Time_points
  List_behavior_saliva[[gene]] <- gene_behavior
}

#### Correlations ####
Z_clinical <- Znorm(Clinical_data[c(9:60,62,63)])
Z_clinical <- as.data.frame(Z_clinical)
colnames(Z_clinical) <- colnames(Clinical_data[c(9:60,62,63)])

#Plasma
 #Imputed
Z_Intensity_plasma <- t(Znorm(t(Intensity_plasma)))
rownames(Z_Intensity_plasma) <- rownames(Intensity_plasma)
Corr_plasma_pre <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_plasma)[which(grps_plasma == "pre"),],10)

Intensity_plasma_delta_T0 <- Znorm(t(Delta_calculator(Intensity_plasma, subj_plasma, grps_plasma, "pre", "T0"))) 
colnames(Intensity_plasma_delta_T0) <- rownames(Intensity_plasma)
Corr_plasma_delta_T0 <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T0,10) #Correlation delta pre vs T0

Intensity_plasma_delta_T30 <- Znorm(t(Delta_calculator(Intensity_plasma, subj_plasma, grps_plasma, "pre", "T30"))) 
colnames(Intensity_plasma_delta_T30) <- rownames(Intensity_plasma)
Corr_plasma_delta_T30 <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T30,10) #Correlation delta pre vs T30

Intensity_plasma_delta_T60 <- Znorm(t(Delta_calculator(Intensity_plasma, subj_plasma, grps_plasma, "pre", "T60"))) 
colnames(Intensity_plasma_delta_T60) <- rownames(Intensity_plasma)
Corr_plasma_delta_T60 <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T60,10) #Correlation delta pre vs T60

Intensity_plasma_delta_T180 <- Znorm(t(Delta_calculator(Intensity_plasma, subj_plasma, grps_plasma, "pre", "T180"))) 
colnames(Intensity_plasma_delta_T180) <- rownames(Intensity_plasma)
Corr_plasma_delta_T180 <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T180,10) #Correlation delta pre vs T180

Intensity_plasma_delta_T24 <- Znorm(t(Delta_calculator(Intensity_plasma, subj_plasma, grps_plasma, "pre", "T24"))) 
colnames(Intensity_plasma_delta_T24) <- rownames(Intensity_plasma)
Corr_plasma_delta_T24 <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T24,10) #Correlation delta pre vs T24

#Not Imputed 
Z_Intensity_plasma_NI <- t(Znorm(t(Intensity_plasma_NI)))
rownames(Z_Intensity_plasma_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_pre_10_NI <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_plasma_NI)[which(grps_plasma == "pre"),],10)

Intensity_plasma_delta_T0_NI <- Znorm(t(Delta_calculator(Intensity_plasma_NI, subj_plasma, grps_plasma, "pre", "T0"))) 
colnames(Intensity_plasma_delta_T0_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_delta_T0_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T0_NI,10) #Correlation delta pre vs T0

Intensity_plasma_delta_T30_NI <- Znorm(t(Delta_calculator(Intensity_plasma_NI, subj_plasma, grps_plasma, "pre", "T30"))) 
colnames(Intensity_plasma_delta_T30_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_delta_T30_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T30_NI,10) #Correlation delta pre vs T30

Intensity_plasma_delta_T60_NI <- Znorm(t(Delta_calculator(Intensity_plasma_NI, subj_plasma, grps_plasma, "pre", "T60"))) 
colnames(Intensity_plasma_delta_T60_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_delta_T60_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T60_NI,10) #Correlation delta pre vs T60

Intensity_plasma_delta_T180_NI <- Znorm(t(Delta_calculator(Intensity_plasma_NI, subj_plasma, grps_plasma, "pre", "T180"))) 
colnames(Intensity_plasma_delta_T180_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_delta_T180_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T180_NI,10) #Correlation delta pre vs T180

Intensity_plasma_delta_T24_NI <- Znorm(t(Delta_calculator(Intensity_plasma_NI, subj_plasma, grps_plasma, "pre", "T24"))) 
colnames(Intensity_plasma_delta_T24_NI) <- rownames(Intensity_plasma_NI)
Corr_plasma_delta_T24_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_plasma_delta_T24_NI,10) #Correlation delta pre vs T24

#Saliva
  #Imputed
Z_Intensity_saliva <- t(Znorm(t(Intensity_saliva)))
rownames(Z_Intensity_saliva) <- rownames(Intensity_saliva)
Corr_saliva_pre <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_saliva)[which(grps_saliva == "pre"),],10)

Intensity_saliva_delta_T0 <- Znorm(t(Delta_calculator(Intensity_saliva, subj_saliva, grps_saliva, "pre", "T0"))) 
colnames(Intensity_saliva_delta_T0) <- rownames(Intensity_saliva)
Corr_saliva_delta_T0 <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T0,10) #Correlation delta pre vs T0

Intensity_saliva_delta_T30 <- Znorm(t(Delta_calculator(Intensity_saliva, subj_saliva, grps_saliva, "pre", "T30"))) 
colnames(Intensity_saliva_delta_T30) <- rownames(Intensity_saliva)
Corr_saliva_delta_T30 <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T30,10) #Correlation delta pre vs T30

Intensity_saliva_delta_T60 <- Znorm(t(Delta_calculator(Intensity_saliva, subj_saliva, grps_saliva, "pre", "T60"))) 
colnames(Intensity_saliva_delta_T60) <- rownames(Intensity_saliva)
Corr_saliva_delta_T60 <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T60,10) #Correlation delta pre vs T60

Intensity_saliva_delta_T180 <- Znorm(t(Delta_calculator(Intensity_saliva, subj_saliva, grps_saliva, "pre", "T180"))) 
colnames(Intensity_saliva_delta_T180) <- rownames(Intensity_saliva)
Corr_saliva_delta_T180 <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T180,10) #Correlation delta pre vs T180

Intensity_saliva_delta_T24 <- Znorm(t(Delta_calculator(Intensity_saliva, subj_saliva, grps_saliva, "pre", "T24"))) 
colnames(Intensity_saliva_delta_T24) <- rownames(Intensity_saliva)
Corr_saliva_delta_T24 <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T24,10) #Correlation delta pre vs T24
  
  #Not imputed
Z_Intensity_saliva_NI <- t(Znorm(t(Intensity_saliva_NI)))
rownames(Z_Intensity_saliva_NI) <- rownames(Intensity_saliva)
Corr_saliva_pre_10_NI <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_saliva_NI)[which(grps_saliva == "pre"),], 10)

Intensity_saliva_delta_T0_NI <- Znorm(t(Delta_calculator(Intensity_saliva_NI, subj_saliva, grps_saliva, "pre", "T0"))) 
colnames(Intensity_saliva_delta_T0_NI) <- rownames(Intensity_saliva_NI)
Corr_saliva_delta_T0_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T0_NI,10) #Correlation delta pre vs T0

Intensity_saliva_delta_T30_NI <- Znorm(t(Delta_calculator(Intensity_saliva_NI, subj_saliva, grps_saliva, "pre", "T30"))) 
colnames(Intensity_saliva_delta_T30_NI) <- rownames(Intensity_saliva_NI)
Corr_saliva_delta_T30_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T30_NI,10) #Correlation delta pre vs T30

Intensity_saliva_delta_T60_NI <- Znorm(t(Delta_calculator(Intensity_saliva_NI, subj_saliva, grps_saliva, "pre", "T60"))) 
colnames(Intensity_saliva_delta_T60_NI) <- rownames(Intensity_saliva_NI)
Corr_saliva_delta_T60_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T60_NI,10) #Correlation delta pre vs T60

Intensity_saliva_delta_T180_NI <- Znorm(t(Delta_calculator(Intensity_saliva_NI, subj_saliva, grps_saliva, "pre", "T180"))) 
colnames(Intensity_saliva_delta_T180_NI) <- rownames(Intensity_saliva_NI)
Corr_saliva_delta_T180_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T180_NI,10) #Correlation delta pre vs T180

Intensity_saliva_delta_T24_NI <- Znorm(t(Delta_calculator(Intensity_saliva_NI, subj_saliva, grps_saliva, "pre", "T24"))) 
colnames(Intensity_saliva_delta_T24_NI) <- rownames(Intensity_saliva_NI)
Corr_saliva_delta_T24_10_NI <- Correlation_TwoMethod(Z_clinical, Intensity_saliva_delta_T24_NI,10) #Correlation delta pre vs T24

#Urine
  #Imputed
Z_Intensity_urine <- t(Znorm(t(Intensity_urine)))
rownames(Z_Intensity_urine) <- rownames(Intensity_urine)
Corr_urine_pre <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_urine)[which(grps_urine == "pre"),],10)

Intensity_urine_delta_T0 <- Znorm(t(Delta_calculator(Intensity_urine, subj_urine, grps_urine, "pre", "T0"))) 
colnames(Intensity_urine_delta_T0) <- rownames(Intensity_urine)
Corr_urine_delta_T0 <- Correlation_TwoMethod(Z_clinical, Intensity_urine_delta_T0,10) #Correlation delta pre vs T0

Intensity_urine_delta_T24 <- Znorm(t(Delta_calculator(Intensity_urine, subj_urine, grps_urine, "pre", "T24"))) 
colnames(Intensity_urine_delta_T24) <- rownames(Intensity_urine)
Corr_urine_delta_T24 <- Correlation_TwoMethod(Z_clinical, Intensity_urine_delta_T24,10) #Correlation delta pre vs T24

  #Not Imputed
Z_Intensity_urine_NI <- t(Znorm(t(Intensity_urine_NI)))
rownames(Z_Intensity_urine_NI) <- rownames(Intensity_urine)
Corr_urine_pre_6_NI <- Correlation_TwoMethod(Z_clinical, t(Z_Intensity_urine_NI)[which(grps_urine == "pre"),],6)

Intensity_urine_delta_T0_NI <- Znorm(t(Delta_calculator(Intensity_urine_NI, subj_urine, grps_urine, "pre", "T0"))) 
colnames(Intensity_urine_delta_T0_NI) <- rownames(Intensity_urine)
Corr_urine_delta_T0_6_NI <- Correlation_TwoMethod(Z_clinical, Intensity_urine_delta_T0_NI,6) #Correlation delta pre vs T0

Intensity_urine_delta_T24_NI <- Znorm(t(Delta_calculator(Intensity_urine_NI, subj_urine, grps_urine, "pre", "T24"))) 
colnames(Intensity_urine_delta_T24_NI) <- rownames(Intensity_urine)
Corr_urine_delta_T24_6_NI <- Correlation_TwoMethod(Z_clinical, Intensity_urine_delta_T24_NI,6) #Correlation delta pre vs T24

#### Plots correlations ####
  #Plasma
#T0
data_plot <- cbind(Clinical_data$`Alanine aminotransferase measurement (34608000)`,t(Delta_calculator(Intensity_plasma_NI,subj_plasma,grps_plasma,"pre","T0"))[,"PI16"])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot1 <- ggplot(data_plot, aes(V1,V2)) + 
            geom_point(aes(size = 4)) + 
            geom_smooth(method = lm) +
            theme_minimal() +
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
            labs(x = "Alanine aminotransferase measurement", y = "log2(PI16) (T0)", title = "Alanine Aminotransferase Measurement vs Plasma PI16 T0") +
            annotate("text", 17.5,5.5,label = (paste("r = 0.85")), fontface = "bold") +
            annotate("text", 17.5,5,label = (paste("pvalue = 0.047")), fontface = "bold")
plot1  
  #Saliva
#Pre
data_plot <- cbind(Clinical_data$`Alanine aminotransferase measurement (34608000)`,Intensity_saliva_NI["P48147",which(grps_saliva == "pre")])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot2 <- ggplot(data_plot, aes(V1,V2)) + 
            geom_point(aes(size = 4)) + 
            geom_smooth(method = lm) +
            theme_minimal() +
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
            labs(x = "Alanine aminotransferase measurement", y = "log2(P48147) (pre)", title = "Alanine Aminotransferase Measurement vs Saliva P48147 pre") +
            annotate("text", 17.5,4.5,label = (paste("r = -0.82")), fontface = "bold") +
            annotate("text", 17.5,4,label = (paste("pvalue = 0.0305")), fontface = "bold")
plot2

#T180
data_plot <- cbind(Clinical_data$`Alanine aminotransferase measurement (34608000)`,t(Delta_calculator(Intensity_saliva_NI,subj_saliva,grps_saliva,"pre","T180"))[,"P07954;P07954-2"])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot3 <- ggplot(data_plot, aes(V1,V2)) + 
            geom_point(aes(size = 4)) + 
            geom_smooth(method = lm) +
            theme_minimal() +
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
            labs(x = "Alanine aminotransferase measurement", y = "log2(P07954;P07954-2) (T180)", title = "Alanine Aminotransferase Measurement vs Saliva P07954;P07954-2 T180") +
            annotate("text", 17.5,4.5,label = (paste("r = 0.83")), fontface = "bold") +
            annotate("text", 17.5,4,label = (paste("pvalue = 0.0203")), fontface = "bold")

plot3

data_plot <- cbind(Clinical_data$`Alanine aminotransferase measurement (34608000)`,t(Delta_calculator(Intensity_saliva_NI,subj_saliva,grps_saliva,"pre","T180"))[,"P15104"])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot4 <- ggplot(data_plot, aes(V1,V2)) + 
            geom_point(aes(size = 4)) + 
            geom_smooth(method = lm) +
            theme_minimal() +
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
            labs(x = "Alanine aminotransferase measurement", y = "log2(P15104) (T180)", title = "Alanine Aminotransferase Measurement vs Saliva P15104 T180") +
            annotate("text", 17.5,4.5,label = (paste("r = 0.79")), fontface = "bold") +
            annotate("text", 17.5,4,label = (paste("pvalue = 0.048")), fontface = "bold")
plot4
  #Urine
#pre
data_plot <- cbind(Clinical_data$`Lean body mass (248362003)`,Intensity_urine_NI["P12273",which(grps_saliva == "pre")])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot5 <- ggplot(data_plot, aes(V1,V2)) + 
            geom_point(aes(size = 4)) + 
            geom_smooth(method = lm) +
            theme_minimal() +
            theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
            labs(x = "Lean Body Mass", y = "log2(P12273) (pre)", title = "Lean Body Mass vs Urine P12273 pre") +
            annotate("text", 43.5,4.5,label = (paste("r = 0.89")), fontface = "bold") +
            annotate("text", 43.5,4,label = (paste("pvalue = 0.0001")), fontface = "bold")

plot5

data_plot <- cbind(Clinical_data$`Lean body mass (248362003)`,Intensity_urine_NI["P01036",which(grps_saliva == "pre")])
data_plot <- na.omit(as.data.frame(data_plot))
data_plot$V1 <- as.numeric(data_plot$V1)
data_plot$V2 <- as.numeric(data_plot$V2)

plot6 <- ggplot(data_plot, aes(V1,V2)) + 
  geom_point(aes(size = 4)) + 
  geom_smooth(method = lm) +
  theme_minimal() +
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  labs(x = "Lean Body Mass", y = "log2(P01036) (pre)", title = "Lean Body Mass vs Urine P01036 pre") +
  annotate("text", 43.5,4.5,label = (paste("r = 0.89")), fontface = "bold") +
  annotate("text", 43.5,4,label = (paste("pvalue = 0.0008")), fontface = "bold")

plot6
ggarrange(plot1,plot2,plot3) #Export PDF 11x11
ggarrange(plot4,plot5,plot6) #Export PDF 11x11

#### Annotations ####
annotation_table <- read.table(gzfile("/Users/Gerard/Desktop/Annotations/mainAnnot.homo_sapiens.txt.gz"), header = T, sep = "\t", fill = T)
matrix_annotations_GOBP_plasma = data.frame()
matrix_annotations_GOMF_plasma = data.frame()
matrix_annotations_GOCC_plasma = data.frame()
matrix_annotations_KEGG_plasma = data.frame()


for(Protein.IDs in rownames(Intensity_plasma)){
  id = strsplit(x = Protein.IDs, split = ";")[[1]][1]
  row_uni = grep(id, annotation_table$UniProt)
  if(length(row_uni) != 0){
    row_to_annotate = annotation_table[row_uni,]
    
    GOBP <- strsplit(x = row_to_annotate$GOBP.name, split = ";")[[1]]
    tojoin <- merge(GOBP,id)
    matrix_annotations_GOBP_plasma = rbind(matrix_annotations_GOBP_plasma,tojoin)

    
    GOCC <- strsplit(x = row_to_annotate$GOCC.name, split = ";")[[1]]
    tojoin <- merge(GOCC,id)
    matrix_annotations_GOCC_plasma = rbind(matrix_annotations_GOCC_plasma,tojoin)
    
    GOMF <- strsplit(x = row_to_annotate$GOMF.name, split = ";")[[1]]
    tojoin <- merge(GOMF,id)
    matrix_annotations_GOMF_plasma = rbind(matrix_annotations_GOMF_plasma,tojoin)
    
    KEGG <- strsplit(x = row_to_annotate$KEGG.name, split = ";")[[1]]
    tojoin <- merge(KEGG,id)
    matrix_annotations_KEGG_plasma = rbind(matrix_annotations_KEGG_plasma,tojoin)
    
  }
}

matrix_annotations_GOBP_plasma <- aggregate(matrix_annotations_GOBP_plasma$y, list(matrix_annotations_GOBP_plasma$x), paste ,collapse=" ")

toGSEA_GOBP_plasma <- list()
for(i in seq_along(rownames(matrix_annotations_GOBP_plasma))){
  toGSEA_GOBP_plasma[[matrix_annotations_GOBP_plasma[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOBP_plasma$x[[i]], split = " "))
}
matrixranks_T0 <- as.vector(T0_Pre_plasma$logFC)
names(matrixranks_T0) <- rownames(T0_Pre_plasma)
fgsea_T0_pre_GOBP_plasma <- fgsea(pathways = toGSEA_GOBP_plasma,
                             stats    = matrixranks_T0)

matrixranks_T30 <- as.vector(T30_Pre_plasma$logFC)
names(matrixranks_T30) <- rownames(T30_Pre_plasma)
fgsea_T30_pre_GOBP_plasma <- fgsea(pathways = toGSEA_GOBP_plasma,
                             stats    = matrixranks_T30)

matrixranks_T60 <- as.vector(T60_Pre_plasma$logFC)
names(matrixranks_T60) <- rownames(T60_Pre_plasma)
fgsea_T60_pre_GOBP_plasma <- fgsea(pathways = toGSEA_GOBP_plasma,
                             stats    = matrixranks_T60)

matrixranks_T180 <- as.vector(T180_Pre_plasma$logFC)
names(matrixranks_T180) <- rownames(T180_Pre_plasma)
fgsea_T180_pre_GOBP_plasma <- fgsea(pathways = toGSEA_GOBP_plasma,
                             stats    = matrixranks_T180)

matrixranks_T24 <- as.vector(T24_Pre_plasma$logFC)
names(matrixranks_T24) <- rownames(T24_Pre_plasma)
fgsea_T24_pre_GOBP_plasma <- fgsea(pathways = toGSEA_GOBP_plasma,
                             stats    = matrixranks_T24)

matrix_annotations_GOCC_plasma <- aggregate(matrix_annotations_GOCC_plasma$y, list(matrix_annotations_GOCC_plasma$x), paste ,collapse=" ")
toGSEA_GOCC_plasma <- list()
for(i in seq_along(rownames(matrix_annotations_GOCC_plasma))){
  toGSEA_GOCC_plasma[[matrix_annotations_GOCC_plasma[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOCC_plasma$x[[i]], split = " "))
}

fgsea_T0_pre_GOCC_plasma <- fgsea(pathways = toGSEA_GOCC_plasma,
                             stats    = matrixranks_T0)

fgsea_T30_pre_GOCC_plasma <- fgsea(pathways = toGSEA_GOCC_plasma,
                              stats    = matrixranks_T30)

fgsea_T60_pre_GOCC_plasma <- fgsea(pathways = toGSEA_GOCC_plasma,
                              stats    = matrixranks_T60)

fgsea_T180_pre_GOCC_plasma <- fgsea(pathways = toGSEA_GOCC_plasma,
                               stats    = matrixranks_T180)

fgsea_T24_pre_GOCC_plasma <- fgsea(pathways = toGSEA_GOCC_plasma,
                              stats    = matrixranks_T24)

matrix_annotations_GOMF_plasma <- aggregate(matrix_annotations_GOMF_plasma$y, list(matrix_annotations_GOMF_plasma$x), paste ,collapse=" ")
toGSEA_GOMF_plasma <- list()
for(i in seq_along(rownames(matrix_annotations_GOMF_plasma))){
  toGSEA_GOMF_plasma[[matrix_annotations_GOMF_plasma[i,1]]] <- unlist(strsplit(x = matrix_annotations_GOMF_plasma$x[[i]], split = " "))
}

fgsea_T0_pre_GOMF_plasma <- fgsea(pathways = toGSEA_GOMF_plasma,
                                  stats    = matrixranks_T0)

fgsea_T30_pre_GOMF_plasma <- fgsea(pathways = toGSEA_GOMF_plasma,
                                   stats    = matrixranks_T30)

fgsea_T60_pre_GOMF_plasma <- fgsea(pathways = toGSEA_GOMF_plasma,
                                   stats    = matrixranks_T60)

fgsea_T180_pre_GOMF_plasma <- fgsea(pathways = toGSEA_GOMF_plasma,
                                    stats    = matrixranks_T180)

fgsea_T24_pre_GOMF_plasma <- fgsea(pathways = toGSEA_GOMF_plasma,
                                   stats    = matrixranks_T24)

matrix_annotations_KEGG_plasma <- aggregate(matrix_annotations_KEGG_plasma$y, list(matrix_annotations_KEGG_plasma$x), paste ,collapse=" ")
toGSEA_KEGG_plasma <- list()
for(i in seq_along(rownames(matrix_annotations_KEGG_plasma))){
  toGSEA_KEGG_plasma[[matrix_annotations_KEGG_plasma[i,1]]] <- unlist(strsplit(x = matrix_annotations_KEGG_plasma$x[[i]], split = " "))
}

fgsea_T0_pre_KEGG_plasma <- fgsea(pathways = toGSEA_KEGG_plasma,
                                  stats    = matrixranks_T0)

fgsea_T30_pre_KEGG_plasma <- fgsea(pathways = toGSEA_KEGG_plasma,
                                   stats    = matrixranks_T30)

fgsea_T60_pre_KEGG_plasma <- fgsea(pathways = toGSEA_KEGG_plasma,
                                   stats    = matrixranks_T60)

fgsea_T180_pre_KEGG_plasma <- fgsea(pathways = toGSEA_KEGG_plasma,
                                    stats    = matrixranks_T180)

fgsea_T24_pre_KEGG_plasma <- fgsea(pathways = toGSEA_KEGG_plasma,
                                   stats    = matrixranks_T24)

#### Area under the curve correlation ####

#Plasma
auc_plasma <- Calculate_AUC(Intensity_plasma_mediascale, as.factor(subj_plasma), grps_plasma)
na_auc <- is.na(auc_plasma[,1])
Z_clinical_red_plasma <- Z_clinical[!na_auc,]
auc_plasma <- na.omit(auc_plasma)

Zauc_plasma <- Znorm(auc_plasma)
colnames(Zauc_plasma) <- colnames(auc_plasma) 

corr_auc_plasma <- Correlation_TwoMethod(Z_clinical_red,Zauc_plasma,1)

#saliva
auc_saliva <- Calculate_AUC(Intensity_saliva_medianscal, as.factor(subj_saliva), grps_saliva)

Zauc_saliva <- Znorm(auc_saliva)
colnames(Zauc_saliva) <- colnames(auc_saliva) 

corr_auc_saliva <- Correlation_TwoMethod(Z_clinical,Zauc_saliva,1)

#### Cluster correlation Plasma ####
clusters_plasma <- read.csv("clusters_plasma.csv",sep = ",", header = F)
which(clusters_plasma$V1 == "NAME")

cluster1_plasma <- clusters_plasma[2:33,1]
cluster2_plasma <- clusters_plasma[35:40,1]
cluster3_plasma <- clusters_plasma[42:72,1]
cluster4_plasma <- clusters_plasma[74:111,1]

#Pre
Corr_plasma_pre_C1 <- Corr_plasma_pre[Corr_plasma_pre$Protein %in% cluster1_plasma,]
Corr_plasma_pre_C1$p.adjust <- p.adjust(Corr_plasma_pre_C1$pVal, method = "BH")
Corr_plasma_pre_C2 <- Corr_plasma_pre[Corr_plasma_pre$Protein %in% cluster2_plasma,]
Corr_plasma_pre_C2$p.adjust <- p.adjust(Corr_plasma_pre_C2$pVal, method = "BH")
Corr_plasma_pre_C3 <- Corr_plasma_pre[Corr_plasma_pre$Protein %in% cluster3_plasma,]
Corr_plasma_pre_C3$p.adjust <- p.adjust(Corr_plasma_pre_C3$pVal, method = "BH")
Corr_plasma_pre_C4 <- Corr_plasma_pre[Corr_plasma_pre$Protein %in% cluster4_plasma,]
Corr_plasma_pre_C4$p.adjust <- p.adjust(Corr_plasma_pre_C4$pVal, method = "BH")

#Delta T0
Corr_plasma_delta_T0_C1 <- Corr_plasma_delta_T0[Corr_plasma_delta_T0$Protein %in% cluster1_plasma,] 
Corr_plasma_delta_T0_C1$p.adjust <- p.adjust(Corr_plasma_delta_T0_C1$pVal, method = "BH")
Corr_plasma_delta_T0_C2 <- Corr_plasma_delta_T0[Corr_plasma_delta_T0$Protein %in% cluster2_plasma,] 
Corr_plasma_delta_T0_C2$p.adjust <- p.adjust(Corr_plasma_delta_T0_C2$pVal, method = "BH")
Corr_plasma_delta_T0_C3 <- Corr_plasma_delta_T0[Corr_plasma_delta_T0$Protein %in% cluster3_plasma,] 
Corr_plasma_delta_T0_C3$p.adjust <- p.adjust(Corr_plasma_delta_T0_C3$pVal, method = "BH")
Corr_plasma_delta_T0_C4 <- Corr_plasma_delta_T0[Corr_plasma_delta_T0$Protein %in% cluster4_plasma,] 
Corr_plasma_delta_T0_C4$p.adjust <- p.adjust(Corr_plasma_delta_T0_C4$pVal, method = "BH")

#Delta T30
Corr_plasma_delta_T30_C1 <- Corr_plasma_delta_T30[Corr_plasma_delta_T30$Protein %in% cluster1_plasma,]
Corr_plasma_delta_T30_C1$p.adjust <- p.adjust(Corr_plasma_delta_T30_C1$pVal, method = "BH")
Corr_plasma_delta_T30_C2 <- Corr_plasma_delta_T30[Corr_plasma_delta_T30$Protein %in% cluster2_plasma,]
Corr_plasma_delta_T30_C2$p.adjust <- p.adjust(Corr_plasma_delta_T30_C2$pVal, method = "BH")
Corr_plasma_delta_T30_C3 <- Corr_plasma_delta_T30[Corr_plasma_delta_T30$Protein %in% cluster3_plasma,]
Corr_plasma_delta_T30_C3$p.adjust <- p.adjust(Corr_plasma_delta_T30_C3$pVal, method = "BH")
Corr_plasma_delta_T30_C4 <- Corr_plasma_delta_T30[Corr_plasma_delta_T30$Protein %in% cluster4_plasma,]
Corr_plasma_delta_T30_C4$p.adjust <- p.adjust(Corr_plasma_delta_T30_C4$pVal, method = "BH")

#Delta T60
Corr_plasma_delta_T60_C1 <- Corr_plasma_delta_T60[Corr_plasma_delta_T60$Protein %in% cluster1_plasma,]
Corr_plasma_delta_T60_C1$p.adjust <- p.adjust(Corr_plasma_delta_T60_C1$pVal, method = "BH")
Corr_plasma_delta_T60_C2 <- Corr_plasma_delta_T60[Corr_plasma_delta_T60$Protein %in% cluster2_plasma,]
Corr_plasma_delta_T60_C2$p.adjust <- p.adjust(Corr_plasma_delta_T60_C2$pVal, method = "BH")
Corr_plasma_delta_T60_C3 <- Corr_plasma_delta_T60[Corr_plasma_delta_T60$Protein %in% cluster3_plasma,]
Corr_plasma_delta_T60_C3$p.adjust <- p.adjust(Corr_plasma_delta_T60_C3$pVal, method = "BH")
Corr_plasma_delta_T60_C4 <- Corr_plasma_delta_T60[Corr_plasma_delta_T60$Protein %in% cluster4_plasma,]
Corr_plasma_delta_T60_C4$p.adjust <- p.adjust(Corr_plasma_delta_T60_C4$pVal, method = "BH")

#Delta T180
Corr_plasma_delta_T180_C1 <- Corr_plasma_delta_T180[Corr_plasma_delta_T180$Protein %in% cluster1_plasma,]
Corr_plasma_delta_T180_C1$p.adjust <- p.adjust(Corr_plasma_delta_T180_C1$pVal, method = "BH")
Corr_plasma_delta_T180_C2 <- Corr_plasma_delta_T180[Corr_plasma_delta_T180$Protein %in% cluster2_plasma,]
Corr_plasma_delta_T180_C2$p.adjust <- p.adjust(Corr_plasma_delta_T180_C2$pVal, method = "BH")
Corr_plasma_delta_T180_C3 <- Corr_plasma_delta_T180[Corr_plasma_delta_T180$Protein %in% cluster3_plasma,]
Corr_plasma_delta_T180_C3$p.adjust <- p.adjust(Corr_plasma_delta_T180_C3$pVal, method = "BH")
Corr_plasma_delta_T180_C4 <- Corr_plasma_delta_T180[Corr_plasma_delta_T180$Protein %in% cluster4_plasma,]
Corr_plasma_delta_T180_C4$p.adjust <- p.adjust(Corr_plasma_delta_T180_C4$pVal, method = "BH")

#Delta T24
Corr_plasma_delta_T24_C1 <- Corr_plasma_delta_T24[Corr_plasma_delta_T24$Protein %in% cluster1_plasma,]
Corr_plasma_delta_T24_C1$p.adjust <- p.adjust(Corr_plasma_delta_T24_C1$pVal, method = "BH")
Corr_plasma_delta_T24_C2 <- Corr_plasma_delta_T24[Corr_plasma_delta_T24$Protein %in% cluster2_plasma,]
Corr_plasma_delta_T24_C2$p.adjust <- p.adjust(Corr_plasma_delta_T24_C2$pVal, method = "BH")
Corr_plasma_delta_T24_C3 <- Corr_plasma_delta_T24[Corr_plasma_delta_T24$Protein %in% cluster3_plasma,]
Corr_plasma_delta_T24_C3$p.adjust <- p.adjust(Corr_plasma_delta_T24_C3$pVal, method = "BH")
Corr_plasma_delta_T24_C4 <- Corr_plasma_delta_T24[Corr_plasma_delta_T24$Protein %in% cluster4_plasma,]
Corr_plasma_delta_T24_C4$p.adjust <- p.adjust(Corr_plasma_delta_T24_C4$pVal, method = "BH")

#### Cluster correlation Saliva ####
clusters_saliva <- read.csv("clusters_saliva.csv",sep = ",", header = F)
which(clusters_saliva$V1 == "NAME")

cluster1_saliva <- clusters_saliva[2:80,1]
cluster2_saliva <- clusters_saliva[82:110,1]
cluster3_saliva <- clusters_saliva[112:183,1]
cluster4_saliva <- clusters_saliva[185:368,1]

#Pre
Corr_saliva_pre_C1 <- Corr_saliva_pre[Corr_saliva_pre$Protein %in% cluster1_saliva,]
Corr_saliva_pre_C1$p.adjust <- p.adjust(Corr_saliva_pre_C1$pVal, method = "BH")
Corr_saliva_pre_C2 <- Corr_saliva_pre[Corr_saliva_pre$Protein %in% cluster2_saliva,]
Corr_saliva_pre_C2$p.adjust <- p.adjust(Corr_saliva_pre_C2$pVal, method = "BH")
Corr_saliva_pre_C3 <- Corr_saliva_pre[Corr_saliva_pre$Protein %in% cluster3_saliva,]
Corr_saliva_pre_C3$p.adjust <- p.adjust(Corr_saliva_pre_C3$pVal, method = "BH")
Corr_saliva_pre_C4 <- Corr_saliva_pre[Corr_saliva_pre$Protein %in% cluster4_saliva,]
Corr_saliva_pre_C4$p.adjust <- p.adjust(Corr_saliva_pre_C4$pVal, method = "BH")

#Delta T0
Corr_saliva_delta_T0_C1 <- Corr_saliva_delta_T0[Corr_saliva_delta_T0$Protein %in% cluster1_saliva,] 
Corr_saliva_delta_T0_C1$p.adjust <- p.adjust(Corr_saliva_delta_T0_C1$pVal, method = "BH")
Corr_saliva_delta_T0_C2 <- Corr_saliva_delta_T0[Corr_saliva_delta_T0$Protein %in% cluster2_saliva,] 
Corr_saliva_delta_T0_C2$p.adjust <- p.adjust(Corr_saliva_delta_T0_C2$pVal, method = "BH")
Corr_saliva_delta_T0_C3 <- Corr_saliva_delta_T0[Corr_saliva_delta_T0$Protein %in% cluster3_saliva,] 
Corr_saliva_delta_T0_C3$p.adjust <- p.adjust(Corr_saliva_delta_T0_C3$pVal, method = "BH")
Corr_saliva_delta_T0_C4 <- Corr_saliva_delta_T0[Corr_saliva_delta_T0$Protein %in% cluster4_saliva,] 
Corr_saliva_delta_T0_C4$p.adjust <- p.adjust(Corr_saliva_delta_T0_C4$pVal, method = "BH")

#Delta T30
Corr_saliva_delta_T30_C1 <- Corr_saliva_delta_T30[Corr_saliva_delta_T30$Protein %in% cluster1_saliva,]
Corr_saliva_delta_T30_C1$p.adjust <- p.adjust(Corr_saliva_delta_T30_C1$pVal, method = "BH")
Corr_saliva_delta_T30_C2 <- Corr_saliva_delta_T30[Corr_saliva_delta_T30$Protein %in% cluster2_saliva,]
Corr_saliva_delta_T30_C2$p.adjust <- p.adjust(Corr_saliva_delta_T30_C2$pVal, method = "BH")
Corr_saliva_delta_T30_C3 <- Corr_saliva_delta_T30[Corr_saliva_delta_T30$Protein %in% cluster3_saliva,]
Corr_saliva_delta_T30_C3$p.adjust <- p.adjust(Corr_saliva_delta_T30_C3$pVal, method = "BH")
Corr_saliva_delta_T30_C4 <- Corr_saliva_delta_T30[Corr_saliva_delta_T30$Protein %in% cluster4_saliva,]
Corr_saliva_delta_T30_C4$p.adjust <- p.adjust(Corr_saliva_delta_T30_C4$pVal, method = "BH")

#Delta T60
Corr_saliva_delta_T60_C1 <- Corr_saliva_delta_T60[Corr_saliva_delta_T60$Protein %in% cluster1_saliva,]
Corr_saliva_delta_T60_C1$p.adjust <- p.adjust(Corr_saliva_delta_T60_C1$pVal, method = "BH")
Corr_saliva_delta_T60_C2 <- Corr_saliva_delta_T60[Corr_saliva_delta_T60$Protein %in% cluster2_saliva,]
Corr_saliva_delta_T60_C2$p.adjust <- p.adjust(Corr_saliva_delta_T60_C2$pVal, method = "BH")
Corr_saliva_delta_T60_C3 <- Corr_saliva_delta_T60[Corr_saliva_delta_T60$Protein %in% cluster3_saliva,]
Corr_saliva_delta_T60_C3$p.adjust <- p.adjust(Corr_saliva_delta_T60_C3$pVal, method = "BH")
Corr_saliva_delta_T60_C4 <- Corr_saliva_delta_T60[Corr_saliva_delta_T60$Protein %in% cluster4_saliva,]
Corr_saliva_delta_T60_C4$p.adjust <- p.adjust(Corr_saliva_delta_T60_C4$pVal, method = "BH")

#Delta T180
Corr_saliva_delta_T180_C1 <- Corr_saliva_delta_T180[Corr_saliva_delta_T180$Protein %in% cluster1_saliva,]
Corr_saliva_delta_T180_C1$p.adjust <- p.adjust(Corr_saliva_delta_T180_C1$pVal, method = "BH")
Corr_saliva_delta_T180_C2 <- Corr_saliva_delta_T180[Corr_saliva_delta_T180$Protein %in% cluster2_saliva,]
Corr_saliva_delta_T180_C2$p.adjust <- p.adjust(Corr_saliva_delta_T180_C2$pVal, method = "BH")
Corr_saliva_delta_T180_C3 <- Corr_saliva_delta_T180[Corr_saliva_delta_T180$Protein %in% cluster3_saliva,]
Corr_saliva_delta_T180_C3$p.adjust <- p.adjust(Corr_saliva_delta_T180_C3$pVal, method = "BH")
Corr_saliva_delta_T180_C4 <- Corr_saliva_delta_T180[Corr_saliva_delta_T180$Protein %in% cluster4_saliva,]
Corr_saliva_delta_T180_C4$p.adjust <- p.adjust(Corr_saliva_delta_T180_C4$pVal, method = "BH")

#Delta T24
Corr_saliva_delta_T24_C1 <- Corr_saliva_delta_T24[Corr_saliva_delta_T24$Protein %in% cluster1_saliva,]
Corr_saliva_delta_T24_C1$p.adjust <- p.adjust(Corr_saliva_delta_T24_C1$pVal, method = "BH")
Corr_saliva_delta_T24_C2 <- Corr_saliva_delta_T24[Corr_saliva_delta_T24$Protein %in% cluster2_saliva,]
Corr_saliva_delta_T24_C2$p.adjust <- p.adjust(Corr_saliva_delta_T24_C2$pVal, method = "BH")
Corr_saliva_delta_T24_C3 <- Corr_saliva_delta_T24[Corr_saliva_delta_T24$Protein %in% cluster3_saliva,]
Corr_saliva_delta_T24_C3$p.adjust <- p.adjust(Corr_saliva_delta_T24_C3$pVal, method = "BH")
Corr_saliva_delta_T24_C4 <- Corr_saliva_delta_T24[Corr_saliva_delta_T24$Protein %in% cluster4_saliva,]
Corr_saliva_delta_T24_C4$p.adjust <- p.adjust(Corr_saliva_delta_T24_C4$pVal, method = "BH")

#### Cluster GSEA Plasma ####
  ##GOBP##
#GOBP T0
fgsea_T0_pre_GOBP_plasma_C1 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster1_plasma]))
fgsea_T0_pre_GOBP_plasma_C2 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster2_plasma]))
fgsea_T0_pre_GOBP_plasma_C3 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster3_plasma]))
fgsea_T0_pre_GOBP_plasma_C4 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster4_plasma]))

#GOBP T30
fgsea_T30_pre_GOBP_plasma_C1 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T30[cluster1_plasma]))
fgsea_T30_pre_GOBP_plasma_C2 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T30[cluster2_plasma]))
fgsea_T30_pre_GOBP_plasma_C3 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T30[cluster3_plasma]))
fgsea_T30_pre_GOBP_plasma_C4 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T30[cluster4_plasma]))

#GOBP T60
fgsea_T60_pre_GOBP_plasma_C1 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T60[cluster1_plasma]))
fgsea_T60_pre_GOBP_plasma_C2 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T60[cluster2_plasma]))
fgsea_T60_pre_GOBP_plasma_C3 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T60[cluster3_plasma]))
fgsea_T60_pre_GOBP_plasma_C4 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T60[cluster4_plasma]))

#GOBP T180
fgsea_T180_pre_GOBP_plasma_C1 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T180[cluster1_plasma]))
fgsea_T180_pre_GOBP_plasma_C2 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T180[cluster2_plasma]))
fgsea_T180_pre_GOBP_plasma_C3 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T180[cluster3_plasma]))
fgsea_T180_pre_GOBP_plasma_C4 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T180[cluster4_plasma]))
#GOBP T24
fgsea_T24_pre_GOBP_plasma_C1 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T24[cluster1_plasma]))
fgsea_T24_pre_GOBP_plasma_C2 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T24[cluster2_plasma]))
fgsea_T24_pre_GOBP_plasma_C3 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T24[cluster3_plasma]))
fgsea_T24_pre_GOBP_plasma_C4 <- fgsea(pathways = toGSEA_GOBP_plasma,
                                     stats    = na.omit(matrixranks_T24[cluster4_plasma]))

  ##GOMF##
#GOMF T0
fgsea_T0_pre_GOMF_plasma_C1 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster1_plasma]))
fgsea_T0_pre_GOMF_plasma_C2 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster2_plasma]))
fgsea_T0_pre_GOMF_plasma_C3 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster3_plasma]))
fgsea_T0_pre_GOMF_plasma_C4 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster4_plasma]))

#GOMF T30
fgsea_T30_pre_GOMF_plasma_C1 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster1_plasma]))
fgsea_T30_pre_GOMF_plasma_C2 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster2_plasma]))
fgsea_T30_pre_GOMF_plasma_C3 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster3_plasma]))
fgsea_T30_pre_GOMF_plasma_C4 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster4_plasma]))

#GOMF T60
fgsea_T60_pre_GOMF_plasma_C1 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster1_plasma]))
fgsea_T60_pre_GOMF_plasma_C2 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster2_plasma]))
fgsea_T60_pre_GOMF_plasma_C3 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster3_plasma]))
fgsea_T60_pre_GOMF_plasma_C4 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster4_plasma]))

#GOMF T180
fgsea_T180_pre_GOMF_plasma_C1 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster1_plasma]))
fgsea_T180_pre_GOMF_plasma_C2 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster2_plasma]))
fgsea_T180_pre_GOMF_plasma_C3 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster3_plasma]))
fgsea_T180_pre_GOMF_plasma_C4 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster4_plasma]))
#GOMF T24
fgsea_T24_pre_GOMF_plasma_C1 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster1_plasma]))
fgsea_T24_pre_GOMF_plasma_C2 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster2_plasma]))
fgsea_T24_pre_GOMF_plasma_C3 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster3_plasma]))
fgsea_T24_pre_GOMF_plasma_C4 <- fgsea(pathways = toGSEA_GOMF_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster4_plasma]))

  ##GOCC##
#GOCC T0
fgsea_T0_pre_GOCC_plasma_C1 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster1_plasma]))
fgsea_T0_pre_GOCC_plasma_C2 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster2_plasma]))
fgsea_T0_pre_GOCC_plasma_C3 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster3_plasma]))
fgsea_T0_pre_GOCC_plasma_C4 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster4_plasma]))

#GOCC T30
fgsea_T30_pre_GOCC_plasma_C1 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster1_plasma]))
fgsea_T30_pre_GOCC_plasma_C2 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster2_plasma]))
fgsea_T30_pre_GOCC_plasma_C3 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster3_plasma]))
fgsea_T30_pre_GOCC_plasma_C4 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster4_plasma]))

#GOCC T60
fgsea_T60_pre_GOCC_plasma_C1 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster1_plasma]))
fgsea_T60_pre_GOCC_plasma_C2 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster2_plasma]))
fgsea_T60_pre_GOCC_plasma_C3 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster3_plasma]))
fgsea_T60_pre_GOCC_plasma_C4 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster4_plasma]))

#GOCC T180
fgsea_T180_pre_GOCC_plasma_C1 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster1_plasma]))
fgsea_T180_pre_GOCC_plasma_C2 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster2_plasma]))
fgsea_T180_pre_GOCC_plasma_C3 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster3_plasma]))
fgsea_T180_pre_GOCC_plasma_C4 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster4_plasma]))
#GOCC T24
fgsea_T24_pre_GOCC_plasma_C1 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster1_plasma]))
fgsea_T24_pre_GOCC_plasma_C2 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster2_plasma]))
fgsea_T24_pre_GOCC_plasma_C3 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster3_plasma]))
fgsea_T24_pre_GOCC_plasma_C4 <- fgsea(pathways = toGSEA_GOCC_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster4_plasma]))

  ## KEGG ##
#KEGG T0
fgsea_T0_pre_KEGG_plasma_C1 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster1_plasma]))
fgsea_T0_pre_KEGG_plasma_C2 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster2_plasma]))
fgsea_T0_pre_KEGG_plasma_C3 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster3_plasma]))
fgsea_T0_pre_KEGG_plasma_C4 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                     stats    = na.omit(matrixranks_T0[cluster4_plasma]))

#KEGG T30
fgsea_T30_pre_KEGG_plasma_C1 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster1_plasma]))
fgsea_T30_pre_KEGG_plasma_C2 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster2_plasma]))
fgsea_T30_pre_KEGG_plasma_C3 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster3_plasma]))
fgsea_T30_pre_KEGG_plasma_C4 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T30[cluster4_plasma]))

#KEGG T60
fgsea_T60_pre_KEGG_plasma_C1 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster1_plasma]))
fgsea_T60_pre_KEGG_plasma_C2 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster2_plasma]))
fgsea_T60_pre_KEGG_plasma_C3 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster3_plasma]))
fgsea_T60_pre_KEGG_plasma_C4 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T60[cluster4_plasma]))

#KEGG T180
fgsea_T180_pre_KEGG_plasma_C1 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster1_plasma]))
fgsea_T180_pre_KEGG_plasma_C2 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster2_plasma]))
fgsea_T180_pre_KEGG_plasma_C3 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster3_plasma]))
fgsea_T180_pre_KEGG_plasma_C4 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                       stats    = na.omit(matrixranks_T180[cluster4_plasma]))
#KEGG T24
fgsea_T24_pre_KEGG_plasma_C1 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster1_plasma]))
fgsea_T24_pre_KEGG_plasma_C2 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster2_plasma]))
fgsea_T24_pre_KEGG_plasma_C3 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster3_plasma]))
fgsea_T24_pre_KEGG_plasma_C4 <- fgsea(pathways = toGSEA_KEGG_plasma,
                                      stats    = na.omit(matrixranks_T24[cluster4_plasma]))

#### Amylase correlation on saliva ####
Corr_saliva_amy_pre_NI <- Corr_saliva_pre_10_NI[which(Corr_saliva_pre_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_pre_NI$p.adjust <- p.adjust(Corr_saliva_amy_pre_NI$pVal, method = "BH")

Corr_saliva_amy_T0_NI <- Corr_saliva_delta_T0_10_NI[which(Corr_saliva_delta_T0_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_T0_NI$p.adjust <- p.adjust(Corr_saliva_amy_T0_NI$pVal, method = "BH")

Corr_saliva_amy_T30_NI <- Corr_saliva_delta_T30_10_NI[which(Corr_saliva_delta_T30_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_T30_NI$p.adjust <- p.adjust(Corr_saliva_amy_T30_NI$pVal, method = "BH")

Corr_saliva_amy_T60_NI <- Corr_saliva_delta_T60_10_NI[which(Corr_saliva_delta_T60_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_T60_NI$p.adjust <- p.adjust(Corr_saliva_amy_T60_NI$pVal, method = "BH")

Corr_saliva_amy_T180_NI <- Corr_saliva_delta_T180_10_NI[which(Corr_saliva_delta_T180_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_T180_NI$p.adjust <- p.adjust(Corr_saliva_amy_T180_NI$pVal, method = "BH")

Corr_saliva_amy_T24_NI <- Corr_saliva_delta_T24_10_NI[which(Corr_saliva_delta_T24_10_NI$Clinical == "Amylase measurement (64435009)"),]
Corr_saliva_amy_T24_NI$p.adjust <- p.adjust(Corr_saliva_amy_T24_NI$pVal, method = "BH")

Corr_saliva_fsgl <-  Corr_saliva_pre_10_NI[which(Corr_saliva_pre_10_NI$Clinical == "Glucose measurement (36048009)"),]
Corr_saliva_fsgl$p.adjust <- p.adjust(Corr_saliva_fsgl$pVal, method = "BH")

Corr_saliva_fsin <-  Corr_saliva_pre_10_NI[which(Corr_saliva_pre_10_NI$Clinical == "Insulin measurement (16890009)"),]
Corr_saliva_fsin$p.adjust <- p.adjust(Corr_saliva_fsin$pVal, method = "BH")

Corr_saliva_fshb <-  Corr_saliva_pre_10_NI[which(Corr_saliva_pre_10_NI$Clinical == "Hemoglobin A1c measurement (43396009)"),]
Corr_saliva_fsihb$p.adjust <- p.adjust(Corr_saliva_fshb$pVal, method = "BH")

#### Padjustment with smaller clinical data ####
padjust_small <- function(corr_dataframe){
  Index_INR <- which(corr_dataframe$Clinical == "INR??")
  Index_cal <- which(corr_dataframe$Clinical == "Calcium measurement (71878006)")
  Index_sod <- which(corr_dataframe$Clinical == "Sodium measurement (25197003)")
  Index_pot <- which(corr_dataframe$Clinical == "Potassium measurement (	59573005)")
  All_index <- c(Index_INR,Index_cal,Index_sod,Index_pot)
  small_corr <- corr_dataframe[-All_index,]
  small_corr$p.adjust <- p.adjust(small_corr$pVal, method = "BH")
  return(small_corr)
}

small_corr_pre_plasma <- padjust_small(Corr_plasma_pre_10_NI)
small_corr_T0_plasma <- padjust_small(Corr_plasma_delta_T0_10_NI)
small_corr_T30_plasma <- padjust_small(Corr_plasma_delta_T30_10_NI)
small_corr_T60_plasma <- padjust_small(Corr_plasma_delta_T60_10_NI)
small_corr_T180_plasma <- padjust_small(Corr_plasma_delta_T180_10_NI)
small_corr_T24_plasma <- padjust_small(Corr_plasma_delta_T24_10_NI)

small_corr_pre_saliva <- padjust_small(Corr_saliva_pre_10_NI)
small_corr_T0_saliva <- padjust_small(Corr_saliva_delta_T0_10_NI)
small_corr_T30_saliva <- padjust_small(Corr_saliva_delta_T30_10_NI)
small_corr_T60_saliva <- padjust_small(Corr_saliva_delta_T60_10_NI)
small_corr_T180_saliva <- padjust_small(Corr_saliva_delta_T180_10_NI)
small_corr_T24_saliva <- padjust_small(Corr_saliva_delta_T24_10_NI)
