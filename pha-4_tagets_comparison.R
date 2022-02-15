################
## Checking pha-4 targets expression in MSxaaX vs MSxapX vs MSxpX
################

library(dplyr)


source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)


setwd(source.path)

resDir = paste0("results/")
if(!dir.exists("results/")){dir.create("results/")}


library(Seurat)
library(ggplot2)
library(dplyr)

ms <- readRDS('~/R_projects/BWMs/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_cleanedBWM_and_Pharynx_iteration_GR15.rds')

Idents(ms) <- ms$manual.annot.ids
all_ids <- c(levels(Idents(ms)))[-74]
ms_new <- subset(ms, idents = all_ids)
ms_new@active.ident <- factor(x = ms_new$manual.annot.ids, levels =  all_ids)

ms_new <- subset(ms, idents =c("MSxpp","MSxapa","MSxap", "MSxp",  "MSxpa", "MSxaa",
                               "MSxaaa","MSxa","MSaaaaaa","MSx","MSxppa","MSxapap",
                               "MSxpaa","MSxppaa","MSaaaaap","MSaaapp", "MSxppp","MSxpapp",
                               "MSxpaaa", "MSxapp","MSaaaap", "MSxppppx","MSxpapa", "MSpaapaa",                    
                               "MSxppap", "MSxpppa", "MSxaap","MSaaapaa","MSxpppp", "MSpappa",
                               "MSxapaa", "MSxappp", "MSxaapapx","MSxpap","MSaappa", "MSaaaaaaa",                 
                               "MSxaapa", "MSpaaaa", "MSxapapa.MSxapapaa","MSpaaap", "MSxapppa","MSxpaaaa",
                               "MSaaaapp","MSxapaap","MSxpapap","MSxapaaa","MSapaapp","MSxpapaa",
                               "MSaaaaapa","MSpaaaaa","MSpaaapa","MSaaaappp/MSxapaapp","MSxppapp","MSxpppap",
                               "MSxppppa","MSxapppp","MSxppaaa","MSxaapap","MSxpappp","MSpaaaap",
                               "MSpaaapp","MSxapapp","MSxpppaa","MSxapppax","MSxppppp","MSxpaaap",
                               "MSaaaapa","MSxppaap","MSxappppx","MSxapappp.like","MSxapaapa","MSxpappa",
                               "MSaaappp","MSpaaappp/MSxapappa","MSxpaap.MSppaapp.like","MSxpppaa_MSxppppa_later_like"))



#####################
# Load the data

## read pha-4 target list FROM:
# Genome-Wide Identification of Binding Sites Defines Distinct Functions for Caenorhabditis elegans PHA-4/FOXA in Development and Environmental Response
# Mei Zhong et al
# Table 1

pha4_targets <- read.delim('./pha-4_genes_list_embryos_snyder_plos_gen_2010.txt',header = FALSE)
pha4_targets <- pha4_targets[,1]


# Load the list of TFs
# the list of genes is from https://genomebiology.biomedcentral.com/articles/10.1186/gb-2005-6-13-r110#MOESM1 
# A compendium of Caenorhabditis elegans regulatory transcription factors: a resource for mapping transcription regulatory networks
tflist <- (read.delim('./all_TFs_celegans.csv', header = FALSE))[,1]


###############
# just sub-setting EARLY MS data, splitting into three groups
ms_new_MSxaax <- subset(ms_new, idents = c('MSxaaa', 'MSxaap', 'MSxapa'))
ms_new_MSxaax$ident <- ms_new_MSxaax@active.ident

ms_new_MSxapp <- subset(ms_new, idents = c('MSxapp'))
ms_new_MSxapp$ident <- ms_new_MSxapp@active.ident

ms_new_MSxpxx <- subset(ms_new, idents = c('MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp'))
ms_new_MSxpxx$ident <- ms_new_MSxpxx@active.ident


# not used later
#ms_new_MSxaax_expr_data <- as.data.frame(AverageExpression(ms_new_MSxaax))
#ms_new_MSxaax_expr_data$mean <- apply(ms_new_MSxaax_expr_data,1, mean)
#ms_new_MSxaax_expr_data


# pha4_targets
# # trying to range genes by they mean expression value and leave only the ones that are present in the target list. Not really using this later.
# ms_new_MSxaax_expr_data <- ms_new_MSxaax_expr_data %>% filter(row.names(ms_new_MSxaax_expr_data) %in% pha4_targets)
# ms_new_MSxaax_expr_data <- arrange(ms_new_MSxaax_expr_data, desc(mean))
# ms_new_MSxaax_expr_data[1:1000,]
# 
# DoHeatmap(ms_new_MSxaax, features = rownames(ms_new_MSxaax_expr_data)[1:50],slot = 'data')



#########################################################
# can try to discard the genes that are expressed in less than 50% of cells per cluster

# For MSxaax_early
MSxaax_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxaax))) {
  ms_tmp<- subset(ms_new_MSxaax, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxaax_gene_name_space <- c(MSxaax_gene_name_space,tmp_list)
}
MSxaax_gene_name_space

# For MSxapp_early
MSxapp_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxapp))) {
  ms_tmp<- subset(ms_new_MSxapp, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxapp_gene_name_space <- c(MSxapp_gene_name_space,tmp_list)
}
MSxapp_gene_name_space

# For MSxpxx_early
MSxpxx_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxpxx))) {
  ms_tmp<- subset(ms_new_MSxpxx, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxpxx_gene_name_space <- c(MSxpxx_gene_name_space,tmp_list)
}
MSxpxx_gene_name_space

# identify genes that are pha-4 targets in all 3 datasets:


# Take all unique genes from each cluster per cohort and select only pha-4 targets
write(intersect(unique(unlist(MSxaax_gene_name_space)),pha4_targets), './MSxaax_pha4_targets.txt')
write(intersect(unique(unlist(MSxapp_gene_name_space)),pha4_targets), './MSxapp_pha4_targets.txt')
write(intersect(unique(unlist(MSxpxx_gene_name_space)),pha4_targets), './MSxpxx_pha4_targets.txt')

write(intersect(intersect(unique(unlist(MSxaax_gene_name_space)),pha4_targets), tflist), './MSxaax_pha4_targets_TFs.txt')
write(intersect(intersect(unique(unlist(MSxapp_gene_name_space)),pha4_targets), tflist), './MSxapp_pha4_targets_TFs.txt')
write(intersect(intersect(unique(unlist(MSxpxx_gene_name_space)),pha4_targets), tflist), './MSxpxx_pha4_targets_TFs.txt')



write(intersect(unique(unlist(MSxaax_gene_name_space))), './MSxaax_all.txt')
write(intersect(unique(unlist(MSxapp_gene_name_space))), './MSxapp_all.txt')
write(intersect(unique(unlist(MSxpxx_gene_name_space))), './MSxpxx_all.txt')
####################################################################################################################################
# Subsetting LATE MSx data
ms_new_MSxaxxx_late <- subset(ms_new, idents = c('MSaaapaa','MSxaapapx','MSaaaappp/MSxapaapp','MSxapaaa','MSpaaappp/MSxapappa'))
ms_new_MSxaxxx_late$ident <- ms_new_MSxaxxx_late@active.ident

ms_new_MSxapp_late <- subset(ms_new, idents = c('MSxappppx','MSxapppax'))
ms_new_MSxapp_late$ident <- ms_new_MSxapp_late@active.ident

ms_new_MSxpxx_late <- subset(ms_new, idents = c('MSxppppp', 'MSxppppa', 'MSxpppap', 'MSxpppaa',
                                                'MSxppapp', 'MSxpappp', 'MSxpappa','MSxpapap',
                                                'MSxpaaap'))
ms_new_MSxpxx_late$ident <- ms_new_MSxpxx_late@active.ident

######
# For MSxaxxx_late
MSxaxxx_late_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxaxxx_late))) {
  ms_tmp<- subset(ms_new_MSxaxxx_late, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxaxxx_late_gene_name_space <- c(MSxaxxx_late_gene_name_space,tmp_list)
}
MSxaxxx_late_gene_name_space

########
# For ms_new_MSxapp_late
MSxapp_late_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxapp_late))) {
  ms_tmp<- subset(ms_new_MSxapp_late, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxapp_late_gene_name_space <- c(MSxapp_late_gene_name_space,tmp_list)
}
MSxapp_late_gene_name_space

########
# For ms_new_MSxpxx_late
MSxpxx_late_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxpxx_late))) {
  ms_tmp<- subset(ms_new_MSxpxx_late, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxpxx_late_gene_name_space <- c(MSxpxx_late_gene_name_space,tmp_list)
}
MSxpxx_late_gene_name_space

# identify genes that are pha-4 targets in all 3 datasets:


# Take all unique genes from each cluster per cohort and select only pha-4 targets

write(intersect(unique(unlist(MSxaxxx_late_gene_name_space)),pha4_targets), './MSxaax_late_pha4_targets.txt')
write(intersect(unique(unlist(MSxapp_late_gene_name_space)),pha4_targets), './MSxapp_late_pha4_targets.txt')
write(intersect(unique(unlist(MSxpxx_late_gene_name_space)),pha4_targets), './MSxpxx_late_pha4_targets.txt')

write(intersect(intersect(unique(unlist(MSxaxxx_late_gene_name_space)),pha4_targets), tflist), './MSxaax_pha4_late_targets_TFs.txt')
write(intersect(intersect(unique(unlist(MSxapp_late_gene_name_space)),pha4_targets), tflist), './MSxapp_pha4_late_targets_TFs.txt')
write(intersect(intersect(unique(unlist(MSxpxx_late_gene_name_space)),pha4_targets), tflist), './MSxpxx_pha4_late_targets_TFs.txt')
