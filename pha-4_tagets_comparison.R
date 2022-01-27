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
## read pha-4 target list FROM:
# Genome-Wide Identification of Binding Sites Defines Distinct Functions for Caenorhabditis elegans PHA-4/FOXA in Development and Environmental Response
# Mei Zhong et al
# Table 1

pha4_targets <- read.delim('./pha-4_genes_list_embryos_snyder_plos_gen_2010.txt',header = FALSE)
pha4_targets <- pha4_targets[,1]


ms_new_MSxaax <- subset(ms_new, idents = c('MSxaaa', 'MSxaap', 'MSxapa'))
ms_new_MSxaax$ident <- ms_new_MSxaax@active.ident

ms_new_MSxapp <- subset(ms_new, idents = c('MSxapp'))
ms_new_MSxapp$ident <- ms_new_MSxapp@active.ident

ms_new_MSxpxx <- subset(ms_new, idents = c('MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp'))
ms_new_MSxpxx$ident <- ms_new_MSxpxx@active.ident


ms_new_MSxaax_expr_data <- as.data.frame(AverageExpression(ms_new_MSxaax))
ms_new_MSxaax_expr_data$mean <- apply(ms_new_MSxaax_expr_data,1, mean)
ms_new_MSxaax_expr_data

######



pha4_targets

ms_new_MSxaax_expr_data <- ms_new_MSxaax_expr_data %>% filter(row.names(ms_new_MSxaax_expr_data) %in% pha4_targets)


ms_new_MSxaax_expr_data <- arrange(ms_new_MSxaax_expr_data, desc(mean))
ms_new_MSxaax_expr_data[1:1000,]

DoHeatmap(ms_new_MSxaax, features = rownames(ms_new_MSxaax_expr_data)[1:50],slot = 'data')

#########################################################
# can try to discard the genes that are expressed in less than 50% of cells per cluster

# For MSxapa
MSxapa_gene_name_space <- list()
for (cluster_ident_name in levels(Idents(ms_new_MSxaax))) {
  ms_tmp<- subset(ms_new_MSxaax, idents = c(cluster_ident_name))
  df_tmp <- as.data.frame(ms_tmp@assays$RNA@data)
  percent_not_zero <- rowSums(df_tmp!=0)/ncol(df_tmp)*100 # percent of cells that expression is not 0
  tmp_list <- list(names(percent_not_zero)[percent_not_zero > 50]) # names of genes that are expressed in more than 50% of cells
  names(tmp_list) <- cluster_ident_name
  MSxapa_gene_name_space <- c(MSxapa_gene_name_space,tmp_list)
}
MSxapa_gene_name_space
library(VennDiagram)
venn.diagram(MSxapa_gene_name_space,output=TRUE, filename = "test.png") # look into the working folder - the image is there



############
# try to filter the TFs
# the list of genes is from https://genomebiology.biomedcentral.com/articles/10.1186/gb-2005-6-13-r110#MOESM1 
# A compendium of Caenorhabditis elegans regulatory transcription factors: a resource for mapping transcription regulatory networks


tflist <- (read.delim('./all_TFs_celegans.csv', header = FALSE))[,1]

# All TFs combined in the lineages.
MSxapa_gene_list <- unique(c(
intersect(MSxapa_gene_name_space[[1]], tflist),
intersect(MSxapa_gene_name_space[[2]], tflist),
intersect(MSxapa_gene_name_space[[3]], tflist)))
DoHeatmap(ms_new_MSxaax, features = MSxapa_gene_list, slot = 'data')

# Intersect of the genes
MSxapa_gene_list <- intersect(intersect(MSxapa_gene_name_space[[1]], MSxapa_gene_name_space[[2]]),MSxapa_gene_name_space[[3]])
MSxapa_gene_list <- intersect(MSxapa_gene_list, tflist)
DoHeatmap(ms_new_MSxaax, features = MSxapa_gene_list, slot = 'data')




#########################################################
# For MSxapp
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
#library(VennDiagram)
# venn.diagram(MSxapp_gene_name_space,output=TRUE, filename = "test.png") # look into the working folder - the image is there


tflist <- (read.delim('./all_TFs_celegans.csv', header = FALSE))[,1]
# All TFs combined in the lineages.
MSxapp_gene_list <- intersect(MSxapp_gene_name_space[[1]], tflist)
  
DoHeatmap(ms_new_MSxapp, features = MSxapp_gene_list, slot = 'data')



#########################################################
# 
# For MSxpxx
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
library(VennDiagram)
venn.diagram(MSxpxx_gene_name_space,output=TRUE, filename = "test.png") # look into the working folder - the image is there



############
# try to filter the TFs
# the list of genes is from https://genomebiology.biomedcentral.com/articles/10.1186/gb-2005-6-13-r110#MOESM1 
# A compendium of Caenorhabditis elegans regulatory transcription factors: a resource for mapping transcription regulatory networks


tflist <- (read.delim('./all_TFs_celegans.csv', header = FALSE))[,1]

# All TFs combined in the lineages.
MSxpxx_gene_list <- unique(c(
  intersect(MSxpxx_gene_name_space[[1]], tflist),
  intersect(MSxpxx_gene_name_space[[2]], tflist),
  intersect(MSxpxx_gene_name_space[[3]], tflist),
  intersect(MSxpxx_gene_name_space[[4]], tflist)))
DoHeatmap(ms_new_MSxpxx, features = MSxpxx_gene_list, slot = 'data')

# Intersect of the genes
MSxpxx_gene_list <- intersect(MSxpxx_gene_name_space[[4]],intersect(intersect(MSxpxx_gene_name_space[[1]], MSxpxx_gene_name_space[[2]]),MSxpxx_gene_name_space[[3]]))
MSxpxx_gene_list <- intersect(MSxpxx_gene_list, tflist)
DoHeatmap(ms_new_MSxpxx, features = MSxpxx_gene_list, slot = 'data')




























###########
#trying CoD
install.packages("cvcqv")
library(cvcqv)
ms_new_MSxapa <- subset(ms_new, idents = c('MSxapa'))

df_tmp <- as.data.frame(ms_new_MSxapa@assays$RNA@data)
df_tmp[row.names(df_tmp) == 'act-4',]

as.numeric(df_tmp[row.names(df_tmp) == 'his-24',])
cqv_versatile(as.numeric(df_tmp[row.names(df_tmp) == 'his-24',]))
cqv_versatile(as.numeric(df_tmp[row.names(df_tmp) == 'mog-4',]))

sqv <- function(x){cqv_versatile(x)[2]}

apply(df_tmp[1:10,], 1, sqv)


cqv_versatile(as.numeric(df_tmp[row.names(df_tmp) == 'ceh-34',]))[2]

sum(df_tmp[row.names(df_tmp) == 'aap-1',] == 0)

sum(rowSums(df_tmp==0)/ncol(df_tmp)*100 < 50)


##############
# 


######
ms_new_MSxaax_expr_data <- arrange(ms_new_MSxaax[["RNA"]]@meta.features, desc(vst.mean))
View(ms_new_MSxaax_expr_data[ms_new_MSxaax_expr_data$vst.mean>10,])

#DoHeatmap(ms_new_MSxaax, features = c('ceh-28', 'nhr-281', 'dyf-19', "pha-4"))

ms_new_MSxaax_expr_data <- ms_new_MSxaax_expr_data %>% filter(row.names(ms_new_MSxaax_expr_data) %in% pha4_targets)
DoHeatmap(ms_new_MSxaax, features = rownames(ms_new_MSxaax_expr_data)[1:50])
DoHeatmap(ms_new_MSxaax, features = 'act-4', slot = 'data')

