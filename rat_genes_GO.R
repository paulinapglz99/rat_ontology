#rat_genes_GO.R
#This script performs Gene ontology for a list of genes

#Install packages --- --- 

install.packages("BiocManager")
install.packages("pacman")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("vroom")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Rn.eg.db")

#Load packages --- ---

pacman::p_load("ggplot2",
               "tidyverse", 
               "clusterProfiler", 
               "org.Rn.eg.db")

#Get data --- --- 

#setwd("tu_carpeta")

prot_lista <- vroom::vroom(file = "UE1034_A_0490.csv")

# Prepare data --- ---

prot_lista$Entry <- gsub("_RAT", "", prot_lista$Entry)

#Enrichment --- ---

enrich_GO.go <- enrichGO(gene = prot_lista$Entry,
                        OrgDb = "org.Rn.eg.db", 
                        keyType = 'UNIPROT',
                        readable = TRUE,
                        ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)

enrich_GO.df <- as.data.frame(enrich_GO.go)

#Gene concept network enrichment

cnet <- cnetplot(enrich_GO.go,
                 showCategory= 20, #cuantos quieres mostrar
                 circular = FALSE, #si quieres que sea circular
                 colorEdge = TRUE) #los edges se coloreen por GO term
cnet

#Dotplot enrichment

dotplot <- dotplot(enrich_GO.go)
dotplot

#Barplot enrichment

barplot <- barplot(enrich_GO.go,
                   showCategory=20) 

#Save data and plots --- ---

#Data
vroom::vroom_write(enrich_GO.df, file = "rat_uniprot_3.csv")

#Plots
ggsave("cnet_rat_uniprot_3.png", 
       plot = cnet, 
       dpi = 300,
       width = 10, #ancho
       height = 10, #altura
       units = "in",
       device = "png")

ggsave("dotplot_rat_uniprot_3.png", 
       plot = cnet, 
       dpi = 300,
       width = 10, #ancho
       height = 10, #altura
       units = "in",
       device = "png")

ggsave("barplot_rat_uniprot_3.png", 
       plot = barplot, 
       dpi = 300,
       width = 10, #ancho
       height = 10, #altura
       units = "in",
       device = "png")

#END
