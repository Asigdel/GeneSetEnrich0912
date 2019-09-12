### August 4 ###################
## First STEP #################
## From SNPs to GENES 
rm(list = ls())
setwd("~/Documents/Documents - Anil’s MacBook Air/CODES_PLOTS/PLOTS_MANHATTAN_GENESET/Geneset_Folder/")

# biomaRT is a package from a bioconductor in R
# This is a collection of biomart databases 
library(biomaRt)
listMarts()

# we try github

# we are using ENSEMBL_... database and btaurus_gene dataset by using the function usemart
# We are trying to query the Ensembl bioMart webservices
# useMart function queries the database Ensembl


# Selecting bioMart database and dataset in one query 

genome = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                 dataset = "btaurus_gene_ensembl", 
                 host = "www.ensembl.org")


# How to build a biomaRt query?
# getBM is the main query function in biomaRt and we are calling it with the attribute names
# We are retrieving gene_symbols and chromosomal_ids.

gene = getBM(c("ensembl_gene_id", "start_position", "end_position", "chromosome_name"), mart = genome)
str(gene)

ndata = read.table("SNP_list_60619.txt", header = T) 
dim(ndata)
names(ndata)[1] = "SNP"; names(ndata)[2] = "Chromosome"; names(ndata)[5] = "Location"

plusi = 15000; plusf = 15000           #upstream (promoter) and downstream (regulation) regions
Name = character()
Chromosome = numeric()
Location = numeric()
Gene = numeric()
m = 1

for(k in 1:29)
{
  
  SNP = subset(ndata, ndata$Chromosome == k)
  genes = subset(gene, gene$chromosome_name == k)
  
  for(i in 1:length(SNP$Chromosome))
  {
    for(j in 1:length(genes$chromosome_name))
    {
      
      if(genes$start_position[j] <= (plusi + SNP$Location[i]) & (SNP$Location[i] - plusf) <= genes$end_position[j])
      {
        Name[m] = as.character(SNP$SNP[i])
        Chromosome[m] = SNP$Chromosome[i]
        Location[m] = SNP$Location[i]
        Gene[m] = genes$ensembl_gene_id[j]
        m = m + 1
      }
    }
  }
}

SNPtoGENES = data.frame(Name,Chromosome,Location,Gene)
save(SNPtoGENES, file = "FINAL_SNPtoGENES_15000and15000_August4.rda")


str(SNPtoGENES)
length(unique(SNPtoGENES$Gene))

# 
x <- load("FINAL_SNPtoGENES_15000and15000_August4.rda")
x
head(SNPtoGENES)

## 1 SNP in the same location, same chromosome position can flag many genes !!!



##################################################
# STEP 2 in GENE SET ANALYSIS IN R
#################################################
## Gene Set Analysis
## FILE freqdata.count.after.clean
filter = read.table("freqdata.count.after.clean", header = F)
colnames(filter) = c("SNP", "Frequency", "Filter")
f = filter$Filter == 0; table(f)


## FILE snp_sol
data = read.table("snp_sol_slope_1", header = F, dec = "."); dim(data)
data = data[f,]; dim(data) ## Clean data through logical vector of TRUE
names(data) = c("Trait","Effect","SNP","Chromosome","Position",
                "SNPSolution", "Weight", "Variance"); str(data)


## Significant SNPs Top 5%
ndata = data; thr = 0.05       #for 1% -> 0.005 * 2 = 0.01 = 1%
summary(abs(ndata$SNPSolution))
filter = abs(ndata$SNPSolution) >= quantile(abs(ndata$SNPSolution), 1 - thr)
table(filter)
sndata = ndata[filter,]; dim(sndata) # logical vector create and create a significant SNP

## load from SNP to GENES
load(file = "FINAL_SNPtoGENES_15000and15000_August4.rda") ## FILE output Rcode SNP to Genes
SNP.GENES = SNPtoGENES; str(SNP.GENES)
dim(SNP.GENES)

## Merge (total genes) SNPtoGENES with SNP.Solution
fdata = merge(ndata, SNP.GENES, 
              by.x = c("Chromosome", "Position"),
              by.y = c("Chromosome", "Location"), sort = F); dim(fdata) 

## Merge (Significant genes ) SNPtoGenes with SNP.Solution
sdata = merge(sndata, SNP.GENES, 
              by.x = c("Chromosome", "Position"),
              by.y = c("Chromosome", "Location"), sort = F); str(sdata); dim(sdata)

## Prepare data for Gene Set Analysis

total.genes = unique(fdata$Gene);length(total.genes)
write.table(total.genes,file="total.genes.slope1.txt",sep=" ",append=F, row.names=T,col.names=T)
sig.genes = unique(sdata$Gene);length(sig.genes)
write.table(sig.genes,file="sig.genes.slope1.txt",sep=" ",append=F, row.names=T,col.names=T)


##################################################################################
# After preparing datasets for GeneSet Analysis,
# We are using different databases to extract genesets  
##################################################################################


###################
## GO/KEGG Analysis
library(goseq)
library(org.Bt.eg.db)
library(GO.db)
library(KEGG.db)
library(GOSemSim)
library(corrplot)

assayed.genes = array(total.genes)
de.genes = array(sig.genes)
gene.vector = as.integer(assayed.genes%in%de.genes)
names(gene.vector) = assayed.genes
table = table(gene.vector)

cat(paste("Significant Genes: ", table[2], 
          " Backgroung Genes:", table[1], "\n"))

# An object containing the gene Names
pwf = nullp(gene.vector, "bosTau4", "ensGene", plot.fit = T)

## KEGG
KEGG.hiper <- goseq(pwf, "bosTau4", "ensGene", 
                    method = "Hypergeometric", test.cats = "KEGG", use_genes_without_cat = TRUE)

nKEGG.hiper = KEGG.hiper[(KEGG.hiper$numInCat <= 500 & KEGG.hiper$numInCat >= 5),] 
enriched.KEGG = nKEGG.hiper$category[nKEGG.hiper$over_represented_pvalue <= 0.05]
length(enriched.KEGG)

kegg = as.list(KEGGPATHID2NAME)
for(j in 1:length(enriched.KEGG)){
  for (i in 1:length(names(kegg))) 
  {
    if(names(kegg[i]) == enriched.KEGG[j]){
      
      cat("#############################################\n")	
      cat(paste("KEGG ID: ", enriched.KEGG[j], "\n"))
      cat(paste("KEGG Term Name: ", kegg[[i]], "\n"))
      pvalue = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$over_represented_pvalue
      pvalue = round(pvalue,4)
      nGenes = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numInCat
      nDEG = nKEGG.hiper[nKEGG.hiper$category == enriched.KEGG[j],]$numDEInCat
      cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
      cat("\n")
    }
  }
}

## GO Biological Processes
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:BP",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}

## GO Molecular Function
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:MF",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}

## GO Celular Component
GO.hiper = goseq(pwf, "bosTau4", "ensGene", 
                 method = "Hypergeometric", test.cats="GO:CC",use_genes_without_cat = TRUE)

nGO.hiper = GO.hiper[(GO.hiper$numInCat <= 500 & GO.hiper$numInCat >= 5),] 
enriched.GO = nGO.hiper$category[nGO.hiper$over_represented_pvalue <= 0.05]
length(enriched.GO)

for (go in enriched.GO[1:length(enriched.GO)]) {	
  
  cat("#############################################\n")
  cat(paste("GOID: ", GOTERM[[go]]@GOID, "\n"))
  cat(paste("GO Term Name: ", GOTERM[[go]]@Term, "\n"))
  cat(paste("Definition: ", GOTERM[[go]]@Definition, "\n"))
  
  pvalue = nGO.hiper[nGO.hiper$category == go,]$over_represented_pvalue
  pvalue = round(pvalue,5)
  nGenes = nGO.hiper[nGO.hiper$category == go,]$numInCat
  nDEG = nGO.hiper[nGO.hiper$category == go,]$numDEInCat
  cat(paste("Total Genes:", nGenes, " Number DEG:", nDEG, " P-value:", pvalue))
  cat("\n")
}


################
## MeSH Analysis

library(org.Bt.eg.db)
library(meshr)
library(MeSH.db)
library(MeSH.Bta.eg.db)

# keys return the keys for the database contained in the MeSHdb object
key.symbol = keys(org.Bt.eg.db,  keytype = c("SYMBOL"))
entrezUniverse = select(org.Bt.eg.db, as.character(key.symbol), 
                        columns = c("ENTREZID", "ENSEMBL"),keytype = "SYMBOL")
entrezUniverse2 <- entrezUniverse[!duplicated(entrezUniverse[,2]),]
entrezUniverse3 <- entrezUniverse2[!duplicated(entrezUniverse2[,1]),]

## FULL GENES
genes.back = data.frame(total.genes)
colnames(genes.back) <- "ENSEMBL"
geneID.back <- merge(genes.back, entrezUniverse3, by ="ENSEMBL")
geneID2.back <- geneID.back[ !duplicated(geneID.back[,2]),]

## SIGNIFICANT GENES
genes.sig = data.frame(sig.genes)
colnames(genes.sig) <- "ENSEMBL"
geneID.sig <- merge(genes.sig, entrezUniverse3, by ="ENSEMBL")
geneID2.sig <- geneID.sig[ !duplicated(geneID.sig[,2]),]

## Total Genes
ns = length(geneID2.sig[,1])
nt = length(geneID2.back[,1])
cat(paste("Significant Genes:", ns, " and Backgroung Genes:", nt - ns), "\n")

### MeSH Phenomena and Processes

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "G", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH_Term_ID", "MeSH_Term_Name", 
                  "Total_Genes", "DE_Genes", "P-value")
print(unique(out), row.names = F)

#write.table(unique(out),file="MESH_Phenomena_process.txt",append = F, quote = F,sep ="\t",na = "NA",row.names=F,col.names=T)



### MeSH Chemicals and Drugs

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "D", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name", 
                  "NT.Genes", "DE.Genes", "P-value")
print(unique(out), row.names = F)



### MeSH Diseases

meshParams <- new("MeSHHyperGParams", geneIds = geneID2.sig[,3], 
                  universeGeneIds = geneID2.back[,3], 
                  annotation = "MeSH.Bta.eg.db", category = "A", database = "gene2pubmed", 
                  pvalueCutoff = 0.05, pAdjust = "none")
meshR <- meshHyperGTest(meshParams)
out = data.frame(meshR@ORA$MESHID, meshR@ORA$MESHTERM, 
                 meshR@ORA$Size, meshR@ORA$Count, signif(meshR@ORA$Pvalue,2))
colnames(out) = c("MeSH Term ID", "MeSH Term Name", 
                  "Total Genes", "DE Genes", "P-value")
print(unique(out), row.names = F)
write.table(unique(out),file="ouput_out.txt")

#####################################################################################


## Revealing Significant Genes in KEGG Terms

KEGGterm = "04670"    #change the term
kegg = getgo(names(gene.vector), "bosTau4", "ensGene", fetch.cats = c("KEGG"))
genes = as.character(); m = 1

for(i in 1:length(kegg)){
  if(length(kegg[[i]]) > 0){
    for(j in 1:length(kegg[[i]])){
      if(kegg[[i]][j] == KEGGterm){
        genes[m] = names(kegg[i])
        m = m + 1
      }
    }
  }
}

table(genes%in%sig.genes)
sort(genes[genes%in%sig.genes])

## Revealing Significant Genes in GO Terms

GOterm = "GO:0070201"
go = getgo(names(gene.vector), "bosTau4", "ensGene")
genes = as.character(); m = 1

for(i in 1:length(go)){
  if(length(go[[i]]) > 0){
    for(j in 1:length(go[[i]])){
      if(go[[i]][j] == GOterm){
        genes[m] = names(go[i])
        m = m + 1
      }
    }
  }
}

table(genes%in%sig.genes)
sort(genes[genes%in%sig.genes])

## Revealing Significant Genes in MeSH Terms

cls = columns(MeSH.Bta.eg.db); cl = cls[c(1,2)]
kts = keytypes(MeSH.Bta.eg.db); kt = kts[2]
ks = head(keys(MeSH.Bta.eg.db, keytype=kts[2]))
d = select(MeSH.Bta.eg.db, keys=ks, columns=cls, keytype=kt)

MESHID = "D005982"
data = d[d$MESHID == MESHID,]
f = geneID2.sig$ENTREZID %in% data$GENEID; table(f);geneID2.sig[f,]















