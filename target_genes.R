## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Autor: Francisco J. Romero-Campero - fran@us.es

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana
if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("HDO.db")

BiocManager::install("enrichplot")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Leer fichero de picos
sob3.peaks <- readPeakFile(peakfile = "SOB3.narrowPeak",header=FALSE)

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

## Anotación de los picos
sob3.peakAnno <- annotatePeak(peak = sob3.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)

plotAnnoPie(sob3.peakAnno)
plotAnnoBar(sob3.peakAnno)
plotDistToTSS(sob3.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
upsetplot(sob3.peakAnno) ##picos que se encuentran en varias zonas simultaneamente##

plotPeakProf2(peak = sob3.peaks, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)

## Convertir la anotación a data frame
sob3.annotation <- as.data.frame(sob3.peakAnno)
head(sob3.annotation)

target.genes <- sob3.annotation$geneId[sob3.annotation$annotation == "Promoter"]

write(x = target.genes,file = "sob3_target_genes.txt")

## Enriquecimiento funcional. Para ver si regulan una funcion biologica comun
library(clusterProfiler)
library(org.At.tair.db)

sob3.enrich.go <- enrichGO(gene = target.genes,
                           OrgDb         = org.At.tair.db,
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable      = FALSE,
                           keyType = "TAIR")

barplot(sob3.enrich.go,showCategory = 20, cex_label_category= 0.5)
dotplot(sob3.enrich.go,showCategory = 20, cex_label_category= 0.5)
emapplot(pairwise_termsim(sob3.enrich.go),showCategory = 20, cex_label_category= 0.5)
cnetplot(sob3.enrich.go,showCategory = 20, cex_label_category= 0.5, cex_label_gene= 0.5)
library(enrichplot)
sob3.enrich.kegg <- enrichKEGG(gene  = target.genes,
                               organism = "ath",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.05)
df.sob3.enrich.kegg <- as.data.frame(sob3.enrich.kegg)
head(df.sob3.enrich.kegg)

## ChIPpeakAnno es un paquete de R de Bioconductor que implementa 
## análisis de los resultados del procesamiento de los datos de ChIP-seq
## alternativo y complementario al presentado anteriormente. 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
