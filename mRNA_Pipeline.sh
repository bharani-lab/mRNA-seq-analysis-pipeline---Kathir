##Quality Pre-Processing :
cd /home/Data/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o File1_R1_filtered.fastq -p File2_R2_filtered.fastq  File1_R1.fastq File2_R2.fastq
##Reference Indexing:
STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /home/STAR/Indexed_Reference/ --genomeFastaFiles /home/Homo_sapiens.GRCh38.77.fa --sjdbGTFfile /home/Homo_sapiens.GRCh38.77.gtf --sjdbOverhang 101
##Mapping:
STAR --genomeDir /home/STAR/Indexed_Reference/ --runThreadN 32 --sjdbGTFfile /home/Homo_sapiens.GRCh38.77.gtf --readFilesIn /home/Data/File1_R1_filtered.fastq  /home/Data/File1_R2_filtered.fastq --outSAMtype BAM SortedByCoordinate --sjdbOverhang 101 --outFileNamePrefix /home/Mapping/File
##Quantification
cd /home/Mapping/File
for prefix in $(ls *.bam | uniq)
do
featureCounts -T 30 -p -t exon -g gene_name -a /home/Homo_sapiens.GRCh38.77.gtf -o ${prefix%.bam}.csv ${prefix}
done
##Move files
mkdir Mapping Quantification DiffExp
mv *.bam *.out *.tab Mapping
mv *.csv Quantification
##post Quantification filter >9 
	#remove unwanted coloumns 
cd Quantification
for prefix in $(ls *sortedByCoord.out.csv | uniq); do awk 'NR>1 {print $1,$7}' ${prefix} > ${prefix}.csv; done
rm *sortedByCoord.out.csv
	#Combine all csv files into one
paste *Aligned.sortedByCoord.out.csv.csv > Counts.csv 
awk '{print $1,$2,$4,$6,$8,$10,$12,$14,$16,$18,$20,$22,$24,$26,$28,$30,$32}' Counts.csv > Counts_1.csv
rm Counts.csv
	#Remove mRNAs <9 counts
awk '$2>9 || $3>9 || $4>9 || $5>9 || $6>9 || $7>9 || $8>9 || $9>9 || $10>9 || $11>9 || $12>9 || $13>9 || $14>9 || $15>9 || $16>9 || $17>9 {print}' Counts_1.csv > Counts.csv
##Differential Expression Analysis
cd DiffExp
R
setwd("/home/Quantification")
library(readr)
Counts <- read_csv("~/Quantification/Counts_R>9.csv")
GeneNames <- Counts$Geneid
Counts <- Counts[, -c(1)]
row.names(Counts) <- GeneNames
#Condition file creation
file <- c('Etoh_Counts1','Etoh_Counts2','Etoh_Counts3','Etoh_Counts4','Dex_Counts1','Dex_Counts2','Dex_Counts3','Dex_Counts4')
condition <- c('1EtoH','1EtoH','1EtoH','1EtoH','2Dex','2Dex','2Dex','2Dex')
Info <- data.frame(file, condition)
library(edgeR)
library(HTSFilter)
d <- Counts
d <- DGEList(counts = d, group = Info$condition)
TMM <- calcNormFactors(d, method="TMM")
d <- estimateCommonDisp(TMM)
d <- estimateTagwiseDisp(d)
et <- exactTest(d)
p <- et$table$PValue
adj.p <- p.adjust(p, method = "BH")
res <- cbind(id=rownames(et$table), et$table, adj.p, threshold=-p)
Down <- res[res$PValue<0.05,]
DEG <- Down[order(Down$logFC), ]
write.csv(DEG, "/home/DiffExp/DEG_R.csv")
#Make a basic volcano plot
# Add colored points: red if padj<0.05, orange of log2FC>2, green if both)
#Make a basic volcano plot
with(res, plot(logFC, -log10(PValue), pch=20, xlim=c(-8,8)))
with(subset(res,PValue<.05 & abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="red"))
library(calibrate)
with(subset(res,adj.p<.05 & abs(logFC)>2), textxy(logFC, -log10(PValue), labs=id, cex=.6))
# Add colored points: red if padj<0.05, orange of log2FC>2, green if both)
library(EnhancedVolcano)
res <- res[, c(-1)]
EnhancedVolcano(res, lab = rownames(res), x = 'logFC', y = 'PValue', xlim = c(-7, 7), ylim = c(0, 25), FCcutoff = 2, pointSize = 2, labSize = 3)
##Pathway Enrichment Analysis
R
library(readr)
Res <- read_csv("/home/DiffExp/DEG_R.csv")
attach(Res)
Res <- na.omit(Res)
Res$fcsign=sign(logFC)
Res$logP= -log10(PValue)
Res$metric=Res$logP/Res$fcsign
y<-Res[,c("Gene", "metric")]
write.table(y, file = "/home/DiffExp/DEG_R_expression.rnk", quote = F, sep = "\t", row.names = F)
#Running GSEA
library(fgsea)
library(tibble)
pathways <- gmtPathways("/home/c2.cp.symbols.gmt")
ranks = deframe(y)
fgseaRes <- fgsea(pathways=pathways, stats=ranks, nperm=1000)
fgseaRes2 = data.frame(lapply(fgseaRes, as.character), stringsAsFactors=FALSE)
write.csv(fgseaRes2, "/home/Pathways_DEG_R.csv")#Plotting
#Plotting
library(ggplot2)
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=pval<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="pathways NES from GSEA") + 
  theme_minimal()
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
