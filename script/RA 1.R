#=============================================Mapping=======================================================
set.seed(42)  # voor reproduceerbaarheid

#laad in de packes
library(BiocManager)
BiocManager::install('Rsubread')
library(Rsubread)

#unzip de map Homo_sapiens.GRCh38.dna.toplevel.fa.gz eerst
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa")
#ONTHOU HET IS TRIAL AND EROR VOOR ZOEKEN NAAR GENOOM. DE GENOOM EN GTF FILE MOET VAN 1 DATABSE AFKOMEN VOOR DIT GEVAL IS HET ENSEMBLE.
buildindex(
  basename = 'ref_RA2',
  reference = 'Homo_sapiens.GRCh38.dna.toplevel.fa',
  memory = 10000,
  indexSplit = TRUE)

#unzip de map Data_RA_raw2 eerst
#setwd("C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Data_RA_raw2/Data_RA_raw")
#controle bam files maken
align.Con19 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785819_1_subset40k.fastq",
                     readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "Con19.BAM")
align.Con20 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785820_1_subset40k.fastq",
                     readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "Con20.BAM")
align.Con28 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785828_1_subset40k.fastq",
                     readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "Con28.BAM")
align.Con31 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785831_1_subset40k.fastq",
                     readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "Con31.BAM")

#RA bam files maken
align.RA79 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785979_1_subset40k.fastq",
                     readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "RA79.BAM")
align.RA80 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785980_1_subset40k.fastq",
                    readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "RA80.BAM")
align.RA86 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785986_1_subset40k.fastq",
                    readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "RA86.BAM")
align.RA88 <- align(index = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.dna.toplevel.fa/ref_RA2", readfile1 = "SRR4785988_1_subset40k.fastq",
                    readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "RA88.BAM")

#installeer en laad in package
BiocManager::install('Rsamtools')
library(Rsamtools)

#sample names op basis van de naam die je hebt gegeven in "output_file"
samples <- c('Con19', 'Con20', 'Con28', 'Con31', 'RA79', 'RA80', 'RA86', 'RA88')

# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

#maak ook sorted bam files
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.bam'))
})

#============================================Count matrix======================================================

# Je definieert een vector met namen van BAM-bestanden. Elke BAM bevat reads van een RNA-seq-experiment (bijv. behandeld vs. controle).
allsamples <- c("con19.BAM", "con20.BAM", "con28.BAM", "con31.BAM", "RA79.BAM", "RA80.BAM", "RA86.BAM", "RA88.BAM")

#setwd("C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Data_RA_raw2/Data_RA_raw")
#maak de count matrix
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "C:/Users/Hp/OneDrive - NHL Stenden/transcriptoom analyse (R)/Homo_sapiens.GRCh38.114.gtf/Homo_sapiens.GRCh38.114.gtf",
  isPairedEnd = TRUE, #van elk sample is een forward en reverse fasta file
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE) #Exons die tot hetzelfde gen behoren worden samengevoegd en het totaal wordt per gen gerapporteerd

#naam van de gen ophalen
head(count_matrix$annotation)
#hoeveel read er zijn op een gen
head(count_matrix$counts)
#de hele structuur van de countmatrix
str(count_matrix)

# Haal alleen de matrix met tellingen eruit
counts <- count_matrix$counts
# voeg kolomnamen aan de counts.
colnames(counts) <- c("Con19", "Con20", "Con28", "Con31", "RA79", "RA80", "RA86", "RA88")
#sla bestand op
write.csv(counts, "bewerkt_countmatrix.csv")

#==========================================================statestiek en analyse=================================================
#een nieuwe count matrix werd gegeven die compleet is.
countstest <- read.table("count_matrix.txt")

#maak een mata data tabel van elke sample en of het contole of RA is
treatment <- c("control", "control", "control", "control", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")

#instaleer en laad in
BiocManager::install('DESeq2')
library(DESeq2)

#instaleer en laad in
BiocManager::install('KEGGREST')
library(KEGGREST)

# Maak DESeqDataSet object aan
dds <- DESeqDataSetFromMatrix(countData = round(countstest),
                              colData = treatment_table,
                              design = ~ treatment)


# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

# Resultaten opslaan in een bestand
#Bij het opslaan van je tabel kan je opnieuw je pad instellen met `setwd()` of het gehele pad waar je de tabel wilt opslaan opgeven in de code.
write.table(resultaten, file = 'ResultatenWC3.csv', row.names = TRUE, col.names = TRUE)

#---------------------------------PCA plot + anova en heatmap-------------------------------------------

#installer package
BiocManager::install('gridGraphics')
BiocManager::install('patchwork')

#instaleer en laad in
library(vegan)
library(ggplot2)
library(grid)
library(patchwork)
library(gridGraphics)

#PCA plot
#Eerst moeten we de ruwe telgegevens transformeren
#De vst-functie voert een variantiestabiliserende transformatie uit. dat sneller de dispersion kan determineren
vsdata <- vst(dds, blind=FALSE) #aantal genen zijn op default (1000)
#plot de PCA
p1 <- plotPCA(vsdata, intgroup="treatment") + stat_ellipse() + geom_text(
  label=dds@colData@rownames, 
  nudge_x = 5.25, nudge_y = 6.25, 
  check_overlap = TRUE) + theme_classic()

#zie plot
p1

# bepaal significantievan de afstand tussen groepen op basis euclidische afstand
# (Bijvoorbeeld berekend uit VST- of rlog-data)
#haalt de matrix met genexpressie op:
#Rijen = genen
#Kolommen = samples
vsd_mat <- assay(vsdata)

# VST matrix
#t() transponeert de matrix.
#Nu: Rijen = samples, kolommen = genen
vsd_t <- t(vsd_mat)                # Transponeer: samples in rijen

#berken de afstand matrix
euc_dist <- vegdist(vsd_t, method = "euclidean") #de pcaplot heeft standaard euclidian afstand

#Zorg dat metadata overeenkomt met samples
treatment_table$treatment <- as.factor(treatment_table$treatment)

#Voer PERMANOVA uit
permanova_result <- adonis2(euc_dist ~ treatment, data = treatment_table, permutations = 999)

#Resultaten bekijken
print(permanova_result)

#controle
# p < 0.05 → ongelijke spreiding (PERMANOVA met voorzichtigheid interpreteren) gebruik dan liever ANOSIM
#als p waarde groter dan 0.05 is dan gebruik PERMANOVA (wat net gedaan is)
dispersion <- betadisper(euc_dist, treatment_table$treatment)
anova(dispersion)   


#heatmap
# Bereken de variantie van elk gen
gene_var <- apply(vsd_mat, 1, var)

# Selecteer de top 30 meest variabele genen
top_genes <- head(order(gene_var, decreasing = TRUE), 40)

# Subset matrix
mat_top <- vsd_mat[top_genes, ]

# Metadata (bijv. behandeling)
annotation_col <- as.data.frame(colData(vsdata)[, "treatment", drop = FALSE])
colnames(annotation_col) <- "Sample"  # Naam voor heatmap-annotatie

# Kies een kleurenschema (optioneel)
ann_colors <- list(Group = c("RA" = "red", "control" = "blue"))  # Pas aan aan jouw groepen


#maak een object van pheatmap anders word het niet herkent als een plot
heatmap_obj <- grid::grid.grabExpr(
  pheatmap(mat_top,
           scale = "row",                     
           annotation_col = annotation_col,    
           annotation_colors = ann_colors,     
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 6,
           fontsize_col = 7,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)),  # jouw pheatmap code hier
)

#combineer de 2 plots
(wrap_elements(p1) + wrap_elements(heatmap_obj)) + 
  plot_annotation(title = "PCA plot en heatmap van 30 variabele genen") &
  theme(plot.title = element_text(hjust = 0.5))


#------------------------------------------------------------------------------------------------------

#hoeveel genen meer zijn tot expressie
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
#hoeveel genen minder zijn in expressie
sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

hoogste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = TRUE), ] # zet de genen met de sterkste opregulatie bovenaan.
laagste_fold_change <- resultaten[order(resultaten$log2FoldChange, decreasing = FALSE), ] #zet de genen met de sterkste neerregulatie bovenaan
laagste_p_waarde <- resultaten[order(resultaten$padj, decreasing = FALSE), ] # zet de genen met de meest significante p-waarden bovenaan.

head(laagste_p_waarde)
head(hoogste_fold_change)
head(laagste_fold_change)

#--------------optioneel:speciefiek zoeken naar genen en vergelijken tussen groepen----------------------------------------------------------

#zoek de genen. HLA-DRA
test <- rownames(resultaten)
write.csv(test, "genen.csv")

#zoek sig. 0.04
padj <- resultaten@listData[["padj"]]
write.csv(padj, "pvalue.csv")

#log. 2.01882992963716
fold <- resultaten@listData[["log2FoldChange"]]
write.csv(fold, "foldlog.csv")

#welke genen activer zijn in controle of treated. baseer intrgroup op treatment.
#check rownames welke genen interesant zijn. hoogste foldchange betekent meer genen actief in RA.
#laagste foldchange werkt andersom.
par(mfrow=c(2,3))
plotCounts(dds, gene="HLA-DRA", intgroup="treatment")

#----------------------------------------volcano plot---------------------------------------------------------------------------

#installeer package
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}

#laad in package
library(EnhancedVolcano)

#volcano plot met sommige genen namen
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj')

# Alternatieve plot zonder p-waarde cutoff (alle genen zichtbaar)
EnhancedVolcano(resultaten,
                lab = rownames(resultaten),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

# uitleg:
#hoe hoger de waarde in de y-as de significanter de gen is.
#de x as de foldchange. de waardes aan de linkerkant zijn downregulated (de - waardes)
# aan de rechterkant zijn de waardes upregulated (de + waarde)

#voor opslaan van plot
dev.copy(png, 'VolcanoplotWC.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#--------------------------------------KEGG bepaalde pathways analyse----------------------------------------------

#installeer package
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}

#laad in package
library(pathview)

# Laad de resultaten opnieuw in als data frame
resultaten <- read.table("ResultatenWC3.csv", header = TRUE, row.names = 1)

#laad in package
library(org.Hs.eg.db)
library(AnnotationDbi)

# verandert genen naar ENTREZ ID-kolom
resultaten$ENTREZID <- mapIds(
  org.Hs.eg.db, #humaan genoom
  keys = rownames(resultaten),   # de SYMBOLS staan in rownames
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

#Verwijdert rijen met NA-waarden
resultaten_clean <- na.omit(resultaten)

#vector met log2foldchange en entrezid
gene_vector <- setNames(resultaten_clean$log2FoldChange, resultaten_clean$ENTREZID)

#laad in package
library(pathview)

#pathway van up of down regulated genen
pathview(
  gene.data = gene_vector,
  pathway.id = "05323",  #een pathway dat je kan kiezen. zie link 1 voor meer keuze
  species = "hsa",
  gene.idtype = "ENTREZ",
  limit = list(gene = 3),
  low = list(gene = "green"), 
  mid = list(gene = "gray"),
  high = list(gene = "red")
)
#link1 : https://www.genome.jp/kegg/pathway.html

#-----------------------------------------GO ENRICHMENT-------------------------------------------

#Installeer benodigde Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#installeer de human genome
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Laad libraries
library(goseq)
library(org.Hs.eg.db)
library(GO.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(dplyr)

# Laad resultaten in
resultaten <- read.table("ResultatenWC3.csv", header = TRUE, row.names = 1)

#Voeg ENTREZ ID's toe
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(resultaten),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

#voeg entrez_id toe aan resultaten
resultaten$ENTREZID <- entrez_ids

#verwijder waarden met NA
resultaten <- na.omit(resultaten)

#Maak de binaire genvector (1 = DE, 0 = niet-DE)
de_genen <- as.integer(resultaten$padj < 0.05 & !is.na(resultaten$padj))

#label de vector met de corresponderende ENTREZ-ID's
names(de_genen) <- resultaten$ENTREZID

#verwijder genen zonder geldige ENTREZ-ID of padj
de_genen <- na.omit(de_genen)

#Haal genlengtes op via TxDb voor hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_lengths <- transcriptLengths(
  txdb,
  with.cds_len = TRUE,
  with.utr5_len = TRUE,
  with.utr3_len = TRUE
)

#Bereken totale lengte
gene_lengths$tot_length <- gene_lengths$cds_len + gene_lengths$utr5_len + gene_lengths$utr3_len

#Gemiddelde lengte per gen (ENTREZ ID)
length_vector <- tapply(gene_lengths$tot_length, gene_lengths$gene_id, mean)

#Zorg dat genen overeenkomen in beide objecten
common_ids <- intersect(names(de_genen), names(length_vector))
de_genen <- de_genen[common_ids]
bias_vector <- length_vector[common_ids]

#Bereken probability weighting function (PWF) oftewel de nul hypothese
pwf <- nullp(de_genen, bias.data = bias_vector)

#uitleg over deze plot
#houd rekenening meer dat genen die langer zijn meer reads hebben en corrigeert dat in de GO enrichment
#xas: is lengte van de genen.
#yas: hoeveel van de lange genen in expressie zijn in %

#Voer GO-enrichment analyse uit
# Maak eigen GO-annotatie
gene2go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("ENTREZID", "GO"),
  keytype = "ENTREZID"
)
gene2go <- na.omit(gene2go)
gene2cat <- split(gene2go$GO, gene2go$ENTREZID)

# GO enrichment analyse met eigen annotatie
GO.wall <- goseq(pwf, gene2cat = gene2cat, use_genes_without_cat = TRUE)


#Filter op significant verrijkte GO-termen (p < 0.05)
GO.sig <- GO.wall[GO.wall$over_represented_pvalue < 0.05, ]

#Voeg GO-term beschrijvingen toe
GO.sig$term <- Term(GO.sig$category)

# Sla op als CSV-bestand
write.csv(GO.sig, "GO_enrichment_significant.csv", row.names = FALSE)

#Bekijk de top GO-termen
head(GO.sig[, c("term", "over_represented_pvalue", "numDEInCat", "numInCat")])

#maak een plot voor BP, CC en MP
GO <- read.table("GO_enrichment_significant.csv", header = TRUE, sep = ";")

GO2 <- head(GO, n = 30)

df <- GO2 %>%
  arrange(ontology, desc(over_represented_pvalue)) %>%
  mutate(term = factor(term, levels = unique(term)))

#maak de plot
ggplot(df, aes(term, -log10(over_represented_pvalue), color = ontology, fill = ontology )) + geom_col() + 
  theme(panel.background = element_rect(fill = NA), axis.text.x = element_text(size = 6)) + coord_flip()

#uitleg
#category:	De GO-ID van de term, bijv. "GO:0006955" voor “immune response”
#term	De beschrijving van de GO-term (toegevoegd via Term())
#over_represented_pvalue:	P-waarde voor verrijking: geeft aan of deze term vaker voorkomt dan verwacht
#under_represented_pvalue:	P-waarde voor ondervertegenwoordiging (zelden belangrijk in DEG-lijsten)
#numDEInCat:	Aantal differentieel tot expressie gekomen genen in deze GO-term
#numInCat:	Totaal aantal genen in deze GO-term (uit het hele genoom)

#dus over_represented_value als die is veel van jouw DE-genen vallen in deze GO-term, met een zeer lage kans op toeval → belangrijk!
#zo laag mogelijke value!!!
#de under_represented_pvalue geeft aan de onderdrukte genen.
#Een p-waarde van 0.0001 wordt dus -log10(0.0001) = 4
#de log10 waarde word gebruikt om het meer visueel te laten zien inplaats van veel nullen.
# DE LAGE P WAARDE BETEKENT MEER UPREGULATED IN RA


#voor een plot. van de top 10 sig over_represented_pvalue.
library(ggplot2)

# Sorteer op p-waarde (laagste bovenaan)
GO.sig.sorted <- GO.sig[order(GO.sig$over_represented_pvalue), ]

# Pak de top 10 meest significante termen
top_terms <- head(GO.sig.sorted, 10)

# Zet termen als factor met juiste volgorde (voor mooie y-as)
top_terms$term <- factor(top_terms$term, levels = rev(top_terms$term))

# Plot
ggplot(top_terms, aes(x = term, y = -log10(over_represented_pvalue))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 verrijkte GO-Enrichment-termen (RA vs. controle)",
    x = "GO-term",
    y = "-log10(p-waarde)"
  ) +
  theme_minimal(base_size = 12)

#-----------------------------------Gene set testing (kijkt naar pathway)------------------------------------------

#installeer de package
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("msigdbr")   # MSigDB gene sets

#laad in de package
library(limma)
library(edgeR)
library(msigdbr)
library(DESeq2)
library(KEGGREST)
library(ggplot2)

#maak de dds weer als je die nog niet hebt
countstest <- read.table("count_matrix.txt")

treatment <- c("control", "control", "control", "control", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c("SRR4785819", "SRR4785820", "SRR4785828", "SRR4785831", "SRR4785979", "SRR4785980", "SRR4785986", "SRR4785988")

#Maak DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(countstest),
                              colData = treatment_table,
                              design = ~ treatment)


# Voer analyse uit
dds <- DESeq(dds)
resultaten <- results(dds)

#Gebruik een catergory naar keuze (kies 1 collection + subcollection en ga verder bij "<-" ). elke verteld iets anders
msig_go_bp <- msigdbr(species = "Homo sapiens", collection = "H") #geeft globaal aan
msig_go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "MF") #vervang BP naar CC of MF voor andere processen
msig_go_bp <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") #over bijv interferon-pathways, antigen-presentatie
msig_go_bp <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS") #over signaling van ziektes of pathways

#zie voor meer combos met collection en subcollection
print(msigdbr_collections(), n = 30)


#Maak een lijst van genensets op basis van SYMBOL <-
gene_sets <- split(msig_go_bp$gene_symbol, msig_go_bp$gs_name)

#Haal ruwe counts
dge <- DGEList(counts = counts(dds), genes = rownames(dds))

#Normaliseer
dge <- calcNormFactors(dge)

# Design matrix op basis van RA/controle
group <- dds$treatment  # bevat metadata van "control" en "RA"
design <- model.matrix(~ group)

# Transformeer met voom naar lineare model
v <- voom(dge, design, plot = TRUE)

#uitleg
#Zwarte stippen	Elk punt is een gen
#X-as	Hoeveel dat gen gemiddeld is uitgedrukt(dus zo laag mogelijk, het is wat je verwacht als genen hoog in expressie zijn)
#Y-as	Hoe variabel dat gen is tussen samples
#Rode lijn	Hoe de variantie gemiddeld afneemt bij hogere expressie (LOESS)
#Vorm van lijn	Als die daalt: alles werkt zoals het hoort

# camera test op genensets door genen te ranken met andere genen
camera_results <- camera(v$E, index = gene_sets, design = design, contrast = 2)

#Toon topverrijkte pathways
head(camera_results[order(camera_results$PValue), ])

#Filter op p < 0.05
camera_sig <- camera_results[camera_results$PValue < 0.05, ]

#uitleg
#gs_name:	Naam van de genenset (hier: GO Biological Process)
#NGenes:	Aantal genen in deze genenset
#Direction:	Of de genenset gemiddeld op- of neer-gereguleerd is
#PValue:	P-waarde voor verrijking
#FDR:	False Discovery Rate (gecorrigeerde p-waarde, dus betrouwbaarder) (de lager de beter)
#focust voornamelijk op de processen.

#sla op bestand
write.csv(camera_sig, "camera_enrichment_RA_vs_control.csv")

#laat de top 10 zien en maak vector
top_terms <- head(camera_sig[order(camera_sig$PValue), ], 10)

#geef rownames
top_terms$gs_name <- rownames(top_terms)

#word als "up" als het DE van 1 heeft
top_terms$Direction_num <- ifelse(top_terms$Direction == "Up", 1, -1)

#de plot. verander de titel en vector naam naar welke subcategory je hebt gebruikt (voor als je 3 plot samenvoegd).
MF <- ggplot(top_terms, aes(x = reorder(gs_name, -Direction_num * log10(PValue)),
                      y = -log10(PValue),
                      fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Genenset (GO term)",
       y = "-log10(p-waarde)",
       title = "GO-MF") +
  theme_minimal() + scale_fill_manual(values = c("#87CEFA")) #vernader de kleur van de proccesen naar wat je wil
MF #zie plot

#maak de combo plot
plot <- (BP / CC / MF) + 
  plot_annotation(title = "3 processen pathways voor reuma geraleerte DE") &
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 6.3))

#opslaan met betere resolutie
ggsave("Rplot_3_paths.png", plot = plot, width = 12, height = 8, dpi = 300)
