setwd("C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/MiSeq_SOP")

#Instalar paquetes y cargarlos
# Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Core packages
BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "Biostrings", "ShortRead"))
install.packages(c("ggplot2", "dplyr", "tidyr", "vegan", "phangorn"))

library(dada2)
library(phyloseq)
library(DECIPHER)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(phangorn)  # still useful e.g. midpoint.root if you want
library(ape)  

#Set path

path <- "C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/MiSeq_SOP"    # where the fastq files live
list.files(path)

#Vectors forward-reverse

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

# Extract sample names from filenames e.g. "F3D0_S188_L001_R1_001.fastq" -> "F3D0"
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#Quality profiles + filtering / trimming (DADA2)

# Pick a couple of samples to show quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#Filter and trim

filt_path <- file.path(path, "filtered")
if(!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Adjust truncLen based on quality plots for your run
out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240, 160),   # typical for MiSeq V4; adjust as needed
  maxN    = 0,
  maxEE   = c(2, 2),
  truncQ  = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
head(out)

#Learn error models + infer ASVs (DADA2 core)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualize error models
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#Dereplication

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Attach sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample-wise inference

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Merge paired data

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

#Build sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)    # samples x ASVs

# Inspect length distribution of ASVs
table(nchar(getSequences(seqtab)))

#Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)  # fraction of non-chimeric reads

#Optional: track how many reads survived each step (nice for teaching):

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input  = out[,1],
  filtered = out[,2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged   = sapply(mergers, getN),
  nonchim  = rowSums(seqtab.nochim)
)
rownames(track) <- sample.names
head(track)

#6. Taxonomy assignment with SILVA
taxa <- assignTaxonomy(seqtab.nochim,
                       "C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/silva_nr_v138_train_set.fa",
                       multithread = TRUE)

taxa <- addSpecies(taxa,
                   "C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/silva_species_assignment_v138.fa.gz")
taxa[1:5, ]


#-------------------------
# 7. Build a phylogenetic tree with IQ-TREE
#    (Model selection by IQ-TREE, ML optimization there)
#-------------------------

# Extract sequences as DNAStringSet
seqs <- DNAStringSet(getSequences(seqtab.nochim))
names(seqs) <- getSequences(seqtab.nochim)  # tip labels = ASV sequences

# Multiple alignment (DECIPHER)
alignment <- AlignSeqs(seqs, anchor = NA, verbose = FALSE)

# Write alignment to FASTA for IQ-TREE
alignment_fasta <- file.path(path, "ASV_alignment.fasta")
Biostrings::writeXStringSet(alignment, filepath = alignment_fasta)

## ===== Locate IQ-TREE 2 binary =====
iqbin <- file.choose()
iqbin <- normalizePath(iqbin)  # cleans up the path

iqbin
file.exists(iqbin)  # MUST be TRUE

# 3) Optional: choose a prefix for output files
iq_prefix <- file.path(path, "ASV_alignment")  # IQ-TREE will use this as prefix

# 4) Build command for IQ-TREE
# -s : alignment
# -nt AUTO : auto threads
# -m MFP : ModelFinder Plus -> best model
# -pre : prefix for output files
cmd <- sprintf('"%s" -s "%s" -nt AUTO -m MFP -pre "%s"',
               iqbin,
               alignment_fasta,
               iq_prefix)

cat("Running command:\n", cmd, "\n")

# 5) Run IQ-TREE
oldwd <- getwd()
setwd(path)
exit_status <- system(cmd)
setwd(oldwd)

cat("IQ-TREE exit status:", exit_status, "\n")

# IQ-TREE will produce ASV_alignment.fasta.treefile in 'path'
treefile <- file.path(path, "ASV_alignment.treefile")

# Read the ML tree
treeIQ <- read.tree(treefile)



#Creacion de arbol filogenetico
library(ape)
library(ggtree)

# 1) Cargar el árbol desde IQ-TREE
arbol <- read.tree("C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/MiSeq_SOP/ASV_alignment.treefile")
# 2) Ver árbol básico
plot(arbol, cex=0.6, no.margin=TRUE)
# 3) Árbol con ggtree (lineal)
ggtree(arbol) +
  geom_tiplab(size = 2)
# 4) Árbol circular
ggtree(arbol, layout = "circular") +
  geom_tiplab(size = 2)

familias <- taxa[, "Family"]
familias[is.na(familias)] <- "Unclassified"
nombres_nuevos <- paste0("ASV", seq_along(familias), "_", familias)
arbol$tip.label <- nombres_nuevos
library(ggtree)
ggtree(arbol, layout = "circular") +
  geom_tiplab(size = 2)

# (Optional but recommended) midpoint rooting
treeIQ <- phangorn::midpoint(treeIQ)
#8. Build a phyloseq object

# OTU/ASV table
otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)

# Taxonomy table
tax <- tax_table(taxa)  # taxa must be a matrix

# Build phyloseq object
ps <- phyloseq(otu, tax, phy_tree(treeIQ))

ps

# For convenience, store SampleID inside sample_data for plotting
sample_data(ps)$SampleID <- sample_names(ps)


#9. Basic data checks

sample_sums(ps) %>% summary()
hist(sample_sums(ps), main = "Reads per sample", xlab = "Read count")


####Ecological analyses

#10. Alpha diversity (richness & evenness)

#Metadata, for testing. Saving in csv.

metadata <- read.csv("C:/Users/YISAO/Desktop/BioInformática/2do parcial/practica/MiSeq_SOP/metadata.csv", stringsAsFactors = FALSE)
head(metadata)

# Sample metadata (ensure rownames match sample names)
metadata$SampleID <- as.character(metadata$SampleID)
row.names(metadata) <- metadata$SampleID
sample_data_ps <- sample_data(metadata)

# Build phyloseq object
ps <- phyloseq(otu, tax, sample_data_ps, phy_tree(treeIQ))

ps

# For convenience, store SampleID inside sample_data for plotting
sample_data(ps)$SampleID <- sample_names(ps)


#9. Basic data checks

sample_sums(ps) %>% summary()
hist(sample_sums(ps), main = "Reads per sample", xlab = "Read count")

#10. Alpha diversity (richness & evenness)

alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
head(alpha_df)


alpha_df$SampleID <- rownames(alpha_df)
alpha_df <- left_join(alpha_df, metadata, by = "SampleID")

# Example: compare alpha diversity by time or group
ggplot(alpha_df, aes(x = Type, y = Shannon)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.6) +
  labs(x = "Sample type", y = "Shannon diversity") +
  theme_bw()


#Using phyloseq built-in:

plot_richness(ps, x = "Type", measures = c("Observed", "Shannon")) +
  geom_boxplot(alpha = 0.6) +
  theme_bw()

#11. Beta diversity + ordination + PERMANOVA

#Distance matrices

# Bray-Curtis
dist_bc <- phyloseq::distance(ps, method = "bray")

# UniFrac (requires tree and rooted tree ideally; we’ll assume ps has a proper rooted tree)
dist_wu <- phyloseq::distance(ps, method = "wunifrac")

dist_bc
dist_wu

#PCoA ordination and plot

# PCoA on Bray-Curtis
ord_bc <- ordinate(ps, method = "PCoA", distance = dist_bc)

p_bc <- plot_ordination(ps, ord_bc, color = "Type") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCoA (Bray-Curtis)", color = "Sample type")

p_bc

#For UniFrac

ord_wu <- ordinate(ps, method = "PCoA", distance = dist_wu)

p_wu <- plot_ordination(ps, ord_wu, color = "Type") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.7) +
  theme_bw() +
  labs(title = "PCoA (weighted UniFrac)", color = "Sample type")

p_wu


#12. Community composition plots (barplots)

# 12.1. Create ps_phylum_rel (this was missing)
ps_phylum <- tax_glom(ps, taxrank = "Phylum", NArm = TRUE)

ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))

# Add SampleID to sample_data for plotting
sample_data(ps_phylum_rel)$SampleID <- sample_names(ps_phylum_rel)


#Barplot
plot_bar(ps_phylum_rel, x = "SampleID", fill = "Phylum") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample")

#Other taxonomy ranks
rank_names(ps)

#Barplot family
# 12.1. Create ps_phylum_rel (this was missing)
ps_fam <- tax_glom(ps, taxrank = "Family", NArm = TRUE)

ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))

# Add SampleID to sample_data for plotting
sample_data(ps_fam_rel)$SampleID <- sample_names(ps_fam_rel)

#Barplot
plot_bar(ps_fam_rel, x = "SampleID", fill = "Family") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Relative abundance", x = "Sample")

rank_names(ps)

ps_fam <- tax_glom(ps, taxrank = "Class", NArm = TRUE)

ps_fam_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))
sample_data(ps_fam_rel)$SampleID <- sample_names(ps_fam_rel)


