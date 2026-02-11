#Script Title: Final singe loci Phylo code
#Authors: Marcus Barela,Ariel Goldsztejn,Maiya Martinez, Gemini thinking model V3
#Date:8.XI.2025


# ==============================================================================
# SINGLE LOCUS PHYLOGENY (CuttyDude -> DECIPHER -> phangorn -> phytools -> SVG)
# ==============================================================================

# --- CONFIGURATION ---
INPUT_DIR      <- "H:/Shared drives/bioinf/Research_Project/Scripts/MB/fasta_localized_dataset"
OUTPUT_DIR     <- "H:/Shared drives/bioinf/Research_Project/Scripts/MB/output_dir"
OUTGROUP_LABEL <- "RHOL_H59_AY764254" # Must match the FASTA header exactly!

# Trimming percentages
TRIM_PCT_START <- 0.01
TRIM_PCT_END   <- 0.01

# --- PACKAGES ---
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(phytools)
library(svglite)

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ==============================================================================
# STEP 1: LOAD & TRIM -CuttyDude
# ==============================================================================
message(">>> Loading and Aligning Sequences...")

# 1. Load FASTA files
filepaths <- list.files(INPUT_DIR, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE)
raw_seqs  <- readDNAStringSet(filepaths)

# 2. Trim Sequences by percentage of sequence length as set by the variables TRIM_PCT_position
seq_widths    <- width(raw_seqs)
start_indices <- floor(seq_widths * TRIM_PCT_START) + 1
end_indices   <- seq_widths - floor(seq_widths * TRIM_PCT_END)
trimmed_seqs  <- subseq(raw_seqs, start=start_indices, end=end_indices)

# ==============================================================================
# STEP 2: ALIGN SEQUENCES -DECIPHER
# ==============================================================================

# 3. Align Sequences
aligned_seqs <- AlignSeqs(trimmed_seqs, verbose = FALSE)
final_dna    <- as.DNAbin(aligned_seqs) 

# ==============================================================================
# STEP 3: BUILD ML TREE & BOOTSTRAP -Phangorn
# ==============================================================================
message(">>> Building Maximum Likelihood Tree...")
#converts alignments to a phangorn readable format
phydat_data <- phyDat(final_dna)

# Generate Starting Tree
  #applies ML algorithm and makes a distance matrix
dist_mat   <- dist.ml(phydat_data)
  #makes a rough neighbor joining tree
tree_start <- NJ(dist_mat)

# Optimize ML (GTR+G+I)
  #
fit    <- pml(tree_start, phydat_data)
fit_ml <- optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE, control=pml.control(trace=0))

# Run Bootstrap (1000 replicates)
message(">>> Running Bootstrap (1000 reps)...")
bs_ml <- bootstrap.pml(fit_ml, bs=1000, optNni=TRUE, control=pml.control(trace=0))

# ==============================================================================
# STEP 4: ROOT & FINALIZE -Phangorn 
# ==============================================================================
message(">>> Processing Tree Data...")
# 1. Root the tree
rooted_tree <- root(fit_ml$tree, outgroup = OUTGROUP_LABEL, resolve.root = TRUE)

# 2. Map Bootstrap values onto the ROOTED tree
tree_with_bs <- plotBS(rooted_tree, bs_ml, type = "none")

# 3. Convert Bootstraps to Percentages (if they are 0.95 instead of 95)
if (!is.null(tree_with_bs$node.label)) {
  lbls <- suppressWarnings(as.numeric(tree_with_bs$node.label))
  if (!all(is.na(lbls)) && max(lbls, na.rm=TRUE) <= 1) {
    tree_with_bs$node.label <- round(lbls * 100)
  }
}

# 4. Turn sample labels to animal labels  (remove ";MVZ15" etc)
tree_with_bs$tip.label <- gsub(";.*$", "", tree_with_bs$tip.label)

# ==============================================================================
# STEP 4: VISUALIZE -Phytools & SVGlite
# ==============================================================================
message(">>> Creating SVG")

svg_path <- file.path(OUTPUT_DIR, "ML_Phylogeny_wbootstrap_Collective Data.svg")

# Open SVG device
svglite(file = svg_path, width = 30, height = 8)

# A. Plot the Tree
# We plot the tree separately so 'lwd' (line width) works correctly
plotTree(tree_with_bs, 
         ftype = "i",       # Italic font
         lwd = 2,           # Line width (This caused the error in plotBS, but works here)
         fsize = 0.8,       # Font size
         offset = 0.5,      # Spacing
         mar = c(1,1,3,1))  # Margins

title("ML Phylogeny (GTR+G+I) with Total Dataset ")

# B. Add Bootstrap Values
# add them manually
bs_values <- as.numeric(tree_with_bs$node.label)
# Filter: Show only values >= 50%
bs_labels <- ifelse(bs_values >= 50, bs_values, "")

nodelabels(text = bs_labels,
           frame = "none",      # No box
           adj = c(1.2, -0.2),  # Position
           cex = 0.7,           # Text size
           col = "blue")        #Color of text

dev.off()

message(paste("DONE! SVG saved to:", svg_path))


