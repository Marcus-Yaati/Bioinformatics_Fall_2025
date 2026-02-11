# ==============================================================================
# MULTI-LOCUS PHYLOGENOMICS PIPELINE (AB1 -> CONCATIPEDE -> TREE)
# ------------------------------------------------------------------------------
# Context: Raw Sanger files (.ab1) -> Cleaned FASTA -> Multi-Locus Tree
# Logic: 
#   1. Conversion: Convert .ab1 chromatograms to FASTA using sangerseqR.
#   2. Aggregation: Split filenames by ";" to separate Animal from Locus.
#   3. Alignment: Align each locus separately (DECIPHER).
#   4. EXPORT: Write aligned loci to 'aligned_loci' folder.
#   5. Concatenation: Run 'concatipede_prepare' then 'concatipede'.
#   6. Tree: Build ML tree + Export.
# ==============================================================================

# --- CONFIGURATION (UPDATED FROM YOUR FILE) ---
setwd("G:/Shared drives/bioinf/Research_Project/Scripts")
# 1. AB1 CONVERSION SETTINGS
# Where your raw .ab1 files are located
AB1_INPUT_DIR  <- "G:/Shared drives/bioinf/Research_Project/Scripts/MB/Localized_dataset_copy"

# Where the converted FASTAs will be saved (and read by the next steps)
FASTA_DIR <- "G:/Shared drives/bioinf/Research_Project/Scripts/MB/fasta_localized_dataset"

# 2. ANALYSIS OUTPUT SETTINGS
# Where final results (Alignments, Trees, Matrices) go
OUTPUT_DIR     <- "G:/Shared drives/bioinf/Research_Project/Scripts/MB/output_dir"

TRIM_PCT_START <- 0.15                  # 15% cut from start
TRIM_PCT_END   <- 0.15                  # 15% cut from end

# ==============================================================================
# STEP 1: PACKAGE INSTALLATION & LOADING 
# ==============================================================================
message(">>> STEP 1: CHECKING DEPENDENCIES...")

# Create directories if they don't exist
if (!dir.exists(FASTA_DIR)) dir.create(FASTA_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

ensure_package <- function(pkg, source="CRAN", git_repo=NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    if (source == "Bioc") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update=FALSE, ask=FALSE)
    } else if (source == "GitHub") {
      if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
      remotes::install_github(git_repo, upgrade="never")
    } else {
      install.packages(pkg)
    }
  }
}

ensure_package("BiocManager", "CRAN")
ensure_package("sangerseqR", "Bioc")
ensure_package("Biostrings", "Bioc")
ensure_package("DECIPHER", "Bioc")
ensure_package("phangorn", "CRAN")
ensure_package("phytools", "CRAN")
ensure_package("ape", "CRAN")
ensure_package("remotes", "CRAN")
ensure_package("concatipede", "GitHub", "tobias-hofmann/concatipede")
ensure_package("svglite", "CRAN") # Added for better SVG rendering

library(sangerseqR)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(phytools)
library(ape)
library(concatipede)
library(svglite)

# ==============================================================================
# STEP 0: AB1 TO FASTA CONVERSION
# ==============================================================================
message(">>> STEP 0: CONVERTING AB1 FILES TO FASTA...")

if (!dir.exists(AB1_INPUT_DIR)) {
  stop(paste("AB1 Input directory not found:", AB1_INPUT_DIR))
}

ab1_files <- list.files(AB1_INPUT_DIR, pattern = "\\.ab1$", full.names = TRUE)

if (length(ab1_files) == 0) {
  message("No .ab1 files found. Checking for existing FASTA files to proceed...")
} else {
  message(paste("Found", length(ab1_files), ".ab1 files. Converting..."))
  
  for (f in ab1_files) {
    # Extract filename without extension
    base_name <- tools::file_path_sans_ext(basename(f))
    out_name <- file.path(FASTA_DIR, paste0(base_name, ".fasta"))
    
    # Skip if already exists to save time (optional, remove check to overwrite)
    if (file.exists(out_name)) next
    
    tryCatch({
      # Read chromatogram
      ab1 <- read.abif(f)
      seq <- sangerseq(ab1)
      
      # Extract primary basecalls
      s <- as.character(primarySeq(seq))
      
      # Write FASTA
      writeLines(c(paste0(">", base_name), s), out_name)
    }, error = function(e) {
      message(paste("  [Error] Failed to convert:", base_name, "-", e$message))
    })
  }
  message("AB1 conversion complete.")
}

# ==============================================================================
# STEP 2: AGGREGATION (Bucket Files by Locus)
# ==============================================================================
message(">>> STEP 2: AGGREGATING FILES BY LOCUS (MVZ Parsing)...")

# We look in FASTA_DIR (where we just saved the converted files)
file_list <- list.files(FASTA_DIR, pattern = "\\.(fasta|fa|fna)$", full.names = TRUE, ignore.case = TRUE)

if (length(file_list) == 0) stop(paste("No FASTA files found in", FASTA_DIR))

loci_buckets <- list()

for (fpath in file_list) {
  fname <- basename(fpath)
  clean_name <- tools::file_path_sans_ext(fname)
  parts <- strsplit(clean_name, ";")[[1]]
  
  if (length(parts) < 2) next # Skip if format doesn't match
  
  animal_id <- parts[1] 
  locus_id  <- parts[2] 
  
  seq <- readDNAStringSet(fpath)
  if (length(seq) == 0) next
  
  # IMPORTANT: Rename sequence to Animal ID for matching later
  names(seq) <- animal_id
  
  # Add to bucket
  if (is.null(loci_buckets[[locus_id]])) {
    loci_buckets[[locus_id]] <- seq
  } else {
    loci_buckets[[locus_id]] <- c(loci_buckets[[locus_id]], seq)
  }
}

found_loci <- names(loci_buckets)
message(paste("Processing complete. Found", length(found_loci), "unique loci."))

# ==============================================================================
# STEP 3: TRIM, ALIGN & WRITE TO DISK
# ==============================================================================
message(">>> STEP 3: ALIGNING & WRITING LOCUS FILES...")

ALIGNED_DIR <- file.path(OUTPUT_DIR, "aligned_loci")
if (!dir.exists(ALIGNED_DIR)) dir.create(ALIGNED_DIR)

# IMPORTANT: Initialize the list in memory so Step 4 can find it!
aligned_loci_list <- list()

for (locus in names(loci_buckets)) {
  message(paste("Processing Locus:", locus))
  raw_seqs <- loci_buckets[[locus]]
  if (length(raw_seqs) < 3) {
    message(paste("  [!] Skipping", locus, "- too few sequences."))
    next
  }
  
  # Trim
  seq_widths <- width(raw_seqs)
  start_indices <- floor(seq_widths * TRIM_PCT_START) + 1
  end_indices   <- seq_widths - floor(seq_widths * TRIM_PCT_END)
  
  valid_mask <- end_indices >= start_indices
  trimmed_seqs <- subseq(raw_seqs[valid_mask], start=start_indices[valid_mask], end=end_indices[valid_mask])
  
  # Align
  aligned <- AlignSeqs(trimmed_seqs, verbose=FALSE)
  
  # 1. Write to disk (Backup)
  out_path <- file.path(ALIGNED_DIR, paste0(locus, "_RHVA", ".fasta"))
  writeXStringSet(aligned, out_path)
  
  # 2. Store in memory (Crucial for Step 4)
  aligned_loci_list[[locus]] <- as.DNAbin(as.matrix(aligned))
}

message(paste("Alignments saved to disk and memory."))

# ==============================================================================
# STEP 4: CONCATIPEDE (Prepare & Stitch)
# ==============================================================================
message(">>> STEP 4: RUNNING CONCATIPEDE...")

# 1. PREPARE: This scans the directory and builds the matching dataframe
#    It finds matching sequence names (taxa) across the different locus files.
message("  -> Preparing correspondence table...")
concat_result <- find_fasta(ALIGNED_DIR) %>% #finds fasta files in the current dir
  concatipede_prepare() %>%                  #prepares the concatenation file
  
  write_xl("template.xlsx") %>%              #sets the excel file to be used to create the finished binned excel file
  
  auto_match_seqs() %>%                      #automatches the sequences in the input excel to create bins based on the genbank names in the cells for the individuals
  
  write_xl("Concat_table.xlsx") %>%  #sets the output excel file
  
   concatipede() %>%                          #does the actual concatenation 
  
  write_fasta("merged-seqs") %>%       #creates a new FASTA
  image(cex=0.3)                             #creates the concatenation as an image in the plot area

#THE HANDOFF!!! 
#use the concatenated file that was written as the output to then input into phangorn

message(paste("Reading back concat file", "merged-seqs"))
final_dna <- read.FASTA("merged-seqs.fasta")

message("Concatenation complete.")

# ==============================================================================
# STEP 5: GENETIC DISTANCE MATRIX
# ==============================================================================
message(">>> STEP 5: CALCULATING GENETIC DISTANCES...")

dist_mat <- dist.dna(final_dna, model="raw", pairwise.deletion = TRUE)
dist_df <- as.matrix(dist_mat)
write.csv(dist_df, file.path(OUTPUT_DIR, "genetic_distance_matrix.csv"))

# ==============================================================================
# STEP 6: TREE CONSTRUCTION
# ==============================================================================
message(">>> STEP 6: TREE CONSTRUCTION...")

phydat_data <- phyDat(final_dna)

# A. Maximum Parsimony
message("  -> Constructing Maximum Parsimony Tree...")
dist_raw <- dist.ml(phydat_data)
if(any(is.na(dist_raw))) dist_raw[is.na(dist_raw)] <- max(dist_raw, na.rm=TRUE) 

tree_start <- NJ(dist_raw)
#parsimony rachet...whatever that is
tree_parsimony <- pratchet(phydat_data, start=tree_start, trace=0, minit=10, k=5)
#assign branch lengths
tree_parsimony <- acctran(tree_parsimony, phydat_data)

# B. Maximum Likelihood
message("  -> Constructing Maximum Likelihood Tree (GTR Model)...")
fit <- pml(tree_parsimony, phydat_data)

fit_ml <- optim.pml(fit, model="GTR", optInv=TRUE, optGamma=TRUE, 
                    rearrangement = "NNI", control = pml.control(trace = 1)) 

message("  -> Running Bootstrap analysis (100 reps)...")
bs_ml <- bootstrap.pml(fit_ml, bs=100, optNni=TRUE, control = pml.control(trace = 0))

# ==============================================================================
# STEP 7: VISUALIZATION & EXPORT
# ==============================================================================

message(">>> STEP 7: VISUALIZATION & EXPORT...")

# Clear devices
graphics.off()

# --- PART A: Process ML Tree Object ---
# Calculate bootstrap values first in memory
# p=0 prevents plotting to screen
final_ML_tree <- plotBS(midpoint(fit_ml$tree), bs_ml, type="phylogram", p=0)

# Convert proportions to percentages
if (!is.null(final_ML_tree$node.label)) {
  lbls <- suppressWarnings(as.numeric(final_ML_tree$node.label))
  if (!all(is.na(lbls)) && max(lbls, na.rm=TRUE) <= 1) {
    message("  -> Converting Bootstrap Proportions to Percentages...")
    final_ML_tree$node.label <- round(lbls * 100)
  }
}

# --- PART C: Save ML Tree (SVG via svglite) ---
# svglite provides much better font rendering and stability than base svg()
svg_ML_filename <- file.path(OUTPUT_DIR, "ML_concatenated_tree.svg")

# Standard Size: 10x12 inches is perfect for ~14-20 animals
svglite(svg_ML_filename, width=8, height=10)

# We plot the PROCESSED tree (final_ML_tree) which has the % labels
# We use p=0 here because the labels are already in the tree object
plotBS(final_ML_tree, type="phylogram", p=0, 
       main="Concatenated ML Tree", 
       bs.col="blue", bs.adj=c(1.5, 1.5), 
       cex=1.0,       # Large text for readability
       edge.width=2)  # Thicker lines for nice visuals

add.scale.bar(cex=1.0) 
dev.off() 
message(paste("ML Tree SVG saved to:", svg_ML_filename))


# --- PART C: Save Parsimony Tree (SVG via svglite) ---
svg_Parsimony_filename <- file.path(OUTPUT_DIR, "Parsimony_concatenated_tree.svg")
svglite(svg_Parsimony_filename, width=8, height=10)

plot(midpoint(tree_parsimony), type="phylogram", 
     main="Maximum Parsimony Tree", 
     cex=1.0, edge.width=2)
add.scale.bar(cex=1.0)

dev.off()
message(paste("Parsimony Tree SVG saved to:", svg_Parsimony_filename))


# --- PART D: Export Newick Files ---
write.tree(final_ML_tree, file=file.path(OUTPUT_DIR, "ML_tree_bootstrapped.tre"))
write.tree(tree_parsimony, file=file.path(OUTPUT_DIR, "Parsimony_tree_bootstrapped.tre"))
message("Newick tree files (.tre) exported.")

message("Pipeline Finished!")

