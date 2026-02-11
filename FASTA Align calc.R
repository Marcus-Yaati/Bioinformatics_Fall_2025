# ==============================================================================
# FASTA SIMILARITY CALCULATOR
# ==============================================================================
# This script calculates the percent similarity between a reference FASTA file
# and every other FASTA file in a specific directory.
#
# Output: An Excel file (.xlsx) with columns: "FileName" and "PercentSimilarity"
#
# PREREQUISITES:
# You need the 'Biostrings' (from Bioconductor) and 'writexl' packages.
# If you have not installed them, uncomment and run the following lines once:
#
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings")
# install.packages("writexl")
# ==============================================================================

# 1. LOAD LIBRARIES
suppressPackageStartupMessages({
  library(Biostrings)
  library(writexl)
  library(tools)
  library(pwalign)
})

# =========================== USER CONFIGURATION ===============================
# PLEASE EDIT THESE 3 VARIABLES BEFORE RUNNING
setwd("H:/Shared drives/bioinf/Research_Project/Scripts/MB")
# 1. Path to your single reference FASTA file
#for MVZ15: "H:/Shared drives/bioinf/Research_Project/Scripts/MB/fasta_localized_dataset/RHVA_13;MVZ15.fasta"
#for MVZ16:"H:/Shared drives/bioinf/Research_Project/Scripts/MB/MVZ16_stuff/RHVA_13;MVZ16.fasta" 
ref_path <- "H:/Shared drives/bioinf/Research_Project/Scripts/MB/fasta_localized_dataset/RHVA_13;MVZ15.fasta"

# 2. Path to the folder containing your query FASTA files
#    (Ensure there are no trailing slashes)
query_dir <- "unchecked_fasta"

# 3. Name of the output Excel file
output_filename <- "fasta_similarity_results_MVZ15.xlsx"

# ==============================================================================

# 2. LOAD REFERENCE SEQUENCE
cat("--------------------------------------------------\n")
cat("Processing...\n")

if (!file.exists(ref_path)) {
  stop(paste("ERROR: Reference file not found at:", ref_path))
}

cat(paste("Reading reference:", basename(ref_path), "\n"))

# NOTE: Using readDNAStringSet. If using Proteins, change to readAAStringSet
tryCatch({
  ref_seq <- readDNAStringSet(ref_path)[1]
}, error = function(e) {
  stop("Error reading reference file. Is it a valid FASTA format?")
})


# 3. GET LIST OF FILES
# Looks for .fasta, .fa, or .fna extensions (case insensitive)
file_list <- list.files(path = query_dir, 
                        pattern = "\\.(fasta|fa|fna)$", 
                        full.names = TRUE,
                        ignore.case = TRUE)

if(length(file_list) == 0) {
  stop(paste("No FASTA files found in directory:", query_dir))
}

cat(paste("Found", length(file_list), "query files. Starting alignment...\n"))


# 4. MAIN LOOP: ALIGNMENT AND CALCULATION
results_list <- list()
counter <- 0

for (file_path in file_list) {
  counter <- counter + 1
  file_name <- basename(file_path)
  
  # Optional: skip if the file is the reference file itself (if it's in the same folder)
  if (normalizePath(file_path) == normalizePath(ref_path)) {
    next
  }
  
  # Read Query Sequence
  # We use tryCatch to ensure one bad file doesn't crash the whole script
  query_seq <- tryCatch({
    readDNAStringSet(file_path)[1]
  }, error = function(e) {
    warning(paste("Could not read:", file_name, "- Skipping."))
    return(NULL)
  })
  
  if (!is.null(query_seq)) {
    # Perform Global Pairwise Alignment
    # type="global" aligns the full length (Needleman-Wunsch)
    alignment <- pairwiseAlignment(pattern = query_seq, subject = ref_seq, type = "global")
    
    # Calculate Percent Identity
    # type="PID1" = (identical matches / length of alignment) * 100
    percent_sim <- pid(alignment, type = "PID1")
    
    # Add to list
    results_list[[length(results_list) + 1]] <- data.frame(
      FileName = file_name, 
      PercentSimilarity = round(percent_sim, 2) # Rounded to 2 decimal places
    )
    
    # Print progress every 10 files
    if (counter %% 10 == 0) {
      cat(paste("Processed", counter, "files...\n"))
    }
  }
}


# 5. EXPORT RESULTS
if (length(results_list) > 0) {
  
  # Combine all rows into one Data Frame
  final_df <- do.call(rbind, results_list)
  
  # Sort by highest similarity
  final_df <- final_df[order(-final_df$PercentSimilarity), ]
  
  # Write to Excel
  write_xlsx(final_df, output_filename)
  
  cat("--------------------------------------------------\n")
  cat(paste("DONE! Successfully analyzed", length(results_list), "files.\n"))
  cat(paste("Results saved to:", file.path(getwd(), output_filename), "\n"))
  
} else {
  cat("No valid alignments were generated. Please check your files.\n")
}