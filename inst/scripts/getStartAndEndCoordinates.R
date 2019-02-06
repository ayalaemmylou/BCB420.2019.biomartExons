# getStartAndEndCoordinates.R
#
# Purpose: For each gene get a set of chromosomal start/end coordinates for the principal isoform.
# Version: 1.1.1
# Date: 2/5/2019
# Author: Emily Ayala
#
# Input: a vector of HGNC IDs
# Output: A vector of chromosomal start/end coordinates for the principal isoform
# Dependencies:
#   BCB420.2019.STRING
#   BiocManager
#   biomaRt
# ToDo:
# Notes:
#
# ==============================================================================


# WARNING: SIDE EFFECTS
# Executing this script will execute code it contains.

# ====  PARAMETERS  ============================================================
# Define and explain all parameters. No "magic numbers" in your code below.

#ensp2sys.RData is produced by Prof. Steipe in his BCB420.2019.STRING package
load(file = file.path("inst", "extdata", "ensp2sym.RData"))




# ====  PACKAGES  ==============================================================
# Load all required packages.

if (requireNamespace("seqinr", quietly=TRUE)) {
  library("seqinr")
} else {
  install.packages("seqinr")
  library(seqinr)
}
# Package information:
#  library(help = seqinr)       # basic information
#  browseVignettes("seqinr")    # available vignettes
#  data(package = "seqinr")     # available datasets
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

library("biomaRt")

if (! requireNamespace("pracma", quietly = TRUE)) {
  install.packages("pracma")
}

library("pracma")
# ====  FUNCTIONS  =============================================================

# Define functions or source external files

getStartAndEndCoodinates <- function(a) {
	# Purpose: For each gene get a set of chromosomal start/end coordinates for the principal isoform.
	#     Describe ...
	# Parameters:
	#     a: A vector of HNGC symbols
	# Value:
	#     result: Start and end coordinates for each gene

	# ensp2sym contains a map of ENSP IDs to HNGC IDs

  # STEP 1: Find the corresponding ENSP ID for each HNGC ID
  # ensp2sym contains a map of ENSP IDs to HNGC IDs
  attributesList <- attributes(ensp2sym)
  typeof(attributesList)
  v <- c()
  for (i in seq(from = 1, to = length(ensp2sym))){
    if(ensp2sym[[i]] %in% a){
      v <- c(v, attributes(strsplit(ensp2sym[i], "\n"))$names)
    }
  }
  # STEP 2: Determine the ENSG ID for each ENSP ID
  ensembl <- useMart("ensembl")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
  ENSG_ids <- getBM( attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
                   filters = "ensembl_peptide_id",
                   values = v,
                   mart = ensembl)
  ENSGToHGNC <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                       filters = "ensembl_gene_id",
                       values = ENSG_ID$ensembl_gene_id,
                       mart = ensembl)
  Transcript_IDs <- getBM( attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "transcript_length"),
                           filters = "ensembl_gene_id",
                           values = ENSG_ids$ensembl_gene_id,
                           mart = ensembl)
  # STEP 3: Dtermine the principal isoform for each ENSP ID
  isoforms <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "transcript_appris", "transcript_length"),
                    filters = c("ensembl_transcript_id"),
                    value = list(Transcript_IDs$ensembl_transcript_id),
                    mart = ensembl)
  isoforms <- isoforms[isoforms$transcript_appris == "principal1",]
  newIsoforms <- isoforms[order(isoforms$ensembl_gene_id),]
  #Choose the longest
  current <- NULL
  selectedIsoforms <- data.frame(ensembl_gene_id = c(),
                                 ensembl_transcript_id = c(),
                                 transcript_appris = c(),
                                 transcript_length = c())
  for (i in seq(from = 1, to = nrow(newIsoforms) + 1)){
    row <- newIsoforms[i,]
    if(!is.na(row$ensembl_gene_id)){
      if(is.null(current)){
        current <- row
      }
      else if (strcmp(row$ensembl_gene_id, current$ensembl_gene_id)){
        if(as.numeric(row$transcript_length) > as.numeric(current$transcript_length)){
          current <- row
        }
      }
      else{
        if(length(selectedIsoforms) == 0){
          selectedIsoforms <- current
          current <- row
        }
        else{
          selectedIsoforms <- rbind(selectedIsoforms, current)
          current <- row
        }
      }
    }
    else{
      selectedIsoforms <- rbind(selectedIsoforms, current)
    }
  }
  # STEP 4: Retrieve the exons for the principal isoform
  exons <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
                 filter = "ensembl_transcript_id",
                 value = selectedIsoforms$ensembl_transcript_id,
                 mart = ensembl)
  # STEP 5: Determine the start and end coordinates
  current <- NULL
  id <- NULL
  start <- NULL
  transcript_id <- NULL
  end <- NULL
  transcriptsStartEnd <- data.frame(ensembl_gene_id = c(),
                                    ensembl_transcript_id = c(),
                                    start = c(),
                                    end = c())
  for (i in seq(from = 1, to = nrow(exons) + 1)){
    #Get the row
    row <- exons[i,]
    #If it's not NA continue on
    if(!is.na(row$ensembl_gene_id)){
      #If id has not been set, set it along with start and end
      if(is.null(id)){
        id <- row$ensembl_gene_id
        transcript_id <- row$ensembl_transcript_id
        start <- row$exon_chrom_start
        end <- row$exon_chrom_end
      }
      else if (strcmp(row$ensembl_gene_id, id)){
        #Find the start and end
        #It's the same gene
        if(row$exon_chrom_end > end){
          #set end
          end <- row$exon_chrom_end
        }
      }
      else{
        if(length(transcriptsStartEnd) == 0){
          transcriptsStartEnd <- data.frame("ensembl_gene_id" = c(id), "ensembl_transcript_id" = c(transcript_id), "start" = c(start), "end" = c(end), stringsAsFactors = FALSE)
          id <- row$ensembl_gene_id
          transcript_id <- row$ensembl_transcript_id
          start <- row$exon_chrom_start
          end <- row$exon_chrom_end
        }
        else{
          transcriptsStartEnd <- rbind(transcriptsStartEnd, c("ensembl_gene_id" = id, "ensembl_transcript_id" = transcript_id, "start" = start, "end" = end))
          id <- row$ensembl_gene_id
          transcript_id <- row$ensembl_transcript_id
          start <- row$exon_chrom_start
          end <- row$exon_chrom_end
        }
      }
    }
    else{
      transcriptsStartEnd <- rbind(transcriptsStartEnd, c(id, transcript_id, start, end))
    }
  }
  merged <- merge(ENSGToHGNC, transcriptsStartEnd, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", sort = TRUE)
  # STEP 6: Output result
	return(merged)
}



# ====  PROCESS  ===============================================================
# Enter the step-by-step process of your project here. Strive to write your
# code so that you can simply run this entire block and re-create all
# intermediate results.
if (FALSE) {

# ...



}

# ====  TESTS  =================================================================
if (FALSE) {
# Enter your function tests here...

}


# [END]
