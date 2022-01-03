### GC_MultiFold_Analysis.R###
#================================#
# Objective: 
# 1. Conduct analysis of all N-fold sites using genomic GC content to generate the probability P(codon) of each nucleotide triplet
# 2. For each Triplet calculated, calculate the weighted probability of each triplet, W(P(codon)) = P(codon)/SUM(P(codon))
# 3. Calculate the expected # of codons Exp(codon) = W(P(codon))*SUM(OBS_codon)
# 4. Chi Square Test: chisq.test(OBS_codon, EXP_codon)

FASTA_names_array <- c()
FASTA_gc_array <- c()
GC_cols <- c('Name', 'GC Content')
GC_output_matrix <- matrix( nrow =0, ncol =2)
colnames(GC_output_matrix) <- GC_cols


directory <- list.files(path_to_FASTA)
setwd(path_to_FASTA)

for(i in 1:length(directory))
{
  FASTA_name <- unlist(strsplit(directory[i], split='.', fixed = TRUE))[1]
  FASTA_names_array <- c(FASTA_names_array, FASTA_name)
  Seq_obj <- getSequence(read.fasta(directory[i]))
  Seq_obj <<- unlist(Seq_obj) 
  Seq_obj <- toupper(Seq_obj)
  GC_obj <- GC(Seq_obj)
  FASTA_gc_array <- c(FASTA_gc_array, GC_obj)
  
  output_row <- c(FASTA_name, GC_obj)
  GC_output_matrix <- rbind(GC_output_matrix, output_row)
}

rownames(GC_output_matrix) <- NULL
GC_output_matrix <- data.frame(GC_output_matrix)
GC_sort_matrix <- GC_output_matrix[order(GC_output_matrix$GC, decreasing = TRUE),]
