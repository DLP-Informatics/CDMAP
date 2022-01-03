#GCcontent.r
#David L Patton
#Sung Lab

# Purpose: The objective of this script is to dynamically generate a sorted list of GC content 
# (Greates to least) for each organism for mapping correlations by GC content.
#path_to_FASTA  <- "/Users/triadge/Desktop/Sung_Lab/testdir/GC_directory/FASTA"

FASTA_names_array <- c() #name vector to hold name values
#FASTA_obj_array <- c() #vector to hold fasta sequences
FASTA_gc_array <- c() #vector to hold gc content

GC_cols <- c('name', 'GC')
GC_output_matrix <- matrix( nrow =0, ncol =2)
colnames(GC_output_matrix) <- GC_cols

directory <- list.files(path_to_FASTA)
setwd(path_to_FASTA)

path_correlate_output <- paste(path_analyze, "/Correlation_Results", sep = "")
if(!(dir.exists(path_correlate_output)))
{
  dir.create(path_correlate_output)
}

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

setwd(path_correlate_output)
write.csv(GC_sort_matrix, "GC_Organism_Content.csv") 

