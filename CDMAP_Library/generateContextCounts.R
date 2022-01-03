##### generateContextCounts.r #######
#=================================#
#Author: David Logan Patton

# Objective: The objective of this Script is to generate the context mutation counts for a given chromosome of a single organism.
# This file parses the MutBaseCalls data frame and then for each location in MutBaseCalls it grabs the 5mer fragment for each 
# base pair substitution and calculates the linear distance each BPS from the given ORI and TERM locations.
#=================================#

#developer note: Add in flanking information for indels from the VCF Parser!

# Setting up output Matrices:
#==============================================
context_col_names <- c("Position", "Left2", "Left1", "Original", "Mutation", "Right1", "Right2", "WRT_ORI", "WRT_TERM") #output columns for the full text file
context_up_names <- c("Position", "Left1", "Original", "Mutation", "Right1", "Right2", "WRT_ORI", "WRT_TERM")
context_down_names <- c("Position", "Left2", "Left1", "Original", "Mutation", "Right1", "WRT_ORI", "WRT_TERM")


output_full_matrix <- matrix( nrow =0, ncol =9)
colnames(output_full_matrix) <- context_col_names

output_upstream_matrix <- matrix( nrow =0, ncol =8)
colnames(output_upstream_matrix) <- context_up_names

output_downstream_matrix <- matrix( nrow =0, ncol =8)
colnames(output_downstream_matrix) <- context_down_names

Path_Basecall_output <- paste(Path_text_output, "/BaseCall_Output", sep = "")
Path_Basecall_output_up <- paste(Path_text_output_up, "/BaseCall_Output", sep = "")
Path_Basecall_output_down <- paste(Path_text_output_down, "/BaseCall_Output", sep = "")
dir.create(Path_Basecall_output)
dir.create(Path_Basecall_output_up)
dir.create(Path_Basecall_output_down)
#==============================================


# Recording the flanking region outputs of the sequence
#==============================================================================



#storage variables and vectors for binding to output CSV and data frames for manipulation
context_output_full_matrix <<- output_full_matrix
context_output_upstream_matrix <<- output_upstream_matrix
context_output_downstream_matrix <<- output_downstream_matrix

#==============================================
# Recording the flanking region outputs of the sequence
#==============================================================================

#context_diag_matrix 
newrow <<- c() #dummy vector used to append context triplets to base output file
newrow_up <<- c() #dummy vector used to append upstream context 4mers to base output file
newrow_down <<- c() #dummy vector used to append downstream context 4mers to base output file

Left1 <<- '' #nucleotide immediately to the left of the context mutation
Right1 <<- '' #nucleotide immediately to the right of the context mutation
Left2 <<- '' #nucleotide two spaces to the left of the context mutation (upstream left)
Right2 <<- '' #nucleotide two spaces to the right of the context mutation (downstream right)


for(i in 1:length(MutBaseCalls$position))
{
  position <- MutBaseCalls$position[i]
  WRT_ORI <- ori_bp - position
  WRT_TERM <- term_bp - position
  source <- as.character(MutBaseCalls$original[i])
  mutant <- as.character(MutBaseCalls$mutant[i])
  Left1 <- toupper(RefSeq_arr[position-1]) # Immediate left neighbhor
  Left2 <- toupper(RefSeq_arr[position-1]) # left neighbhor 2 spaces
  Right1 <- toupper(RefSeq_arr[position+1]) # Immediate right neighbhor
  Right2 <- toupper(RefSeq_arr[position+2]) # Right neighbhor 2 spaces
  
  newrow <- c()
  newrow_up <- c()
  newrow_down <- c()
  
  newrow <- c(position, Left2, Left1, source, mutant, Right1, Right2, WRT_ORI, WRT_TERM)
  newrow_up <- c(position, Left1, source, mutant, Right1, Right2, WRT_ORI, WRT_TERM)
  newrow_down <- c(position, Left2, Left1, source, mutant, Right1, WRT_ORI, WRT_TERM)
  #print (newrow)
  #print(newrow_abbrev)
  context_output_full_matrix <- rbind(context_output_full_matrix, newrow) #appends new mutant to full output text
  context_output_upstream_matrix <- rbind(context_output_upstream_matrix, newrow_up) #appends new mutant to abbreviated text
  context_output_downstream_matrix <- rbind(context_output_downstream_matrix, newrow_down) #appends new mutant to abbreviated text
  
  i <- i+1
}


rownames(context_output_upstream_matrix) <- NULL
rownames(context_output_downstream_matrix) <- NULL
rownames(context_output_full_matrix) <- NULL
setwd(Path_Basecall_output)
write.csv(context_output_full_matrix, "Context_Output_full.csv") #full 5mer output of a given context mutation
setwd(Path_Basecall_output_up)
write.csv(context_output_upstream_matrix, "Context_Output_upstream.csv") #full 4mer upstream output of a given mutation
setwd(Path_Basecall_output_down)
write.csv(context_output_downstream_matrix, "Context_Output_downstream.csv") #full 4mer downstream output of a given mutation
