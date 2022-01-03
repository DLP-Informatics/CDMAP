#Config Dump File
# Dump file of information to allow users to allocate information regarding genome and replichore size, ORI/TERM positions
# and other relevant information for downstream visualization usage.

setwd(Path_output_organism)

report_suffix <- "DataSummary.txt"
Report_Name <- paste(organism, report_suffix, sep = "_")
fasta_file_name <- paste(organism, ".fasta", sep = "" )
detach("package:vcfR", unload = TRUE)
write.fasta(RefFile, organism, fasta_file_name)


#Generating a Data report of all pertinent Information about the Organism post-pipeline
print(paste("Generating Post Analysis Report of", organism, sep = " "))

sink(file = Report_Name, type = "output")

print("ORI, TERM, and Replichore Information: ")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print(paste("Name of organism analyzed: ", organism, sep = ''))
print(paste("Number of Mutation Accumulation Generations carried out: ", generations, sep = ''))
print(paste("Number of Mutation Accumulation Lines Ran: ", malines, sep = ''))
  
print(paste("Number of mutation accumulation lines during experiment: ", malines, sep = ''))
print(paste("Number of experimental generations during MA expiriment: ", generations, sep = ''))
print(paste("Reference Sequence length (Genome Size): ", len_refseq, sep = ''))
print(paste("Left Replichore Sequence length: ", sum(GWTC_Left)-1, sep = '')) # n-1 since last triplet is double counted in each replichore
print(paste("Right Replichore Sequence length: ", sum(GWTC_Right)-1, sep = '')) #see
cat("\n")
print(paste("ORI position value in oriloc: ", ori_pos, sep = ''))
print(paste("ORI KB adjusted position value in oriloc: ", ori_bp, sep = ''))
cat("\n")
print(paste("Terminus position value in oriloc: ", term_pos, sep = ''))
print(paste("Terminus adjusted position value in oriloc: ", term_bp, sep = ''))
cat("\n")
print(paste("Sum total length of all Coding regions: ", CodingRegionSize, sep = ''))
print(paste("Percent of Genome that's a Coding region: ", CodingRegionPercent, sep = ''))
print(paste("Coding Region GC Content: ", CodonGCcontent, sep = ''))
print(paste("Coding Region AT Content: ", CodonATcontent, sep = ''))
cat("\n")
print(paste("Total number of Mutants detected in BaseCall file: ", length(MutBaseCalls[,1]), sep = ''))
print(paste("Number of Mutants Occuring in the Left Replichore: ", length(query_Lcore_full[,1]), sep = ''))
print(paste("Number of Mutants Occuring in the Right Replichore: ", length(query_Rcore_full[,1]), sep = ''))
sum_mut <- length(query_Lcore_full[,1])+length(query_Rcore_full[,1])
print(paste("Total number of processed mutants after analysis: ", sum_mut, sep = ''))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")



sink()
Path_correlate_repoGC <- paste(Path_home, correlate_repoGC, sep = "")
if(!(dir.exists(Path_correlate_repoGC)))
{
  dir.create(Path_correlate_repoGC)
}
setwd(Path_correlate_repoGC)
write.fasta(RefFile, organism, fasta_file_name)
