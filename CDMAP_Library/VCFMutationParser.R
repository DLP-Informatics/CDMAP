#VCF file CDMAP import
#written by David L Patton

library("vcfR")
library("pracma")
library("beepr")
library("seqinr")
library("BiocManager")
library("tidyverse")
library("genbankr")
#Objective: this script is designed to read in non-modified VCF mutation call files, and to extract the 
# position of point mutations for processing context mutation rates in CDMAP. In addition, this script will
# record indel mutations (but does not compute their mutation rates) for ease of use downstream.

#debug paths
VCFfile <- "/Users/triadge/Documents/BINF_Coursework/BINF6203_PROJECT_IMPORTANT/Genome_Mapping/Bmmr1/Bmmr59_1.vcf"
VCFOutputDir <- "/Users/triadge/Documents/BINF_Coursework/BINF6203_PROJECT_IMPORTANT/Genome_Mapping/Bmmr1"
setwd(VCFOutputDir)

#VCFfile <- Path_InFile
#VCFOutputDir <- Path_output_organism


VCFobj <- read.vcfR(VCFfile)
VCFobjConv <- vcfR2tidy(VCFobj)
VCFfix <- VCFobjConv$fix #gives us variant call data
RefData <-VCFobjConv$fix$REF
AltData <- VCFobjConv$fix$ALT
PosData <<- which(AltData != is.na(AltData))

context_col_names <- c("postion", "original", "mutant")
output_vcf_matrix <- matrix( nrow =0, ncol = 3)
colnames(output_vcf_matrix) <- context_col_names



context_col_names <- c("postion", "original", "mutant")
output_indel_matrix <- matrix( nrow =0, ncol = 3)
colnames(output_indel_matrix) <- context_col_names

#conditions we need to check for in our loop
#loop length? length(AltData)
# is the variant nucleotide NA? is.na(AltData[i])
# is the variant or reference an indel? (nchar(RefData[i]) == 1 && nchar(AltData[i]) == 1)
# if the conditions are met, then we record the position, original, and variant nucleotide.
i <-1


for(i in 1:length(PosData))
{
  if((nchar(RefData[PosData[i]]) == 1 && nchar(AltData[PosData[i]]) == 1))
  {
    position <- PosData[i]
    original <- RefData[PosData[i]]
    mutant <- AltData[PosData[i]]
    newrow <- c(position, original, mutant)
    output_vcf_matrix <- rbind(output_vcf_matrix, newrow)
    next
  }
  if(nchar(RefData[PosData[i]]) > 1 | nchar(AltData[PosData[i]]) > 1) #correct this line VCF works, indel does not
  {
    position <- PosData[i]
    original <- RefData[PosData[i]]
    mutant <- AltData[PosData[i]]
    newrow <- c(position, original, mutant)
    output_indel_matrix <- rbind(output_indel_matrix, newrow)
  }
  #i <- i+1
}

#Path_Basecall_output <- paste(Path_text_output, "/BaseCall_Output", sep = "")
#dir.create(Path_Basecall_output)
#setwd(Path_Basecall_output)

rownames(output_vcf_matrix) <- NULL
write.csv(output_vcf_matrix, "VCF_mutations.csv")

rownames(output_indel_matrix) <- NULL
output_indel_matrix <- data.frame(output_indel_matrix)
indelbyoriginal <- output_indel_matrix[order(output_indel_matrix$original),]
indelbymutant <- output_indel_matrix[order(output_indel_matrix$mutant),]

write.csv(output_indel_matrix, "Indel_Mutations_Position.csv")
write.csv(indelbyoriginal, "Indel_Mutations_Original.csv")
write.csv(indelbymutant, "Indel_Mutations_Mutant.csv")

MutBaseCalls <- output_vcf_matrix  #read in the basecall file
colnames(MutBaseCalls) <- c("position", "original", "mutant")
MutBaseCalls <- data.frame(MutBaseCalls)
