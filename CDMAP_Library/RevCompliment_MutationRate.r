#Mutation rate Calculations

# This script generates the mutation rates for each replichore, and the chromosome of a given organim using the raw
# Base substitution count, the GWTC, and the number of generations of an organism

Comp_MutationRatio_Left <- matrix(0L, nrow =16, ncol =4)
Comp_MutationRatio_Right <- matrix(0L, nrow =16, ncol =4)
Comp_MutationRatio_LeftRep <- matrix(0L, nrow =16, ncol =4)
Comp_MutationRatio_RightRep <- matrix(0L, nrow =16, ncol =4)
Comp_MutationRatio_Chromosome <- matrix(0L, nrow =16, ncol =4)


Comp_Sung_MutLeft <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_MutRight <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_MutLeftRep <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_MutRightRep <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_MutChrome <- matrix(0L, nrow =16, ncol =4)

Comp_Sung_ContextLeft <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_ContextRight <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_ContextLeftRep <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_ContextRightRep <- matrix(0L, nrow =16, ncol =4)
Comp_Sung_ContextChrome <- matrix(0L, nrow =16, ncol =4)

cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
sung_cols <- c("T", "G", "C", "A")

rownames(Comp_MutationRatio_Left) <- rows
rownames(Comp_MutationRatio_Right) <- rows
rownames(Comp_MutationRatio_LeftRep) <- rows
rownames(Comp_MutationRatio_RightRep) <- rows
rownames(Comp_MutationRatio_Chromosome) <- rows

rownames(Comp_Sung_MutLeft) <- rows
rownames(Comp_Sung_MutRight) <- rows
rownames(Comp_Sung_MutLeftRep) <- rows
rownames(Comp_Sung_MutRightRep) <- rows
rownames(Comp_Sung_MutChrome) <- rows

rownames(Comp_Sung_ContextLeft) <- rows
rownames(Comp_Sung_ContextRight) <- rows
rownames(Comp_Sung_ContextLeftRep) <- rows
rownames(Comp_Sung_ContextRightRep) <- rows
rownames(Comp_Sung_ContextChrome) <- rows

colnames(Comp_MutationRatio_Left) <- sung_cols
colnames(Comp_MutationRatio_Right) <- sung_cols
colnames(Comp_MutationRatio_LeftRep) <- sung_cols
colnames(Comp_MutationRatio_RightRep) <- sung_cols
colnames(Comp_MutationRatio_Chromosome) <- sung_cols

colnames(Comp_Sung_MutChrome) <- sung_cols
colnames(Comp_Sung_MutRight) <- sung_cols
colnames(Comp_Sung_MutLeft) <- sung_cols
colnames(Comp_Sung_MutRightRep) <- sung_cols
colnames(Comp_Sung_MutLeftRep) <- sung_cols

colnames(Comp_Sung_ContextLeft) <- sung_cols
colnames(Comp_Sung_ContextRight) <- sung_cols
colnames(Comp_Sung_ContextLeftRep) <- sung_cols
colnames(Comp_Sung_ContextRightRep) <- sung_cols
colnames(Comp_Sung_ContextChrome) <- sung_cols



generations <- as.numeric(generations)
num_triplets <- 64

#Calculate the Base Substitution Mutation Rate and Context Mutation Rate for the Entire Chromosome
for(row in 1:nrow(Comp_Sung_Matrix)) 
{
  for(col in 1:ncol(Comp_Sung_Matrix)) 
  {
    Rev_GWTC_count <- as.numeric(Rev_GWTC_Current[row, col]) #Chromosome GWTC for a given Codon
    Comp_MutationRatio_Chromosome[row, col] <- Comp_Sung_Matrix[row, col]/Rev_GWTC_count #mutant/codon site ratio
    Comp_Sung_MutChrome[row, col] <- Comp_Sung_Matrix[row, col]/(generations*Rev_GWTC_count) #Base Substitution Mutation rate
  }
}

#Calculate the Mutation Rate and Context Mutation Rate for the Left Replichore
for(row in 1:nrow(Comp_Sung_Matrix_Left)) 
{
  for(col in 1:ncol(Comp_Sung_Matrix_Left)) 
  {
    #print(Mut_Matrix_Total[row,])
    Rev_GWTC_count <- as.numeric(Rev_GWTC_Current[row, col]) #Chromosome wide GWTC
    Rev_GWTC_Left_Count <- as.numeric(Rev_GWTC_Left_Current[row, col]) #Left Replichore GWTC
    Comp_MutationRatio_Left[row, col] <- Comp_Sung_Matrix_Left[row, col]/Rev_GWTC_count #mutant/codon site ratio (Chromosome-Wide)
    Comp_MutationRatio_LeftRep[row, col] <- Comp_Sung_Matrix_Left[row, col]/Rev_GWTC_Left_Count #mutant/codon site ratio (Left Replichore Specific)
    Comp_Sung_MutLeft[row, col] <- Comp_Sung_Matrix_Left[row, col]/(generations*Rev_GWTC_count) #Base Substitution Mutation Rate (Left/chromosome codon GWTC)
    Comp_Sung_MutLeftRep[row, col] <- Comp_Sung_Matrix_Left[row, col]/(generations*Rev_GWTC_Left_Count) #Base Substitution Mutation Rate (Left/Left Rep codon GWTC)
  }
  
}

#Calculate the Mutation Rate and Context Mutation Rate for the Left Replichore
for(row in 1:nrow(Comp_Sung_Matrix_Right)) 
{
  for(col in 1:ncol(Comp_Sung_Matrix_Right)) 
  {
    #print(Mut_Matrix_Total[row,])
    Rev_GWTC_count <- as.numeric(Rev_GWTC_Current[row, col])
    Rev_GWTC_Right_Count <- as.numeric(Rev_GWTC_Right_Current[row, col]) #Right Replichore GWTC
    Comp_MutationRatio_Right[row, col] <- Comp_Sung_Matrix_Right[row, col]/GWTC_count #mutant/codon site ratio (Right,Chromosome-Wide)
    Comp_MutationRatio_RightRep[row, col] <- Comp_Sung_Matrix_Right[row, col]/Rev_GWTC_Right_Count #mutant/codon site ratio (Left Replichore Specific)
    Comp_Sung_MutRight[row, col] <- Comp_Sung_Matrix_Right[row, col]/(generations*Rev_GWTC_count) #Base Substitution Mutation Rate (Right/Chromosome codon GWTC)
    Comp_Sung_MutRightRep[row, col] <- Comp_Sung_Matrix_Right[row, col]/(generations*Rev_GWTC_Right_Count) #Base Substitution Mutation Rate (Right/Right Rep codon GWTC)
  }
}
setwd(Path_MutationRates_output)
write.csv(Comp_MutationRatio_Left, "RevCompliment_MutationRatio_Left.csv")
write.csv(Comp_MutationRatio_Right, "RevCompliment_MutationRatio_Right.csv")
write.csv(Comp_MutationRatio_Chromosome, "RevCompliment_MutationRatio_Chromosome.csv")
write.csv(Comp_MutationRatio_LeftRep, "RevCompliment_MutationRatio_LeftRep.csv")
write.csv(Comp_MutationRatio_RightRep, "RevCompliment_MutationRatio_RightRep.csv")

write.csv(Comp_Sung_MutChrome, "RevCompliment_MutationRate_SungChrome.csv")
write.csv(Comp_Sung_MutLeft, "RevCompliment_MutationRate_SungLeft.csv")
write.csv(Comp_Sung_MutRight, "RevCompliment_MutationRate_SungRight.csv")
write.csv(Comp_Sung_MutLeftRep, "RevCompliment_MutationRate_SungLeft_Replichore.csv")
write.csv(Comp_Sung_MutRightRep, "RevCompliment_MutationRate_SungRight_Replichore.csv")


setwd(Path_correlate_RevComp_Chromosome_current)
output_name <- paste(organism, "RevCompliment_Context_MutChrome.csv", sep = "_")
write.csv(Comp_Sung_MutChrome, output_name)


setwd(Path_correlate_RevComp_Lcore_current)
output_name <- paste(organism, "RevCompliment_Context_MutLeft.csv", sep = "_")
write.csv(Comp_Sung_MutLeft, output_name)
setwd(Path_correlate_RevComp_Lcore_Replichore_current) 
output_name <- paste(organism, "RevCompliment_Context_MutLeft_Replichore.csv", sep = "_")
write.csv(Comp_Sung_MutLeftRep, output_name)

setwd(Path_correlate_RevComp_Rcore_current)
output_name <- paste(organism, "RevCompliment_Context_MutRight.csv", sep = "_")
write.csv(Comp_Sung_MutRight, output_name)
setwd(Path_correlate_RevComp_Rcore_Replichore_current) 
output_name <- paste(organism, "RevCompliment_Context_MutRight_Replichore.csv", sep = "_")
write.csv(Comp_Sung_MutRightRep, output_name)

