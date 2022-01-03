#Mutation rate Calculations

# This script generates the mutation rates for each replichore, and the chromosome of a given organim using the raw
# Base substitution count, the GWTC, and the number of generations of an organism


#mutation ratio calculation for each triplet site chromosome wide and replichore specific
MutationRatio_Left <- matrix(0L, nrow =16, ncol =4)
MutationRatio_Right <- matrix(0L, nrow =16, ncol =4)
MutationRatio_LeftRep <- matrix(0L, nrow =16, ncol =4)
MutationRatio_RightRep <- matrix(0L, nrow =16, ncol =4)
MutationRatio_Chromosome <- matrix(0L, nrow =16, ncol =4)

#Base Substitution mutation rate matrices for the 16x12 mutation rates at each nucleotide site
Base_MutLeft <- matrix(0L, nrow =16, ncol =12)
Base_MutRight <- matrix(0L, nrow =16, ncol =12)
Base_MutChrome <- matrix(0L, nrow =16, ncol =12)

#Base Substitution mutation Rate matrices for the 16x4 mutation matrices
Sung_MutLeft <- matrix(0L, nrow =16, ncol =4)
Sung_MutRight <- matrix(0L, nrow =16, ncol =4)
Sung_MutLeftRep <- matrix(0L, nrow =16, ncol =4)
Sung_MutRightRep <- matrix(0L, nrow =16, ncol =4)
Sung_MutChrome <- matrix(0L, nrow =16, ncol =4)

#Context Mutation rates for the Chromosome, and each replichore with respect to chromosome wide and replichore specific rates
Sung_ContextLeft <- matrix(0L, nrow =16, ncol =4)
Sung_ContextRight <- matrix(0L, nrow =16, ncol =4)
Sung_ContextLeftRep <- matrix(0L, nrow =16, ncol =4)
Sung_ContextRightRep <- matrix(0L, nrow =16, ncol =4)
Sung_ContextChrome <- matrix(0L, nrow =16, ncol =4)

cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
sung_cols <- c("T", "G", "C", "A")

rownames(Base_MutLeft) <- rows
rownames(Base_MutRight) <- rows
rownames(Base_MutChrome) <- rows

rownames(MutationRatio_Left) <- rows
rownames(MutationRatio_Right) <- rows
rownames(MutationRatio_LeftRep) <- rows
rownames(MutationRatio_RightRep) <- rows
rownames(MutationRatio_Chromosome) <- rows

rownames(Sung_MutLeft) <- rows
rownames(Sung_MutRight) <- rows
rownames(Sung_MutLeftRep) <- rows
rownames(Sung_MutRightRep) <- rows
rownames(Sung_MutChrome) <- rows


colnames(Base_MutLeft) <- cols
colnames(Base_MutRight) <- cols
colnames(Base_MutChrome) <- cols

colnames(MutationRatio_Left) <- sung_cols
colnames(MutationRatio_Right) <- sung_cols
colnames(MutationRatio_LeftRep) <- sung_cols
colnames(MutationRatio_RightRep) <- sung_cols
colnames(MutationRatio_Chromosome) <- sung_cols

colnames(Sung_MutChrome) <- sung_cols
colnames(Sung_MutRight) <- sung_cols
colnames(Sung_MutLeft) <- sung_cols
colnames(Sung_MutRightRep) <- sung_cols
colnames(Sung_MutLeftRep) <- sung_cols





num_triplets <- 64
full_col <- 1
#Calculate the Base Substitution Mutation Rate and Context Mutation Rate for the Entire Chromosome
for(row in 1:nrow(Sung_Matrix)) 
{
  for(col in 1:ncol(Sung_Matrix)) 
  {
    GWTC_count <- as.numeric(GWTC_Current[row, col]) #Chromosome GWTC for a given Codon CHECK THIS VALUE NEXT!!!
    MutationRatio_Chromosome[row, col] <- Sung_Matrix[row, col]/GWTC_count #mutant/codon site ratio
    Sung_MutChrome[row, col] <- Sung_Matrix[row, col]/(generations*GWTC_count*malines) #Base Substitution Mutation rate
    
    Base_MutChrome[row, full_col] <- Mut_Matrix_Total_Current[row, full_col]/(generations*GWTC_count*malines)
    #print(full_col+1)
    Base_MutChrome[row, full_col+1] <- Mut_Matrix_Total_Current[row, full_col+1]/(generations*GWTC_count*malines)
    #print(full_col+2)
    Base_MutChrome[row, full_col+2] <- Mut_Matrix_Total_Current[row, full_col+2]/(generations*GWTC_count*malines)
    full_col <- full_col+3
  }
  full_col <- 1 #reset for the next row calculations until all rows calculated.
}

full_col <-1
#Calculate the Mutation Rate and Context Mutation Rate for the Left Replichore
for(row in 1:nrow(Sung_Matrix_Left)) 
{
  for(col in 1:ncol(Sung_Matrix_Left)) 
  {
    #print(Mut_Matrix_Total_Current[row,])
    GWTC_count <- as.numeric(GWTC_Current[row, col]) #Chromosome wide GWTC
    GWTC_Left_Count <- as.numeric(GWTC_Left_Current[row, col]) #Left Replichore GWTC
    MutationRatio_Left[row, col] <- Sung_Matrix_Left[row, col]/GWTC_count #mutant/codon site ratio (Chromosome-Wide)
    MutationRatio_LeftRep[row, col] <- Sung_Matrix_Left[row, col]/GWTC_Left_Count #mutant/codon site ratio (Left Replichore Specific)
    Sung_MutLeft[row, col] <- Sung_Matrix_Left[row, col]/(generations*GWTC_count*malines) #Base Substitution Mutation Rate (Left/chromosome codon GWTC)
    Sung_MutLeftRep[row, col] <- Sung_Matrix_Left[row, col]/(generations*GWTC_Left_Count*malines) #Base Substitution Mutation Rate (Left/Left Rep codon GWTC)

    Base_MutLeft[row, full_col] <- Mut_Matrix_Left_Current[row, full_col]/(generations*GWTC_Left_Count*malines)
    Base_MutLeft[row, full_col+1] <- Mut_Matrix_Left_Current[row, full_col+1]/(generations*GWTC_Left_Count*malines)
    Base_MutLeft[row, full_col+2] <- Mut_Matrix_Left_Current[row, full_col+2]/(generations*GWTC_Left_Count*malines)
    full_col <- full_col+3
    }
  full_col <- 1
}

full_col <- 1
#Calculate the Mutation Rate and Context Mutation Rate for the Left Replichore
for(row in 1:nrow(Sung_Matrix_Right)) 
{
  for(col in 1:ncol(Sung_Matrix_Right)) 
  {
    #print(Mut_Matrix_Total_Current[row,])
    GWTC_count <- as.numeric(GWTC_Current[row, col])
    GWTC_Right_Count <- as.numeric(GWTC_Right_Current[row, col]) #Right Replichore GWTC
    MutationRatio_Right[row, col] <- Sung_Matrix_Right[row, col]/GWTC_count #mutant/codon site ratio (Right,Chromosome-Wide)
    MutationRatio_RightRep[row, col] <- Sung_Matrix_Right[row, col]/GWTC_Right_Count #mutant/codon site ratio (Left Replichore Specific)
    Sung_MutRight[row, col] <- Sung_Matrix_Right[row, col]/(generations*GWTC_count*malines) #Base Substitution Mutation Rate (Right/Chromosome codon GWTC)
    Sung_MutRightRep[row, col] <- Sung_Matrix_Right[row, col]/(generations*GWTC_Right_Count*malines) #Base Substitution Mutation Rate (Right/Right Rep codon GWTC)
   
    Base_MutRight[row, full_col] <- Mut_Matrix_Right_Current[row, full_col]/(generations*GWTC_Right_Count*malines)
    Base_MutRight[row, full_col+1] <- Mut_Matrix_Right_Current[row, full_col+1]/(generations*GWTC_Right_Count*malines)
    Base_MutRight[row, full_col+2] <- Mut_Matrix_Right_Current[row, full_col+2]/(generations*GWTC_Right_Count*malines)
    full_col <- full_col +3
  }
  full_col <- 1
}
setwd(Path_output)
#Path_MutationRates_output <- paste(Path_text_output_current, "/MutationRates_Output", sep = "")
#dir.create(Path_MutationRates_output)
setwd(Path_MutationRates_output)

write.csv(MutationRatio_Left, "MutationRatio_Left.csv")
write.csv(MutationRatio_Right, "MutationRatio_Right.csv")
write.csv(MutationRatio_Chromosome, "MutationRatio_Chromosome.csv")

write.csv(MutationRatio_LeftRep, "MutationRatio_LeftRep.csv")
write.csv(MutationRatio_RightRep, "MutationRatio_RightRep.csv")

write.csv(Base_MutChrome, "BaseMutationRate_Chrome.csv")
write.csv(Base_MutLeft, "BaseMutationRate_Left.csv")
write.csv(Base_MutRight, "BaseMutationRate_Right.csv")

write.csv(Sung_MutChrome, "MutationRate_SungChrome.csv")
write.csv(Sung_MutLeft, "MutationRate_SungLeft.csv")
write.csv(Sung_MutRight, "MutationRate_SungRight.csv")
write.csv(Sung_MutLeftRep, "MutationRate_SungLeft_Replichore.csv")
write.csv(Sung_MutRightRep, "MutationRate_SungRight_Replichore.csv")


setwd(Path_correlate_Chromosome_current)
output_name <- paste(organism, "Context_MutChrome.csv", sep = "_")
write.csv(Sung_MutChrome, output_name)

setwd(Path_correlate_Lcore_current)
output_name <- paste(organism, "Context_MutLeft.csv", sep = "_")
write.csv(Sung_MutLeft, output_name)
setwd(Path_correlate_Lcore_Replichore_current) 
output_name <- paste(organism, "Context_MutLeft_Replichore.csv", sep = "_")
write.csv(Sung_MutLeftRep, output_name)

setwd(Path_correlate_Rcore_current)
output_name <- paste(organism, "Context_MutRight.csv", sep = "_")
write.csv(Sung_MutRight, output_name)
setwd(Path_correlate_Rcore_Replichore_current) 
output_name <- paste(organism, "Context_MutRight_Replichore.csv", sep = "_")
write.csv(Sung_MutRightRep, output_name)
