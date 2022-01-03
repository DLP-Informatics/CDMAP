#Sung Style output table

rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")
cols <- c("T", "G", "C", "A")

Sung_Matrix <- matrix(0L, nrow =16, ncol =4)
rownames(Sung_Matrix) <- rows
colnames(Sung_Matrix) <- cols

i <- 1
while(i <= 16 )
{
  Sung_Matrix[i,1] <- Mut_Matrix_Total[i,1] + Mut_Matrix_Total[i,2] + Mut_Matrix_Total[i,3]
  Sung_Matrix[i,2] <- Mut_Matrix_Total[i,4] + Mut_Matrix_Total[i,5] + Mut_Matrix_Total[i,6]
  Sung_Matrix[i,3] <- Mut_Matrix_Total[i,7] + Mut_Matrix_Total[i,8] + Mut_Matrix_Total[i,9]
  Sung_Matrix[i,4] <- Mut_Matrix_Total[i,10] + Mut_Matrix_Total[i,11] + Mut_Matrix_Total[i,12]
  
  i <- i+1
}

Comp_Sung_Matrix <- matrix(0L, nrow =16, ncol =4)
rownames(Comp_Sung_Matrix) <- rows
colnames(Comp_Sung_Matrix) <- cols

i <- 1
while(i <= 16 )
{
  Comp_Sung_Matrix[i,1] <- Comp_Mut_Matrix_Total[i,1] + Comp_Mut_Matrix_Total[i,2] + Comp_Mut_Matrix_Total[i,3]
  Comp_Sung_Matrix[i,2] <- Comp_Mut_Matrix_Total[i,4] + Comp_Mut_Matrix_Total[i,5] + Comp_Mut_Matrix_Total[i,6]
  Comp_Sung_Matrix[i,3] <- Comp_Mut_Matrix_Total[i,7] + Comp_Mut_Matrix_Total[i,8] + Comp_Mut_Matrix_Total[i,9]
  Comp_Sung_Matrix[i,4] <- Comp_Mut_Matrix_Total[i,10] + Comp_Mut_Matrix_Total[i,11] + Comp_Mut_Matrix_Total[i,12]
  
  i <- i+1
}

setwd(Path_BaseSubstitution_output_current)
write.csv(Sung_Matrix, "Sung_Chromosome.csv")
write.csv(Comp_Sung_Matrix, "RevCompliment_Sung_Chromosome.csv")

#setwd(Path_correlate_Chromosome_current)
#output_name <- paste(organism, "_Sung_Chromosome.csv", sep = "")
#write.csv(Sung_Matrix, output_name)

#setwd(Path_correlate_Chromosome)
#output_name <- paste(organism, "_RevCompliment_Sung_Chromosome.csv", sep = "")
#write.csv(Comp_Sung_Matrix, output_name)