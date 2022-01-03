cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")
Path_BaseSubstitution_output <- paste(Path_text_output, "/BaseSubstitution_Output", sep = "")
dir.create(Path_BaseSubstitution_output)

Mut_Matrix <- matrix(0L, nrow =16, ncol =12)
rownames(Mut_Matrix) <- rows
colnames(Mut_Matrix) <- cols

Mut_Matrix_Left <- Mut_Matrix #Left replichore
Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome

Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome


#Populate the output matrix with left replichore mutants
#================
i <- 2

while(i <= length(query_Lcore_full[,1]))
{
  #print(i)
  #print(length(query_Lcore_full))
  #print(query_Lcore_full[i,2])
  #print(query_Lcore_full[i,5])
  row_switch1 <- query_Lcore_full[i,2]
  row_switch2 <- query_Lcore_full[i,5]
  rowswitch <- paste(row_switch1, row_switch2, sep = '')
  
  #print(query_Lcore_full[i,3])
  #print(query_Lcore_full[i,4])
  col_switch1 <- query_Lcore_full[i,3]
  col_switch2 <- query_Lcore_full[i,4]
  colswitch <- paste(col_switch1, col_switch2, sep = '')
  
  rowiter <- switch(rowswitch,
                    "TT" = 1,
                    "TG" = 2,
                    "TC" = 3,
                    "TA" = 4,
                    "GT" = 5,
                    "GG" = 6,
                    "GC" = 7,
                    "GA" = 8,
                    "CT" = 9,
                    "CG" = 10,
                    "CC" = 11,
                    "CA" = 12,
                    "AT" = 13,
                    "AG" = 14,
                    "AC" = 15,
                    "AA" = 16
  )
  
  row_comp <- switch(rowswitch,
                     "TT" = 'AA',
                     "TG" = 'CA',
                     "TC" = 'GA',
                     "TA" = 'TA',
                     "GT" = 'AC',
                     "GG" = 'CC',
                     "GC" = 'GC',
                     "GA" = 'TC',
                     "CT" = 'AG',
                     "CG" = 'CG',
                     "CC" = 'GG',
                     "CA" = 'TG',
                     "AT" = 'AT',
                     "AG" = 'CT',
                     "AC" = 'GT',
                     "AA" = 'TT'
  )
  
  comp_rowiter <- switch(row_comp,
                         "TT" = 1,
                         "TG" = 2,
                         "TC" = 3,
                         "TA" = 4,
                         "GT" = 5,
                         "GG" = 6,
                         "GC" = 7,
                         "GA" = 8,
                         "CT" = 9,
                         "CG" = 10,
                         "CC" = 11,
                         "CA" = 12,
                         "AT" = 13,
                         "AG" = 14,
                         "AC" = 15,
                         "AA" = 16
  )
  
  coliter <- switch(colswitch,
                    "TG" = 1,
                    "TC" = 2,
                    "TA" = 3,
                    "GT" = 4,
                    "GC" = 5,
                    "GA" = 6,
                    "CT" = 7,
                    "CG" = 8,
                    "CA" = 9,
                    "AT" = 10,
                    "AG" = 11,
                    "AC" = 12
  )
  
  col_comp <- switch(colswitch,
                     "TG" = 'AC',
                     "TC" = 'AG',
                     "TA" = 'AT',
                     "GT" = 'CA',
                     "GC" = 'CG',
                     "GA" = 'CT',
                     "CT" = 'GA',
                     "CG" = 'GC',
                     "CA" = 'GT',
                     "AT" = 'TA',
                     "AG" = 'TC',
                     "AC" = 'TG'
  )
  
  comp_coliter <- switch(col_comp,
                         "TG" = 1,
                         "TC" = 2,
                         "TA" = 3,
                         "GT" = 4,
                         "GC" = 5,
                         "GA" = 6,
                         "CT" = 7,
                         "CG" = 8,
                         "CA" = 9,
                         "AT" = 10,
                         "AG" = 11,
                         "AC" = 12
  )
  
  increment <- as.integer(Mut_Matrix_Left[rowiter,coliter])
  increment <- increment + 1
  Mut_Matrix_Left[rowiter, coliter] <- increment       #record in Left Replichore Matrix
  
  comp_increment <- as.integer(Comp_Mut_Matrix_Left[comp_rowiter,comp_coliter])
  comp_increment <- comp_increment + 1
  Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment       #record in Left Replichore Reverse Compliment Matrix
  
  
  i <- i+1
  
}
#================

## Populate the output matrix with right replichore mutants
#================
i <- 1

while(i <= length(query_Rcore_full[,2]))
{
  row_switch1 <- query_Rcore_full[i,2]
  row_switch2 <- query_Rcore_full[i,5]
  rowswitch <- paste(row_switch1, row_switch2, sep = '')
  
  col_switch1 <- query_Rcore_full[i,3]
  col_switch2 <- query_Rcore_full[i,4]
  colswitch <- paste(col_switch1, col_switch2, sep = '')
  
  row_comp <- switch(rowswitch,
                     "TT" = 'AA',
                     "TG" = 'CA',
                     "TC" = 'GA',
                     "TA" = 'TA',
                     "GT" = 'AC',
                     "GG" = 'CC',
                     "GC" = 'GC',
                     "GA" = 'TC',
                     "CT" = 'AG',
                     "CG" = 'CG',
                     "CC" = 'GG',
                     "CA" = 'TG',
                     "AT" = 'AT',
                     "AG" = 'CT',
                     "AC" = 'GT',
                     "AA" = 'TT'
  )
  
  rowiter <- switch(rowswitch,
                    "TT" = 1,
                    "TG" = 2,
                    "TC" = 3,
                    "TA" = 4,
                    "GT" = 5,
                    "GG" = 6,
                    "GC" = 7,
                    "GA" = 8,
                    "CT" = 9,
                    "CG" = 10,
                    "CC" = 11,
                    "CA" = 12,
                    "AT" = 13,
                    "AG" = 14,
                    "AC" = 15,
                    "AA" = 16
  )
  
  comp_rowiter <- switch(row_comp,
                         "TT" = 1,
                         "TG" = 2,
                         "TC" = 3,
                         "TA" = 4,
                         "GT" = 5,
                         "GG" = 6,
                         "GC" = 7,
                         "GA" = 8,
                         "CT" = 9,
                         "CG" = 10,
                         "CC" = 11,
                         "CA" = 12,
                         "AT" = 13,
                         "AG" = 14,
                         "AC" = 15,
                         "AA" = 16
  )
  
  col_comp <- switch(colswitch,
                     "TG" = 'AC',
                     "TC" = 'AG',
                     "TA" = 'AT',
                     "GT" = 'CA',
                     "GC" = 'CG',
                     "GA" = 'CT',
                     "CT" = 'GA',
                     "CG" = 'GC',
                     "CA" = 'GT',
                     "AT" = 'TA',
                     "AG" = 'TC',
                     "AC" = 'TG'
  )
  
  coliter <- switch(colswitch,
                    "TG" = 1,
                    "TC" = 2,
                    "TA" = 3,
                    "GT" = 4,
                    "GC" = 5,
                    "GA" = 6,
                    "CT" = 7,
                    "CG" = 8,
                    "CA" = 9,
                    "AT" = 10,
                    "AG" = 11,
                    "AC" = 12
  )
  
  comp_coliter <- switch(col_comp,
                         "TG" = 1,
                         "TC" = 2,
                         "TA" = 3,
                         "GT" = 4,
                         "GC" = 5,
                         "GA" = 6,
                         "CT" = 7,
                         "CG" = 8,
                         "CA" = 9,
                         "AT" = 10,
                         "AG" = 11,
                         "AC" = 12
  )
  
  increment <- as.integer(Mut_Matrix_Right[rowiter,coliter])
  increment <- increment + 1
  Mut_Matrix_Right[rowiter, coliter] <- increment #record on right replichore matrix
  
  comp_increment <- as.integer(Comp_Mut_Matrix_Right[comp_rowiter,comp_coliter])
  comp_increment <- comp_increment + 1
  Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment       #record in Left Replichore Reverse Compliment Matrix
  
  
  
  i=i+1
}
#================

Mut_Matrix_Total <- Comp_Mut_Matrix_Left + Mut_Matrix_Right
Comp_Mut_Matrix_Total <- Mut_Matrix_Left + Comp_Mut_Matrix_Right


setwd(Path_BaseSubstitution_output)
write.csv(Mut_Matrix_Total, "Context_Mut_Clockwise_Chromosome.csv")
write.csv(Mut_Matrix_Left, "Context_Mutations_Left.csv")
write.csv(Mut_Matrix_Right, "Context_Mutations_Right.csv")


write.csv(Comp_Mut_Matrix_Total, "RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Comp_Mut_Matrix_Left, "RevCompliment_Context_Mutations_Left.csv")
write.csv(Comp_Mut_Matrix_Right, "RevCompliment_Context_Mutations_Right.csv")


#setwd(Path_correlate_Chromosome)
#output_name <- paste(organism, "_Mut_Clockwise_Chromosome.csv", sep = "")
#write.csv(Mut_Matrix_Total, output_name)
#output_name <- paste(organism, "_Mut_Counterclockwise_Chromosome.csv", sep = "")
#write.csv(Mut_Matrix_Total, output_name)

#setwd(Path_correlate_Lcore)
#output_name <- paste(organism, "_Mutations_Lcore.csv", sep = "")
#write.csv(Mut_Matrix_Left, output_name)

#setwd(Path_correlate_Rcore)
#output_name <- paste(organism, "_Mutations_Rcore.csv", sep = "")
#write.csv(Mut_Matrix_Right, output_name)
