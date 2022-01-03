#DevNote: Figure out and debug rate calcuation for covid19!!!


cols <- c("T→G", "T→C", "T→A",	"G→T",	"G→C",	"G→A",	"C→T",	"C→G",	"C→A",	"A→T",	"A→G",	"A→C")
rows <- c("T[X-->Y]T", "T[X-->Y]G","T[X-->Y]C","T[X-->Y]A","G[X-->Y]T","G[X-->Y]G","G[X-->Y]C","G[X-->Y]A",
          "C[X-->Y]T","C[X-->Y]G","C[X-->Y]C","C[X-->Y]A", "A[X-->Y]T","A[X-->Y]G", "A[X-->Y]C","A[X-->Y]A")

Path_BaseSubstitution_output <- paste(Path_text_output, "/BaseSubstitution_Output", sep = "")
dir.create(Path_BaseSubstitution_output)
Path_BaseSubstitution_output <- paste(Path_BaseSubstitution_output, "/Triplet", sep = "")
dir.create(Path_BaseSubstitution_output)
Path_BaseSubstitution_outputGC3C <- paste(Path_text_output, "/BaseSubstitution_Output/GC3C", sep = "")
dir.create(Path_BaseSubstitution_outputGC3C)

Path_BaseSubstitution_output_up <- paste(Path_text_output_up, "/BaseSubstitution_Output", sep = "")
dir.create(Path_BaseSubstitution_output_up)
Path_BaseSubstitution_output_upA <- paste(Path_BaseSubstitution_output_up, "/A", sep = "")
dir.create(Path_BaseSubstitution_output_upA)
Path_BaseSubstitution_output_upC <- paste(Path_BaseSubstitution_output_up, "/C", sep = "")
dir.create(Path_BaseSubstitution_output_upC)
Path_BaseSubstitution_output_upT <- paste(Path_BaseSubstitution_output_up, "/T", sep = "")
dir.create(Path_BaseSubstitution_output_upT)
Path_BaseSubstitution_output_upG <- paste(Path_BaseSubstitution_output_up, "/G", sep = "")
dir.create(Path_BaseSubstitution_output_upG)

Path_BaseSubstitution_output_down <- paste(Path_text_output_down, "/BaseSubstitution_Output", sep = "")
dir.create(Path_BaseSubstitution_output_down)
Path_BaseSubstitution_output_downA <- paste(Path_BaseSubstitution_output_down, "/A", sep = "")
dir.create(Path_BaseSubstitution_output_downA)
Path_BaseSubstitution_output_downC <- paste(Path_BaseSubstitution_output_down, "/C", sep = "")
dir.create(Path_BaseSubstitution_output_downC)
Path_BaseSubstitution_output_downT <- paste(Path_BaseSubstitution_output_down, "/T", sep = "")
dir.create(Path_BaseSubstitution_output_downT)
Path_BaseSubstitution_output_downG <- paste(Path_BaseSubstitution_output_down, "/G", sep = "")
dir.create(Path_BaseSubstitution_output_downG)

Mut_Matrix <- matrix(0L, nrow =16, ncol =12)
rownames(Mut_Matrix) <- rows
colnames(Mut_Matrix) <- cols

#=== Base GC3C mutation recordings and complimentary matrix == #
Mut_Matrix_Left <- Mut_Matrix #Left replichore
Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome

#Genome coding 3mer matrices#
GC3C_Mut_Matrix_Left <- Mut_Matrix #Left replichore
GC3C_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
GC3C_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome

Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
#=========================#

#========Upstream left, right, and chrome recording matrices with compliments ==== #
Aup_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Aup_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Aup_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Cup_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Cup_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Cup_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Gup_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Gup_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Gup_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Tup_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Tup_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Tup_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome

Aup_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Aup_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Aup_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Cup_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Cup_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Cup_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Gup_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Gup_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Gup_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Tup_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Tup_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Tup_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
#=====================================#

#======== Downstream Left, Right, and chrome matrices with compliments ====#
Adown_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Adown_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Adown_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Cdown_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Cdown_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Cdown_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Gdown_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Gdown_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Gdown_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Tdown_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Tdown_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Tdown_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome

Adown_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Adown_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Adown_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Cdown_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Cdown_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Cdown_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Gdown_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Gdown_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Gdown_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
Tdown_Comp_Mut_Matrix_Left <- Mut_Matrix #Left replichore
Tdown_Comp_Mut_Matrix_Right <- Mut_Matrix #Right Replichore
Tdown_Comp_Mut_Matrix_Total <- Mut_Matrix #Mutations accross entire Chromosome
#=====================#



Lcore_Frame <- data.frame(query_Lcore_full)
Rcore_Frame <- data.frame(query_Rcore_full)


#Populate the output matrix with left replichore mutants
#================
i <<- 1

while(i < length(Lcore_Frame$Position))
{
  codeRegion <- FALSE
  #print(paste("iteratior: ", i, sep = ''))
  #print(paste(" mutation row: ",query_Lcore_full[i,], sep = ""))
  #print(paste("Left neighbor: ", query_Lcore_full[i,3], sep = ''))
  #print(paste("Right neighbor: ", query_Lcore_full[i,6], sep = ''))
  position <- as.integer(Lcore_Frame$Position[i])
  row_switch1 <- Lcore_Frame$Left1[i]
  row_switch2 <- Lcore_Frame$Right1[i]
  rowswitch <- paste(row_switch1, row_switch2, sep = '')
  
  #print(query_Lcore_full[i,4])
  #print(query_Lcore_full[i,5])
  col_switch1 <- Lcore_Frame$Original[i]
  col_switch2 <- Lcore_Frame$Mutation[i]
  colswitch <- paste(col_switch1, col_switch2, sep = '')
  
  upstream <- Lcore_Frame$Left2[i]
  downstream <- Lcore_Frame$Right2[i]
  
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
  

    #print(i)
  j <- 1
  for(j in 1:length(startindex))
  {
    start <- as.numeric(startindex[j])
    end <- as.numeric(endindex[j])
    if(start < position && position < end)
    {
      #print("Inisde a coding region!")
      codeRegion <- TRUE
      increment <- as.integer(GC3C_Mut_Matrix_Left[rowiter,coliter])
      increment <- increment + 1
      GC3C_Mut_Matrix_Left[rowiter, coliter] <- increment
    }
  }
    #print(as.character(codeRegion))
    
  
  if(upstream == 'A')
  {
    increment <- as.integer(Aup_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Aup_Mut_Matrix_Left[rowiter, coliter] <- increment 
    
    comp_increment <- as.integer(Aup_Comp_Mut_Matrix_Left[comp_rowiter,comp_coliter])
    comp_increment <- comp_increment + 1
    Aup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment 
  }
  if(upstream == 'C')
  {
    increment <- as.integer(Cup_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Cup_Mut_Matrix_Left[rowiter, coliter] <- increment 
    
    comp_increment <- as.integer(Cup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Cup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment 
  }
  if(upstream == 'T')
  {
    increment <- as.integer(Tup_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Tup_Mut_Matrix_Left[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Tup_Comp_Mut_Matrix_Left[comp_rowiter,comp_coliter])
    comp_increment <- comp_increment + 1
    Tup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(upstream == 'G')
  {
    increment <- as.integer(Gup_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Gup_Mut_Matrix_Left[rowiter, coliter] <- increment 
    
    comp_increment <- as.integer(Gup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Gup_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment 
  }
  
  if(downstream == 'A')
  {
    increment <- as.integer(Adown_Mut_Matrix_Left[rowiter,coliter]) 
    increment <- increment + 1
    Adown_Mut_Matrix_Left[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Adown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter]) 
    comp_increment <- comp_increment + 1
    Adown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(downstream == 'C')
  {
    increment <- as.integer(Cdown_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Cdown_Mut_Matrix_Left[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Cdown_Comp_Mut_Matrix_Left[comp_rowiter,comp_coliter])
    comp_increment <- comp_increment + 1
    Cdown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(downstream == 'T')
  {
    increment <- as.integer(Tdown_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Tdown_Mut_Matrix_Left[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Tdown_Comp_Mut_Matrix_Left[comp_rowiter,comp_coliter])
    comp_increment <- comp_increment + 1
    Tdown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(downstream == 'G')
  {
    increment <- as.integer(Gdown_Mut_Matrix_Left[rowiter,coliter])
    increment <- increment + 1
    Gdown_Mut_Matrix_Left[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Gdown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Gdown_Comp_Mut_Matrix_Left[comp_rowiter, comp_coliter] <- comp_increment
  } 
  
 i <- i+1 
}
#================

## Populate the output matrix with right replichore mutants
#================
i <<- 1

while(i < length(Rcore_Frame$Position))
{
  #print(i)
  #print(length(query_Lcore_full))
  #print(query_Rcore_full[i,3])
  #print(query_Rcore_full[i,6])
  position <- as.numeric(Rcore_Frame$Position[i])
  row_switch1 <- Rcore_Frame$Left1[i]
  row_switch2 <- Rcore_Frame$Right1[i]
  rowswitch <- paste(row_switch1, row_switch2, sep = '')
  
  #print(query_Lcore_full[i,4])
  #print(query_Lcore_full[i,5])
  col_switch1 <- Rcore_Frame$Original[i]
  col_switch2 <- Rcore_Frame$Mutation[i]
  colswitch <- paste(col_switch1, col_switch2, sep = '')
  
  upstream <- Rcore_Frame$Left2[i]
  downstream <- Rcore_Frame$Right2[i]
  
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
  
  ##check to see if its in a coding region#
  j <- 1
  for(j in 1:length(startindex))
  {
    start <- as.numeric(startindex[j])
    end <- as.numeric(endindex[j])
    if(start < position && position < end)
    {
      #print("Inisde a coding region!")
      codeRegion <- TRUE
      increment <- as.integer(GC3C_Mut_Matrix_Right[rowiter,coliter])
      increment <- increment + 1
      GC3C_Mut_Matrix_Right[rowiter, coliter] <- increment
    }
  }
  #print(as.character(codeRegion))
  
  
  if(upstream == 'A')
  {
    increment <- as.integer(Aup_Mut_Matrix_Right[rowiter,coliter]) 
    increment <- increment + 1
    Aup_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Aup_Comp_Mut_Matrix_Right[comp_rowiter,comp_coliter]) 
    comp_increment <- comp_increment + 1
    Aup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(upstream == 'C')
  {
    increment <- as.integer(Cup_Mut_Matrix_Right[rowiter,coliter])
    increment <- increment + 1
    Cup_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Cup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Cup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(upstream == 'T')
  {
    increment <- as.integer(Tup_Mut_Matrix_Right[rowiter, coliter])
    increment <- increment + 1
    Tup_Mut_Matrix_Right[rowiter, coliter] <- increment 
    
    comp_increment <- as.integer(Tup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Tup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment 
  }
  if(upstream == 'G')
  {
    increment <- as.integer(Gup_Mut_Matrix_Right[rowiter, coliter])
    increment <- increment + 1
    Gup_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Gup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Gup_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment 
  }
  
  if(downstream == 'A')
  {
    increment <- as.integer(Adown_Mut_Matrix_Right[rowiter,coliter]) 
    increment <- increment + 1
    Adown_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Adown_Comp_Mut_Matrix_Right[comp_rowiter,comp_coliter]) 
    comp_increment <- comp_increment + 1
    Adown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment 
  }
  if(downstream == 'C')
  {
    increment <- as.integer(Cdown_Mut_Matrix_Right[rowiter,coliter])
    increment <- increment + 1
    Cdown_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Cdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Cdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(downstream == 'T')
  {
    increment <- as.integer(Tdown_Mut_Matrix_Right[rowiter,coliter])
    increment <- increment + 1
    Tdown_Mut_Matrix_Right[rowiter, coliter] <- increment
    
    comp_increment <- as.integer(Tdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Tdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment
  }
  if(downstream == 'G')
  {
    increment <- as.integer(Gdown_Mut_Matrix_Right[rowiter,coliter])
    increment <- increment + 1
    Gdown_Mut_Matrix_Right[rowiter, coliter] <- increment 
    
    comp_increment <- as.integer(Gdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter])
    comp_increment <- comp_increment + 1
    Gdown_Comp_Mut_Matrix_Right[comp_rowiter, comp_coliter] <- comp_increment 
  } 
  
  i=i+1
}
#================

Mut_Matrix_Total <- Comp_Mut_Matrix_Left + Mut_Matrix_Right
Comp_Mut_Matrix_Total <- Mut_Matrix_Left + Comp_Mut_Matrix_Right

#Clockwise upstream
Aup_Mut_Matrix_Total <- Aup_Comp_Mut_Matrix_Left + Aup_Mut_Matrix_Right
Cup_Mut_Matrix_Total <- Cup_Comp_Mut_Matrix_Left + Cup_Mut_Matrix_Right
Gup_Mut_Matrix_Total <- Gup_Comp_Mut_Matrix_Left + Gup_Mut_Matrix_Right
Tup_Mut_Matrix_Total <- Tup_Comp_Mut_Matrix_Left + Tup_Mut_Matrix_Right

#Counterclockwise upstream
Aup_Comp_Mut_Matrix_Total <- Aup_Mut_Matrix_Left + Aup_Comp_Mut_Matrix_Right
Cup_Comp_Mut_Matrix_Total <- Cup_Mut_Matrix_Left + Cup_Comp_Mut_Matrix_Right
Gup_Comp_Mut_Matrix_Total <- Gup_Mut_Matrix_Left + Gup_Comp_Mut_Matrix_Right
Tup_Comp_Mut_Matrix_Total <- Tup_Mut_Matrix_Left + Tup_Comp_Mut_Matrix_Right

#Clockwise downstream
Adown_Mut_Matrix_Total <- Adown_Comp_Mut_Matrix_Left + Adown_Mut_Matrix_Right
Cdown_Mut_Matrix_Total <- Cdown_Comp_Mut_Matrix_Left + Cdown_Mut_Matrix_Right
Gdown_Mut_Matrix_Total <- Gdown_Comp_Mut_Matrix_Left + Gdown_Mut_Matrix_Right
Tdown_Mut_Matrix_Total <- Tdown_Comp_Mut_Matrix_Left + Tdown_Mut_Matrix_Right

#Counterclockwise downstream
Adown_Comp_Mut_Matrix_Total <- Adown_Mut_Matrix_Left + Adown_Comp_Mut_Matrix_Right
Cdown_Comp_Mut_Matrix_Total <- Cdown_Mut_Matrix_Left + Cdown_Comp_Mut_Matrix_Right
Gdown_Comp_Mut_Matrix_Total <- Gdown_Mut_Matrix_Left + Gdown_Comp_Mut_Matrix_Right
Tdown_Comp_Mut_Matrix_Total <- Tdown_Mut_Matrix_Left + Tdown_Comp_Mut_Matrix_Right



setwd(Path_BaseSubstitution_output)
write.csv(Mut_Matrix_Total, "Context_Mut_Clockwise_Chromosome.csv")
write.csv(Mut_Matrix_Left, "Context_Mutations_Left.csv")
write.csv(Mut_Matrix_Right, "Context_Mutations_Right.csv")
write.csv(Comp_Mut_Matrix_Total, "RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Comp_Mut_Matrix_Left, "RevCompliment_Context_Mutations_Left.csv")
write.csv(Comp_Mut_Matrix_Right, "RevCompliment_Context_Mutations_Right.csv")

#===== Recording all upstream output for replichores, chromosomes, and the compliments ====#
setwd(Path_BaseSubstitution_output_upA)
write.csv(Aup_Mut_Matrix_Total, "Aup_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Aup_Mut_Matrix_Left, "Aup_Context_Mutations_Left.csv")
write.csv(Aup_Mut_Matrix_Right, "Aup_Context_Mutations_Right.csv")
write.csv(Aup_Comp_Mut_Matrix_Total, "Aup_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Aup_Comp_Mut_Matrix_Left, "Aup_RevCompliment_Context_Mutations_Left.csv")
write.csv(Aup_Comp_Mut_Matrix_Right, "Aup_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_upC)
write.csv(Cup_Mut_Matrix_Total, "Cup_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Cup_Mut_Matrix_Left, "Cup_Context_Mutations_Left.csv")
write.csv(Cup_Mut_Matrix_Right, "Cup_Context_Mutations_Right.csv")
write.csv(Cup_Comp_Mut_Matrix_Total, "Cup_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Cup_Comp_Mut_Matrix_Left, "Cup_RevCompliment_Context_Mutations_Left.csv")
write.csv(Cup_Comp_Mut_Matrix_Right, "Cup_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_upT)
write.csv(Tup_Mut_Matrix_Total, "Tup_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Tup_Mut_Matrix_Left, "Tup_Context_Mutations_Left.csv")
write.csv(Tup_Mut_Matrix_Right, "Tup_Context_Mutations_Right.csv")
write.csv(Tup_Comp_Mut_Matrix_Total, "Tup_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Tup_Comp_Mut_Matrix_Left, "Tup_RevCompliment_Context_Mutations_Left.csv")
write.csv(Tup_Comp_Mut_Matrix_Right, "Tup_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_upG)
write.csv(Gup_Mut_Matrix_Total, "Gup_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Gup_Mut_Matrix_Left, "Gup_Context_Mutations_Left.csv")
write.csv(Gup_Mut_Matrix_Right, "Gup_Context_Mutations_Right.csv")
write.csv(Gup_Comp_Mut_Matrix_Total, "Gup_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Gup_Comp_Mut_Matrix_Left, "Gup_RevCompliment_Context_Mutations_Left.csv")
write.csv(Gup_Comp_Mut_Matrix_Right, "Gup_RevCompliment_Context_Mutations_Right.csv")
#=========================================#

#========= Recording Downstream replichores, chromosomes and compliments output =====#
setwd(Path_BaseSubstitution_output_downA)
write.csv(Adown_Mut_Matrix_Total, "Adown_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Adown_Mut_Matrix_Left, "Adown_Context_Mutations_Left.csv")
write.csv(Adown_Mut_Matrix_Right, "Adown_Context_Mutations_Right.csv")
write.csv(Adown_Comp_Mut_Matrix_Total, "Adown_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Adown_Comp_Mut_Matrix_Left, "Adown_RevCompliment_Context_Mutations_Left.csv")
write.csv(Adown_Comp_Mut_Matrix_Right, "Adown_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_downC)
write.csv(Cdown_Mut_Matrix_Total, "Cdown_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Cdown_Mut_Matrix_Left, "Cdown_Context_Mutations_Left.csv")
write.csv(Cdown_Mut_Matrix_Right, "Cdown_Context_Mutations_Right.csv")
write.csv(Cdown_Comp_Mut_Matrix_Total, "Cdown_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Cdown_Comp_Mut_Matrix_Left, "Cdown_RevCompliment_Context_Mutations_Left.csv")
write.csv(Cdown_Comp_Mut_Matrix_Right, "Cdown_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_downT)
write.csv(Tdown_Mut_Matrix_Total, "Tdown_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Tdown_Mut_Matrix_Left, "Tdown_Context_Mutations_Left.csv")
write.csv(Tdown_Mut_Matrix_Right, "Tdown_Context_Mutations_Right.csv")
write.csv(Tdown_Comp_Mut_Matrix_Total, "Tdown_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Tdown_Comp_Mut_Matrix_Left, "Tdown_RevCompliment_Context_Mutations_Left.csv")
write.csv(Tdown_Comp_Mut_Matrix_Right, "Tdown_RevCompliment_Context_Mutations_Right.csv")

setwd(Path_BaseSubstitution_output_downG)
write.csv(Gdown_Mut_Matrix_Total, "Gdown_Context_Mut_Clockwise_Chromosome.csv")
write.csv(Gdown_Mut_Matrix_Left, "Gdown_Context_Mutations_Left.csv")
write.csv(Gdown_Mut_Matrix_Right, "Gdown_Context_Mutations_Right.csv")
write.csv(Gdown_Comp_Mut_Matrix_Total, "Gdown_RevCompliment_Context_Mut_CounterClockwise_Chromosome.csv")
write.csv(Gdown_Comp_Mut_Matrix_Left, "Gdown_RevCompliment_Context_Mutations_Left.csv")
write.csv(Gdown_Comp_Mut_Matrix_Right, "Gdown_RevCompliment_Context_Mutations_Right.csv")
#=================================#


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
