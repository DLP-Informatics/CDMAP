## Genome Wide Triplet Count (GWTC.r)
# David L Patton
# 800728881

# Purpose: Take the Input Reference fasta file (Path_RefFile) and compute the Reverse Compliment Genome Wide triplet count for each
# amino acid triplet.

TripEnd <- length(RefSeq_arr)-2 #Last left oriented position of nucleotide triplet

cols <- c("T", "G", "C",	"A")
rows <- c("T[X]T", "T[X]G","T[X]C","T[X]A","G[X]T","G[X]G","G[X]C","G[X]A",
          "C[X]T","C[X]G","C[X]C","C[X]A", "A[X]T","A[X]G", "A[X]C","A[X]A")

Rev_GWTC <- matrix(0L, nrow =16, ncol = 4)#Matrix to store Genome wide Triplet count
Rev_GWTC_Left <- matrix(0L, nrow =16, ncol = 4)
Rev_GWTC_Right <- matrix(0L, nrow =16, ncol = 4)

rownames(Rev_GWTC) <- rows
colnames(Rev_GWTC) <- cols
rownames(Rev_GWTC_Left) <- rows
colnames(Rev_GWTC_Left) <- cols
rownames(Rev_GWTC_Right) <- rows
colnames(Rev_GWTC_Right) <- cols


RefFile_new <- getSequence(read.fasta(Path_RefFile))
SourceSeq <- unlist(RefFile_new)
SourceSeq <- toupper(SourceSeq)
Scenario <- (ori_bp < term_bp)
#setwd(path_output)
#write.table(SourceSeq, "SourceSeq.txt")


i <- 1
while(i <= TripEnd)
{
  #Grab each nucleotide in the triplet
  LeftNuc <- RefSeq_arr[i]
  MiddleNuc <- RefSeq_arr[i+1]
  RightNuc <- RefSeq_arr[i+2]
  #print(paste("Triplet: ", LeftNuc, MiddleNuc, RightNuc, sep = ""))
  pos <- i+1 #position of center nucleotide in codon
  rowswitch <- paste(LeftNuc, RightNuc, sep = "")
  colswitch <- MiddleNuc
  
  if(LeftNuc != 'A' & LeftNuc != 'T' & LeftNuc != 'C' & LeftNuc != 'G')
  {
    i <- i+1
    next
  }
  if(MiddleNuc != 'A' & MiddleNuc != 'T' & MiddleNuc != 'C' & MiddleNuc != 'G')
  {
    i <- i+1
    next
  }
  if(RightNuc != 'A' & RightNuc != 'T' & RightNuc != 'C' & RightNuc != 'G')
  {
    i <- i+1
    next
  }
  #Determine the row of the GWTC entry
  #print(paste("Left: ", LeftNuc,"Right: ", RightNuc, sep = " "))
  #print(paste("Left Position: ", i, "Consensus Position: ", i+1, "Right Position: ", i+2, sep = " "))
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

  rowiter <- switch(row_comp,
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
 #Determine the column of the GWTC entry
  col_comp <- switch(colswitch,
                    "T" = 'A',
                    "G" = 'C',
                    "C" = 'G',
                    "A" = 'T'
  )
    coliter <- switch(col_comp,
                    "T" = 1,
                    "G" = 2,
                    "C" = 3,
                    "A" = 4
  )
  #Note: the order is done T -> A, rather than A -> Alphabetically
  #to be consistent with other results generated in the pipeline.
  increment <- as.integer(Rev_GWTC[rowiter,coliter])
  increment <- increment + 1
  Rev_GWTC[rowiter, coliter] <- increment 
  
  if(Scenario == "TRUE")
  {
    Bool1 <- isTRUE(pos >= ori_bp && pos <= term_bp) #Right Core
    Bool2 <- isTRUE(pos <= ori_bp || pos >= term_bp) #Left
  }
  if(Scenario == "FALSE")
  {
    Bool1 <- isTRUE(pos <= term_bp || pos >= ori_bp) #Right Core
    Bool2 <- isTRUE(pos >= term_bp && pos <= ori_bp) #Left Core
  }
  
  if(Bool1)
  {
    increment <- as.integer(Rev_GWTC_Right[rowiter,coliter])
    increment <- increment + 1
    Rev_GWTC_Right[rowiter, coliter] <- increment 
  }
  if(Bool2) 
  {
    increment <- as.integer(Rev_GWTC_Left[rowiter,coliter])
    increment <- increment + 1
    Rev_GWTC_Left[rowiter, coliter] <- increment 
  }
  
  
  i <- i+1
}

#For N-1 position in the sequence
#====================================
LeftNuc <- RefSeq_arr[length(RefSeq_arr)-1]
MiddleNuc <- RefSeq_arr[length(RefSeq_arr)]
RightNuc <- RefSeq_arr[1]
rowswitch <- paste(LeftNuc, RightNuc, sep = "")
colswitch <- MiddleNuc

row_comp <- switch(rowswitch,
                   "TT" = 'AA',
                   "TG" = 'CA',
                   "TC" = 'GA',
                   "TA" = 'AT',
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

rowiter <- switch(row_comp,
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


#Determine the column of the GWTC entry
col_comp <- switch(colswitch,
                   "T" = 'A',
                   "G" = 'C',
                   "C" = 'G',
                   "A" = 'T'
)
coliter <- switch(col_comp,
                  "T" = 1,
                  "G" = 2,
                  "C" = 3,
                  "A" = 4
)

#Note: the order is done T -> A, rather than A -> Alphabetically
#to be consistent with other results generated in the pipeline.
increment <- as.integer(Rev_GWTC[rowiter,coliter])
increment <- increment + 1
Rev_GWTC[rowiter, coliter] <- increment
pos <- length(RefSeq_arr)

if(Scenario == "TRUE")
{
  Bool1 <- isTRUE(pos >= ori_bp && pos <= term_bp) #Right Core
  Bool2 <- isTRUE(pos <= ori_bp || pos >= term_bp) #Left
}
if(Scenario == "FALSE")
{
  Bool1 <- isTRUE(pos <= term_bp || pos >= ori_bp) #Right Core
  Bool2 <- isTRUE(pos >= term_bp && pos <= ori_bp) #Left Core
}

if(Bool1)
{
  increment <- as.integer(Rev_GWTC_Right[rowiter,coliter])
  increment <- increment + 1
  Rev_GWTC_Right[rowiter, coliter] <- increment 
}
if(Bool2) 
{
  increment <- as.integer(Rev_GWTC_Left[rowiter,coliter])
  increment <- increment + 1
  Rev_GWTC_Left[rowiter, coliter] <- increment 
}


#Nth position in the sequence
#=====================================================
LeftNuc <- RefSeq_arr[length(RefSeq_arr)]
MiddleNuc <- RefSeq_arr[1]
RightNuc <- RefSeq_arr[2]
rowswitch <- paste(LeftNuc, RightNuc, sep = "")
colswitch <- MiddleNuc

row_comp <- switch(rowswitch,
                   "TT" = 'AA',
                   "TG" = 'CA',
                   "TC" = 'GA',
                   "TA" = 'AT',
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

rowiter <- switch(row_comp,
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
#Determine the column of the GWTC entry
col_comp <- switch(colswitch,
                   "T" = 'A',
                   "G" = 'C',
                   "C" = 'G',
                   "A" = 'T'
)
coliter <- switch(col_comp,
                  "T" = 1,
                  "G" = 2,
                  "C" = 3,
                  "A" = 4
)
#Note: the order is done T -> A, rather than A -> Alphabetically
#to be consistent with other results generated in the pipeline.
increment <- as.integer(Rev_GWTC[rowiter,coliter])
increment <- increment + 1
Rev_GWTC[rowiter, coliter] <- increment 
pos <- 1

if(Scenario == "TRUE")
{
  Bool1 <- isTRUE(pos >= ori_bp && pos <= term_bp) #Right Core
  Bool2 <- isTRUE(pos <= ori_bp || pos >= term_bp) #Left
}
if(Scenario == "FALSE")
{
  Bool1 <- isTRUE(pos <= term_bp || pos >= ori_bp) #Right Core
  Bool2 <- isTRUE(pos >= term_bp && pos <= ori_bp) #Left Core
}

if(Bool1)
{
  increment <- as.integer(Rev_GWTC_Right[rowiter,coliter])
  increment <- increment + 1
  Rev_GWTC_Right[rowiter, coliter] <- increment 
}
if(Bool2) 
{
  increment <- as.integer(Rev_GWTC_Left[rowiter,coliter])
  increment <- increment + 1
  Rev_GWTC_Left[rowiter, coliter] <- increment 
}

output_name <- paste("Rev_GWTC_", organism, ".csv", sep = "")
output_Left <- paste("Rev_GWTC_Left_", organism, ".csv", sep = "")
output_Right <- paste("Rev_GWTC_Right_", organism, ".csv", sep = "")
setwd(path_GWTC_output)
write.csv(Rev_GWTC, output_name)
write.csv(Rev_GWTC_Left, output_Left)
write.csv(Rev_GWTC_Right, output_Right)