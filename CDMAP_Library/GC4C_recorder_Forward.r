

if(iter == 2) #assigns upNuc to up if its the first 4mer in analysis
{
  up <- upNuc
}else 
{
  up <- gene_arr[iter-2] #Assigns upstream 4mer wrt triplet
}

if(iter == gene_end) #assigns downNuc to down if its the last 4mer in analysis
{
  down <- downNuc 
}else
{
  down <- gene_arr[iter+2]  #Assigns downstream 4mer wrt triplet
}

LeftNuc <- gene_arr[iter-1] #left neighbor nucleotide
MiddleNuc <- gene_arr[iter] #center nucleotide
RightNuc <- gene_arr[iter+1] #right neighbor nucelotide

# DiagTriplet <- paste(LeftNuc, MiddleNuc, RightNuc, sep = "")
# DiagArr <- c(DiagArr, DiagTriplet)

#GC content checker
GC_Check <-  c("C", "c", "G", "g")

if(any(grepl(LeftNuc, GC_Check, ignore.case = TRUE)))
{
  GC_counter <- GC_counter + 1
}
if(any(grepl(MiddleNuc, GC_Check, ignore.case = TRUE)))
{
  GC_counter <- GC_counter + 1
}
if(any(grepl(RightNuc, GC_Check, ignore.case = TRUE)))
{
  GC_counter <- GC_counter + 1
}



rowswitch <- paste(LeftNuc, RightNuc, sep = "") #row matrix assignment operator
colswitch <- MiddleNuc #column matrix assignment operator

#null reference checker. checking the 5mer block for null references
if(is.na(LeftNuc))
{
  j <- j+3
  next
}
if(is.na(MiddleNuc))
{
  j <- j+3
  next
}
if(is.na(RightNuc))
{
  j <- j+3
  next
}
if(is.na(up))
{
  j <- j+3
  next
}
if(is.na(down))
{
  j <- j+3
  next
}


#checks the nucleotide triplets for non-ambiguous entries, skips if not atcg
if(LeftNuc != 'A' & LeftNuc != 'T' & LeftNuc != 'C' & LeftNuc != 'G')
{
  j <- j+3
  next
}
if(MiddleNuc != 'A' & MiddleNuc != 'T' & MiddleNuc != 'C' & MiddleNuc != 'G')
{
  j <- j+3
  next
}
if(RightNuc != 'A' & RightNuc != 'T' & RightNuc != 'C' & RightNuc != 'G')
{
  j <- j+3
  next
}


#Row matrix assignment operator for 3mer
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
#column assignment operator
coliter <- switch(colswitch,
                  "T" = 1,
                  "G" = 2,
                  "C" = 3,
                  "A" = 4
)
#row and column complement converters
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

col_comp <- switch(colswitch,
                   "T" = 'A',
                   "G" = 'C',
                   "C" = 'G',
                   "A" = 'T'
)
comp_coliter <- switch(col_comp,
                  "T" = 1,
                  "G" = 2,
                  "C" = 3,
                  "A" = 4
)

#manually extracts and increments the entry in the matrix, then reassigns it to i,j
# increment <- as.integer(GC3C[rowiter,coliter])
# increment <- increment + 1
# GC3C[rowiter, coliter] <- increment 

increment <- as.integer(TestTotalMatrix[rowiter,coliter])
increment <- increment + 1
TestTotalMatrix[rowiter, coliter] <- increment


if(orientation == '+' && flagcheck == 'Right')
{
  increment <- as.integer(GC3C_Right[rowiter,coliter])
  increment <- increment + 1
  GC3C_Right[rowiter, coliter] <- increment 
  
}
if(orientation == '-' && flagcheck == 'Left')
{
  increment <- as.integer(Rev_GC3C_Left[rowiter,coliter])
  increment <- increment + 1
  Rev_GC3C_Left[rowiter, coliter] <- increment 
  
}


#checks the upstream and downstream nucleotide, and whether
#it sits on the right or left replichore and updates the corresponding matrix and updates it.

#Retrieves to temporary matrix based on left/right, and corresponding up/down nucleotide

#Upstream
#===========================================
if(flagcheck == 'Left')
{
  if(up == "T")
  {
    upstream_matrix <- GC4C_Tup_Left
    comp_upstream_matrix <- GC4C_Tup_Left_Comp
  }
  if(up == "G")
  {
    upstream_matrix <- GC4C_Gup_Left
    comp_upstream_matrix <- GC4C_Gup_Left_Comp
  }
  if(up == "C")
  {
    upstream_matrix <- GC4C_Cup_Left
    comp_upstream_matrix <- GC4C_Cup_Left_Comp
  }
  if(up == "A")
  {
    upstream_matrix <- GC4C_Aup_Left
    comp_upstream_matrix <- GC4C_Aup_Left_Comp
  }
}
if(flagcheck == 'Right')
{
  if(up == "T")
  {
    upstream_matrix <- GC4C_Tup_Right
    comp_upstream_matrix <- GC4C_Tup_Right_Comp
  }
  if(up == "G")
  {
    upstream_matrix <- GC4C_Gup_Right
    comp_upstream_matrix <- GC4C_Gup_Right_Comp
  }
  if(up == "C")
  {
    upstream_matrix <- GC4C_Cup_Right
    comp_upstream_matrix <- GC4C_Cup_Right_Comp
  }
  if(up == "A")
  {
    upstream_matrix <- GC4C_Aup_Right
    comp_upstream_matrix <- GC4C_Aup_Right_Comp
  }
}
#Extracts and updates the corresponding i,j value of each matrix
increment <- as.integer(upstream_matrix[rowiter,coliter])
increment <- increment + 1
upstream_matrix[rowiter, coliter] <- increment 

comp_increment <- as.integer(comp_upstream_matrix[comp_rowiter, comp_coliter])
comp_increment <- comp_increment + 1
comp_upstream_matrix[comp_rowiter, comp_coliter] <- comp_increment 

#reassigns it to original matrix in updated form
if(flagcheck == 'Left')
{
  if(up == "T")
  {
    GC4C_Tup_Left <- upstream_matrix
    GC4C_Tup_Left_Comp <- comp_upstream_matrix
  }
  if(up == "G")
  {
    GC4C_Gup_Left <- upstream_matrix
    GC4C_Gup_Left_Comp <- comp_upstream_matrix
  }
  if(up == "C")
  {
    GC4C_Cup_Left <- upstream_matrix
    GC4C_Cup_Left_Comp <- comp_upstream_matrix
  }
  if(up == "A")
  {
    GC4C_Aup_Left <- upstream_matrix
    GC4C_Aup_Left_Comp <- comp_upstream_matrix
  }
}
#===========================================

#Downstream
#===========================================
if(flagcheck == 'Right')
{
  if(up == "T")
  {
    GC4C_Tup_Right <- upstream_matrix
    GC4C_Tup_Right_Comp <- comp_upstream_matrix
  }
  if(up == "G")
  {
    GC4C_Gup_Right <- upstream_matrix
    GC4C_Gup_Right_Comp <- comp_upstream_matrix
  }
  if(up == "C")
  {
    GC4C_Cup_Right <- upstream_matrix
    GC4C_Cup_Right_Comp <- comp_upstream_matrix
  }
  if(up == "A")
  {
    GC4C_Aup_Right <- upstream_matrix
    GC4C_Aup_Right_Comp <- comp_upstream_matrix
  }
}

if(flagcheck == 'Left')
{
  if(down == "T")
  {
    downstream_matrix <- GC4C_Tdown_Left
    comp_downstream_matrix <- GC4C_Tdown_Left_Comp
  }
  if(down == "G")
  {
    downstream_matrix <- GC4C_Gdown_Left
    comp_downstream_matrix <- GC4C_Gdown_Left_Comp
  }
  if(down == "C")
  {
    downstream_matrix <- GC4C_Cdown_Left
    comp_downstream_matrix <- GC4C_Cdown_Left_Comp
  }
  if(down == "A")
  {
    downstream_matrix <- GC4C_Adown_Left
    comp_downstream_matrix <- GC4C_Adown_Left_Comp
  }
}
if(flagcheck == 'Right')
{
  if(down == "T")
  {
    downstream_matrix <- GC4C_Tdown_Right
    comp_downstream_matrix <- GC4C_Tdown_Right_Comp
  }
  if(down == "G")
  {
    downstream_matrix <- GC4C_Gdown_Right
    comp_downstream_matrix <- GC4C_Gdown_Right_Comp
  }
  if(down == "C")
  {
    downstream_matrix <- GC4C_Cdown_Right
    comp_downstream_matrix <- GC4C_Cdown_Right_Comp
  }
  if(down == "A")
  {
    downstream_matrix <- GC4C_Adown_Right
    comp_downstream_matrix <- GC4C_Adown_Right_Comp
  }
}

increment <- as.integer(downstream_matrix[rowiter,coliter])
increment <- increment + 1
downstream_matrix[rowiter, coliter] <- increment 

comp_increment <- as.integer(comp_downstream_matrix[comp_rowiter,comp_coliter])
comp_increment <- comp_increment + 1
comp_downstream_matrix[comp_rowiter, comp_coliter] <- comp_increment 

if(flagcheck == 'Left')
{
  if(down == "T")
  {
    GC4C_Tdown_Left <-  downstream_matrix
    GC4C_Tdown_Left_Comp <-  comp_downstream_matrix
  }
  if(down == "G")
  {
    GC4C_Gdown_Left <-  downstream_matrix
    GC4C_Gdown_Left_Comp <-  comp_downstream_matrix
  }
  if(down == "C")
  {
    GC4C_Cdown_Left <-  downstream_matrix
    GC4C_Cdown_Left_Comp <-  comp_downstream_matrix
  }
  if(down == "A")
  {
    GC4C_Adown_Left <-  downstream_matrix
    GC4C_Adown_Left_Comp <-  comp_downstream_matrix
  }
}

if(flagcheck == 'Right')
{
  if(down == "T")
  {
    GC4C_Tdown_Right <-  downstream_matrix
    GC4C_Tdown_Right_Comp <-  comp_downstream_matrix
  }
  if(down == "G")
  {
    GC4C_Gdown_Right <-  downstream_matrix
    GC4C_Gdown_Right_Comp <-  comp_downstream_matrix
  }
  if(down == "C")
  {
    GC4C_Cdown_Right <-  downstream_matrix
    GC4C_Cdown_Right_Comp <-  comp_downstream_matrix
  }
  if(down == "A")
  {
    GC4C_Adown_Right <-  downstream_matrix
    GC4C_Adown_Right_Comp <-  comp_downstream_matrix
  }
} 