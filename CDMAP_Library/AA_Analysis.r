#AA_Analysis.R

#    Things to do
# 2. generate a mutation probability matrix for each amino acid (if applicable)
    #2a. Implement length checker to ensure that we do skip sites with zero observed mutation
#3. generate the Chi Square analysis for each amino acid.

#Vars#
#FoldObj 
#FoldObjName
#########

AA_Pval <- c()
AA_Dir <- paste(FoldDir, "/AA", sep = "")

if(!(dir.exists(AA_Dir)))
{
  dir.create(AA_Dir)
}

# 1. Output each respective amino acid input into the corresponding fold specific directory
#=================================#
setwd(FoldDir)
for(p in 1: length(FoldObjName))
{
  AA_MutRateProbMatrix <- matrix(nrow = 0, ncol = 6)
  colnames(AA_MutRateProbMatrix) <- MutRateColnames
  # 1. Output each respective amino acid input into the corresponding fold specific directory
  #=================================#
  AA_ObjName <- FoldObjName[p]
  AA_matrix <- get(AA_ObjName)
  filetitle <- paste(AA_ObjName, "_input.csv", sep = "")
  
  #remove square brackets
  AA_matrix[,2] <- gsub("\\[|\\]", "", AA_matrix[,2])
  setwd(AA_Dir)
  write.csv(AA_matrix, filetitle)
  setwd(FoldDir)
  #=================================#
  
  #remove nonzero entries
  for(s in length(AA_matrix[,1]):1) #reverse iterate through matrix to avoid subscript out of bound errors
  {
    if(AA_matrix[s,6] == 0) #detect if there is a 0 mutation rate
    {
      AA_matrix <- AA_matrix[-s,] 
    }
  } 
  
  
  if(dim(AA_matrix)[1] == 0)
  {
  AA_ChiSq_Obj <- c(AA_ObjName, "NA", "NA")
  AA_Pval <- rbind(AA_Pval, AA_ChiSq_Obj)
  next
  }
  
  # 2. Collect all codon triplets that encode a given amino acid
  #=================================#
  for(m in 1:length(AA_matrix[,2]))
  {
    TripObj <- AA_matrix[m,2]
    Tripind <- which(grepl(TripObj, MutRateProbMatrix[,1]))
    if(isempty(Tripind))
    {
      next
    }
    AA_MutRateProbMatrix <- rbind(AA_MutRateProbMatrix, MutRateProbMatrix[Tripind, ])
  }
  
  #manual recalculation of expected codon usage w.r.t amino acid
  AA_MutRateProbMatrix[,5:6] <- 0
  AA_MutRateProbMatrix <- data.frame(AA_MutRateProbMatrix)
  SumInvMutRate <- sum(as.numeric(AA_MutRateProbMatrix$Inverse.Mutation.Rate))
  SumObs <- sum(as.numeric(AA_MutRateProbMatrix$Observed.Counts))
  AA_MutRateProbMatrix <- as.matrix(AA_MutRateProbMatrix)
  AA_MutRateProbMatrix[,5] <- as.numeric(AA_MutRateProbMatrix[,4])/SumInvMutRate
  AA_MutRateProbMatrix[,6] <- as.numeric(AA_MutRateProbMatrix[,5])*SumObs
  
  setwd(AA_Dir)
  AA_MutProbMatName <- paste(AA_ObjName, "Mutation_Probability_Matrix.csv")
  write.csv(AA_MutRateProbMatrix, AA_MutProbMatName)
  setwd(FoldDir)
  #=================================#
  
  AA_DF <- data.frame(as.numeric(AA_MutRateProbMatrix[,2]), as.numeric(AA_MutRateProbMatrix[,6]))
  AA_ChiSq <- chisq.test(AA_DF)
  
  AA_ChiSq_Obj <- c(AA_ObjName, AA_ChiSq$statistic, AA_ChiSq$p.value)
  AA_Pval <- rbind(AA_Pval, AA_ChiSq_Obj)
}
rownames(AA_Pval) <- NULL
AAColNames <- c("AminoAcid", "MutChiStat", "MutChiPval")
colnames(AA_Pval) <- AAColNames


filetitle <- paste(fold, "AminoAcid_ChiSquare.csv", sep = "_")
write.csv(AA_Pval, filetitle)