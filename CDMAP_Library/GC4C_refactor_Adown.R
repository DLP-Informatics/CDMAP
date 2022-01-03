
#====== Data Refactoring template for single organism analysis ======#

#===== refactoring and creating directories for downstream A storage =====#
Path_BaseSubstitution_output_current <- Path_BaseSubstitution_output_downA
#Path_text_output_current <- paste(Path_text_output_down, "/A", sep = "")
#dir.create(Path_text_output_current)
Path_MutationRates_output <- paste(Path_text_output_down, "/MutationRates_Output", sep = "")
dir.create(Path_MutationRates_output)
Path_MutationRates_output <- paste(Path_MutationRates_output, "/A", sep = "")
dir.create(Path_MutationRates_output)

Path_graphical_output_down_current <- paste(Path_graphical_output_down, "/A", sep = "")
dir.create(Path_graphical_output_down_current)
Path_graphical_output_down_current_counts <- paste(Path_graphical_output_down_current, "/MutationCount_Images", sep = "")
dir.create(Path_graphical_output_down_current_counts)
Path_graphical_output_down_current_rates <- paste(Path_graphical_output_down_current, "/MutationRate_Images", sep = "")
dir.create(Path_graphical_output_down_current_rates)

Path_rates_output <- Path_graphical_output_down_current_rates
Path_counts_output <- Path_graphical_output_down_current_counts

#======== Top strand correlation directories ==========#
Path_correlate_Chromosome_current <- paste(Path_correlate_Chromosome_down, "/A", sep = "")
dir.create(Path_correlate_Chromosome_current)

Path_correlate_Lcore_current <- paste(Path_correlate_Lcore_down, "/A", sep = "")
dir.create(Path_correlate_Lcore_current)
Path_correlate_Rcore_current <- paste(Path_correlate_Rcore_down, "/A", sep = "")
dir.create(Path_correlate_Rcore_current)

Path_correlate_Lcore_Replichore_current <- paste(Path_correlate_Lcore_replichore_down, "/A", sep = "")
dir.create(Path_correlate_Lcore_Replichore_current)
Path_correlate_Rcore_Replichore_current <- paste(Path_correlate_Rcore_replichore_down, "/A", sep = "")
dir.create(Path_correlate_Rcore_Replichore_current)
#======================================================#

#======== Bottom strand correlation directories ==========#
Path_correlate_RevComp_Chromosome_current <- paste(Path_correlate_RevComp_Chromosome_down, "/A", sep = "")
dir.create(Path_correlate_RevComp_Chromosome_current)

Path_correlate_RevComp_Lcore_current <- paste(Path_correlate_RevComp_Lcore_down, "/A", sep = "")
dir.create(Path_correlate_RevComp_Lcore_current)
Path_correlate_RevComp_Lcore_Replichore_current <- paste(Path_correlate_RevComp_Lcore_Replichore_down, "/A", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore_current)

Path_correlate_RevComp_Rcore_current <- paste(Path_correlate_RevComp_Rcore_down, "/A", sep = "")
dir.create(Path_correlate_RevComp_Rcore_current)
Path_correlate_RevComp_Rcore_Replichore_current <- paste(Path_correlate_RevComp_Rcore_Replichore_down, "/A", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore_current)
#======================================================#



#===== Refactoring Variables for operations on the 4mers ====$
GWTC_Current <- GC4C_Adown
GWTC_Left_Count_Current <- GC4C_Adown_Left
GWTC_Right_Count_Current <- GC4C_Adown_Right
Rev_GWTC_Current <- GC4C_Adown_Comp
Rev_GWTC_Left_Current <- GC4C_Adown_Left_Comp
Rev_GWTC_Right_Current <- GC4C_Adown_Right_Comp

Mut_Matrix_Total <- Adown_Mut_Matrix_Total
Mut_Matrix_Left <- Adown_Mut_Matrix_Left
Mut_Matrix_Right <- Adown_Mut_Matrix_Right
Comp_Mut_Matrix_Total <- Adown_Comp_Mut_Matrix_Total
Comp_Mut_Matrix_Left <- Adown_Comp_Mut_Matrix_Left
Comp_Mut_Matrix_Right <- Adown_Comp_Mut_Matrix_Right

len_row <- length(GWTC_Current[,1])
len_col <- length(GWTC_Current[1,])

for(i in 1:len_row)
{
  for(j in 1:len_col)
  {
    if(GWTC_Current[i,j] == 0)
    {
      GWTC_Current[i,j] <- 1
    }
    if(GWTC_Left_Current[i,j] == 0)
    {
      GWTC_Left_Current[i,j] <- 1
    }
    if(GWTC_Right_Current[i,j] == 0)
    {
      GWTC_Right_Current[i,j] <- 1
    }
    if(Rev_GWTC_Current[i,j] == 0)
    {
      Rev_GWTC_Current[i,j] <- 1
    }
    if(Rev_GWTC_Left_Current[i,j] == 0)
    {
      Rev_GWTC_Left_Current[i,j] <- 1
    }
    if(Rev_GWTC_Right_Current[i,j] == 0)
    {
      Rev_GWTC_Right_Current[i,j] <- 1
    }
  }
}