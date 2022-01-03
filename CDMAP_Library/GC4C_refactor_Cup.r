
#====== Data Refactoring template for single organism analysis ======#

#===== refactoring and creating directories for upstream A storage =====#
Path_BaseSubstitution_output_current <- Path_BaseSubstitution_output_upC
#Path_text_output_current <- paste(Path_text_output_up, "/C", sep = "")
#dir.create(Path_text_output_current)
Path_MutationRates_output <- paste(Path_text_output_up, "/MutationRates_Output", sep = "")
dir.create(Path_MutationRates_output)
Path_MutationRates_output <- paste(Path_MutationRates_output, "/C", sep = "")
dir.create(Path_MutationRates_output)

Path_graphical_output_up_current <- paste(Path_graphical_output_up, "/C", sep = "")
dir.create(Path_graphical_output_up_current)
Path_graphical_output_up_current_counts <- paste(Path_graphical_output_up_current, "/MutationCount_Images", sep = "")
dir.create(Path_graphical_output_up_current_counts)
Path_graphical_output_up_current_rates <- paste(Path_graphical_output_up_current, "/MutationRate_Images", sep = "")
dir.create(Path_graphical_output_up_current_rates)

Path_rates_output <- Path_graphical_output_up_current_rates
Path_counts_output <- Path_graphical_output_up_current_counts

#======== Top strand correlation directories ==========#
Path_correlate_Chromosome_current <- paste(Path_correlate_Chromosome_up, "/C", sep = "")
dir.create(Path_correlate_Chromosome_current)

Path_correlate_Lcore_current <- paste(Path_correlate_Lcore_up, "/C", sep = "")
dir.create(Path_correlate_Lcore_current)
Path_correlate_Rcore_current <- paste(Path_correlate_Rcore_up, "/C", sep = "")
dir.create(Path_correlate_Rcore_current)

Path_correlate_Lcore_Replichore_current <- paste(Path_correlate_Lcore_replichore_up, "/C", sep = "")
dir.create(Path_correlate_Lcore_Replichore_current)
Path_correlate_Rcore_Replichore_current <- paste(Path_correlate_Rcore_replichore_up, "/C", sep = "")
dir.create(Path_correlate_Rcore_Replichore_current)
#======================================================#

#======== Bottom strand correlation directories ==========#
Path_correlate_RevComp_Chromosome_current <- paste(Path_correlate_RevComp_Chromosome_up, "/C", sep = "")
dir.create(Path_correlate_RevComp_Chromosome_current)

Path_correlate_RevComp_Lcore_current <- paste(Path_correlate_RevComp_Lcore_up, "/C", sep = "")
dir.create(Path_correlate_RevComp_Lcore_current)
Path_correlate_RevComp_Lcore_Replichore_current <- paste(Path_correlate_RevComp_Lcore_Replichore_up, "/C", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore_current)

Path_correlate_RevComp_Rcore_current <- paste(Path_correlate_RevComp_Rcore_up, "/C", sep = "")
dir.create(Path_correlate_RevComp_Rcore_current)
Path_correlate_RevComp_Rcore_Replichore_current <- paste(Path_correlate_RevComp_Rcore_Replichore_up, "/C", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore_current)
#======================================================#



#===== Refactoring Variables for operations on the 4mers ====$
GWTC_Current <- GC4C_Cup
GWTC_Left_Current <- GC4C_Cup_Left
GWTC_Right_Current <- GC4C_Cup_Right
Rev_GWTC_Current <- GC4C_Cup_Comp
Rev_GWTC_Left_Current <- GC4C_Cup_Left_Comp
Rev_GWTC_Right_Current <- GC4C_Cup_Right_Comp

Mut_Matrix_Total <- Cup_Mut_Matrix_Total
Mut_Matrix_Left <- Cup_Mut_Matrix_Left
Mut_Matrix_Right <- Cup_Mut_Matrix_Right
Comp_Mut_Matrix_Total <- Cup_Comp_Mut_Matrix_Total
Comp_Mut_Matrix_Left <- Cup_Comp_Mut_Matrix_Left
Comp_Mut_Matrix_Right <- Cup_Comp_Mut_Matrix_Right

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