
#====== Data Refactoring template for single organism analysis ======#

#===== refactoring and creating directories for upstream A storage =====#
Path_BaseSubstitution_output_current <- Path_BaseSubstitution_output_upA
#Path_text_output_current <- paste(Path_text_output_up, "/A", sep = "")
#dir.create(Path_text_output_current)
Path_MutationRates_output <- paste(Path_text_output_up, "/MutationRates_Output", sep = "")
dir.create(Path_MutationRates_output)
Path_MutationRates_output <- paste(Path_MutationRates_output, "/A", sep = "")
dir.create(Path_MutationRates_output)

Path_graphical_output_up_current <- paste(Path_graphical_output_up, "/A", sep = "")
dir.create(Path_graphical_output_up_current)
Path_graphical_output_up_current_counts <- paste(Path_graphical_output_up_current, "/MutationCount_Images", sep = "")
dir.create(Path_graphical_output_up_current_counts)
Path_graphical_output_up_current_rates <- paste(Path_graphical_output_up_current, "/MutationRate_Images", sep = "")
dir.create(Path_graphical_output_up_current_rates)

Path_rates_output <- Path_graphical_output_up_current_rates
Path_counts_output <- Path_graphical_output_up_current_counts

#======== Top strand correlation directories ==========#
Path_correlate_Chromosome_current <- paste(Path_correlate_Chromosome_up, "/A", sep = "")
dir.create(Path_correlate_Chromosome_current)

Path_correlate_Lcore_current <- paste(Path_correlate_Lcore_up, "/A", sep = "")
dir.create(Path_correlate_Lcore_current)
Path_correlate_Rcore_current <- paste(Path_correlate_Rcore_up, "/A", sep = "")
dir.create(Path_correlate_Rcore_current)

Path_correlate_Lcore_Replichore_current <- paste(Path_correlate_Lcore_replichore_up, "/A", sep = "")
dir.create(Path_correlate_Lcore_Replichore_current)
Path_correlate_Rcore_Replichore_current <- paste(Path_correlate_Rcore_replichore_up, "/A", sep = "")
dir.create(Path_correlate_Rcore_Replichore_current)
#======================================================#

#======== Bottom strand correlation directories ==========#
Path_correlate_RevComp_Chromosome_current <- paste(Path_correlate_RevComp_Chromosome_up, "/A", sep = "")
dir.create(Path_correlate_RevComp_Chromosome_current)

Path_correlate_RevComp_Lcore_current <- paste(Path_correlate_RevComp_Lcore_up, "/A", sep = "")
dir.create(Path_correlate_RevComp_Lcore_current)
Path_correlate_RevComp_Lcore_Replichore_current <- paste(Path_correlate_RevComp_Lcore_Replichore_up, "/A", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore_current)

Path_correlate_RevComp_Rcore_current <- paste(Path_correlate_RevComp_Rcore_up, "/A", sep = "")
dir.create(Path_correlate_RevComp_Rcore_current)
Path_correlate_RevComp_Rcore_Replichore_current <- paste(Path_correlate_RevComp_Rcore_Replichore_up, "/A", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore_current)
#======================================================#



#===== Refactoring Variables for operations on the 4mers ====$
GWTC <- GC4C_Aup
Mut_Matrix_Total <- Aup_Mut_Matrix_Total
Mut_Matrix_Left <- Aup_Mut_Matrix_Left
Mut_Matrix_Right <- Aup_Mut_Matrix_Right
Comp_Mut_Matrix_Total <- Aup_Comp_Mut_Matrix_Total
Comp_Mut_Matrix_Left <- Aup_Comp_Mut_Matrix_Left
Comp_Mut_Matrix_Right <- Aup_Comp_Mut_Matrix_Right
#=====================================#
