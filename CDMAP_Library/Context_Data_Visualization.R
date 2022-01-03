#Path_test_output <- "/Users/triadge/Desktop/Sung_Lab/testdir/test_output"




data_rates <- c(
                "Sung_MutChrome",
                "Sung_MutLeft",
                "Sung_MutLeftRep",
                "Sung_MutRight",
                "Sung_MutRightRep")
#if(Flag_3mer)
#{
  data_rates_comp <- c(
                     "Comp_Sung_MutChrome",
                     "Comp_Sung_MutLeft",
                     "Comp_Sung_MutLeftRep",
                     "Comp_Sung_MutRight",
                     "Comp_Sung_MutRightRep")
  
  data_rates <- c(data_rates, data_rates_comp)
#}  
# data_ratios <- c("MutationRatio_Chromosome",
#                 "MutationRatio_Left",
#                 "MutationRatio_LeftRep",
#                 "MutationRatio_Right",
#                 "MutationRatio_RightRep",
#                 "Comp_MutationRatio_Chromosome",
#                 "Comp_MutationRatio_Left",
#                 "Comp_MutationRatio_LeftRep",
#                 "Comp_MutationRatio_Right",
#                 "Comp_MutationRatio_RightRep")

data_counts <- c("Mut_Matrix_Left", 
               "Mut_Matrix_Right",
               "Mut_Matrix_Total",
               "Sung_Matrix",
               "Sung_Matrix_Left",
               "Sung_Matrix_Right",
               "GWTC",
               "GWTC_Left",
               "GWTC_Right",
               "Comp_Mut_Matrix_Left",
               "Comp_Mut_Matrix_Right",
               "Comp_Mut_Matrix_Total",
               "Comp_Sung_Matrix",
               "Comp_Sung_Matrix_Left",
               "Comp_Sung_Matrix_Right",
               "Rev_GWTC",
               "Rev_GWTC_Left",
               "Rev_GWTC_Right"
               ) 
flag_correlate <- FALSE

#param_flag <- readline("How would you like to scale your output data for the mutation counts, triplets, and mutation rates? 
#         (0 for default, 1 for scaling to the average mean, or 2, for custom parameters): ")

if(param_flag == 0)
{
  mut_count <- 1
  #triplet_count <- 1e3 #DEPRECIATE
  mutrate_count <- 1e8
}
if(param_flag == 1)
{
  mut_count <- 1
  triplet_count <- mean(MutationRatio_Chromosome)
  mutrate_count <- mean(Sung_ContextChrome)
}
if(param_flag == 2)
{
  mut_count <- readline("please input the numeric scalar for raw mutant counts: ")
  triplet_count <- readline("please input the numeric scalar for triplet counts: ")
  mutrate_count <- readline("please input the numeric scalar for mutation rate calculations: ")
}

if(param_flag == 4)
{
  mut_count <- 1
  triplet_count <- 1e2
  mutrate_count <- 1e2
}

#BaseSubstitution and Raw Triplet Counts
for(i in 1:length(data_counts))
{
  setwd(Path_to_scripts)
  matrix_name <- as.character(data_counts[i])
  output_data_matrix <- get(matrix_name)
  Path_image_output <-  Path_counts_output
  output_data_matrix <- output_data_matrix*mut_count
  matrix_name_graph <- matrix_name
  source("lattice_visualizer.R")
}

#MutationRatio Counts
# for(i in 1:length(data_ratios))
# {
#  setwd(Path_to_scripts)
#  matrix_name <- as.character(data_ratios[i])
#  output_data_matrix <- get(matrix_name)
#  Path_image_output <-  Path_ratios_output
#  matrix_name_graph <- paste(as.character(data_ratios[i]), "(Scaled to ", triplet_count, ")", sep =" ")
#  output_data_matrix <- output_data_matrix*triplet_count
#  source("lattice_visualizer.R")
# 
# }

#MutationRate Counts
for(i in 1:length(data_rates))
{
  setwd(Path_to_scripts)
  matrix_name <- as.character(data_rates[i])
  output_data_matrix <- get(matrix_name)
  Path_image_output <-  Path_rates_output
  matrix_name_graph <- paste(as.character(data_rates[i]), "(Scaled to ", mutrate_count, ")", sep =" ")
  output_data_matrix <- output_data_matrix*mutrate_count #Rate is multiplied by scalar so it's visible on heatmap
  source("lattice_visualizer.R")
}