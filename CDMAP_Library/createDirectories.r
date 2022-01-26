##### createDirectories.r #######
#=================================
#Author: David Logan Patton
# Objective: The primary objective of this R script is to generate all of the ordered output directories dynamically
#for a given organisms triplet counts and context dependent mutation rates with respect to the leading and lagging strand.
#For reach directory, CDMAP checks if the directory exists, then proceeds to create it if it does not exist.

#create the macro output directory for all output media if it does not exist
if(!(dir.exists(Path_MainRepo)))
{
  dir.create(Path_MainRepo)
}

if(!(dir.exists(Path_output)))
{
  dir.create(Path_output) #creates the base level output directory for all organisms
}

if(!(dir.exists(Path_repoGC)))
{
  dir.create(Path_repoGC)
}

if(!(dir.exists(Path_repoFASTA)))
{
  dir.create(Path_repoFASTA)
}


if(!(dir.exists(Path_output_organism)))
{
  dir.create(Path_output_organism) #creates the organism specific directory in the base directory
}
if(!(dir.exists(Path_output_downstream)))
{
  dir.create(Path_output_downstream) #creates the downstream output directory within the organism specific directory
}

if(!(dir.exists(Path_output_upstream)))
{
  dir.create(Path_output_upstream) #creates the upstream output directory within the organism specific directory
}
if(!(dir.exists(Path_output_triplet)))
{
  dir.create(Path_output_triplet) #creates the base triplet output directory within the organism specific directory
}
#create the macro level repository for replichore and chromosome output for correlation analysis
if(!(dir.exists(Path_correlate_repo)))
{
  dir.create(Path_correlate_repo) #creates the base level correlation dump repository for all organisms for 3mer and 4mer
}

if(!(dir.exists(Path_correlate_triplet)))
{
  dir.create(Path_correlate_triplet) #creates the triplet correlation dump repository
}

if(!(dir.exists(Path_correlate_repo_down)))
{
  dir.create(Path_correlate_repo_down) #creates the downstream 4mer correlation dump repository
}

if(!(dir.exists(Path_correlate_repo_up)))
{
  dir.create(Path_correlate_repo_up) #Creates the upstream 4mer correlation dump repository
}


#This creates all of the text and graphical output directories for both the organism specific output, and the correlation
#repositories later on for statistical analysis of the base triplets
#===================================

#macro repositories for text and graphical output
Path_graphical_output <- paste(Path_output_triplet, "/Output_Images", sep = "")
dir.create(Path_graphical_output)
Path_graphical_output_counts <- paste(Path_graphical_output, "/MutationCount_Images", sep = "")
dir.create(Path_graphical_output_counts)
Path_graphical_output_rates <- paste(Path_graphical_output, "/MutationRate_Images", sep = "")
dir.create(Path_graphical_output_rates)

Path_text_output <- paste(Path_output_triplet, "/Output_Text", sep = "")
dir.create(Path_text_output)

#chromosome correlation repository
Path_correlate_Chromosome <- paste(Path_correlate_triplet, "/Chromosome", sep = "")
dir.create(Path_correlate_Chromosome)
Path_correlate_RevComp_Chromosome <- paste(Path_correlate_triplet, "/RevComp_Chromosome", sep = "")
dir.create(Path_correlate_RevComp_Chromosome)

#left replichore correlation repository
Path_correlate_Lcore <- paste(Path_correlate_triplet, "/Lcore", sep = "")
dir.create(Path_correlate_Lcore)
Path_correlate_Lcore_Replichore <- paste(Path_correlate_triplet, "/Lcore_Replichore", sep = "")
dir.create(Path_correlate_Lcore_Replichore)
Path_correlate_RevComp_Lcore <- paste(Path_correlate_triplet, "/RevComp_Lcore", sep = "")
dir.create(Path_correlate_RevComp_Lcore)
Path_correlate_RevComp_Lcore_Replichore <- paste(Path_correlate_triplet, "/RevComp_Lcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore)

#right replichore correlation repository
Path_correlate_Rcore <- paste(Path_correlate_triplet, "/Rcore", sep = "")
dir.create(Path_correlate_Rcore)
Path_correlate_Rcore_Replichore <- paste(Path_correlate_triplet, "/Rcore_Replichore", sep = "")
dir.create(Path_correlate_Rcore_Replichore)
Path_correlate_RevComp_Rcore <- paste(Path_correlate_triplet, "/RevComp_Rcore", sep = "")
dir.create(Path_correlate_RevComp_Rcore)
Path_correlate_RevComp_Rcore_Replichore <- paste(Path_correlate_triplet, "/RevComp_Rcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore)


#===================================


#This creates all of the text and graphical output directories for both the organism specific output, and the correlation
#repositories later on for statistical analysis of the Upstream 4mers
#===================================

#macro repositories for text and graphical output
Path_graphical_output_up <- paste(Path_output_upstream, "/Output_Images", sep = "")
dir.create(Path_graphical_output_up)
Path_text_output_up <- paste(Path_output_upstream, "/Output_Text", sep = "")
dir.create(Path_text_output_up)

#chromosome correlation repository
Path_correlate_Chromosome_up <- paste(Path_correlate_repo_up, "/Chromosome", sep = "")
dir.create(Path_correlate_Chromosome_up)
Path_correlate_RevComp_Chromosome_up <- paste(Path_correlate_repo_up, "/RevComp_Chromosome", sep = "")
dir.create(Path_correlate_RevComp_Chromosome_up)

#left replichore correlation repository
Path_correlate_Lcore_up <- paste(Path_correlate_repo_up, "/Lcore", sep = "")
dir.create(Path_correlate_Lcore_up)
Path_correlate_Lcore_replichore_up <- paste(Path_correlate_repo_up, "/Lcore_Replichore", sep = "")
dir.create(Path_correlate_Lcore_replichore_up)
Path_correlate_RevComp_Lcore_up <- paste(Path_correlate_repo_up, "/RevComp_Lcore", sep = "")
dir.create(Path_correlate_RevComp_Lcore_up)
Path_correlate_RevComp_Lcore_Replichore_up <- paste(Path_correlate_repo_up, "/RevComp_Lcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore_up)

#right replichore correlation repository
Path_correlate_Rcore_up <- paste(Path_correlate_repo_up, "/Rcore", sep = "")
dir.create(Path_correlate_Rcore_up)
Path_correlate_Rcore_replichore_up <- paste(Path_correlate_repo_up, "/Rcore_Replichore", sep = "")
dir.create(Path_correlate_Rcore_replichore_up)
Path_correlate_RevComp_Rcore_up <- paste(Path_correlate_repo_up, "/RevComp_Rcore", sep = "")
dir.create(Path_correlate_RevComp_Rcore_up)
Path_correlate_RevComp_Rcore_Replichore_up <- paste(Path_correlate_repo_up, "/RevComp_Rcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore_up)


#===================================

#This creates all of the text and graphical output directories for both the organism specific output, and the correlation
#repositories later on for statistical analysis of the Downstream 4mers
#===================================

#macro repositories for text and graphical output
Path_graphical_output_down <- paste(Path_output_downstream, "/Output_Images", sep = "")
dir.create(Path_graphical_output_down)
Path_text_output_down <- paste(Path_output_downstream, "/Output_Text", sep = "")
dir.create(Path_text_output_down)

#chromosome correlation repository
Path_correlate_Chromosome_down <- paste(Path_correlate_repo_down, "/Chromosome", sep = "")
dir.create(Path_correlate_Chromosome_down)
Path_correlate_RevComp_Chromosome_down <- paste(Path_correlate_repo_down, "/RevComp_Chromosome", sep = "")
dir.create(Path_correlate_RevComp_Chromosome_down)

#left replichore correlation repository
Path_correlate_Lcore_down <- paste(Path_correlate_repo_down, "/Lcore", sep = "")
dir.create(Path_correlate_Lcore_down)
Path_correlate_Lcore_replichore_down <- paste(Path_correlate_repo_down, "/Lcore_Replichore", sep = "")
dir.create(Path_correlate_Lcore_replichore_down)
Path_correlate_RevComp_Lcore_down <- paste(Path_correlate_repo_down, "/RevComp_Lcore", sep = "")
dir.create(Path_correlate_RevComp_Lcore_down)
Path_correlate_RevComp_Lcore_Replichore_down <- paste(Path_correlate_repo_down, "/RevComp_Lcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Lcore_Replichore_down)

#right replichore correlation repository
Path_correlate_Rcore_down <- paste(Path_correlate_repo_down, "/Rcore", sep = "")
dir.create(Path_correlate_Rcore_down)
Path_correlate_Rcore_replichore_down <- paste(Path_correlate_repo_down, "/Rcore_Replichore", sep = "")
dir.create(Path_correlate_Rcore_replichore_down)
Path_correlate_RevComp_Rcore_down <- paste(Path_correlate_repo_down, "/RevComp_Rcore", sep = "")
dir.create(Path_correlate_RevComp_Rcore_down)
Path_correlate_RevComp_Rcore_Replichore_down <- paste(Path_correlate_repo_down, "/RevComp_Rcore_Replichore", sep = "")
dir.create(Path_correlate_RevComp_Rcore_Replichore_down)
#===================================

