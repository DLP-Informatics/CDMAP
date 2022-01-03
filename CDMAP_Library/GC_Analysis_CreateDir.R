flag <- readline("Home or guest?")

path_to_scripts <- "/Users/triadge/Desktop/CDMP_Output/R_workspace/Context_Project/Context_Pipeline_0.1.0"

if(flag == "Home")
{
  path_to_FASTA  <- "/Users/triadge/Desktop/CDMP_Output/GC_directory/FASTA"
  path_to_scripts <- "/Users/triadge/Desktop/CDMP_Output/R_workspace/Context_Project/Context_Pipeline_0.1.0"
  
  path_to_chrome  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Chromosome"
  path_to_Lcore  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Lcore"
  path_to_Rcore <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Rcore"
  path_to_Lcore_Rep  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Lcore_Replichore"
  path_to_Rcore_Rep <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Rcore_Replichore"
  
  path_to_RevComp_chrome  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/RevComp_Chromosome"
  path_to_RevComp_Lcore  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/RevComp_Lcore"
  path_to_RevComp_Rcore <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/RevComp_Rcore"
  path_to_RevComp_Lcore_Rep  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/RevComp_Lcore_Replichore"
  path_to_RevComp_Rcore_Rep <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/RevComp_Rcore_Replichore"
  
  
  #path_to_chrome_rate  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Chromosome_Rates"
  #path_to_Lcore_rate  <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Lcore_Rates"
  #path_to_Rcore_rate <- "/Users/triadge/Desktop/CDMP_Output/Correlation_Repository/Rcore_Rates"
}

if(flag == "Guest")
{
  path_to_FASTA  <- readline("please specify the path to your FASTA files: ")
  path_to_chrome  <- readline("please specify the path to your Chromosome output: ")
  path_to_Lcore  <- readline("please specify the path to your Left Replichore output: ")
  path_to_Rcore <- readline("please specify the path to your Right Replichore output: ")
  path_to_chrome_rate  <- readline("please specify the path to your Clockwise (5 to 3) Chromosome mutation rate output: ")
  path_to_Lcore_rate  <- readline("please specify the path to your Left (5 to 3) Replichore mutation rate output: ")
  path_to_Rcore_rate <- readline("please specify the path to your Right (5 to 3) Replichore mutation rate output: ")
  path_to_Lcore_rate  <- readline("please specify the path to your replichore specific Left (5 to 3) Replichore Mutation Rate output: ")
  path_to_Rcore_rate <- readline("please specify the path to your replichore specific Right (5 to 3) Replichore Mutation Rate output: ")
  
  path_to_RevComp_chrome  <- readline("please specify the path to your Counterclockwise (3 to 5) Chromosome mutation rate output: ")
  path_to_RevComp_Lcore  <- readline("please specify the path to your Left (3 to 5) Replichore mutation rate output: ")
  path_to_RevComp_Rcore <- readline("please specify the path to your Right (3 to 5) Replichore mutation rate output: ")
  path_to_RevComp_Lcore_Rep  <- readline("please specify the path to your replichore specific Left (3 to 5) Replichore Mutation Rate output: ")
  path_to_RevComp_Rcore_Rep <- readline("please specify the path to your replichore specific Right (3 to 5) Replichore Mutation Rate output: ")
  
  
  #path_to_chrome_rate  <- readline("please specify the path to your Chromosome mutation rate output: ")
  #path_to_Lcore_rate  <- readline("please specify the path to your Left (5 to 3) Replichore mutation rate output: ")
  #path_to_Rcore_rate <- readline("please specify the path to your Right (5 to 3) Replichore mutation rate output: ")
  
}