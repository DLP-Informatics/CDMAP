##### partitionReplichores.r #######
#=================================
#Author: David Logan Patton
# Objective: Purpose: This script is responsible for taking the output context_output_full, 
#context_output_upstream, and context_output_downstream matrices and then to partition them 
#accordingly into left and right replichore outputs for each matrix. 



#========================================================================
# CONVERSION INTO LEFT AND RIGHT REPLICHORES
#========================================================

# to convert them into left and right replichores, we have to
# handle it in a multicase if statement

query_Rcore_full <<- matrix( nrow =0, ncol =9)
colnames(query_Rcore_full) <- context_col_names


query_Lcore_full <<- matrix( nrow =0, ncol =9)
colnames(query_Lcore_full) <- context_col_names

#CASES
#Scenario 1: Origin < Terminus
i <-1
if(ori_bp < term_bp)
{
  setwd(Path_to_scripts)
  source("Core_OriTerm.r")
}


#Scenario 2: Origin > Terminus
i <- 1
if(ori_bp > term_bp)
{
  setwd(Path_to_scripts)
  source("Core_TermOri.r")
} 

setwd(Path_Basecall_output)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")


setwd(Path_Basecall_output_up)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")


setwd(Path_Basecall_output_down)
write.csv(query_Lcore_full, "Left_Context_Core_Full.csv")
write.csv(query_Rcore_full, "Right_Context_Core_Full.csv")

