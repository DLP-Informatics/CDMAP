#GC 1:N Triplet analysis




#path_analyze_rates <- path_to_chrome_rate
setwd(LibDir)
name_addon <- "Chromosome"
source("Correlation_Script.r")
setwd(LibDir)
#source("ChiSquare_Script.r")

setwd(LibDir)
path_analyze <- path_to_RevComp_chrome
setwd(LibDir)
name_addon <- "Counterclockwise_Chromosome"
source("Correlation_Script.r")

# #----- insert mutation simulation, and proportion graphing scripts here
# #==========#
# 
#==Left Replichore Analysis==#
setwd(LibDir)
path_analyze <- path_to_Lcore
name_addon <- "Left_Replichore"
source("Correlation_Script.r")
# #source("ChiSquare_Script.r")
# 
setwd(LibDir)
path_analyze <- path_to_RevComp_Lcore
name_addon <- "Left_Replichore_ReverseCompliment"
source("Correlation_Script.r")
# source("ChiSquare_Script.r")
# 
setwd(LibDir)
path_analyze <- path_to_Lcore_Rep
name_addon <- "Left_Replichore_RepSpecific"
source("Correlation_Script.r")
# 
setwd(LibDir)
path_analyze <- path_to_RevComp_Lcore_Rep
name_addon <- "Left_Replichore_RepSpecific_ReverseCompliment"
source("Correlation_Script.r")
# 
# #----- insert mutation simulation, and proportion graphing scripts here
# #==========#
# 
# #==Right Replichore Analysis==#
setwd(LibDir)
path_analyze <- path_to_Rcore
name_addon <- "Right Replichore"
source("Correlation_Script.r")
# source("ChiSquare_Script.r")
# #----- insert mutation simulation, and proportion graphing scripts here
# 
setwd(LibDir)
path_analyze <- path_to_RevComp_Rcore
name_addon <- "Right_Replichore_ReverseCompliment"
source("Correlation_Script.r")
# source("ChiSquare_Script.r")
# 
setwd(LibDir)
path_analyze <- path_to_Rcore_Rep
name_addon <- "Right_Replichore_RepSpecific"
source("Correlation_Script.r")
# 
setwd(LibDir)
path_analyze <- path_to_RevComp_Rcore_Rep
name_addon <- "Right_Replichore_RepSpecific_ReverseCompliment"
source("Correlation_Script.r")
# #==========#