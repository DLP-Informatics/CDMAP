#Heatmap Correlations in R using dynamic directories

packages <- c("seqinr", "BiocManager", "pracma", "beepr", "lattice", "tidyverse", "vcfR", "stringr", "reshape2", "scales")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#BiocManager::install("genbankr", force = TRUE)
#library("genbankr")

username <- Sys.info()[7]

MainDir <- paste("/Users/", username, "/Desktop/CDMAP", sep = "")
LibDir <- paste("/Users/", username, "/Desktop/CDMAP/CDMAP_Library", sep = "")
UserOutput <- paste("/Users/", username, sep ="")
path_to_FASTA  <- paste(UserOutput,"/Desktop/CDMAP_Output/GC_directory/FASTA", sep = "")



setwd(LibDir)
source("GC_Analysis_InitDirTriplet.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GCcontent.r")
setwd(LibDir)
source("GC_Analysis.R")

######## Upstream Analysis ##################
setwd(LibDir)
source("GC_Analysis_InitDirUpA.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirUpC.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirUpG.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirUpT.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")


######## End Upstream Analysis ##################

######## Downstream Analysis ##################
setwd(LibDir)
source("GC_Analysis_InitDirDownA.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirDownC.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirDownG.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")

setwd(LibDir)
source("GC_Analysis_InitDirDownT.r")
path_analyze <- path_to_chrome
setwd(LibDir)
source("GC_Analysis.R")


######## END Downstream Analysis ##################




##Scratch space======================##

# check if a string contains a variable
# grep('Cat', testvar) where testvar <- "This is a Cat."
# note: grep will return 1 (true) or integer(0) for false

#list files in a directory#
# list.files(testdir), where testdir is a path to a given directory

#access the ith object in list.files as a string
#list.files(testdir)[2], where testdir is a directory path, and 2 is the 2nd element of the array

#round function for the heatmaps (for the titles)
# round(varname)

