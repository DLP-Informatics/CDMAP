#GC4C.r
# Purpose: This R script is an augmentation of GWTC.r and Rev_GWTC.r. In this script, for a given
# Chromosome we are obtaining both the Genome Coding 3mer (triplet) Count i.e GC3C and the 
# Genomic Intergenic 4mer Count (GC4C) to compare and calculate mutation rates of 3 fold and 4fold
# sites upstream and downstream of each gene. We take the Genbank file, and for each Gene in the genbank
# file, we compute the 3mer, 4mer (up), and 4mer(down) triplet counts for each chromosome.

#DevNotes
#=========================

#testobj1 <- Cleaned reference fasta sequence
#testobj2 <- reference amino acid sequence

#command to split and unlist them
#testsplit1 <- unlist(strsplit(testobj1, ""))

#It turns out that we need to seperate codon usage out into '+' (forward) and '-' (reverse) strand orientation features, and recombine them
#with respect to a given 5->3 or 3->5 orientation.

#we retrieve this information with the command;
# orientation <- OrganismGB$FEATURES[[i]]$strand

#then implement the following logic

#if(orientation == '+')
#{
  #run forward pipeline
#}
#if(orientation == '-')
#{
  #run reverse orientation pipeline
#}

#Then after each strand has been tablulated, They both need to be converted into their opposite matrices (5-3/ 3-5) and then
#merged with the other half.

#=========================

starttime <- Sys.time()

cols <- c("T", "G", "C",	"A")
rows <- c("T[X]T", "T[X]G","T[X]C","T[X]A","G[X]T","G[X]G","G[X]C","G[X]A",
          "C[X]T","C[X]G","C[X]C","C[X]A", "A[X]T","A[X]G", "A[X]C","A[X]A")

GC3C <- matrix(0L, nrow =16, ncol = 4) #Matrix to store Genome Coding Triplet count
rownames(GC3C) <- rows
colnames(GC3C) <- cols

Rev_GC3C <- GC3C
GC3C_Left <- GC3C
GC3C_Right <- GC3C
Rev_GC3C_Left <- GC3C
Rev_GC3C_Right <- GC3C

TestTotalMatrix <- GC3C



GC4C_Base <- matrix(0L, nrow =16, ncol = 4) #Base matrix for storing the Genome Coding 4-mer counts
rownames(GC4C_Base) <- rows
colnames(GC4C_Base) <- cols

upstream_matrix <- matrix(0L, nrow =16, ncol = 4) #Base matrix for storing the Genome Coding 4-mer counts
rownames(upstream_matrix) <- rows
colnames(upstream_matrix) <- cols

downstream_matrix <- matrix(0L, nrow =16, ncol = 4) #Base matrix for storing the Genome Coding 4-mer counts
rownames(downstream_matrix) <- rows
colnames(downstream_matrix) <- cols

comp_upstream_matrix <- matrix(0L, nrow =16, ncol = 4) #Base matrix for storing the Complement of Genome Coding 4-mer counts
rownames(upstream_matrix) <- rows
colnames(upstream_matrix) <- cols

comp_downstream_matrix <- matrix(0L, nrow =16, ncol = 4) #Base matrix for storing the Complement of Genome Coding 4-mer counts
rownames(downstream_matrix) <- rows
colnames(downstream_matrix) <- cols

#The Genome Wide 4-mer Count upstream matrices
#====================
GC4C_Tup_Left <- GC4C_Base
GC4C_Gup_Left <- GC4C_Base
GC4C_Cup_Left <- GC4C_Base
GC4C_Aup_Left <- GC4C_Base
GC4C_Tup_Left_Comp <- GC4C_Base
GC4C_Gup_Left_Comp <- GC4C_Base
GC4C_Cup_Left_Comp <- GC4C_Base
GC4C_Aup_Left_Comp <- GC4C_Base
GC4C_Tup_Right <- GC4C_Base
GC4C_Gup_Right <- GC4C_Base
GC4C_Cup_Right <- GC4C_Base
GC4C_Aup_Right <- GC4C_Base
GC4C_Tup_Right_Comp <- GC4C_Base
GC4C_Gup_Right_Comp <- GC4C_Base
GC4C_Cup_Right_Comp <- GC4C_Base
GC4C_Aup_Right_Comp <- GC4C_Base
#====================

#The Genome Wide 4-mer Count downstream matrices
#====================
GC4C_Tdown_Left <- GC4C_Base
GC4C_Gdown_Left <- GC4C_Base
GC4C_Cdown_Left <- GC4C_Base
GC4C_Adown_Left <- GC4C_Base
GC4C_Tdown_Right <- GC4C_Base
GC4C_Gdown_Right <- GC4C_Base
GC4C_Cdown_Right <- GC4C_Base
GC4C_Adown_Right <- GC4C_Base

GC4C_Tdown_Left_Comp <- GC4C_Base
GC4C_Gdown_Left_Comp <- GC4C_Base
GC4C_Cdown_Left_Comp <- GC4C_Base
GC4C_Adown_Left_Comp <- GC4C_Base
GC4C_Tdown_Right_Comp <- GC4C_Base
GC4C_Gdown_Right_Comp <- GC4C_Base
GC4C_Cdown_Right_Comp <- GC4C_Base
GC4C_Adown_Right_Comp <- GC4C_Base
#====================


OrganismGB <- parseGenBank(Path_GBFile) #parsing the genbank file and converting to a dataframe

setwd(Path_to_scripts)
print("Starting the timer!")

#i <<- 2 #ignore initial entry start/end is the entire sequence
featlength <- length(OrganismGB$FEATURES) #number of gene features contained in the genbank file

startindex <- c() #array of all all gene starting positions 
endindex <- c() #array of all gene ending positions
featlengthindex <- c() #array of the length of all gene features

geneFeature_matrix <- matrix(0L, nrow = featlength, ncol = 3) #Base matrix for storing the Complement of Genome Coding 4-mer counts
gene_cols <- c("Coding Region Start Position", "Coding Region End Position", "Coding Region Length")
colnames(geneFeature_matrix) <- gene_cols

#======== DEV NOTES ==========
#use this feature to extract the gene coding features 
# OrganismGB$FEATURES[[i]]$type
#if 'gene', then execute, if not 'gene', then i+1
#=============================


flagcheck <<- ''
#ticker <- 0
for(i in  1:featlength)
{
  #ticker <- ticker + 1
  if(i == 1)
  {
    #print("skipping first entry")
    next()
  }
  featuretype <- OrganismGB$FEATURES[[i]]$type
  print(paste("feature #: ", i, " feature type: ", featuretype, sep = ""))
  
  #Index the start and end of each gene feature in the Genbank file, skipping the first entry
  #Checks for duplicate entry in the startindex array, if true skips to next entry in gbk file
  if(featuretype != "gene") #is.element(featend, endindex)
  {
    print(paste("Feature", i, "skipped. Is not a gene", sep = " "))
    next()
  }
  dataobj <- OrganismGB[["FEATURES"]][[i]]
  matobj <- data.matrix(OrganismGB[["FEATURES"]][[i]]) #extracting an individual gene feature from the data frame
  #matobj <- testlist[[1]]$type
  
  #Grab the start and end position of the coding region
  featstart <- as.numeric(matobj[1,2]) #gene nucleotide start position
  featend <- as.numeric(matobj[1,3]) #gene nucelotide end position

  
  #Boolean Structure to assign to left and right replichore
  if(ori_bp < term_bp)
  {
    Bool1 <- isTRUE(featstart >= ori_bp && featstart <= term_bp) #Right Core
    Bool2 <- isTRUE(featstart <= ori_bp || featstart >= term_bp) #Left Core
    if(Bool1)
    {
      flagcheck <- "Right"
    }
    if(Bool2) 
    {
      flagcheck <- "Left"
    }
  }
  if(ori_bp > term_bp)
  {
    Bool1 <- isTRUE(featstart <= term_bp || featstart >= ori_bp) #Right Core
    Bool2 <- isTRUE(featstart >= term_bp && featstart <= ori_bp) #Left Core
    if(Bool1) 
    {
      flagcheck <- "Right"
    }
    if(Bool2)
    {
      flagcheck <- "Left"
    }
  }


  startindex <- c(startindex, featstart) 
  endindex <- c(endindex, featend) 
  
  if(featstart == 1)
  {
    upNuc <- RefSeq_arr[length(RefSeq_arr)] 
    downNuc <- RefSeq_arr[featend+1]
  } else if(featend == length(RefSeq_arr))
  {
    upNuc <- RefSeq_arr[featstart-1] 
    downNuc <- RefSeq_arr[1]
  }else
  {
  upNuc <- RefSeq_arr[featstart-1] #leftmost upstream nucleotide 4mer position in gene coding region
  downNuc <- RefSeq_arr[featend+1] # rightmost downstream nucleotide 4mer position in gene coding region
  }
  
  
  #Checks for null reference, skips if a null reference is induced
  if(identical(upNuc, character(0)) | identical(downNuc, character(0)) ) #is.element(featend, endindex)
  {
    next()
  }

  
  gene_arr <- RefSeq_arr[featstart:featend] #isolation of a specific gene feature
  gene_end <- length(gene_arr)-1 
  featlengthindex <- c(featlengthindex, length(gene_arr)) #length of the gene feature
  print(paste("Replichore: ", flagcheck, sep = ""))
  print(paste("analyzing feature ", as.character(i), " with length: ", as.character(length(gene_arr)), sep = ""))

  setwd(Path_text_output)
  options(max.print=999999)
  gc3cOrgobj <- paste(organism, "Nucleotide_Coding_Regions", sep = "_")
  #zz <- file(gc3cOrgobj,"w")
  sink(gc3cOrgobj, append = TRUE)
  print(paste("Replichore: ", flagcheck, sep = ""))
  print(paste("start position: ", featstart, " End Position: ", featend, sep = ""))
  print(paste("analyzing feature ", as.character(i), " with length: ", as.character(length(gene_arr)), sep = ""))
  print(gene_arr)
  sink()
  #close(gc3cOrgobj)
  
  
  #checks for null reference, skips if null detected
  if(is.na(upNuc))
  {
    i <- i+1
    next
  }
  if(is.na(downNuc))
  {
    i <- i+1
    next
  }  
  
  #checks if an ambiguous nucleotide 'N' assigned, if true, it skips it
  
  if(upNuc != 'A' & upNuc != 'T' & upNuc != 'C' & upNuc != 'G')
  {
    i <- i+1
    next
  }
  if(downNuc != 'A' & downNuc != 'T' & downNuc != 'C' & downNuc != 'G')
  {
    i <- i+1
    next
  }
  
  orientation <- OrganismGB$FEATURES[[i]]$strand
  
  if(orientation == '+')
  {
    setwd(Path_to_scripts)
    print("Strand in the forward '+' Orientation")
  }
  if(orientation == '-')
  {
    setwd(Path_to_scripts)
    print("Strand in the reverse '-' Orientation")
  }
  
  # GC3C <- matrix(0L, nrow =16, ncol = 4) #Matrix to store Genome Coding Triplet count
  # rownames(GC3C) <- rows
  # colnames(GC3C) <- cols
  # 
  # GC3C_Left <- GC3C
  # GC3C_Right <- GC3C
  # Rev_GC3C_Left <- GC3C
  # Rev_GC3C_Right <- GC3C
  # DiagArr <- c()
  

  j <- 2 #iterator instantiated with respect to center nucleotide, grabs j-1 and j+1 for left and right nucleotide.
  while(j <= gene_end)
  {
    iter <- j
    
    #Leading and Lagging# ========
    if(orientation == '+' && flagcheck == 'Right')
    {
      setwd(Path_to_scripts)
      #print("Strand in the forward '+' Orientation")
      source("GC4C_recorder_Forward.r")
    }
    if(orientation == '-' && flagcheck == 'Right')
    {
      setwd(Path_to_scripts)
      #print("Strand in the reverse '-' Orientation")
      source("GC4C_recorder_Reverse.r")
    }
    
    if(orientation == '+' && flagcheck == 'Left')
    {
      setwd(Path_to_scripts)
      #print("Strand in the forward '+' Orientation")
      source("GC4C_recorder_Reverse.r")
    }
    if(orientation == '-' && flagcheck == 'Left')
    {
      setwd(Path_to_scripts)
      #print("Strand in the reverse '-' Orientation")
      source("GC4C_recorder_Forward.r")
    }

    #4fold sites===========
    #setwd(Path_to_scripts)
    #source("GC4C_recorder.r")
    #======================
  

    j <- j+3 #iterates by 3, so it checks with respect to each nucleotide triplet.
  }
  #DiagArr <- sort(DiagArr)
  
}

GC3C <- GC3C_Right + Rev_GC3C_Left
Rev_GC3C <- GC3C_Left + Rev_GC3C_Right
  
GC4C_Aup <- GC4C_Aup_Left + GC4C_Aup_Right
GC4C_Cup <- GC4C_Cup_Left + GC4C_Cup_Right
GC4C_Gup <- GC4C_Gup_Left + GC4C_Gup_Right
GC4C_Tup <- GC4C_Tup_Left + GC4C_Tup_Right

GC4C_Aup_Comp <- GC4C_Aup_Left_Comp + GC4C_Aup_Right_Comp
GC4C_Cup_Comp <- GC4C_Cup_Left_Comp + GC4C_Cup_Right_Comp
GC4C_Gup_Comp <- GC4C_Gup_Left_Comp + GC4C_Gup_Right_Comp
GC4C_Tup_Comp <- GC4C_Tup_Left_Comp + GC4C_Tup_Right_Comp

GC4C_Adown <- GC4C_Adown_Left + GC4C_Adown_Right
GC4C_Cdown <- GC4C_Cdown_Left + GC4C_Cdown_Right
GC4C_Gdown <- GC4C_Gdown_Left + GC4C_Gdown_Right
GC4C_Tdown <- GC4C_Tdown_Left + GC4C_Tdown_Right

GC4C_Adown_Comp <- GC4C_Adown_Left_Comp + GC4C_Adown_Right_Comp
GC4C_Cdown_Comp <- GC4C_Cdown_Left_Comp + GC4C_Cdown_Right_Comp
GC4C_Gdown_Comp <- GC4C_Gdown_Left_Comp + GC4C_Gdown_Right_Comp
GC4C_Tdown_Comp <- GC4C_Tdown_Left_Comp + GC4C_Tdown_Right_Comp

GC3C_up <- GC4C_Aup+GC4C_Cup+GC4C_Gup+GC4C_Tup
GC3C_down <- GC4C_Adown+GC4C_Cdown+GC4C_Gdown+GC4C_Tdown

GC3C_doesNotEqual <- function(x){
                        nucdiff <- sum(GC3C_up)-sum(GC3C_down)
                        warning("Warning: there is a ", nucdiff, "nucleotide difference between upstream and downstream coding matrices.")
                        FALSE
}

if(sum(GC3C_up) != sum(GC3C_down))
{
  GC3C_doesNotEqual()
}



#Record all of the Output GC3C and GC4C matrices here
#==============================
output_GC3C <- paste("GC3C_", organism, ".csv", sep = "")
output_Rev_GC3C <- paste("Rev_GC3C_", organism, ".csv", sep = "")

output_GC3C_Left <- paste("GC3C_Left_", organism, ".csv", sep = "")
output_GC3C_Right <- paste("GC3C_Right_", organism, ".csv", sep = "")

output_Rev_GC3C_Left <- paste("Rev_GC3C_Left_", organism, ".csv", sep = "")
output_Rev_GC3C_Right <- paste("Rev_GC3C_Right_", organism, ".csv", sep = "")

output_GC3C_CUB <- paste("GC3C_", organism,"_CUB", ".txt", sep = "")
output_genePos <- paste("GenePositions_", organism, ".csv", sep = "")

#=========== UPSTREAM =======#
output_Aup_Left <- paste("GC4C_Aup_Left_", organism, ".csv", sep = "")
output_Cup_Left <- paste("GC4C_Cup_Left_", organism, ".csv", sep = "")
output_Gup_Left <- paste("GC4C_Gup_Left_", organism, ".csv", sep = "")
output_Tup_Left <- paste("GC4C_Tup_Left_", organism, ".csv", sep = "")
output_Aup_Right <- paste("GC4C_Aup_Right_", organism, ".csv", sep = "")
output_Cup_Right <- paste("GC4C_Cup_Right_", organism, ".csv", sep = "")
output_Gup_Right <- paste("GC4C_Gup_Right_", organism, ".csv", sep = "")
output_Tup_Right <- paste("GC4C_Tup_Right_", organism, ".csv", sep = "")
output_Aup <- paste("GC4C_Aup", organism, ".csv", sep = "")
output_Cup <- paste("GC4C_Cup", organism, ".csv", sep = "")
output_Gup <- paste("GC4C_Gup", organism, ".csv", sep = "")
output_Tup <- paste("GC4C_Tup", organism, ".csv", sep = "")

comp_output_Aup_Left <- paste("GC4C_Aup_Comp_Left_", organism, ".csv", sep = "")
comp_output_Cup_Left <- paste("GC4C_Cup_Comp_Left_", organism, ".csv", sep = "")
comp_output_Gup_Left <- paste("GC4C_Gup_Comp_Left_", organism, ".csv", sep = "")
comp_output_Tup_Left <- paste("GC4C_Tup_Comp_Left_", organism, ".csv", sep = "")
comp_output_Aup_Right <- paste("GC4C_Aup_Comp_Right_", organism, ".csv", sep = "")
comp_output_Cup_Right <- paste("GC4C_Cup_Comp_Right_", organism, ".csv", sep = "")
comp_output_Gup_Right <- paste("GC4C_Gup_Comp_Right_", organism, ".csv", sep = "")
comp_output_Tup_Right <- paste("GC4C_Tup_Comp_Right_", organism, ".csv", sep = "")
comp_output_Aup <- paste("GC4C_Aup_Comp_", organism, ".csv", sep = "")
comp_output_Cup <- paste("GC4C_Cup_Comp_", organism, ".csv", sep = "")
comp_output_Gup <- paste("GC4C_Gup_Comp_", organism, ".csv", sep = "")
comp_output_Tup <- paste("GC4C_Tup_Comp_", organism, ".csv", sep = "")


output_Aup_CUB <- paste("GC4C_Aup_", organism, "_CUB", ".txt", sep = "")
output_Cup_CUB <- paste("GC4C_Cup_", organism, "_CUB", ".txt", sep = "")
output_Gup_CUB <- paste("GC4C_Gup_", organism, "_CUB", ".txt", sep = "")
output_Tup_CUB <- paste("GC4C_Tup_", organism, "_CUB", ".txt", sep = "")
#======================================#

#======DOWNSTREAM=====================#

output_Adown_Left <- paste("GC4C_Adown_Left_", organism, ".csv", sep = "")
output_Cdown_Left <- paste("GC4C_Cdown_Left_", organism, ".csv", sep = "")
output_Gdown_Left <- paste("GC4C_Gdown_Left_", organism, ".csv", sep = "")
output_Tdown_Left <- paste("GC4C_Tdown_Left_", organism, ".csv", sep = "")
output_Adown_Right <- paste("GC4C_Adown_Right_", organism, ".csv", sep = "")
output_Cdown_Right <- paste("GC4C_Cdown_Right_", organism, ".csv", sep = "")
output_Gdown_Right <- paste("GC4C_Gdown_Right_", organism, ".csv", sep = "")
output_Tdown_Right <- paste("GC4C_Tdown_Right_", organism, ".csv", sep = "")
output_Adown <- paste("GC4C_Adown_", organism, ".csv", sep = "")
output_Cdown <- paste("GC4C_Cdown_", organism, ".csv", sep = "")
output_Gdown <- paste("GC4C_Gdown_", organism, ".csv", sep = "")
output_Tdown <- paste("GC4C_Tdown_", organism, ".csv", sep = "")

comp_output_Adown_Left <- paste("GC4C_Adown_Comp_Left_", organism, ".csv", sep = "")
comp_output_Cdown_Left <- paste("GC4C_Cdown_Comp_Left_", organism, ".csv", sep = "")
comp_output_Gdown_Left <- paste("GC4C_Gdown_Comp_Left_", organism, ".csv", sep = "")
comp_output_Tdown_Left <- paste("GC4C_Tdown_Comp_Left_", organism, ".csv", sep = "")
comp_output_Adown_Right <- paste("GC4C_Adown_Comp_Right_", organism, ".csv", sep = "")
comp_output_Cdown_Right <- paste("GC4C_Cdown_Comp_Right_", organism, ".csv", sep = "")
comp_output_Gdown_Right <- paste("GC4C_Gdown_Comp_Right_", organism, ".csv", sep = "")
comp_output_Tdown_Right <- paste("GC4C_Tdown_Comp_Right_", organism, ".csv", sep = "")
comp_output_Adown <- paste("GC4C_Adown_Comp_", organism, ".csv", sep = "")
comp_output_Cdown <- paste("GC4C_Cdown_Comp_", organism, ".csv", sep = "")
comp_output_Gdown <- paste("GC4C_Gdown_Comp_", organism, ".csv", sep = "")
comp_output_Tdown <- paste("GC4C_Tdown_Comp_", organism, ".csv", sep = "")


output_Adown_CUB <- paste("GC4C_Adown_", organism, "_CUB", ".txt", sep = "")
output_Cdown_CUB <- paste("GC4C_Cdown_", organism, "_CUB", ".txt", sep = "")
output_Gdown_CUB <- paste("GC4C_Gdown_", organism, "_CUB", ".txt", sep = "")
output_Tdown_CUB <- paste("GC4C_Tdown_", organism, "_CUB", ".txt", sep = "")
#=================================#


path_GC3C <- paste(Path_text_output, "/GC3C_Triplet", sep = "")
dir.create(path_GC3C)
path_GC3C_CUB <- paste(Path_text_output, "/GC3C_CUB", sep = "")
dir.create(path_GC3C_CUB)

path_GC4C_output <- paste(Path_text_output_up, "/GC4C_Output", sep = "")
path_GC4C_upCUB <- paste(path_GC4C_output, "/GC4C_Upstream_CUB", sep = "")
path_GC4C_Chrome <- paste(path_GC4C_output, "/Chromosome", sep = "")
path_GC4C_output_Left <- paste(path_GC4C_output, "/Left", sep = "")
path_GC4C_output_Right <- paste(path_GC4C_output, "/Right", sep = "")
path_GC4C_up_Chrome <- path_GC4C_Chrome
path_GC4C_up_Left <- path_GC4C_output_Left
path_GC4C_up_Right <- path_GC4C_output_Right
dir.create(path_GC4C_output)
dir.create(path_GC4C_Chrome)
dir.create(path_GC4C_output_Left)
dir.create(path_GC4C_output_Right)
dir.create(path_GC4C_upCUB)

path_GC4C_output <- paste(Path_text_output_down, "/GC4C_Output", sep = "")
path_GC4C_downCUB <- paste(path_GC4C_output, "/GC4C_Downstream_CUB", sep = "")
path_GC4C_Chrome <- paste(path_GC4C_output, "/Chromosome", sep = "")
path_GC4C_output_Left <- paste(path_GC4C_output, "/Left", sep = "")
path_GC4C_output_Right <- paste(path_GC4C_output, "/Right", sep = "")
path_GC4C_down_Chrome <- path_GC4C_Chrome
path_GC4C_down_Left <- path_GC4C_output_Left
path_GC4C_down_Right <- path_GC4C_output_Right
dir.create(path_GC4C_output)
dir.create(path_GC4C_downCUB)

dir.create(path_GC4C_Chrome)
dir.create(path_GC4C_output_Left)
dir.create(path_GC4C_output_Right)

#Grab the Gene start and end locations, and the Gene feature length

for(i in 1:featlength)
{
  geneFeature_matrix[i,1] <- startindex[i] #i-1 is used since we skip the first feature in the GBK file, but still
  # have to store gene feature information in the first row position.
  geneFeature_matrix[i,2] <- endindex[i]   
  geneFeature_matrix[i,3] <- featlengthindex[i]

}
CodingRegionSize <- sum(featlengthindex)
print(paste("Final Size: ", CodingRegionSize, sep = ""))




setwd(path_GC3C)
write.csv(GC3C, output_GC3C)
write.csv(Rev_GC3C, output_Rev_GC3C)
write.csv(GC3C_Left, output_GC3C_Left)
write.csv(GC3C_Right, output_GC3C_Right)
write.csv(Rev_GC3C_Left, output_Rev_GC3C_Left)
write.csv(Rev_GC3C_Right, output_Rev_GC3C_Right)


write.csv(geneFeature_matrix, output_genePos)

setwd(path_GC4C_up_Chrome)
write.csv(GC4C_Aup, output_Aup)
write.csv(GC4C_Cup, output_Cup)
write.csv(GC4C_Gup, output_Gup)
write.csv(GC4C_Tup_Comp, comp_output_Tup)
write.csv(GC4C_Aup_Comp, comp_output_Aup)
write.csv(GC4C_Cup_Comp, comp_output_Cup)
write.csv(GC4C_Gup_Comp, comp_output_Gup)


setwd(path_GC4C_up_Left)
write.csv(GC4C_Aup_Left, output_Aup)
write.csv(GC4C_Cup_Left, output_Cup)
write.csv(GC4C_Gup_Left, output_Gup)
write.csv(GC4C_Tup_Left, output_Tup)
write.csv(GC4C_Tup_Left_Comp, comp_output_Tup_Left)
write.csv(GC4C_Aup_Left_Comp, comp_output_Aup_Left)
write.csv(GC4C_Cup_Left_Comp, comp_output_Cup_Left)
write.csv(GC4C_Gup_Left_Comp, comp_output_Gup_Left)


setwd(path_GC4C_up_Right)
write.csv(GC4C_Aup_Right, output_Aup)
write.csv(GC4C_Cup_Right, output_Cup)
write.csv(GC4C_Gup_Right, output_Gup)
write.csv(GC4C_Tup_Right, output_Tup)
write.csv(GC4C_Aup_Right_Comp, comp_output_Aup_Right)
write.csv(GC4C_Cup_Right_Comp, comp_output_Cup_Right)
write.csv(GC4C_Gup_Right_Comp, comp_output_Gup_Right)
write.csv(GC4C_Tup_Right_Comp, comp_output_Tup_Right)

setwd(path_GC4C_down_Chrome)
write.csv(GC4C_Adown, output_Adown)
write.csv(GC4C_Cdown, output_Cdown)
write.csv(GC4C_Gdown, output_Gdown)
write.csv(GC4C_Tdown, output_Tdown)
write.csv(GC4C_Tdown_Comp, comp_output_Tdown)
write.csv(GC4C_Adown_Comp, comp_output_Adown)
write.csv(GC4C_Cdown_Comp, comp_output_Cdown)
write.csv(GC4C_Gdown_Comp, comp_output_Gdown)
write.csv(GC4C_Tdown_Comp, comp_output_Tdown)

setwd(path_GC4C_down_Left)
write.csv(GC4C_Adown_Left, output_Adown)
write.csv(GC4C_Cdown_Left, output_Cdown)
write.csv(GC4C_Gdown_Left, output_Gdown)
write.csv(GC4C_Tdown_Left, output_Tdown)
write.csv(GC4C_Tdown_Left_Comp, comp_output_Tdown_Left)
write.csv(GC4C_Adown_Left_Comp, comp_output_Adown_Left)
write.csv(GC4C_Cdown_Left_Comp, comp_output_Cdown_Left)
write.csv(GC4C_Gdown_Left_Comp, comp_output_Gdown_Left)

setwd(path_GC4C_down_Right)
write.csv(GC4C_Adown_Right, output_Adown)
write.csv(GC4C_Cdown_Right, output_Cdown)
write.csv(GC4C_Gdown_Right, output_Gdown)
write.csv(GC4C_Tdown_Right, output_Tdown)
write.csv(GC4C_Tdown_Right_Comp, comp_output_Tdown_Right)
write.csv(GC4C_Adown_Right_Comp, comp_output_Adown_Right)
write.csv(GC4C_Cdown_Right_Comp, comp_output_Cdown_Right)
write.csv(GC4C_Gdown_Right_Comp, comp_output_Gdown_Right)
#==============================

#Record the Codon Usage Bias for each organism both up and downstream along with triplet
#=============================

# input_matrix <- GC3C
# output_path <- path_GC3C_CUB
# output_name <- output_GC3C_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# output_path <- path_GC4C_upCUB
# 
# input_matrix <- GC4C_Aup
# output_name <- output_Aup_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Cup
# output_name <- output_Cup_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Gup
# output_name <- output_Gup_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Tup
# output_name <- output_Tup_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# output_path <- path_GC4C_downCUB
# 
# input_matrix <- GC4C_Adown
# output_name <- output_Adown_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Cdown
# output_name <- output_Cdown_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Gdown
# output_name <- output_Gdown_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
# 
# input_matrix <- GC4C_Tdown
# output_name <- output_Tdown_CUB
# setwd(Path_to_scripts)
# source("CUB_Output.r")
#=================================

#endtime <- Sys.time()
#runtime <- endtime - starttime

#print(paste("The total runtime for GC3C and GC4C calculation: ", runtime, sep = ""))
beep(3)
