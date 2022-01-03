#testcompare (ori_bp < term_bp)
#query_Rcore <<- matrix( nrow =1, ncol =3)
#query_Lcore <<- matrix( nrow =1, ncol =3)

print("enter loop(SCENARIO 1)")
i <- 1
flagcheck <- ""

#Full format output replichore matrices
#=======================================
for(i in 1:length(context_output_full_matrix[,1]))  
{
  
  #Explicit casting of BaseCall position, ORI, and Terminus as numeric for Comparison
  #test_pos <- context_output_full_matrix[i,1]
  mut_position <- as.numeric(context_output_full_matrix[i,1])
  terminus <- as.numeric(term_bp)
  origin <- as.numeric(ori_bp)
  
  #boolean operators to check for left or right core placement 
  Bool1 <- isTRUE(mut_position >= origin && mut_position <= terminus) #Right Core
  Bool2 <- isTRUE(mut_position <= origin || mut_position >= terminus) #Left Core
  
  if(Bool1)
  {
    flagcheck <- "1"
  }
  if(Bool2) 
  {
    flagcheck <- "2"
  }
  
  switch(flagcheck,
         "1" = query_Rcore_full <- rbind(query_Rcore_full, context_output_full_matrix[i,]),
         "2" = query_Lcore_full <- rbind(query_Lcore_full, context_output_full_matrix[i,])
  )
  i=i+1
  
}

return()

#return(query_Lcore)
#return(query_Rcore)