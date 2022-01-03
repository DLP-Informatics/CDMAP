#testcompare (ori_bp > term_bp)


  print("enter loop(SCENARIO 2) ")
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
    Bool1 <- isTRUE(mut_position <= terminus || mut_position >= origin) #Right Core
    Bool2 <- isTRUE(mut_position >= terminus && mut_position <= origin) #Left Core
    
    #print(i)
    #print(context_output_full_matrix[i,1] > term_bp & context_output_full_matrix[i,1] < ori_bp)
    #print(context_output_full_matrix[i,1] > ori_bp | context_output_full_matrix[i,1] < term_bp)
    
    # Check if its on the right core (ORI to Terminus)
    if(Bool1) 
    {
      flagcheck <- "1"
    }
    
    #Check condition for the Left core (Terminus to ORI)
    if(Bool2)
    {
      flagcheck <- "2"
    }
    #else
    #{
    #  flagcheck <- "1"
    #}
    
    
    switch(flagcheck,
           "1" = query_Rcore_full <- rbind(query_Rcore_full, context_output_full_matrix[i,]), #Assign to the Right Core
           "2" = query_Lcore_full <- rbind(query_Lcore_full, context_output_full_matrix[i,])  #Assign to the Left Core
    )
    #i=i+1
    
  }
  #=======================================
  

return()  
  
#return(query_Lcore)
#return(query_Rcore)