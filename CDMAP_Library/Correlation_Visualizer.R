#Output Conversion script
#David L Patton / 800728881
#Sung Lab
library("lattice")


output_pic_name <- paste(matrix_name, ".png", sep = "")
#output_pic_name <- "test.png"
GC_flag <- grepl("GC", matrix_name_graph)
dlong_sort <- matrix(0L, nrow = 0, ncol = 3)
sort_array <- matrix(0L, nrow = 0, ncol = 3)


d <- data.frame(output_data_matrix)
#d <- data.frame(GC_corr_matrix)
d <- t(output_data_matrix)
dlin <- linspace(-1, 1, n = 30)
dlin <- round(dlin, digits = 3)

setwd(path_correlate_output)

png(output_pic_name, width = 1280, height = 1080, units = "px", res = 150)



#Current Correlation Visualizer Script (needs to be tested and debugged)

rgb.palette <- colorRampPalette(c("red","white", "blue"), space = "Lab")#((length(dlin)*2)+1)
breaks=c(seq(-1,-0,length=length(dlin)),0,seq(0,1,length=length(dlin)))

heatplot <- levelplot(GC_corr_matrix, col.regions=rgb.palette, breaks= breaks, scales=list(x=list(rot=45), cex=1.5),
                      xlab = "", ylab = "", main = matrix_name_graph,
                      aspect="fill",
                      colorkey=list(at=as.numeric(factor(c(seq(dlin)))),
                                    labels=as.character(dlin),
                                    col=(rgb.palette)
                                    ))
#Current Working version
# heatplot <- levelplot(GC_corr_matrix, col.regions=rgb.palette, scales=list(x=list(rot=45)),
#                       xlab = "", ylab = "", main = matrix_name_graph,
#                       aspect="fill",
#                       colorkey=list(at=as.numeric(factor(c(seq(dlin)))),
#                                     labels=as.character(dlin),
#                                     col=(rgb.palette)))

print(heatplot)
  dev.off()