#Output Conversion script
#David L Patton / 800728881
#Sung Lab
library("lattice")


output_pic_name <- paste(matrix_name, ".png", sep = "")
#output_pic_name <- "lattice_test_image.png"

d <- t(output_data_matrix)
dlin <- linspace(min(d), max(d), n = 20) #this serves as the key scale for the heatmap label
dlin <- round(dlin, digits = 3)
## FIGURE OUT HOW TO MANUALLY CREATE THE NUMERIC SCALE IN LATTICE##
setwd(Path_image_output)
png(output_pic_name, width = 800, height = 600, units = "px")
rgb.palette <- colorRampPalette(c("white", "blue"), space = "Lab")
heatplot <- levelplot(d, col.regions=rgb.palette, scales=list(x=list(rot=90)),
                      xlab = "", ylab = "", main = matrix_name_graph,
                      aspect="fill",
                      colorkey=list(at=as.numeric(factor(c(seq(dlin)))),
                                    labels=as.character(dlin),
                                    col=(rgb.palette)))

print(heatplot)
dev.off()



