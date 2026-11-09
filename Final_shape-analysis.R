
install.packages("vegan")
install.packages("shapeR")
install.packages("readr")
install.packages("ggplot2")
library("vegan")
library("shapeR")
library("readr")
library("ggplot2")


#### 1) Extracting otolith outlines ####-------------------------------------------------------------

#read in metadata
Metadata <- read_csv("Final_shape-analysis.csv")

#create the shape object
shape <- shapeR(project.path = getwd(), info.file = "Final_shape-analysis.csv")

#extract otolith outlines from images
shape <- detect.outline(shape, threshold = 0.2, write.outline.w.org = TRUE)

warnings()

#smooth otolith outlines
shape <- smoothout(shape, n = 100)

#check success of outline extraction

#create extracted_names object containing successfully extracted picnames
#these are stored in the slot @outline.list
extracted_names <- unlist(lapply(shape@outline.list, names))
length(extracted_names) 
#number of successfully extracted outlines matches number of rows in metadata
#one successfully extracted otolith for each row of otolith metadata



#### 2) Extracting Fourier Descriptors (FDs) ####-----------------------------------------------------------

#update shape object with shape coefficients (wavelet and fourier)
shape <- generateShapeCoefficients(shape)

#connect shape information to metadata
shape <- enrich.master.list(shape)

#standardize coefficients for fish total length with Bonferroni adjustment
shape <- stdCoefs(shape, classes = "location", "total_length_cm", bonferroni = TRUE)

#the shape object stores the otoliths in the order they were extracted 
#this is not the same order as the metadata rows
#to combine the FDs with the correct otoliths and their metadata, merge by picname

#make dataframe for the standardised FDs
fd_std <- as.data.frame(shape@fourier.coef.std)

#add picname column to attach picnames from extracted_names 
#this is the order they are stored in the shape object
fd_std$picname <- extracted_names

#merge with metadata by picname
rostrata_FDs <- merge(fd_std, Metadata, by = "picname")

#write FDs csv
write.csv(rostrata_FDs, file ="A-rostrata_FDs.csv", row.names = FALSE)



#### 3) Extracting morphometrics ####------------------------------------------------

#get morphometrics from shape object and assign to new morphs object
morphs <- getMeasurements(shape)

#make morphs dataframe and add picname column from extracted_names 
#same process as for FDs to save in correct order
morphs_df <- as.data.frame(morphs)
morphs_df$picname <- extracted_names

#merge with metadata by picname
rostrata_morphs <- merge(morphs_df, Metadata, by = "picname")

#write morphometrics csv
write.csv(rostrata_morphs, file = "A-rostrata_morphs.csv", row.names = FALSE)

#double check correct morphometrics have been assigned to correct fish
#morph should increase with fish total length

#assign morphometrics
morphometrics <- c("otolith.area", "otolith.length", "otolith.width", "otolith.perimeter")

#for loop to print a plot for each morph - TL relationship
for (morph in morphometrics) {
  morph_plot <- ggplot(rostrata_morphs, aes(x = total_length_cm, y = .data[[morph]])) +
    geom_point() +
    labs(x = "Fish total length (cm)", y = morph, title = paste(morph, "- fish total length")) +
    theme_classic()
  print(morph_plot)
}



#### 4) Mean otolith outlines plot ####----------------------------------------------

#want to change order of locations from alphabetical to latitudinal

#get master.list
ml <- slot(shape, "master.list")
#set the desired geographic order of locations
ml$location <- factor(ml$location, levels = c("GL", "MAR", "RT", "FI", "SG"))
#put the modified master.list back into the shape object
slot(shape, "master.list") <- ml

#set location colours
location_colours <- c(
  "GL"   = "#FFCC00",  #yellow
  "MAR"  = "#33CC00",  #green
  "RT"   = "#33CCFF",  #blue
  "FI"   = "#9933FF",  #purple
  "SG"   = "#FF33CC"   #pink
)

#edit plot fourier shape function to change how the legend looks
#otherwise legend will overlap with otolith outlines

plotFourierShape2 <- plotFourierShape #make new plot fourier shape
fix(plotFourierShape2)
#change the legend line to: 
#legend("bottomleft", lev.classes, cex = 1.1, lty = lty, col = col, bty = "n", ...)
#changes cex from 0.9 to 1.1 for larger legend text to match size of angle numbers
#add bty = "n" to remove legend box 

#save plot as pdf

#create pdf with best dimensions for otolith proportions (10x6), name it Mean_otolith_outlines
pdf("Mean_otolith_outlines.pdf", width = 10, height = 6)

#draw the plot into the pdf
plotFourierShape2(shape, "location", show.angle = TRUE, lwd = 1.8, lty = 1, col = location_colours)
#line width 1.8 for better discrimination among lines

#close and save the pdf
dev.off()
