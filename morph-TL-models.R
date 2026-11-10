#script to test the relationships between Antimora rostrata otolith morphometrics and TL with sampling location
#uses 3 nested linear models for each morphometric and compares them using anova
#plots the residuals vs fitted and QQ for each model for supplementary materials
#plots the morph-TL relationships for first results figure

install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")



#### 1) Prepare data ####--------------------------------------------------------

#read in morphometrics calculated during shapeR shape analysis
rostrata_morphs <- read_csv("A-rostrata_morphs.csv")

#define morphometric columns
morphometrics <- c("otolith.area", "otolith.length", "otolith.width", "otolith.perimeter")

#make location a factor
rostrata_morphs$location <- as.factor(rostrata_morphs$location)



#### 2) Morph~TL relationships with location ####-----------------------------------

#nested linear models
#1. morphometric ~ TL
#2. morphometric ~ TL + location 
#3. morphometric ~ TL * location

#anova used to compare the models
#m1 and m2 comparison - does location affect the intercept of the morph~TL relationship?
#m2 and m3 comparison - does location affect the slope of the morph~TL relationship?

#test simplest lm (m1) for all morphometrics
length_m1 <- lm(otolith.length ~ total_length_cm, data = rostrata_morphs)
plot(length_m1) #needs log adjustment
width_m1 <- lm(otolith.width ~ total_length_cm, data = rostrata_morphs)
plot(width_m1) #needs log adjustment
area_m1 <- lm(otolith.area ~ total_length_cm, data = rostrata_morphs)
plot(area_m1) #needs log adjustment
perimeter_m1 <- lm(otolith.perimeter ~ total_length_cm, data = rostrata_morphs)
plot(perimeter_m1) #needs log adjustment

#log-transform TL
rostrata_morphs$log_TL <- log(rostrata_morphs$total_length_cm)

#log-transform all morphs
for (morph in morphometrics) {
  rostrata_morphs[[paste0("log_", morph)]] <- log(rostrata_morphs[[morph]])
}

#otolith length m1
length_m1 <- lm(log_otolith.length ~ log_TL, data = rostrata_morphs)
plot(length_m1) #assumptions fine
summary(length_m1) #significant 

#otolith length m2 
length_m2 <- lm(log_otolith.length ~ log_TL + location, data = rostrata_morphs)
plot(length_m2) #assumptions fine
summary(length_m2) #significant

#otolith length m3
length_m3 <- lm(log_otolith.length ~ log_TL * location, data = rostrata_morphs)
plot(length_m3) #assumptions fine
summary(length_m3) #significant

#anova of otolith length m1, m2, m3
anova_length <- anova(length_m1, length_m2, length_m3)
anova_length
#each model has less residual error than the last (explains more variation)
#location impacts both intercept and slope 

#otolith width m1
width_m1 <- lm(log_otolith.width ~ log_TL, data = rostrata_morphs)
plot(width_m1) #assumptions fine
summary(width_m1) #significant

#otolith width m2 
width_m2 <- lm(log_otolith.width ~ log_TL + location, data = rostrata_morphs)
plot(width_m2) #assumptions fine
summary(width_m2) #significant

#otolith width m3
width_m3 <- lm(log_otolith.width ~ log_TL * location, data = rostrata_morphs)
plot(width_m3) #assumptions fine
summary(width_m3) #significant

#anova of otolith width m1, m2, m3
anova_width <- anova(width_m1, width_m2, width_m3)
anova_width #location only affects intercept, not slope

#otolith area m1
area_m1 <- lm(log_otolith.area ~ log_TL, data = rostrata_morphs)
plot(area_m1) #assumptions fine
summary(area_m1) #significant

#otolith area m2 
area_m2 <- lm(log_otolith.area ~ log_TL + location, data = rostrata_morphs)
plot(area_m2) #assumptions fine
summary(area_m2) #significant

#otolith area m3
area_m3 <- lm(log_otolith.area ~ log_TL * location, data = rostrata_morphs)
plot(area_m3) #assumptions fine
summary(area_m3) #significant

#anova of otolith area m1, m2, m3
anova_area <- anova(area_m1, area_m2, area_m3)
anova_area #location affects both intercept and slope

#otolith perimeter m1
perimeter_m1 <- lm(log_otolith.perimeter ~ log_TL, data = rostrata_morphs)
plot(perimeter_m1) #assumptions fine
summary(perimeter_m1) #significant

#otolith perimeter m2 
perimeter_m2 <- lm(log_otolith.perimeter ~ log_TL + location, data = rostrata_morphs)
plot(perimeter_m2) #assumptions fine
summary(perimeter_m2) #significant

#otolith perimeter m3
perimeter_m3 <- lm(log_otolith.perimeter ~ log_TL * location, data = rostrata_morphs)
plot(perimeter_m3) #assumptions fine
summary(perimeter_m3) #significant

#anova of otolith perimeter m1, m2, m3
anova_perimeter <- anova(perimeter_m1, perimeter_m2, perimeter_m3)
anova_perimeter #location affects both intercept and slope



#### 3) Supplementary plots####------------------------------------------------------

#4x figures, one figure for each morphometric 
#within each figure, 6x faceted plots - 2 diagnostic plots for each model
#3x residuals vs fitted, 3x QQ

#use patchwork package to stitch plots together
install.packages("patchwork")
library("patchwork")

#collect the models for each morphometric
model_list <- list(
  "Otolith length"    = list(length_m1,    length_m2,    length_m3),
  "Otolith width"     = list(width_m1,     width_m2,     width_m3),
  "Otolith area"      = list(area_m1,      area_m2,      area_m3),
  "Otolith perimeter" = list(perimeter_m1, perimeter_m2, perimeter_m3)
)

model_labels <- c("M1: log TL",
                  "M2: log TL + location",
                  "M3: log TL x location")

#create function called plot_diagnostics
#needs 2 inputs: a morph model and its corresponding label
#residuals vs fitted (top) and Q-Q (bottom) for one model
plot_diagnostics <- function(model, label) {
  
  #build a dataframe (d) with 3 columns: 
  #fitted (predicted values), residuals (observed - predicted), standardised residuals (for QQ)
  d <- data.frame(
    fitted    = fitted(model),
    resid     = resid(model),
    std_resid = rstandard(model)   
  )
  
  #residuals vs fitted plot
  p_rvf <- ggplot(d, aes(x = fitted, y = resid)) +
    geom_point(shape = 1, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + #grey horizontal line at 0
    geom_line(data = as.data.frame(lowess(d$fitted, d$resid, f = 2/3, iter = 3)),
              aes(x = x, y = y), colour = "red", linewidth = 0.6,
              inherit.aes = FALSE) +
    labs(title = label, x = "Fitted values", y = "Residuals") +
    theme_classic(base_size = 10)
  
  #QQ plot
  p_qq <- ggplot(d, aes(sample = std_resid)) + #QQ plot uses sample not x and y (ggplot works out theoretical quantiles itself)
    stat_qq(shape = 1, alpha = 0.6) +
    stat_qq_line(colour = "red", linewidth = 0.6) + #red reference line
    labs(x = "Theoretical quantiles", y = "Standardised residuals") +
    theme_classic(base_size = 10)
  
  #/ from patchwork package stacks the first plot above the second
  p_rvf / p_qq
}

#create a folder in the directory called supplementary-figures
dir.create("supplementary-figures", showWarnings = FALSE)

#for loop to repeat for each morphometric
#each time round it builds one figure and saves it as a png and pdf
for (morph_name in names(model_list)) {
  
  #make panels - list of three stacked plot pairs, one per model
  #Map() runs the plot_diagnostics function on each model and label 
  panels <- Map(plot_diagnostics, model_list[[morph_name]], model_labels)
  
  #wrap_plots() from patchwork arranges the list of 3 plots 
  figure <- wrap_plots(panels, nrow = 1) + #nrow = 1 puts them side by side in a single row (gives three columns)
    plot_annotation(title = morph_name) #adds title clarifying which morphometric
  
  #work out file name for the figure
  file_stub <- paste0("supplementary-figures/diagnostic-plots_",
                      gsub(" ", "_", tolower(morph_name)))#lower case, underscores not spaces
  
  #save as png and pdf
  ggsave(paste0(file_stub, ".png"), figure, width = 10, height = 6, dpi = 300)
  ggsave(paste0(file_stub, ".pdf"), figure, width = 10, height = 6)
}



#### 4) Morph - TL relationship plots ####----------------------------------------------

#line graphs faceted per morphometric

#define location order used for factor levels and legend - latitudinal not alphabetical
legend_order <- c("GL", "MAR", "RT", "FI", "SG")

#make location a factor in this order
rostrata_morphs$location <- factor(rostrata_morphs$location, levels = legend_order)

#set location colours
location_colours <- c(
  "GL"   = "#FFCC00",  #yellow
  "MAR"  = "#33CC00",  #green
  "RT"   = "#33CCFF",  #blue
  "FI"   = "#9933FF",  #purple
  "SG"   = "#FF33CC"   #pink
)

#reshape to long format - morphometric (name) and value (of morphometric) are 2 separate columns
morphs_long <- rostrata_morphs %>%
  pivot_longer(
    cols = all_of(c("otolith.length", "otolith.width",
                    "otolith.area", "otolith.perimeter")),
    names_to = "morphometric",
    values_to = "value"
  ) %>%
  
  #set order of morphometrics (logical not alphabetical)
  mutate(morphometric = factor(morphometric, levels = c(
    "otolith.length", "otolith.width", "otolith.area", "otolith.perimeter")))

#facet labels, single quotes inside double quotes to allow normal text
#^2* is not inside single quotes, so R will read it as code not text (parsed)
morphometric_labels <- c(
  "otolith.length"    = "'Otolith length (mm)'",
  "otolith.width"     = "'Otolith width (mm)'",
  "otolith.area"      = "'Otolith area (mm'^2*')'",
  "otolith.perimeter" = "'Otolith perimeter (mm)'"
)

#create plot figure_morph-TL using morphs_long object
figure_morph_TL <- ggplot(morphs_long, 
                          
  aes(x = total_length_cm, y = value, fill = location, colour = location)) +
  #plot total length on x axis, morphometric value on y axis
  
  geom_point(shape = 21, size = 1, stroke = 0.4, colour = "black", alpha = 0.4) +
  #shape 21 is a circle with different coloured fill and outline (outline black, fill location)
  #alpha 0.4 so points are 40% transparent and can overlap
  
  geom_smooth(aes(group = location), method = "lm", 
              formula = y ~ x, se = FALSE, colour = "black", linewidth = 1.6) +  
  #black line slightly thicker drawn underneath line filled by location colour
  #allows location colours to be seen more clearly (acts like a black outline)
  
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.1) +  
  #location-coloured line on top
  
  scale_fill_manual(values = location_colours, breaks = legend_order) +
  scale_colour_manual(values = location_colours, breaks = legend_order) +
  #applies chosen location colours, puts legend in correct order
  #does this for both points (fill) and lines (colour)
      
  facet_wrap(~morphometric, scales = "free_y", ncol = 2, strip.position = "left",
            labeller = as_labeller(morphometric_labels, default = label_parsed)) +
  #splits the plot into one panel per morphometric
  #scales = "free_y" - each panel gets its own scale for the y axis
  #ncol = 2 - 2x2 grid
  #labeller swaps the column names for the nicer labels
  #label_parsed converts text to expressions so superscript works for area
  
  labs(x = "Total length (cm)", y = NULL, colour = "Location", fill = "Location") +
  #y labels null because panel labels will act as y axis labels
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)),
        fill   = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  #make legend points opaque and larger than points on the plot
  
  theme_classic(base_size = 16) +
    
  theme(
    legend.position  = "right",
    axis.line        = element_line(linewidth = 0.5),
    axis.ticks       = element_line(linewidth = 0.5),
    panel.spacing    = unit(1, "lines"),
    legend.title     = element_text(face = "plain"),
    axis.title.x     = element_text(face = "plain"),
    strip.text       = element_text(face = "plain", angle = 0, size = 16),
    strip.background = element_blank(),
    strip.placement  = "outside"
  )

figure_morph_TL

#create a folder called morph-TL_figures in the current directory
dir.create("morph-TL_figures", showWarnings = FALSE)

#save the figure as a png and pdf in the folder
ggsave("figure_morph-TL_by_location.png", figure_morph_TL, path = "morph-TL_figures",
       width = 11, height = 8, dpi = 300)
ggsave("figure_morph-TL_by_location.pdf", figure_morph_TL, path = "morph-TL_figures",
       width = 11, height = 8)
