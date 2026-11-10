
#Comparison of Antimora rostrata otolith shape indices among 5 sampling locations in the Atlantic
#shape indices have been tested for correlation with fish total length (TL), no relationship
#(previously size standardised) - no confounding effect of fish size
#input: 'A-rostrata_shape_indices.csv' (calculated in 'A-rostrata_shape_indices.R')

install.packages("readr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")



#### 1) Prepare data ####-------------------------------------------------------------

#colours for each location (also sets location order to latitudinal not alphabetical)
location_colours <- c(
  "GL"   = "#FFCC00",  #yellow
  "MAR"  = "#33CC00",  #green
  "RT"   = "#33CCFF",  #blue
  "FI"   = "#9933FF",  #purple
  "SG"   = "#FF33CC"   #pink
)

#read in data, set location as factor with correct colours and order for plots
rostrata_indices <- read_csv("A-rostrata_shape_indices.csv") %>%
  mutate(location = factor(location, levels = names(location_colours)))

#name indices for ease
index_cols <- c("aspect_ratio", "circularity", "ellipticity",
                "form_factor", "rectangularity", "roundness")

#plot labels
index_labels <- c(aspect_ratio   = "Aspect ratio",
                  circularity    = "Circularity",
                  ellipticity    = "Ellipticity",
                  form_factor    = "Form factor",
                  rectangularity = "Rectangularity",
                  roundness      = "Roundness")



#### 2) Testing indices among locations ####------------------------------------

#make an anova for each index, check assumptions for all
anova_AR <- aov(aspect_ratio ~ location, data = rostrata_indices)
plot(anova_AR) #assumptions fine
summary(anova_AR) #location significant
TukeyHSD(anova_AR) #significant:
                   #FI-GL
                   #FI-MAR

anova_C <- aov(circularity ~ location, data = rostrata_indices)
plot(anova_C) #assumptions fine
summary(anova_C) #location significant
TukeyHSD(anova_C) #significant:
                  #GL-MAR
                  #FI-MAR

anova_E <- aov(ellipticity ~ location, data = rostrata_indices)
plot(anova_E) #assumptions fine
summary(anova_E) #location significant
TukeyHSD(anova_E) #significant:
                  #FI-GL
                  #FI-MAR

anova_FF <- aov(form_factor ~ location, data = rostrata_indices)
plot(anova_FF) #assumptions fine
summary(anova_FF) #location significant
TukeyHSD(anova_FF) #significant:
                   #MAR-GL
                   #FI-MAR

anova_R <- aov(rectangularity ~ location, data = rostrata_indices)
plot(anova_R) #assumptions fine
summary(anova_R) #location significant
TukeyHSD(anova_R) #significant:
                  #FI-GL
                  #SG-GL
                  #FI-MAR
                  #SG-MAR

anova_Ro <- aov(roundness ~ location, data = rostrata_indices)
plot(anova_Ro) #assumptions fine
summary(anova_Ro) #location not significant



#### 3) Faceted boxplots ####--------------------------------------------------------

#long format
locations_long <- rostrata_indices %>%
  select(location, all_of(index_cols)) %>%
  pivot_longer(cols = -location, names_to = "Shape_Index", values_to = "value") %>%
  mutate(Shape_Index = factor(Shape_Index, levels = index_cols))

#create dataframe to turn ANOVA p-value for each index into significance stars for the plot
anova_stars <- data.frame(
  Shape_Index = factor(index_cols, levels = index_cols),
  p = sapply(index_cols, function(idx) {
    summary(aov(reformulate("location", response = idx),
                data = rostrata_indices))[[1]][["Pr(>F)"]][1]
  })
) %>%
  mutate(label = case_when(p < 0.001 ~ "***",
                           p < 0.01  ~ "**",
                           p < 0.05  ~ "*",
                           TRUE      ~ ""))

#create faceted boxplot
index_plot <- ggplot(locations_long, aes(x = location, y = value, fill = location)) +
  geom_boxplot(notch = TRUE) +
  geom_text(data = anova_stars,
            aes(x = 3, y = Inf, label = label),   #x = 3 is the middle of the 5 locations
            inherit.aes = FALSE, vjust = 1.2, size = 6) +
  facet_wrap(~ Shape_Index, scales = "free_y", ncol = 3,
             labeller = labeller(Shape_Index = index_labels)) +
  scale_fill_manual(values = location_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  #headroom for the stars
  theme_classic() +
  labs(x = "Location", y = "Value", fill = "Location") +
  theme(legend.position = "right",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

print(index_plot)

#save as pdf and png
ggsave("shape_indices_boxplots.pdf", index_plot, width = 9, height = 6.5, dpi = 300)
ggsave("shape_indices_boxplots.png", index_plot, width = 9, height = 6.5, dpi = 300)



#### 4) Combining faceted boxplot with mean outlines plot ####-----------------------------------

#'Mean_otolith_outlines.pdf' was the output of shape analysis
#here combine the outlines plot with the shape indices plot for a faceted overall otolith shape plot
#can visually compare shape data from both shape indices and FDs (outlines)

install.packages("pdftools")
install.packages("magick")
library("pdftools")
library("magick")

#read the saved outlines PDF as a high-resolution image, trim the white margins
outline_img <- image_read_pdf("Mean_otolith_outlines.pdf", density = 600) %>%
  image_trim()

#height/width ratio of the trimmed image, so the panel keeps the right proportions
outline_info  <- image_info(outline_img)
outline_ratio <- outline_info$height / outline_info$width

#turn the image into a panel patchwork can handle
outline_panel <- wrap_elements(
  full = grid::rasterGrob(as.raster(outline_img), interpolate = TRUE)
)

#remove the legend from the boxplots (the outlines panel already has the location key)
index_plot_nl <- index_plot + theme(legend.position = "none")

#panel sizes in inches
fig_width      <- 8
box_height     <- 6
outline_height <- fig_width * outline_ratio

#A = boxplots, B = mean outlines
figure_shape <- index_plot_nl / outline_panel +
  plot_layout(heights = c(box_height, outline_height)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

#create folder for the plots
dir.create("indices-outlines_figures", showWarnings = FALSE)

#save as png and pdf
ggsave("indices-outlines_figures/shape_indices_and_outlines.png", figure_shape,
       width = fig_width, height = box_height + outline_height, dpi = 300)
ggsave("indices-outlines_figures/shape_indices_and_outlines.pdf", figure_shape,
       width = fig_width, height = box_height + outline_height)



#### 5) Supplementary ANOVA diagnostic plots ####------------------------------------------------

#one figure: 6 columns (one per shape index), 2 rows (residuals vs fitted on top, QQ below)

#use patchwork package to stitch plots together
install.packages("patchwork")
library("patchwork")

#create function called plot_diagnostics
#needs 2 inputs: an aov model and its corresponding label
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
  #in a one-way ANOVA the fitted values are just the 5 location means, so points fall in
  #5 vertical columns - the dashed line at 0 and the spread of each column is what to look at 
  p_rvf <- ggplot(d, aes(x = fitted, y = resid)) +
    geom_point(shape = 1, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") + #grey horizontal line at 0
    labs(title = label, x = "Fitted values", y = "Residuals") +
    theme_classic(base_size = 10)
  
  #QQ plot
  p_qq <- ggplot(d, aes(sample = std_resid)) + #QQ plot uses sample not x and y
    stat_qq(shape = 1, alpha = 0.6) +
    stat_qq_line(colour = "red", linewidth = 0.6) + #red reference line
    labs(x = "Theoretical quantiles", y = "Standardised residuals") +
    theme_classic(base_size = 10)
  
  #/ from patchwork package stacks the first plot above the second
  p_rvf / p_qq
}

#fit one ANOVA per index and store in a named list
#reformulate() builds the formula "index ~ location" from the column name
anova_list <- lapply(index_cols, function(idx) {
  aov(reformulate("location", response = idx), data = rostrata_indices)
})
names(anova_list) <- index_cols

#make panels - list of six stacked plot pairs, one per index
#Map() runs plot_diagnostics on each model and its label
panels <- Map(plot_diagnostics, anova_list, index_labels[index_cols])

#wrap_plots() arranges the six stacked pairs side by side in a single row
figure <- wrap_plots(panels, ncol = 3) +
  plot_annotation(title = "ANOVA diagnostic plots: otolith shape indices by location")

plot(figure)

#save as png and pdf in supplementary figures folder created during morph-TL-models.R script
ggsave("supplementary-figures/diagnostic-plots_shape_indices.png",
       figure, width = 8, height = 10, dpi = 300)
ggsave("supplementary-figures/diagnostic-plots_shape_indices.pdf",
       figure, width = 8, height = 10)

