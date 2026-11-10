#script to calculate shape indices from TL-standardised otolith morphometrics (see 'Elliot-adj.R')
#formulae from Moore et al. 2022
#input is 'A-rostrata_adjusted_morphs.csv' from 'Elliot-adj.R'

#   Aspect ratio (AR)      = OL / OW                overall elongation
#   Circularity (C)        = OP^2 / OA              closeness to a perfect circle
#   Ellipticity (E)        = (OL - OW) / (OL + OW)  closeness to an ellipse
#   Form factor (FF)       = (4 * pi * OA) / OP^2   surface area irregularity
#   Rectangularity (Re)    = OA / (OL * OW)         closeness to a perfect rectangle
#   Roundness (Ro)         = (4 * OA) / (pi * OL^2) closeness to a perfect circle

#OL = otolith length, OW = otolith width, OP = otolith perimeter, OA = otolith area
#(all are the TL-adjusted versions, i.e. the *_adj columns)

install.packages("readr")
install.packages("dplyr")
library(readr)
library(dplyr)



#### 1) Prepare data ####-------------------------------------------------------------

#read in the adjusted morphometrics, make location a factor
rostrata_adj <- read_csv("A-rostrata_adjusted_morphs.csv") %>%
  mutate(location = as.factor(location))

#check the adjusted columns are present
stopifnot(all(c("otolith.area_adj", "otolith.length_adj",
                "otolith.width_adj", "otolith.perimeter_adj") %in% names(rostrata_adj)))

#name length column for ease
length_col <- "total_length_cm"



#### 2) Calculate shape indices ####--------------------------------------------------

#dplyr pipeline to add new shape indices columns to a dataframe called rostrata_indices
rostrata_indices <- rostrata_adj %>%
  
  #add temporary columns with names matching the formulae components
  mutate(OA = otolith.area_adj,
         OL = otolith.length_adj,
         OW = otolith.width_adj,
         OP = otolith.perimeter_adj,
         
         #formulae
         #each line creates a new column that includes all the fish
         aspect_ratio     = OL / OW,
         circularity      = OP^2 / OA,
         ellipticity      = (OL - OW) / (OL + OW),
         form_factor      = (4 * pi * OA) / OP^2,
         rectangularity   = OA / (OL * OW),
         roundness        = (4 * OA) / (pi * OL^2)) %>%
  
  #drop the temporary short-name columns
  select(-OA, -OL, -OW, -OP)

#names of the indices for later use
index_cols <- c("aspect_ratio", "circularity", "ellipticity",
                "form_factor", "rectangularity", "roundness")

#quick look at the values (check for NA, Inf or odd ranges)
summary(rostrata_indices[, index_cols])



#### 3) Save ####---------------------------------------------------------------------

write_csv(rostrata_indices, "A-rostrata_shape_indices.csv")



#### 4) Checks: are the indices free of TL effects? ####------------------------------

#within-location tests: for each index (6) in every location (5), 
#calculate a Spearman correlation between TL and the index (30 correlations)
within_tests <- bind_rows(lapply(index_cols, function(i) {
  rostrata_indices %>%
    group_by(location) %>%
    summarise(index = i,
              rho   = cor(.data[[length_col]], .data[[i]], method = "spearman", use = "complete.obs"),
              p     = cor.test(.data[[length_col]], .data[[i]], method = "spearman", exact = FALSE)$p.value,
              .groups = "drop")
})) %>%
  
  #holm adjustment for multiple testing
  mutate(p_holm = p.adjust(p, method = "holm"))

print(within_tests, n = 30)

#result: spearman rank correlations between TL and each shape index, 
#calculated separately for each location,
#showed no relationship after Holm adjustment for multiple testing (all adjusted p > 0.7)
#effect of TL successfully removed, shape indices ready for among-location difference testing
