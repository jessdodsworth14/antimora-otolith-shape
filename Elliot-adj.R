#morphometrics in A-rostrata_morphs.csv have been tested in 'morph-TL-models.R' script
#location affects the slope and intercept of the morphometric-total length (TL) relationships
#therefore differences in shape among locations could be confounded by fish size 
#and location-specific scaling relationships
#morphometrics therefore need to be standardised for mean TL 
#done in this script using method of Elliot et al. 1995
#the output (standardised morphs) can be used to calculate shape indices

install.packages("readr")
install.packages("dplyr")
library("readr")
library("dplyr")

#Elliott et al. 1995 formula:
#   Ms = Mo * (Ls / Lo)^b

#where Ms is the standardised morphometric
#Mo is the original (observed) morphometric of the specimen
#Ls is the overall mean TL of all samples across locations (one reference length for all fish)
#Lo is the TL of that specimen
#b is the slope of the regression of log(Mo) on log(Lo)
#b is calculated separately per location for otolith length, area and perimeter
#because location affected the slope of those morph - TL relationships (see 'morph-TL-models.R')
#otolith width - TL relationship intercept was affected by location, not the slope
#so a common slope (b) is used across all locations for otolith width



#### 1) Prepare data ####-------------------------------------------------------------

#read in morphometric data, make location a factor
rostrata_morphs <- read_csv("A-rostrata_morphs.csv") %>%
  mutate(location = as.factor(location))

#store column names for ease
morphometrics <- c("otolith.area", "otolith.length", "otolith.width", "otolith.perimeter")
length_col    <- "total_length_cm"
group_col     <- "location"

#which morphometrics use a common within-location slope?
#TRUE  = location only affects the intercept (TL x location interaction not significant)
#FALSE = location affects the slope, so b is estimated separately per location
common_slope <- c(otolith.area      = FALSE,
                  otolith.length    = FALSE,
                  otolith.width     = TRUE,
                  otolith.perimeter = FALSE)

#### 2) Elliott adjustment (common or per-location b, global Ls) ####---------------------------

#create function with 5 inputs: dataframe, morphometric, TL, location, common slope (TRUE/FALSE)
#returns a list: adjusted values (adj) and the slope for each location (b)
adjust_morphometric <- function(data, morph, length_col, group_col, common = FALSE) {
  
  #Ls: mean TL of all fish across all locations (reference length)
  Ls <- mean(data[[length_col]], na.rm = TRUE)
  
  #dataframe with all the parts needed for the elliot equation
  d <- data.frame(Mo    = data[[morph]],
                  Lo    = data[[length_col]],
                  logM  = log(data[[morph]]),
                  logTL = log(data[[length_col]]),
                  grp   = data[[group_col]])
  
  if (common) {
    
    #common slope - one regression for all fish on log transformed data
    #each location gets its own intercept (+ grp) but all share one slope (b)
    fit <- lm(logM ~ logTL + grp, data = d)
    
    #return the intercept and slope in vector b, remove the 'logTL' label
    b   <- unname(coef(fit)["logTL"])
    
    #Elliott equation: Ms = Mo * (Ls / Lo)^b
    adj <- d$Mo * (Ls / d$Lo)^b
    
    #create output table for value of b
    #b is repeated for locations in this case
    b_by_loc <- setNames(rep(b, nlevels(d$grp)), levels(d$grp))
    
  } else {
    
    #per-location slopes (b) - common = FALSE
    #create empty results vector, one numeric NA per fish to start
    adj <- rep(NA_real_, nrow(d))
    
    #empty vector to collect the slope for each location
    b_by_loc <- numeric(0)
    
    #loop to fit a separate regression for each location
    for (g in levels(d$grp)) {
      
      rows <- which(d$grp == g)
      fit <- lm(logM ~ logTL, data = d[rows, ])
      b_by_loc[g] <- unname(coef(fit)["logTL"])
      
      #Elliott equation using this location's slope, global Ls
      adj[rows] <- d$Mo[rows] * (Ls / d$Lo[rows])^b_by_loc[g]
    }
  }
  
  #return both adjusted value and b
  list(adj = adj, b = b_by_loc)
}


#table to store the slope for each location and morphometric
b_table <- data.frame(location = levels(rostrata_morphs[[group_col]]))


#loop to apply adjust_morphometric function for each morph
for (m in morphometrics) {
  
  adjustment <- adjust_morphometric(rostrata_morphs, m, length_col, group_col,
                                    common = common_slope[[m]])
  
  #create new column on rostrata morphs called eg. otolith.length_adj
  rostrata_morphs[[paste0(m, "_adj")]] <- adjustment$adj
  
  #store the b slopes for each morphometric (one row per location)
  b_table[[m]] <- unname(adjustment$b[b_table$location])
}


#b slopes rounded to 3 decimal places
print(b_table %>% mutate(across(where(is.numeric), ~ round(.x, 3))))



#### 3) Save ####---------------------------------------------------------------------

#create new dataframe rostrata_morphs_adj where original morphometrics are replaced with adjusted
#adjusted have already been added during the function, select(-) drops original morphs
rostrata_morphs_adj <- rostrata_morphs %>% select(-all_of(morphometrics))

#save in a csv
write_csv(rostrata_morphs_adj, "A-rostrata_adjusted_morphs.csv")

#save slopes too
write_csv(b_table, "A-rostrata_b_slopes.csv")



#### 4) Checks ####-------------------------------------------------------------------

#is the size effect removed?
#loop to apply spearman's correlation for each adjusted morph and TL

adj_cols <- paste0(morphometrics, "_adj")

for (a in adj_cols) {
  
  #prints a heading for output table (\n starts a new line)
  cat("\n", a, "- within location:\n")
  
  #treat each location as its own mini dataset
  rostrata_morphs %>%
    group_by(location) %>%
    
    #summarise each group into one row of results
    #n() is number of fish in that location
    summarise(n   = n(),
              
              #correlation test using spearman's rank
              #.data[[length_col]] is TL column, .data[[a]] is adjusted morph 
              #extract rho (spearman's rank correlation)
              rho = cor(.data[[length_col]], .data[[a]], method = "spearman", use = "complete.obs"),
              
              #extract p value
              p   = cor.test(.data[[length_col]], .data[[a]], method = "spearman", exact = FALSE)$p.value) %>%
    print()
}

#adjustment has worked, rho values are all close to 0, p values are all >0.05
#now can use adjusted morphs to calculate shape indices
