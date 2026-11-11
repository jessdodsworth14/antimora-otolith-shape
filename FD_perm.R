
#Testing Antimora rostrata otolith Fourier Descriptors (FDs) for differences among 5 sampling locations
#using permanova to test differences among group FDs
#input is 'A-rostrata_FDs.csv' created during shape analysis (see 'Final_shape-analysis.R')

install.packages("readr")
install.packages("vegan")
library("readr")
library("vegan")



#### 1) Prepare data #### ---------------------------------------------------------

#read in metadata
rostrata_FDs <- read_csv("A-rostrata_FDs.csv")

#define FD columns
FD_columns <- paste0("V", 1:27)

#define location as a factor
rostrata_FDs$location <- factor(rostrata_FDs$location)

#PERMANOVA needs data to be split into measurements and group labels
#create matrix of FDs, drop other metadata
FD_matrix <- as.matrix(rostrata_FDs[, FD_columns])
#create grouping variable
location  <- rostrata_FDs$location

#check location has been stored correctly and has correct number of fish per location
table(location) 



#### 2) PERMANOVA on FDs #### ------------------------------------------------------

#distance matrix (Euclidean distance) on FDs
#this contains all the distances between every pair of otoliths
FD_dist <- dist(FD_matrix, method = "euclidean")

#for reproducibility
set.seed(123)

#run permanova using adonis2
#are the distances (differences in shape) explained by location
perm_result <- adonis2(FD_dist ~ location, 
                       data = data.frame(location = location), #give grouping variable
                       permutations = 9999,
                       method = "euclidean")

print(perm_result)

#test dispersion 
disp <- betadisper(FD_dist, location)
anova(disp)
permutest(disp, permutations = 9999)
plot(disp)
disp$group.distances



#### 3) Pairwise PERMANOVA between locations #### ----------------------------------

#create pairwise permanova function
#needs inputs: FD matrix, location groups, number of permutations, p value adjustment method
pairwise_permanova <- function(mat, groups, nperm = 9999, adjust = "holm") {
  
  #get location names
  lv    <- levels(groups)
  
  #list every possible pair of locations
  pairs <- combn(lv, 2, simplify = FALSE)
  
  #loop over every pair and stack results into a table using rbind
  res <- do.call(rbind, lapply(pairs, function(p) {
    
    #mark which otoliths belong to either location in this pair
    keep <- groups %in% p
    
    #keep only those otoliths and calculate distances between them
    d    <- dist(mat[keep, , drop = FALSE], method = "euclidean")
    
    #keep only the labels for those otoliths, drop unused locations
    g    <- droplevels(groups[keep])
    
    #run permanova on just these 2 locations
    a    <- adonis2(d ~ g, permutations = nperm)
    
    #get F value, R2 and p value for this pair, store in one row of table with the 2 location names
    data.frame(group1 = p[1], group2 = p[2],
               F = a$F[1], R2 = a$R2[1], p = a$`Pr(>F)`[1])
  }))
  
  #once all pairs are done, adjust p values for running 10 tests at once
  res$p_adj <- p.adjust(res$p, method = adjust)
  
  #return final table
  res
}

set.seed(123)

#run function on the data
pairwise_results <- pairwise_permanova(FD_matrix, location)
print(pairwise_results)

#run pairwise dispersion test (same seed approach for reproducibility)
set.seed(123)
pw_disp <- permutest(disp, pairwise = TRUE, permutations = 9999)

#extract permuted p values for each pair
disp_p <- pw_disp$pairwise$permuted

#build a table and apply the same Holm adjustment across the 10 pairs
pairwise_disp_results <- data.frame(
  pair  = names(disp_p),
  p     = disp_p,
  p_adj = p.adjust(disp_p, method = "holm"),
  row.names = NULL
)
print(pairwise_disp_results)
