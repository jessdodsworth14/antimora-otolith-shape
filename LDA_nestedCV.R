
#linear discriminant analysis (LDA) with nested cross validation (CV) on otolith fourier descriptors (FDs)
#testing discriminatory power to distinguish between locations using FDs
#input is 'A-rostrata_FDs.csv' (output of 'Final_shape-analysis.R')

#nested cross validation has 2 loops:

#OUTER loop (10-fold, stratified by location) - estimates predictive accuracy
#one fold of otoliths set aside as unseen test data
#these otoliths do not affect variable selection or model fitting

#     |
#     v

#INNER loop (10-fold, run inside klaR stepclass) - variable selection using LDA and CV
#forward-backward stepwise selection on training folds, scoring each FD subset with 10-fold CV

#     |
#     v

#back in OUTER loop: predict held-out fold with LDA model
#fit LDA on training folds using selected FDs
#the selected FD subset may differ between folds - this is expected
#(aim is to evaluate method, not one particular subset)
#repeat so every fold has a turn as the test fold, then pool all predictions
#accuracy = total correct predictions / total predictions 
#(the subsets are NOT combined into a single model)

install.packages("readr")
install.packages("MASS")
install.packages("klaR")
install.packages("dplyr")
install.packages("ggplot2")
library("readr")
library("MASS")
library("klaR")
library("dplyr")
library("ggplot2")



#### 1) Settings for nested CV ------------------------------------------------------------

seed    <- 123
k_outer <- 10 #number of outer folds
k_inner <- 10 #number of inner folds 



#### 2) Prepare data --------------------------------------------------------------------

#read in metadata
rostrata_FDs <- read_csv("A-rostrata_FDs.csv")

#define FD columns
FD_columns <- paste0("V", 1:27)

#create LDA dataframe with FDs and location 
LDA_df <- as.data.frame(rostrata_FDs[, c(FD_columns, "location")])

#define location as a factor
LDA_df$location <- factor(LDA_df$location)

#check sample size per location
table(LDA_df$location) 

#number of otoliths
n <- nrow(LDA_df)



#### 3) Assign stratified outer folds -----------------------------------------------

set.seed(seed)

#each location is split as evenly as possible across the folds,
#so every fold contains some otoliths from every location (incl. the small RT group)

#create a vector of n (569) zeros to store fold labels
fold_id <- integer(n)

#loop for every location (loc)
for (loc in levels(LDA_df$location)) {
  
  #find all the otoliths in the current location and their row numbers
  #sample() shuffles the rows into a random order
  idx <- sample(which(LDA_df$location == loc))    
  
  #shuffle folds 1 to 10 (seq_len(k_outer) is 10 for 10 outer folds)
  #repeat that random sequence and assign to the otoliths in the current location
  #until all the otoliths in the current location have a fold number
  fold_id[idx] <- rep_len(sample(seq_len(k_outer)), length(idx))  
}

#check: each location should be spread across all folds
table(LDA_df$location, fold_id)   



#### 4) Nested CV loop ----------------------------------------------------------------

#object to store the predicted location for each otolith
pred_nested <- factor(rep(NA, n), levels = levels(LDA_df$location))

#per-fold records: number of FDs selected, which FDs, and whether the fallback was used
n_sel         <- integer(k_outer)
sel_list      <- vector("list", k_outer)
used_fallback <- logical(k_outer)

#record clock time so message can report time elapsed
start_time <- Sys.time()

#repeat 10 times so every otolith has a go in the test fold, f is current fold
for (f in seq_len(k_outer)) { 
  
  #OUTER LOOP
  
  #find row numbers of otoliths assigned to current fold f 
  test_idx <- which(fold_id == f)
  
  #training data is everything except fold f
  train <- LDA_df[-test_idx, ]                
  
  #test data is fold f only
  test  <- LDA_df[test_idx, , drop = FALSE] #drop = FALSE keeps as dataframe  
  
    #INNER LOOP
    #variable selection using only training data
    #(stepclass runs its own k_inner-fold CV = the inner loop)
    #tryCatch() means if it crashes, store NULL and carry on
    sel <- tryCatch(
    
    #run stepwise selection
    stepclass(
      x = train[, FD_columns],
      grouping = train$location,
      method = "lda",
      improvement = 0.01, #stop when no change improves accuracy by at least 1%
      direction = "both", #forward-backward selection
      fold = k_inner,
      output = FALSE #silence per-step printing
    ),
    error = function(e) {   #if stepclass fails, say why and fall back below instead of crashing
      message("stepclass failed in fold ", f, ": ", conditionMessage(e))
      NULL
    }
  )
  
  #REST OF OUTER LOOP
      
  #character vector of selected FD names (chosen by inner loop)
  vars <- if (is.null(sel)) character(0) else sel$model$name
  
  #fallback: use all FDs if nothing was selected
  used_fallback[f] <- length(vars) == 0
  if (length(vars) == 0) vars <- FD_columns
  
  #store the names and the count of selected FDs for fold f
  sel_list[[f]] <- vars
  n_sel[f] <- length(vars)
  
  #fit LDA on training data with the selected FDs, predict the held-out fold
  fit <- MASS::lda(as.formula(paste("location ~", paste(vars, collapse = "+"))), data = train)
  pred_nested[test_idx] <- predict(fit, newdata = test)$class
  
  #message saying when each fold is complete
  message("fold ", f, "/", k_outer, " done; selected ", length(vars),
          " FDs; elapsed: ", round(difftime(Sys.time(), start_time, units = "mins"), 1), " min")
}



#### 5) Results -------------------------------------------------------------------

#confusion matrix (rows = observed, columns = predicted)
conf_nested <- table(Observed = LDA_df$location, Predicted = pred_nested)
conf_nested

#overall accuracy, pooled across all folds
correct <- pred_nested == LDA_df$location
overall_acc <- mean(correct)
overall_acc

#95% confidence interval for overall accuracy
binom.test(sum(correct), n)$conf.int

#accuracy within each outer fold (shows how much it varies between folds)
round(tapply(correct, fold_id, mean), 3)

#per-location recall (proportion of each location's otoliths correctly assigned)
recall <- diag(conf_nested) / rowSums(conf_nested)
round(recall, 3)

#balanced accuracy = mean of per-location recall (not inflated by large groups)
balanced_acc <- mean(recall)
balanced_acc

#what you'd get by chance
baseline_majority <- max(table(LDA_df$location)) / n   #always guessing the biggest location
baseline_balanced <- 1 / nlevels(LDA_df$location)      #expected balanced accuracy by guessing
baseline_majority
baseline_balanced

#is overall accuracy better than always guessing the biggest location?
binom.test(sum(correct), n, p = baseline_majority, alternative = "greater")



#### 6) Plots ------------------------------------------------------------------------

#set colour palette
location_colours <- c(
  "GL"   = "#FFCC00",  #yellow
  "MAR"  = "#33CC00",  #green
  "RT"   = "#33CCFF",  #blue
  "FI"   = "#9933FF",  #purple
  "SG"   = "#FF33CC"   #pink
)

#fit final LDA model on all otoliths
#variable selection on the full dataset (same settings as the nested CV)
set.seed(seed)
sel_final <- stepclass(
  x = LDA_df[, FD_columns],
  grouping = LDA_df$location,
  method = "lda",
  improvement = 0.01,
  direction = "both",
  fold = k_inner,
  output = FALSE
)
vars_final <- sel_final$model$name

#a 2D plot needs at least 2 FDs; stop if less than 2 are selected
stopifnot(length(vars_final) >= 2)

#the FDs used in the final model
vars_final

#fit the final LDA
lda_final_model <- MASS::lda(as.formula(paste("location ~", paste(vars_final, collapse = "+"))),
                             data = LDA_df)

#get LDA scores for each otolith
lda_proj <- predict(lda_final_model)$x

#order locations the same as the palette (controls legend and boxplot order)
lda_plot_df <- data.frame(LD1 = lda_proj[, 1],
                          LD2 = lda_proj[, 2],
                          location = factor(LDA_df$location, levels = names(location_colours)))

#proportion of between-group variance explained by each LD
prop_var <- lda_final_model$svd^2 / sum(lda_final_model$svd^2)
x_lab <- paste0("LD1 (", round(100 * prop_var[1], 2), "%)")
y_lab <- paste0("LD2 (", round(100 * prop_var[2], 2), "%)")

#plot otoliths + shaded 95% ellipses
LD_plot <- ggplot(lda_plot_df, aes(x = LD1, y = LD2, color = location, fill = location)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95,
               alpha = 0.1, linewidth = 1) +
  geom_point(shape = 21, size = 1.2, alpha = 0.8, stroke = 0.5, color = "black") +
  scale_color_manual(values = location_colours) +
  scale_fill_manual(values = location_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = x_lab,
       y = y_lab,
       color = "Location",
       fill = "Location")
LD_plot

#LD1 boxplot
LD1_boxplot <- ggplot(lda_plot_df, aes(x = location, y = LD1, fill = location)) +
  geom_boxplot(alpha = 0.8, notch = TRUE) +                          
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = location_colours) +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Location",
       y = "LD1")
LD1_boxplot

#combine LD1-LD2 plot and LD1 boxplot into one panelled figure
#A = LD1-LD2 plot, B = LD1 boxplot

#patchwork stitches plots together
install.packages("patchwork")
library("patchwork")

#panel sizes in inches
fig_width  <- 7
A_height   <- 5
B_height   <- 4

# / stacks A above B, heights set the relative panel sizes
figure_LDA <- LD_plot / LD1_boxplot +
  plot_layout(heights = c(A_height, B_height)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

figure_LDA

#save as png and pdf
ggsave("LDA_panel_plot.png", figure_LDA, width = fig_width, height = A_height + B_height, dpi = 300)
ggsave("LDA_panel_plot.pdf", figure_LDA, width = fig_width, height = A_height + B_height)

