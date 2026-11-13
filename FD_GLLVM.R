
#Generalised Latent Linear Variable Models (GLLVMs) on Antimora rostrata otolith Fourier Descriptors (FDs)
#FDs are the output of 'Final_shape-analysis.R' and are standardised for fish total length
#will test location's ability to explain variation in FD data as a predictor, 
#and investigate if adding depth and sex improves model fit
#and/or explains more variation in FDs than location alone
#sampling station will be a random effect in all models
#will try models with 1-3 latent variables and choose the one with lowest AIC to fit the final models

install.packages("readr")
install.packages("gllvm")
install.packages("tidyr")
install.packages("dplyr")
library("readr")
library("gllvm")
library("tidyr")
library("dplyr")


#### 1) Prepare the data ------------------------------------------------------------

rostrata_FDs <- read_csv("A-rostrata_FDs.csv")

#define FDs
FDs <- as.matrix(rostrata_FDs[, paste0("V", 1:27)])

#create new station ID column - MAR and GL share some station numbers
rostrata_FDs$station_id <- paste(rostrata_FDs$location, rostrata_FDs$station, sep = "_")

#create scaled depth column - GLLVM handles continuous data better this way
#scale() standardises to a mean of 0 and standard deviation of 1
rostrata_FDs$depth.s <- as.numeric(scale(rostrata_FDs$depth))

#scale FDs
FDs.s <- scale(FDs)

#set N to total sample size
N <- nrow(rostrata_FDs)

#set N.FDs to number of FD columns
N.FDs <- ncol(FDs)



#### 2) GLLVM with location and depth predictors -----------------------------------------------

#define covariates
CovData <- as.data.frame(rostrata_FDs[, c("location", "depth.s")])

#make sure FDs and CovData have matching row names
rownames(FDs.s) <- paste("r", 1:N, sep = "")
rownames(CovData) <- rownames(FDs.s)

#define studyDesign argument for gllvm function - used for random effects 
#random effect is sampling station
studyDesignData <- data.frame(stationEffect = factor(rostrata_FDs$station_id))

#give it same row numbers
rownames(studyDesignData) <- rownames(FDs.s)

#fit GLLVM for location and depth covariates, station random effect
#start with 1 latent variable then 2 and 3, compare with AIC to see best fit

#1 latent variable
locdepth_gllvm1 <- gllvm(y = FDs.s, X = CovData,
                         studyDesign = studyDesignData,
                         row.eff = ~ (1 | stationEffect),
                         formula = ~ location + depth.s, #location and depth
                         num.lv = 1, #1 latent variable
                         control.start = list(n.init = 5, jitter.var = 0.1),
                         family = "gaussian")

#2 latent variables
locdepth_gllvm2 <- gllvm(y = FDs.s, X = CovData,
                         studyDesign = studyDesignData,
                         row.eff = ~ (1 | stationEffect),
                         formula = ~ location + depth.s,
                         num.lv = 2, #2 latent variables
                         control.start = list(n.init = 5, jitter.var = 0.1),
                         family = "gaussian")

#3 latent variables
locdepth_gllvm3 <- gllvm(y = FDs.s, X = CovData,
                      studyDesign = studyDesignData,
                      row.eff = ~ (1 | stationEffect),
                      formula = ~ location + depth.s,
                      num.lv = 3, #3 latent variables
                      control.start = list(n.init = 5, jitter.var = 0.1),
                      family = "gaussian")

#compare AIC values among the 3 models
AIC(locdepth_gllvm1, locdepth_gllvm2, locdepth_gllvm3)
#3 LV gives best model fit - locdepth_gllvm3 is final model for location + depth predictors

#check assumptions of location + depth model at 3LV
plot(locdepth_gllvm3)

#test significance
coefplot(locdepth_gllvm3)

#see how much variation in FDs is explained by covariates

#repeat model at 3 LV with no covariates for comparison
NoCov1 <- gllvm(y = FDs.s, #no X = CovData
                studyDesign = studyDesignData,
                row.eff = ~ (1 | stationEffect),
                formula = ~ 1, #instead of covariates
                num.lv = 3,
                family = "gaussian")

#get trace of residual covariance matrix for locdepth_gllvm3 and NoCov1
rcov_locdepth_gllvm3 <- getResidualCov(locdepth_gllvm3, adjust = 0)
rcov_NoCov1 <- getResidualCov(NoCov1, adjust = 0)

#how much of variation is explained by location + depth
1 - rcov_locdepth_gllvm3$trace / rcov_NoCov1$trace
#2.7 %



#### 3) GLLVM with location as only predictor -------------------------------------------

#fit new set of models with location as only predictor (remove depth)
#compare AIC among LV1, 2 and 3 models
#compare AIC with location + depth model at same LV number

#1 LV
loc_gllvm1 <- gllvm(y = FDs.s, X = CovData,
                    studyDesign = studyDesignData,
                    row.eff = ~ (1 | stationEffect),
                    formula = ~ location, #location only
                    num.lv = 1,
                    control.start = list(n.init = 5, jitter.var = 0.1),
                    family = "gaussian")

#compare AIC with location + depth model at same LV number
AIC(loc_gllvm1, locdepth_gllvm1)
#location only gives slightly better model fit

#2 LV
loc_gllvm2 <- gllvm(y = FDs.s, X = CovData,
                    studyDesign = studyDesignData,
                    row.eff = ~ (1 | stationEffect),
                    formula = ~ location, 
                    num.lv = 2,
                    control.start = list(n.init = 5, jitter.var = 0.1),
                    family = "gaussian")

#compare AIC with location + depth model at same LV number
AIC(loc_gllvm2, locdepth_gllvm2)
#location only gives slightly better model fit

#3 LV
loc_gllvm3 <- gllvm(y = FDs.s, X = CovData,
                    studyDesign = studyDesignData,
                    row.eff = ~ (1 | stationEffect),
                    formula = ~ location, 
                    num.lv = 3,
                    control.start = list(n.init = 5, jitter.var = 0.1),
                    family = "gaussian")

#compare AIC values among the 3 models
AIC(loc_gllvm1, loc_gllvm2, loc_gllvm3)
#3 LV gives best model fit - loc_gllvm3 is final model for location as only predictor

#compare AIC with location + depth model at same LV number
AIC(loc_gllvm3, locdepth_gllvm3)
#location only gives slightly better model fit

#check assumptions of final location only model
plot(loc_gllvm3)

#test significance
coefplot(loc_gllvm3)

##see how much variation in FDs is explained by covariates

#get trace of residual covariance matrix for loc_gllvm3 
rcov_loc_gllvm3 <- getResidualCov(loc_gllvm3, adjust = 0)

#how much of variation is explained by location only
1 - rcov_loc_gllvm3$trace / rcov_NoCov1$trace
#2.5 %

#RESULT: location + depth as predictors do not improve model fit compared to location only
#however combining location + depth increases explained variation from 2.5 % to 2.7 %



#### 4) Prepare reduced dataset for models including sex -----------------------------------------

#see how many males and females were caught per location
table(rostrata_FDs$location, rostrata_FDs$sex)     
#only MAR has no males

#make smaller dataset corresponding to only M and F fish (no U), and no MAR
keep <- rostrata_FDs$location != "MAR" & rostrata_FDs$sex %in% c("F", "M")
rostrata_MandF <- rostrata_FDs[keep, ]

#new CovData for reduced dataset
CovData_reduced <- data.frame(location = rostrata_MandF$location, sex = rostrata_MandF$sex)

#new FDs for reduced dataset
FDs_reduced   <- as.matrix(rostrata_MandF[, paste0("V", 1:27)])
FDs_reduced.s <- scale(FDs_reduced)

#new studyDesignData for reduced dataset
studyDesignData_reduced <- data.frame(stationEffect = factor(rostrata_MandF$station_id))

#make sure new gllvm input datasets have matching row names
N <- nrow(rostrata_MandF)
rownames(FDs_reduced.s) <- paste("r", 1:N, sep = "")
rownames(CovData_reduced) <- rownames(FDs_reduced.s)
rownames(studyDesignData_reduced) <- rownames(FDs_reduced.s)



#### 5) Reduced daataset GLLVM with location and sex predictors -----------------------------

#LV1 
locsex_gllvm1 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                                  studyDesign = studyDesignData_reduced,
                                  row.eff = ~ (1 | stationEffect),
                                  formula = ~ location + sex, #location and sex
                                  num.lv = 1,
                                  control.start = list(n.init = 5, jitter.var = 0.1),
                                  family = "gaussian")

#LV2 
locsex_gllvm2 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                       studyDesign = studyDesignData_reduced,
                       row.eff = ~ (1 | stationEffect),
                       formula = ~ location + sex, 
                       num.lv = 2,
                       control.start = list(n.init = 5, jitter.var = 0.1),
                       family = "gaussian")

#LV3 
locsex_gllvm3 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                       studyDesign = studyDesignData_reduced,
                       row.eff = ~ (1 | stationEffect),
                       formula = ~ location + sex, 
                       num.lv = 3,
                       control.start = list(n.init = 5, jitter.var = 0.1),
                       family = "gaussian")

#compare AIC values among the 3 models
AIC(locsex_gllvm1, locsex_gllvm2, locsex_gllvm3)
#3 LV gives best model fit - locsex_gllvm3 is final model for location + sex as predictors

#check assumptions of location + sex model at 3LV
plot(locsex_gllvm3)

#test significance
coefplot(locsex_gllvm3)

#see how much variation in FDs is explained by covariates

#repeat model at 3 LV with no covariates for comparison (reduced dataset)
NoCov2 <- gllvm(y = FDs_reduced.s, #no X = CovData_reduced
                studyDesign = studyDesignData_reduced,
                row.eff = ~ (1 | stationEffect),
                formula = ~ 1, #instead of covariates
                num.lv = 3,
                family = "gaussian")

#get trace of residual covariance matrix for locsex_gllvm3 and NoCov2
rcov_locsex_gllvm3 <- getResidualCov(locsex_gllvm3, adjust = 0)
rcov_NoCov2 <- getResidualCov(NoCov2, adjust = 0)

#how much of variation is explained by location + sex
1 - rcov_locsex_gllvm3$trace / rcov_NoCov2$trace
#1.8 %



#### 6) Reduced dataset GLLVM with location as only predictor --------------------------------

#LV1 
locReduced_gllvm1 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                       studyDesign = studyDesignData_reduced,
                       row.eff = ~ (1 | stationEffect),
                       formula = ~ location, #location only
                       num.lv = 1,
                       control.start = list(n.init = 5, jitter.var = 0.1),
                       family = "gaussian")

#compare AIC with location + sex model at same LV number
AIC(locsex_gllvm1, locReduced_gllvm1)
#location + sex gives slightly better model fit

#LV2 
locReduced_gllvm2 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                           studyDesign = studyDesignData_reduced,
                           row.eff = ~ (1 | stationEffect),
                           formula = ~ location, 
                           num.lv = 2,
                           control.start = list(n.init = 5, jitter.var = 0.1),
                           family = "gaussian")

#compare AIC with location + sex model at same LV number
AIC(locsex_gllvm2, locReduced_gllvm2)
#location + sex gives slightly better model fit

#LV3 
locReduced_gllvm3 <- gllvm(y = FDs_reduced.s, X = CovData_reduced,
                           studyDesign = studyDesignData_reduced,
                           row.eff = ~ (1 | stationEffect),
                           formula = ~ location, 
                           num.lv = 3,
                           control.start = list(n.init = 5, jitter.var = 0.1),
                           family = "gaussian")

#compare AIC with location + sex model at same LV number
AIC(locsex_gllvm3, locReduced_gllvm3)
#location + sex gives slightly better model fit

#compare AIC values among the 3 models
AIC(locReduced_gllvm1, locReduced_gllvm2, locReduced_gllvm3)
#3 LV gives best model fit - locReduced_gllvm3 is final model for location as only predictor (reduced dataset)

#check assumptions of final location only (reduced dataset) model
plot(locReduced_gllvm3)

#test significance
coefplot(loc_gllvm3)

##see how much variation in FDs is explained by covariates

#get trace of residual covariance matrix for locReduced_gllvm3 
rcov_locReduced_gllvm3 <- getResidualCov(locReduced_gllvm3, adjust = 0)

#how much of variation is explained by location only
1 - rcov_locReduced_gllvm3$trace / rcov_NoCov2$trace
#1.1 %

#RESULT: location + sex as predictors improves model fit compared to location alone
#also explains 1.8 % of variation in FDs compared to 1.1 %



#### 7) Save diagnostic plots for supplementary -----------------------------------------------

#create new directory
dir.create("gllvm_diagnostic_plots", showWarnings = FALSE)

#full dataset location + depth predictors (LV 3)
pdf("gllvm_diagnostic_plots/diagnostics_locations-depth.pdf", width = 10, height = 8)
plot(locdepth_gllvm3)
dev.off()

#full dataset location only predictor (LV 3)
pdf("gllvm_diagnostic_plots/diagnostics_locations.pdf", width = 10, height = 8)
plot(loc_gllvm3)
dev.off()

#reduced dataset location + sex predictors (LV 3)
pdf("gllvm_diagnostic_plots/diagnostics_locations-sex.pdf", width = 10, height = 8)
plot(locsex_gllvm3)
dev.off()

#reduced dataset location only predictor (LV 3)
pdf("gllvm_diagnostic_plots/diagnostics_locations-reduced.pdf", width = 10, height = 8)
plot(locReduced_gllvm3)
dev.off()

