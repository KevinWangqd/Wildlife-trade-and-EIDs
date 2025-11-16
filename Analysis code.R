require(INLA)
library(car)
#library(dplyr)
# Load the Matrix package
library(Matrix)
# Functions #
slcpo <- function(m) {
  -sum(log(m$cpo$cpo),na.rm=TRUE)
}


RMSE <- function(set, outcome, data, fit) {
  res <- data[set, outcome] - fit[set,1]
  RMSE_val <- sqrt(mean(res^2, na.rm = TRUE))
  return(RMSE_val)
}


# Data #
data_total <- read.csv("data_total.csv")
data_zoo <- read.csv("data_zoo.csv")
data_nonzoo <- read.csv("data_nonzoo.csv")


# adjusted zoonotic EID

ad_Y_zoo = matrix(NA, 8466, 2) 
ad_Y_zoo[1:8346, 1] = data_total$zoo_EID_counts_binom # Occurrence (8346 rows) in the first column, followed by 120 NAs for frequency
ad_Y_zoo[8347:8466, 2] = data_zoo$Ad_zoo_EID_counts_pois # Frequency data (8347 to 8466) in the second col

# nonzoonotic EID

ad_Y_notzoo = matrix(NA, 8426, 2) 
ad_Y_notzoo[1:8346, 1] = data_total$notzoo_EID_counts_binom # Occurrence (8346 rows) in the first column, followed by 80 NAs for frequency
ad_Y_notzoo[8347:8426, 2] = data_nonzoo$Ad_notzoo_EID_counts_pois # Frequency data (8347 to 8426) in the second col



# PC priors for random intercept

prior.prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01))) 
beta.hyper <- list(beta = list(prior = 'normal', param = c(1,10))) # default choice
rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.95)))

#str(data_total)

#m <- lm(zoo_EID_counts_binom~ Code+ Im_counts_scaled +  Ex_counts_scaled + Mammalian_richness_scaled +
#          Population_density_scaled + Popolation_growth + GDP_scaled + Precipitation_scaled+
 #         Agricultural_land_scaled + Forest_coverage_scaled + Temperature_scaled, data= data_total)

#car::vif(m)


# Model Analysis #

# ZINB for adj zoonotic EIDs

ldat_zoo_adEID_binom = list(
  Y = ad_Y_zoo,
  mu0 = rep(1:0, c(8346, 120)), 
  mu1 = rep(0:1, c(8346, 120)), 
  Country0 = c(as.factor(data_total$Code), rep(NA,120)),
  Imp0 = c(data_total$Im_counts_scaled, rep(NA,120)), 
  EX0 = c(data_total$Ex_counts_scaled, rep(NA,120)),
  Mam0 = c(data_total$Mammalian_richness_scaled, rep(NA,120)), 
  Pop_d0 = c(data_total$Population_density_scaled, rep(NA,120)),
  Pop_g0 = c(data_total$Popolation_growth, rep(NA,120)),
  GDP0 = c(data_total$GDP_scaled, rep(NA,120)),
  Pre0 = c(data_total$Precipitation_scaled, rep(NA,120)),
  Agr0 = c(data_total$Agricultural_land_scaled, rep(NA,120)),
  For0 = c(data_total$Forest_coverage_scaled, rep(NA,120)),
  Tem0 = c(data_total$Temperature_scaled, rep(NA,120)), 
  Country1 = c(rep(NA, 8346), data_zoo$Code), 
  Imp1 = c(rep(NA, 8346), data_zoo$Im_counts_scaled), 
  EX1 = c(rep(NA, 8346), data_zoo$Ex_counts_scaled),
  Mam1 = c(rep(NA, 8346), data_zoo$Mammalian_richness_scaled), 
  Pop_d1 = c(rep(NA, 8346), data_zoo$Population_density_scaled),
  Pop_g1 = c(rep(NA, 8346), data_zoo$Popolation_growth),
  GDP1 = c(rep(NA, 8346), data_zoo$GDP_scaled),
  Pre1 = c(rep(NA, 8346), data_zoo$Precipitation_scaled),
  Agr1 = c(rep(NA, 8346), data_zoo$Agricultural_land_scaled),
  For1 = c(rep(NA, 8346), data_zoo$Forest_coverage_scaled),
  Tem1 = c(rep(NA, 8346), data_zoo$Temperature_scaled), 
  Year0 = c(data_total$Year,rep(NA,120)), 
  Year1 = c(rep(NA,8346),data_zoo$Year))

W <- read.csv("scaled_matrix.csv",header = F)
W <- as.matrix(W)
D <- diag(colSums(W))
# Custom precision matrix incorporating asymmetry
Q <- D - W
diag(Q) <- 1
# Convert to a sparse matrix
sparse_matrix <- Matrix(Q, sparse = TRUE)
# Optional: Symmetrize the matrix if needed (depends on the context)
sparse_matrix <- (sparse_matrix + t(sparse_matrix)) / 2
# Ensure the matrix is symmetric and suitable for INLA (if required)
sparse_matrix <- forceSymmetric(sparse_matrix)



form_zoo_adEID_binom = Y ~ 0 + mu0 + mu1 + 
   EX0 + Mam0 + Pop_d0 + Pop_g0 + GDP0 + Pre0 + Agr0 + For0 + Tem0 +
   EX1 + Mam1 + Pop_d1 + Pop_g1 + GDP1 + Pre1 + Agr1 + For1 + Tem1 + 
  #f(Country0, model = "generic0", Cmatrix = sparse_matrix, rankdef=1, constr=TRUE, diagonal=1e-05)+
  #f(Country1,copy = "Country0", fixed = FALSE, hyper = beta.hyper)+
  f(Year0, model = "rw1", hyper = prior.prec) + f(Year1,copy = "Year0", fixed = FALSE, hyper = beta.hyper)

# ZINB model
res_zoo_zinb <- inla(form_zoo_adEID_binom, data = ldat_zoo_adEID_binom,
                      family=c('binomial', 'zeroinflatednbinomial0'), 
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                      control.family=list(list(), list(hyper = list(prob = list(initial = -20,fixed = TRUE))))) 
  
summary(res_zoo_zinb)

options(digits=8)



index <- 8347:8466
RMSE(index, 2, ad_Y_zoo, res_zoo_zinb$summary.fitted.values)

#res_zoo_zinb$summary.random$Country0[,c(1,3)]


#write.csv(res_zoo_zinb$summary.random$Country0[,c(1,3)],"spatial_sd.csv")


fitted_value <- cbind(res_zoo_zinb$summary.fitted.values[1:8346,c(1,3,5)],data_total$Im_counts)
fitted_value <- cbind(res_zoo_zinb$summary.fitted.values[1:8346,c(1,3,5)],data_total$Im_counts_scaled,data_total$zoo_EID_counts_binom)
fitted_value$residual <- fitted_value$`data_total$zoo_EID_counts_binom`-fitted_value$mean
plot(fitted_value$`data_total$Im_counts_scaled`, fitted_value$residual)


fitted_value <- cbind(res_zoo_zinb$summary.fitted.values[8347:8466,c(1,3,5)],data_zoo$Im_counts_scaled,data_zoo$Ad_zoo_EID_counts_pois)
fitted_value$residual <- fitted_value$`data_zoo$Ad_zoo_EID_counts_pois`-fitted_value$mean
fitted_value$residual

plot(fitted_value$`data_zoo$Im_counts_scaled`, fitted_value$residual,xlim = c(-0.5,2))





ggplot(fitted_value, aes(x = `data_total$Im_counts`, y = mean)) +
  geom_line(color = "#fb6a4b", size = 1) +
  geom_ribbon(aes(ymin ='0.025quant', ymax = '0.975quant'), alpha = 0.3, fill = "#fb6a4b") +
  labs(x = "X", y = "Predicted Probability", title = "Partial Effect of X") +
  theme_minimal(base_size = 14)




# Load the necessary libraries
library(pROC)
library(ggplot2)

# Assuming 'y_true' is the column of 0/1 actual event incidence
# and 'y_pred_prob' is the column of predicted probabilities from logistic regression

# Example:
y_true <- fitted_value[,2]
y_pred_prob <- fitted_value[,1]

# Compute ROC and AUC
roc_obj <- roc(y_true, y_pred_prob)
auc_value <- auc(roc_obj)

# Print AUC value
print(paste("AUC:", auc_value))

# Plot ROC curve
plot(roc_obj, col = "blue", main = "ROC Curve", print.auc = TRUE)

optimal_coords <- coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "closest.topleft")

# Print optimal threshold, sensitivity, and specificity
print(optimal_coords)




#Sensitivity non-zoonotic


ldat_notzoo_adEID_binom = list(
  Y = ad_Y_notzoo,
  mu0 = rep(1:0, c(8346, 80)), 
  mu1 = rep(0:1, c(8346, 80)), 
  Country0 = c(as.factor(data_total$Code), rep(NA,80)),
  Imp0 = c(data_total$Im_counts_scaled, rep(NA,80)), 
  EX0 = c(data_total$Ex_counts_scaled, rep(NA,80)), 
  Mam0 = c(data_total$Mammalian_richness_scaled, rep(NA,80)), 
  Pop_d0 = c(data_total$Population_density_scaled, rep(NA,80)),
  Pop_g0 = c(data_total$Popolation_growth, rep(NA,80)),
  GDP0 = c(data_total$GDP_scaled, rep(NA,80)),
  Pre0 = c(data_total$Precipitation_scaled, rep(NA,80)),
  Agr0 = c(data_total$Agricultural_land_scaled, rep(NA,80)),
  For0 = c(data_total$Forest_coverage_scaled, rep(NA,80)),
  Tem0 = c(data_total$Temperature_scaled, rep(NA,80)), 
  Country1 = c(rep(NA, 8346), as.factor(data_nonzoo$Code)), 
  Imp1 = c(rep(NA, 8346), data_nonzoo$Im_counts_scaled), 
  EX1 = c(rep(NA, 8346), data_nonzoo$Ex_counts_scaled),
  Mam1 = c(rep(NA, 8346), data_nonzoo$Mammalian_richness_scaled), 
  Pop_d1 = c(rep(NA, 8346), data_nonzoo$Population_density_scaled),
  Pop_g1 = c(rep(NA, 8346), data_nonzoo$Popolation_growth),
  GDP1 = c(rep(NA, 8346), data_nonzoo$GDP_scaled),
  Pre1 = c(rep(NA, 8346), data_nonzoo$Precipitation_scaled),
  Agr1 = c(rep(NA, 8346), data_nonzoo$Agricultural_land_scaled),
  For1 = c(rep(NA, 8346), data_nonzoo$Forest_coverage_scaled),
  Tem1 = c(rep(NA, 8346), data_nonzoo$Temperature_scaled), 
  Year0 = c(data_total$Year,rep(NA,80)), 
  Year1 = c(rep(NA,8346),data_nonzoo$Year))


prior.prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01))) 
beta.hyper <- list(beta = list(prior = 'normal', param = c(1,10))) # default choice


form_notzoo_adEID_binom = Y ~ 0 + mu0 + mu1 + 
  Imp0 + EX0 + Mam0 + Pop_d0 + Pop_g0 + GDP0 + Pre0 + Agr0 + For0 + Tem0 +
  Imp1 + EX1 + Mam1 + Pop_d1 + Pop_g1 + GDP1 + Pre1 + Agr1 + For1 + Tem1 + 
  f(Country1, model = "generic0", Cmatrix = sparse_matrix, rankdef=1, constr=TRUE, diagonal=1e-05)+
  f(Country0,copy = "Country1", fixed = FALSE, hyper = beta.hyper)+
  f(Year0, model = "rw1", hyper = prior.prec) + f(Year1,copy = "Year0", fixed = FALSE, hyper = beta.hyper)


res_notzoo_zinb2 <- inla(form_notzoo_adEID_binom, data=ldat_notzoo_adEID_binom,
                         family=c('binomial', 'zeroinflatednbinomial0'), 
                         #control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                         control.family=list(list(), list(hyper = list(prob = list(initial = -20,fixed = TRUE))))) 

summary(res_notzoo_zinb2)



