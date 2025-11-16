# packages
#library(ggplot2)
library(pacman)
p_load(dplyr, extRemes, fExtremes, MASS, scales, nortest, evd, POT, gridExtra, evmix)
# 导入数据
data <- read.csv("EVT.csv", header = T)


# extreme analysis

# 1: LB test, p value > alpha 
ts_eid <- ts(data$total_adjusted_zooEID_r)
Box.test(ts_eid, lag = 5, type = "Ljung-Box") # p-value = 
ts_eid <- ts(data$total_adjusted_nonzooEID_r)
Box.test(ts_eid, lag = 5, type = "Ljung-Box") # p-value = 
ts_im <- ts(data$total_Im_counts)
Box.test(ts_im, lag = 5, type = "Ljung-Box") # p-value = 
ts_ex <- ts(data$total_Ex_counts)
Box.test(ts_ex, lag = 5, type = "Ljung-Box") # p-value = 



## zoonotic VS import
Adzoo_im <- data[,c(3,5,6,9)]


Adzoo_im_model <- Adzoo_im[,c(3,4)]

model1_x <- quantile(Adzoo_im_model[,1], 0.8)
model1_y <- quantile(Adzoo_im_model[,2], 0.8)

model1 <- fitbvgpd(Adzoo_im_model, threshold = c(9.2, 81892.32), model = "nlog")
model1$chi  #Get the dependence coefficient Chi
AIC(model1)
print(model1)

model1_x <- quantile(Adzoo_im_model[,1], 0.81)
model1_y <- quantile(Adzoo_im_model[,2], 0.81)

model1 <- fitbvgpd(Adzoo_im_model, threshold = c(model1_x, model1_y), model = "nlog")
model1$chi  #Get the dependence coefficient Chi
AIC(model1)
print(model1)



## Export
Adzoo_ex <- data[,c(3,5,6,10)]
quantile(Adzoo_ex[,2],0.8)
quantile(Adzoo_ex[,3],0.8)
quantile(Adzoo_ex[,4],0.8)

Adzoo_ex_model <- Adzoo_ex[,c(2,4)]

model1_x <- quantile(Adzoo_ex_model[,1], 0.8)
model1_y <- quantile(Adzoo_ex_model[,2], 0.8)
model1 <- fitbvgpd(Adzoo_ex_model, threshold = c(9.2, 217621.7), model = "amix")
model1$chi  #Get the dependence coefficient Chi
AIC(model1)
print(model1)
model1$param
model1$std.err




Adnzoo_im <- data[,c(4,7,8,9)]
Adnzoo_im_model <- Adnzoo_im[,c(2,4)]
model2_x <- 2
model2_y <- quantile(Adnzoo_im_model[,2], 0.8)

quantile(Adnzoo_im_model[,1], 0.85)
quantile(Adnzoo_im_model[,2], 0.85)

model2 <- fitbvgpd(Adnzoo_im_model, threshold = c(12.1, 81892.32), model = "nlog")
model2$chi  #Get the dependence coefficient Chi
AIC(model2)
print(model2)


Adnzoo_Ex <- data[,c(4,7,8,10)]
Adnzoo_Ex_model <- Adnzoo_Ex[,c(2,4)]
quantile(Adnzoo_Ex_model[,2], 0.8)

model2 <- fitbvgpd(Adnzoo_Ex_model, threshold = c(12.1, 217621.7), model = "nlog")
model2$chi  #Get the dependence coefficient Chi
AIC(model2)
print(model2)
















Adzoo_ex <- data[,c(8,11)]
Adnzoo_ex <- data[,c(9,11)]





m1 <- fevd(Adzoo_im[,3],threshold = 9, type = "GP")
m2 <- fevd(Adzoo_im[,2],threshold = 81892.32, type = "GP")

m3 <- fevd(Adnzoo_ex[,1],threshold = 5, type = "GP")
m4 <- fevd(Adnzoo_ex[,2],threshold = 208286, type = "GP")

plot(m1)
plot(m2)
plot(m3)
plot(m4)

m3
m4


#Adjusted zoonotic EID Counts versus import WOE

threshold_Adzoo_im <- bvtcplot(Adzoo_im)
threshold_Adzoo_im$k0 # 45

par(mfrow=c(1,1))






#Adjusted nonzoonotic EID Counts versus import WOE

threshold_Adnzoo_im <- bvtcplot(Adnzoo_im)
threshold_Adnzoo_im$k0 # 41




#Adjusted zoonotic EID Counts versus export WOE

threshold_Adzoo_ex <- bvtcplot(Adzoo_ex)
threshold_Adzoo_ex$k0 # 45

model3_x <- 9
model3_y <- quantile(Adzoo_ex[,2], 0.8)
model3 <- fitbvgpd(Adzoo_ex, threshold = c(model3_x, model3_y), model = "nlog")
model3$chi  #Get the dependence coefficient Chi
AIC(model3)
print(model3)



jpeg("threshold_Adzoo_ex.jpeg", quality = 100, units = "in", width = 6, height = 6, res = 300)
#Print the plot or visualization
bvtcplot(Adzoo_ex)
# Close the jpeg device
dev.off()



#Adjusted nonzoonotic EID Counts versus export WOE
threshold_Adnzoo_ex <- bvtcplot(Adnzoo_ex)
threshold_Adnzoo_ex$k0 # 25

model4_x <- 5
model4_y <- quantile(Adnzoo_ex[,2], 0.8)
model4 <- fitbvgpd(Adnzoo_ex, threshold = c(model4_x, model4_y), model = "nlog")
model4$chi  #Get the dependence coefficient Chi
AIC(model4)
print(model4)



jpeg("threshold_Adnzoo_ex.jpeg", quality = 100, units = "in", width = 6, height = 6, res = 300)
#Print the plot or visualization
bvtcplot(Adnzoo_ex)
# Close the jpeg device
dev.off()



















