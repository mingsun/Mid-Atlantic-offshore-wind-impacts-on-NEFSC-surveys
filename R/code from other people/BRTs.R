library(dismo)
library(gbm)

##Load calibrated data
LRS_1982 <- read.csv("results/stratified.mean.indices/summer.flounder/bay.anchovy.pysl.calibrated.csv")
LRS_1982=LRS_1982[LRS_1982$year>=1988,]

#remove uneeded columns
LRS_1982=LRS_1982[c(18,19,21:24)]
#removed DO
LRS_1982=LRS_1982[-c(3)]

#scale habitat data- this step is important, allows for direct comparison of covariate importance
##but not necessarily needed, I've come across papers who don't scale data first
LRS_1982[c(2:5)]<- scale(LRS_1982[c(2:5)])

#create matrix where the results from the BRT runs will be save, ncol = 10 (or however many simulations), nrow = number of habitat variables
BRT_res=matrix(data=0, ncol=10, nrow=3)

#names for how the variables are in the dataframe, don't need to be exactly the same
rownames(BRT_res)=c("Temp", "Conductivity", "Depth")
#gbm.x are the column numbers of the habitat variables you are testing, gbm.y is the column number for bay anchovy abundance. 
for (i in 1:10){          
  LRS_brt = gbm.step(data = LRS_1982, gbm.x = c(2:4), gbm.y = 1, family = "gaussian", tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.75)
  q = summary(LRS_brt)
  reorder = c( "WATER.TEMPERATURE", "CONDUCTIVITY", "RIVER.DEPTH")#reorder them so that they are in the same order as you made them 4 lines up. In this case they do have to match the exact column names
  q = q[match(reorder, q$var),]
  BRT_res[,i] = q$rel.inf
}
LRS_inf=as.data.frame(rowMeans(BRT_res))
summary(LRS_inf)
# LRS_Inf = as.data.frame(rowMeans(LRS_brt))

#Code for just running one simulation
#the tree complexity I tested are 1,3,5 and the learning rates I tested with 0.005, 0.01, 0.02, and 0.03. The bag fraction was always kept at 0.75.
#You optimize the BRTs by messing around with the TC and LR to get the model that has the lowest estimated cv deviance
#Only use models that create at least 1000 trees
LRS_brt = gbm.step(data = LRS_1982, gbm.x = c(2:4), gbm.y = 1, family = "gaussian", tree.complexity = 5, learning.rate = 0.03, bag.fraction = 0.75)
summary(LRS_brt)

#Load uncalibrated data
LRS_1982=read.csv("~/Desktop/HRBMP/BayAnchovy/Data/uncalibrated0118.csv")
LRS_1982=LRS_1982[-c(1)]
LRS_1982=LRS_1982[c(9,15,16,17,22)]
#scale habitat data
LRS_1982[c(1:4)]<- scale(LRS_1982[c(1:4)])
LRS_brt = gbm.step(data = LRS_1982, gbm.x = c(1,2,4), gbm.y = 6, family = "gaussian", tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.75)

write.csv(LRS_1982, "uncalibratedBRTdata.csv")
