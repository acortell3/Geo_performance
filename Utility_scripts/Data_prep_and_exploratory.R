

### LOAD DATASETS AND MATRICES

################################################################################
###########################       BASIC DATA       #############################
################################################################################

### Load metric info geos and throw results

## Add datasets with the characteristics of the geos
D1 <- read.csv("./Data/Geo_throws_D1.csv")
D1$Day <- rep("D1", nrow(D1))

D2 <- read.csv("./Data/Geo_throws_D2.csv")
D2$Day <- rep("D2", nrow(D2))

D3 <- read.csv("./Data/Geo_throws_D3.csv")
D3$Day <- rep("D3", nrow(D3))

Geo_dat <- rbind(D1,D2,D3)

## Substitute commas by points
Geo_dat[,c(6:10)] <- apply(Geo_dat[,c(6:10)], 2, function(x) {as.numeric(gsub(",", ".", x))})

## Dataset with the results of the throws
Throws <- read.csv("./Data/Results_clean.csv", header = TRUE)
Throws <- Throws[,-c(1,ncol(Throws))] ## Remove columns of day and machine ID (already in the other df)

Exp_data <- cbind(Geo_dat,Throws)

## Include Original type ID (for comparison)
Rep.type <- c(rep(2622,2),rep(2321,2),rep(658,2),rep(240,2),rep(2831,2),rep(126,2),rep(154,2),
              rep(109,2),rep(2836,2),rep(187,2),rep(101,2),rep(2181,2),rep(23,2),rep(1810,2),2321)
Exp_data <- cbind(Exp_data,Rep.type)
write.csv(Exp_data,"./Data/Exp_data.csv")

################################################################################
###########################          GMM           #############################
################################################################################


## Load outlines for GMM analysis and perform PCA to extract morphometric values
library(Momocs)
set.seed(12345)

## Extract outlines
dir <- getwd()
setwd("./Outlines/")
geo_list <- Exp_data$ID.geo[1:29]
geo_list <- paste0(geo_list,".jpg")
geo_out <- Out(import_jpg(geo_list, auto.notcentered = T), Exp_data$ID.geo[1:29])
setwd(dir)

## Elliptical fourier analysis at n=7 harmonics
geo_out_F <- efourier(coo_center(coo_scale(geo_out)), nb.h = 7, norm=F, start=T)

## PCA
geo_out_pca <- PCA(geo_out_F)

## Figure 

png("./Figures/Figure_GMM.png", width = 1600, height = 1200, res = 200)
plot(geo_out_pca, col = "brown4", cex = 1.2, pch = 20, col.shp = "lightblue", 
     border.shp = "blue3") ## 2PCs justify 88.2% variance
dev.off()

Full_data <- cbind(Exp_data,
                   "PC1" = rep(geo_out_pca$x[,1],3),
                   "PC2" = rep(geo_out_pca$x[,2],3))

colnames(Full_data) <- c("ID_machine","ID_geo","Symmetry","Inclination","Break",
                         "Parallel_width_shaft","Perpendicular_width_shaft",
                         "Width_before_shaft","Max_width","Distance_tip_to_shaft",
                         "Day","Area_under_curve","Max_force","Displacement_at_max_force",
                         "Modulus_automatic_young","Break_location","Compressive_displacement_max_force",
                         "Compressive_strain_displacement_max_force","Compressive_stress_at_max_force",
                         "Energy_at_max_force","Energy_to_X_intercept_at_modulus_automatic_young",
                         "Rep_type","PC1","PC2")

write.csv(Full_data,"./Data/Full_data.csv")

################################################################################
#########       EXPLORATION AND PRELIMINARY VARIABLE SELECTION       ###########
################################################################################

## Remove unnecessary variables
## Remove ID machine, "Symmetry" (as it is captures by GMM),
## "break" (because neither geometric broke), "day" and "unnecesary 
## machine info"
Sel_data <- Full_data[,-c(1,3,5,11,15:22)]

## Explore variable collinearity (of hand metrics)
shaft <- data.frame("Par" = Full_data$Parallel_width_shaft,
                    "Per" = Full_data$Perpendicular_width_shaft)

Colli_data <- cbind(Sel_data[,c(5:7)],"Shaft" = apply(shaft,1,mean))
psych::pairs.panels(Colli_data) ## Not enought collinearity to discard anything

## Let's try a PCA, to see if we can reduce the dimensionality of our dataset
## while still having MEANINGFUL co-variates
met_var <- Sel_data[,c(3:7)]

# Because there are sinificant differences in scale, we use a correlation matrix
#apply(met_var,2,max)-apply(met_var,2,min)
pca_met <- princomp(met_var, cor = TRUE)

## See pca values
summary(pca_met) ## Two components explain ~70% and three ~85%

## scree plot
plot(pca_met, type = "l")

## loadings
pca_met$loadings

## loading plot
plot(pca_met$loadings[, 1:2], type = "n", xlim = c(-1, 1), ylim = c(-1,1))
for (i in 1:nrow(pca_met$loadings)) {
  arrows(0, 0, pca_met$loadings[i, 1], pca_met$loadings[i, 2])
}
text(pca_met$loadings[, 1:2], dimnames(pca_met$loadings)[[1]])
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

## Max width seems to follow its own path, as does distance tip to shaft.
## Shaft-related width seem to pull in the same direction
## Two clear MEANINGFUL components cannot be extracted (e.g. loadings are 
## too spread to give a single meaning to any individual component), but 
## this can be taken into account in the next analysis

## We CANNOT discard co-variates based on linear simmilarity. Thus, the
## data set goes like this to the model, and it will be further studied there

## Outlier detection, based on Mahalanobis distance
out_de <- df[,-1]
par(mfrow = c(1,1))
ma <- mahalanobis(out_de, apply(out_de, 2, mean, na.rm = TRUE), cov(out_de, use = "na.or.complete"))
k <- dim(out_de)[2] #Covariate number
Lim <- k + 3 * sqrt(k * 2) #Mahalanobis distance limit
plot(ma, pch = 20, ylim = c(0, max(ma, Lim, na.rm = TRUE)))
text(ma, rownames(out_de), pos = 2)
abline(h = Lim, col= "red")
title("Mahalanobis distance on covariates")

## Check response variable "displacement" for censoring
disp_check <- Sel_data$Displacement_at_max_force
disp_check_29.5 <- disp_check[disp_check<29.5]
disp_check_30 <- disp_check[disp_check<30]

## Censoring will be at 29.5
hist(disp_check, breaks = 30)
hist(disp_check_29.5, breaks = 30)
hist(disp_check_30, breaks = 30)


Sel_data$Inclination <- as.factor(Sel_data$Inclination)

write.csv(Sel_data,"./Data/Sel_data.csv")


