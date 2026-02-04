rm(list = ls())

setwd("~/Documents/lab/pop_size/EcolEvol_23sep25/rev1/analyses/")

library(ape)
library(funspace)
library(dplyr)
library(doParallel)
library(foreach)
library(ggplot2)

dat <- read.csv("Anura_data_complete_imputed.csv")

tree <- read.tree("amph_shl_new_Consensus_7238 (1).tre")
trees <- read.nexus("Anura_tree_Jetz_Pyron_2018.nex")

subdat <- subset(dat, select = c(Species, Family, Threat_category_IUCN,
                                 Population_trend_IUCN, SVL_Mean_databases,
                                 LatitudeModule, Range_size, Range_size_log,
                                 bio1, bio7, Moisture))

## excluding extinct species
subdat <- subset(subdat, !is.na(Population_trend_IUCN))
## final N=3271 spp (non-imputed); N=3845 spp (imputed).

## transforming unknown population trend into NA
subdat$Population_trend_IUCN <- replace(
  subdat$Population_trend_IUCN, subdat$Population_trend_IUCN == "Unknown", NA
) 

unknown <- sum(is.na(subdat$Population_trend_IUCN))

subdat$Population_trend_IUCN <- as.factor(subdat$Population_trend_IUCN)

## transforming variable
subdat$SVL_Mean_databases <- log10(subdat$SVL_Mean_databases)
scale_stats <- scale(subdat$Moisture)
subdat$Moisture <- scale_stats

tree <- keep.tip(tree, subdat$Species)
for (i in 1:length(trees)) trees[[i]] <- keep.tip(trees[[i]], subdat$Species)

subdat <- subdat[match(tree$tip.label, subdat$Species), ]
rownames(subdat) <- subdat$Species

imp_consensus <- impute(traits = subdat[, c(4:6, 8:11)], phylo = tree,
                        addingSpecies = FALSE)
subdat_imp <- imp_consensus$imputed

## validation with the consensus tree
subdat_without_NA <- subdat[complete.cases(subdat), ]
tree_without_NA <- keep.tip(tree, subdat_without_NA$Species)
subdat_without_NA <- subdat_without_NA[match(tree_without_NA$tip.label,
                                             subdat_without_NA$Species), ]

cross_val <- function(dat, tr) {

  n_species <- nrow(dat)
  predictions <- dat[, c(4:6, 8:11)]
  predictions[, c(1:2, 5:7)] <- NA
  
  for (i in 145:1) {
    print(i)
    train_data <- dat[, c(4:6, 8:11)]
    train_data[i, c(1:2, 5:7)] <- NA
    
    model_impute <- impute(traits = train_data, phylo = tr, addingSpecies = F)
    
    sp_name <- dat$Species[i]
    
    imputed_values <- model_impute$imputed[sp_name, ]
    predictions[i, c(1:2, 5:7)] <- imputed_values[, c(1:2, 5:7)]
    saveRDS(predictions, file = paste0("imputation/predictions_", i, ".rds"))
  }
  
  return(predictions)
}

cross_val_dat <- cross_val(subdat_without_NA, tree_without_NA)
cross_val_dat$Population_trend_IUCN <- as.character(cross_val_dat$Population_trend_IUCN)
cross_val_dat$Population_trend_IUCN[cross_val_dat$Population_trend_IUCN == "1"] <- "Decreasing" 
cross_val_dat$Population_trend_IUCN[cross_val_dat$Population_trend_IUCN == "2"] <- "Increasing" 
cross_val_dat$Population_trend_IUCN[cross_val_dat$Population_trend_IUCN == "3"] <- "Stable" 

#saveRDS(cross_val_dat, file = "cross_val_dat.rds")
#cross_val_dat <- readRDS(file = "cross_val_dat.rds")

sp_mean_abs_error <- colMeans(cross_val_dat[, -c(1, 3:4)] - subdat_without_NA[, -c(1:4, 6:8)])
sp_mean_sq_error <- colMeans((cross_val_dat[, -c(1, 3:4)] - subdat_without_NA[, -c(1:4, 6:8)])^2)

sp_mean_abs_error_data <- sp_mean_abs_error %>%
  as.data.frame() %>%
  dplyr::rename(mean_abs_error = ".")

sp_mean_sq_error_data <- sp_mean_sq_error %>%
  as.data.frame() %>%
  dplyr::rename(sp_mean_sq_error = ".")

results <- cbind(round(sp_mean_abs_error_data, 3),
                 round(sp_mean_sq_error_data, 3),
                 pred_coeff = NA)

R2_pred <- function(obs, pred) {
  res <- 1 - mean((obs - pred)^2, na.rm = TRUE) / var(obs, na.rm = TRUE)
  res <- round(res, 3)
  return(res)
}

results["SVL_Mean_databases", "pred_coeff"] <- 
  R2_pred(subdat_without_NA$SVL_Mean_databases, cross_val_dat$SVL_Mean_databases)
results["bio1", "pred_coeff"] <- 
  R2_pred(subdat_without_NA$bio1, cross_val_dat$bio1)
results["bio7", "pred_coeff"] <- 
  R2_pred(subdat_without_NA$bio7, cross_val_dat$bio7)
results["Moisture", "pred_coeff"] <- 
  R2_pred(subdat_without_NA$Moisture, cross_val_dat$Moisture)

library(caret) 

obs_val <- subdat_without_NA$Population_trend_IUCN
pred_val <- cross_val_dat$Population_trend_IUCN

levels_order <- c("Decreasing", "Stable", "Increasing")

obs_factor <- factor(obs_val, levels = levels_order)
pred_factor <- factor(pred_val, levels = levels_order)

cm <- confusionMatrix(pred_factor, obs_factor) 

results <- as.data.frame(rbind(results, NA))
rownames(results)[5] <- "Population_trend_IUCN"

results$accuracy <- c(NA, NA, NA, NA, round(cm$overall["Accuracy"], 3))
results$kappa <- c(NA, NA, NA, NA, round(cm$overall["Kappa"], 3))

write.csv(results, "cross_val.csv")

## repeating for 100 topologies

cl <- makeCluster(2)
registerDoParallel(cl)

foreach(i = 1:length(trees), .packages = c("funspace", "ape")) %dopar% {
  
  trees_i <- trees[[i]]
  
  subdat_i <- subdat[match(trees_i$tip.label, subdat$Species), ]
  subdat_i <- subdat_i[, c(4:6, 8:11)]
  
  imp_i <- impute(traits = subdat_i, phylo = trees_i, addingSpecies = FALSE)
  
  write.csv(imp_i$imputed, file = paste0("imputation/dat_imp_", i, ".csv"))
  
  NULL
}

stopCluster(cl)

for (i in 1:length(trees)) {
  print(i)
	trees_i <- trees[[i]]
	subdat_i <- subdat[match(trees_i$tip.label, subdat$Species), ]

	subdat_i <- subdat_i[, c(4:6, 8:11)]

	imp_i <- impute(traits = subdat_i, phylo = trees_i, addingSpecies = FALSE)

	write.csv(imp_i$imputed, file = paste0("imputation/dat_imp_", i, ".csv"))
}

subdat <- subdat %>% 
  mutate(SVL_Mean_databases = 10^(SVL_Mean_databases),
         Moisture = (Moisture * attr(scale_stats, "scaled:scale")) + attr(scale_stats, "scaled:center"))

subdat_imp <- subdat_imp %>% 
    tibble::rownames_to_column(var = "Species") %>%
    left_join(subdat[, c("Species", "Family", "Threat_category_IUCN", "Range_size")], by = "Species") %>%
  mutate(SVL_Mean_databases = 10^(SVL_Mean_databases),
         Moisture = (Moisture * attr(scale_stats, "scaled:scale")) + attr(scale_stats, "scaled:center"))

files_imp <- list.files("imputation/", full.names = T)

imp_dat <- purrr::map(files_imp, \(file_path) {
  dat <- read.csv(file_path, row.names = 1) %>% 
    tibble::rownames_to_column(var = "Species") %>%
    left_join(subdat[, c("Species", "Family", "Threat_category_IUCN", "Range_size")], by = "Species") %>%
    mutate(SVL_Mean_databases = 10^(SVL_Mean_databases),
           Moisture = (Moisture * attr(scale_stats, "scaled:scale")) + attr(scale_stats, "scaled:center"))
})

## calculating prevalence

## searching for the prevalence after imputation (because it is better than imput the prevalence - this way, we can have a sistematic prevalence value)

prev_lut <- read.csv("prevalence_lut.csv")
prev_lut[is.na(prev_lut)] <- 0

prev_lut2 <- read.csv("prevalence_lut2.csv")
prev_lut2[is.na(prev_lut2)] <- 0

## bio1 (temperature)
labelstemplutcut <- c(prev_lut$bio1-0.5, 32.5) # breaks (I added 32.5 to be the 
## upper limit of the label "32", that is present in bio1 column
labelstemplut <- c(prev_lut$bio1) # labels
prevalencetemplut <- c(prev_lut$prevBio1) # backup variable

## the same, but with LUT2 (double spaced intervals)
labelstemplut2cut <- c(prev_lut2$bio1lut2-0.5, 32.5)
labelstemplut2 <- c(prev_lut2$bio1lut2)
prevalencetemplut2 <- c(prev_lut2$prevBio1lut2)

## transforming the continuous variable in intervals
subdat$bio1_int <- cut(subdat$bio1, breaks = labelstemplutcut, 
                       labels = labelstemplut, 
                       include.lowest = FALSE)

subdat_imp$bio1_int <- cut(subdat_imp$bio1, breaks = labelstemplutcut, 
                           labels = labelstemplut, 
                           include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$bio1_int <- cut(x$bio1, breaks = labelstemplutcut, labels = labelstemplut, include.lowest = FALSE)
  return(x)
})

## doing the same as above, but including the prevalence as label, to assess  
## the prevalence
subdat$prev_bio1 <- cut(subdat$bio1, breaks = labelstemplutcut, 
                        labels = prevalencetemplut, 
                        include.lowest = FALSE)

subdat_imp$prev_bio1 <- cut(subdat_imp$bio1, breaks = labelstemplutcut, 
                            labels = prevalencetemplut, 
                            include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$prev_bio1 <- cut(x$bio1, breaks = labelstemplutcut, 
                     labels = prevalencetemplut, include.lowest = FALSE)
  return(x)
})

subdat$bio1_int2 <- cut(subdat$bio1, breaks = labelstemplut2cut, 
                        labels = labelstemplut2, 
                        include.lowest = FALSE)

subdat_imp$bio1_int2 <- cut(subdat_imp$bio1, breaks = labelstemplut2cut, 
                            labels = labelstemplut2, 
                            include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$bio1_int2 <- cut(x$bio1, breaks = labelstemplut2cut, 
                     labels = labelstemplut2, include.lowest = FALSE)
  return(x)
})

subdat$prev2_bio1 <- cut(subdat$bio1, breaks = labelstemplut2cut, 
                         labels = prevalencetemplut2, 
                         include.lowest = FALSE)

subdat_imp$prev2_bio1 <- cut(subdat_imp$bio1, breaks = labelstemplut2cut, 
                             labels = prevalencetemplut2, 
                             include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$prev2_bio1 <- cut(x$bio1, breaks = labelstemplut2cut, 
                      labels = prevalencetemplut2, 
                      include.lowest = FALSE)
  return(x)
})

## moisture CHELSA
labelsmoislutcut <- c(-200, prev_lut$Moisture - 11.70667/2, 14000)
labelsmoislut <- c(min(subdat$Moisture, na.rm = T), prev_lut$Moisture)
prevalencemmoislut <- c(0, prev_lut$prevMoisture)

## the same, but with LUT2 (double spaced intervals)
labelsmoislut2cut <- sort(c(-200, prev_lut2$Moisture - 200/2, 14000))
labelsmoislut2 <- c(min(subdat$Moisture, na.rm = T),prev_lut2$Moisture)
prevalencemmoislut2 <- c(0, prev_lut2$prevMoisture)

subdat$Moisture_int <- cut(subdat$Moisture, breaks = labelsmoislutcut, 
                           labels = labelsmoislut,
                           include.lowest = FALSE)

subdat_imp$Moisture_int <- cut(subdat_imp$Moisture, breaks = labelsmoislutcut, 
                               labels = labelsmoislut,
                               include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$Moisture_int <- cut(x$Moisture, breaks = labelsmoislutcut, 
                        labels = labelsmoislut,
                        include.lowest = FALSE)
  return(x)
})

subdat$prev_Moisture <- cut(subdat$Moisture, breaks = labelsmoislutcut, 
                            labels = prevalencemmoislut, 
                            include.lowest = FALSE)

subdat_imp$prev_Moisture <- cut(subdat_imp$Moisture, breaks = labelsmoislutcut, 
                                labels = prevalencemmoislut, 
                                include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$prev_Moisture <- cut(x$Moisture, breaks = labelsmoislutcut, 
                         labels = prevalencemmoislut, 
                         include.lowest = FALSE)
  return(x)
})

subdat$Moisture_int2 <- cut(subdat$Moisture, breaks = labelsmoislut2cut, 
                            labels = labelsmoislut2, 
                            include.lowest = FALSE)

subdat_imp$Moisture_int2 <- cut(subdat_imp$Moisture, breaks = labelsmoislut2cut,
                                labels = labelsmoislut2, 
                                include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$Moisture_int2 <- cut(x$Moisture, breaks = labelsmoislut2cut, 
                         labels = labelsmoislut2, 
                         include.lowest = FALSE)
  return(x)
})

subdat$prev2_Moisture <- cut(subdat$Moisture, breaks = labelsmoislut2cut, 
                             labels = prevalencemmoislut2, 
                             include.lowest = FALSE)

subdat_imp$prev2_Moisture <- cut(subdat_imp$Moisture, 
                                 breaks = labelsmoislut2cut, 
                                 labels = prevalencemmoislut2, 
                                 include.lowest = FALSE)

imp_dat <- lapply(imp_dat, function(x) {
  x$prev2_Moisture <- cut(x$Moisture, breaks = labelsmoislut2cut, 
                          labels = prevalencemmoislut2, 
                          include.lowest = FALSE)
  return(x)
})

write.csv(subdat, file = "subdat.csv", row.names = FALSE)
write.csv(subdat_imp, file = "subdat_imp.csv", row.names = FALSE)
saveRDS(imp_dat, file = "imp_dat.rds")

## Figure 1

family <- unique(subdat_imp$Family)

keep <- subset(subdat_imp, Family==family[1], select=c("Family", "Species"))
keep <- keep[sample(nrow(keep), 1), ]

for (i in 1:length(family)) {
  keep[i, ] <- subset(subdat_imp, Family == family[i], 
                      select = c("Family", "Species"))[sample(nrow(keep[i, ]), 1), ]
}

orders.tr <- keep.tip(tree, as.character(keep[, 2]))

for (k in 1:length(family)) {
  orders.tr$tip.label[k] <- as.character(keep[, 1][keep[, 2]==orders.tr$tip.label[k]])
}

orders <- orders.tr$tip.label
aa <- table(subdat_imp$Family, subdat_imp$Population_trend_IUCN)

GG <- as.data.frame(aa)
GG$Var1 <- factor(GG$Var1, levels = orders)
GG$Var2 <- as.character(GG$Var2)
GG$Var2 <- ifelse(GG$Var2 == "Increasing", "c", GG$Var2)
GG$Var2 <- ifelse(GG$Var2 == "Stable", "b", GG$Var2)
GG$Var2 <- ifelse(GG$Var2 == "Decreasing", "a", GG$Var2)

plot1 <- ggplot(GG, aes(fill = Var2, y = Var1, x = Freq)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#8c4733", "#ecc05f", "#93afa7")) +
  theme_classic() +
  xlab("Number of species") +
  ylab("Family")

pdf("Fig1_bars_rev1.pdf")
print(plot1)
dev.off()

pdf("Fig1_tree_rev1.pdf")
plot(orders.tr, show.node.label = TRUE, cex = 0.8, align.tip.label = TRUE, direction = "right")
dev.off()
