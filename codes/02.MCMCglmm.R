rm(list = ls())

setwd("~/Documents/lab/pop_size/EcolEvol_23sep25/rev1/analyses")

library(phytools)
library(MCMCglmm)
library(foreach)
library(doParallel)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringi)
library(phylolm)

dat <- read.csv("subdat.csv")
dat_imp <- read.csv("subdat_imp.csv")

dat_rep_imp <- readRDS("imp_dat.rds")

tree <- read.tree("amph_shl_new_Consensus_7238 (1).tre")
trees <- read.nexus("Anura_tree_Jetz_Pyron_2018.nex")

subdat <- subset(dat, select = c(Species, Family, Threat_category_IUCN,
                                 Population_trend_IUCN, 
                                 SVL_Mean_databases, LatitudeModule, 
                                 Range_size_log, bio1, bio7, 
                                 Moisture, prev_bio1, prev2_bio1,
                                 prev_Moisture, prev2_Moisture))

subdat <- subdat[complete.cases(subdat), ]

subdat_imp <- subset(dat_imp, select = c(Species, Family, Threat_category_IUCN,
                                         Population_trend_IUCN, 
                                         SVL_Mean_databases, LatitudeModule, 
                                         Range_size_log, bio1, bio7, 
                                         Moisture, prev_bio1, prev2_bio1,
                                         prev_Moisture, prev2_Moisture))

subdat_rep_imp <- lapply(dat_rep_imp, function(x) {
  subset(x, select = c(Species, Family, Population_trend_IUCN, 
                       SVL_Mean_databases, LatitudeModule, 
                       Range_size_log, bio1, bio7, 
                       Moisture, prev_bio1, prev2_bio1,
                       prev_Moisture, prev2_Moisture))
})

## transforming the population trend into a binary data decreasing //// stable 
## (binary incomplete), removing the increasing species considering 
## "increasing" as NA because we will remove NA values, so they will be removed

subdat$Population_trend_binary <- subdat$Population_trend_IUCN
subdat$Population_trend_binary <- as.character(subdat$Population_trend_binary)
subdat$Population_trend_binary <- ifelse(subdat$Population_trend_binary == "Increasing", NA, subdat$Population_trend_binary)
subdat$Population_trend_binary <- ifelse(subdat$Population_trend_binary == "Stable", 1, subdat$Population_trend_binary)
subdat$Population_trend_binary <- ifelse(subdat$Population_trend_binary == "Decreasing", 0, subdat$Population_trend_binary)
subdat$Population_trend_binary <- factor(subdat$Population_trend_binary, levels = c(1, 0))

subdat_imp$Population_trend_binary <- subdat_imp$Population_trend_IUCN
subdat_imp$Population_trend_binary <- as.character(subdat_imp$Population_trend_binary)
subdat_imp$Population_trend_binary <- ifelse(subdat_imp$Population_trend_binary == "Increasing", NA, subdat_imp$Population_trend_binary)
subdat_imp$Population_trend_binary <- ifelse(subdat_imp$Population_trend_binary == "Stable", 1, subdat_imp$Population_trend_binary)
subdat_imp$Population_trend_binary <- ifelse(subdat_imp$Population_trend_binary == "Decreasing", 0, subdat_imp$Population_trend_binary)
subdat_imp$Population_trend_binary <- factor(subdat_imp$Population_trend_binary, levels = c(1, 0))

subdat_rep_imp <- lapply(subdat_rep_imp, function(x) {
  x$Population_trend_binary <- x$Population_trend_IUCN
  x$Population_trend_binary <- as.character(x$Population_trend_binary)
  x$Population_trend_binary <- ifelse(x$Population_trend_binary == "Increasing", NA, x$Population_trend_binary)
  x$Population_trend_binary <- ifelse(x$Population_trend_binary == "Stable", 1, x$Population_trend_binary)
  x$Population_trend_binary <- ifelse(x$Population_trend_binary == "Decreasing", 0, x$Population_trend_binary)
  x$Population_trend_binary <- factor(x$Population_trend_binary, levels = c(1, 0))
  return(x)
})

## transforming variables
subdat$SVL_log <- log10(subdat$SVL_Mean_databases)
subdat$bio1_scale <- scale(subdat$bio1, center = F)
subdat$Moisture_scale <- scale(subdat$Moisture, center = F)
subdat$prev_bio1_scale <- scale(subdat$prev_bio1, center = F)
subdat$prev_Moisture_scale <- scale(subdat$prev_Moisture, center = F)
subdat$Latitude_scale <- scale(subdat$LatitudeModule, center = F)

subdat_imp$SVL_log <- log10(subdat_imp$SVL_Mean_databases)
subdat_imp$bio1_scale <- scale(subdat_imp$bio1, center = F)
subdat_imp$Moisture_scale <- scale(subdat_imp$Moisture, center = F)
subdat_imp$prev_bio1_scale <- scale(subdat_imp$prev_bio1, center = F)
subdat_imp$prev_Moisture_scale <- scale(subdat_imp$prev_Moisture, center = F)
subdat_imp$Latitude_scale <- scale(subdat_imp$LatitudeModule, center = F)

subdat_rep_imp <- lapply(subdat_rep_imp, function(x) {
  x$SVL_log <- log10(x$SVL_Mean_databases)
  x$bio1_scale <- scale(x$bio1, center = F)
  x$Moisture <- as.numeric(x$Moisture)
  x$Moisture_scale <- scale(x$Moisture, center = F)
  x$prev_bio1 <- as.numeric(x$prev_bio1)
  x$prev_bio1_scale <- scale(x$prev_bio1, center = F)
  x$prev_Moisture <- as.numeric(x$prev_Moisture)
  x$prev_Moisture_scale <- scale(x$prev_Moisture, center = F)
  x$Latitude_scale <- scale(x$LatitudeModule, center = F)
  return(x)
})

## 

tab_cont <- subdat[complete.cases(subdat$Threat_category_IUCN), ]
tab_cont_imp <- subdat_imp[complete.cases(subdat_imp$Threat_category_IUCN), ]

rownames(tab_cont) <- tab_cont$Species
rownames(tab_cont_imp) <- tab_cont_imp$Species

m_phylo <- phyloglm(Population_trend_binary ~ Threat_category_IUCN, 
                    data = tab_cont,  
                    phy = tree, method = "logistic_MPLE")
m_phylo_imp <- phyloglm(Population_trend_binary ~ Threat_category_IUCN, 
                        data = tab_cont_imp,  
                        phy = tree, method = "logistic_MPLE")

summary(m_phylo)
summary(m_phylo_imp)

## testing and model selection wth the consensus tree

## manual model selection
predictors <- c("SVL_log", "Range_size_log", "bio1_scale", "bio7", 
                "Moisture_scale", "prev_bio1_scale", "prev_Moisture_scale", 
                "Latitude_scale")

all_formulas <- list()
formula_index <- 1
for (i in 0:length(predictors)) {
  combinations <- combn(predictors, i, simplify = FALSE)
  
  for (vars in combinations) {
    if (length(vars) == 0) {
      formula_str <- "Population_trend_binary ~ 1"
    } else {
      formula_str <- paste("Population_trend_binary ~", paste(vars, collapse = " + "))
    }
    
    all_formulas[[formula_index]] <- as.formula(formula_str)
    formula_index <- formula_index + 1
  }
}

nitt_val <- 1000000
burnin_val <- 100000
thin_val <- 500

num_cores <- 5

prior_mcmc <- list(G = list(G1 = list(V = 1, nu = 0.02)),
                   R = list (V = 1, fix = 1))

## non-imputed

cl <- makeCluster(num_cores)
registerDoParallel(cl)

dir_output <- "partials_MCMCglmm"
if (!dir.exists(dir_output)) dir.create(dir_output)

Ainv <- inverseA(force.ultrametric(keep.tip(tree, subdat$Species)))$Ainv

foreach(i = 1:length(all_formulas), .packages = c("MCMCglmm")) %dopar% {
  
  current_formula <- all_formulas[[i]]

  model_mcmc <- MCMCglmm(
    fixed = current_formula,
    random = ~ Species,
    family = "categorical",
    data = subdat,
    ginverse = list(Species = Ainv),
    prior = prior_mcmc,
    nitt = nitt_val,
    burnin = burnin_val,
    thin = thin_val,
    verbose = FALSE
  )
  
  current_dic <- model_mcmc$DIC
  
  formula_txt <- format(current_formula) 
  
  res_row <- data.frame(
    Model_ID = i,
    Formula = formula_txt,
    DIC = current_dic,
    stringsAsFactors = FALSE
  )
  
  file_name <- paste0(dir_output, "/model_", i, "_dic.rds")
  
  saveRDS(res_row, file = file_name)

  return(NULL)
}

stopCluster(cl)

print("the end")

## imputed

cl <- makeCluster(num_cores)
registerDoParallel(cl)

dir_output <- "partials_imputed_MCMCglmm"
if (!dir.exists(dir_output)) dir.create(dir_output)

Ainv_imp <- inverseA(phytools::force.ultrametric(keep.tip(tree, subdat_imp$Species)))$Ainv

foreach(i = 1:length(all_formulas), .packages = c("MCMCglmm")) %dopar% {
  
  current_formula <- all_formulas[[i]]

  model_mcmc <- MCMCglmm(
    fixed = current_formula,
    random = ~ Species,
    family = "categorical",
    data = subdat_imp,
    ginverse = list(Species = Ainv_imp),
    prior = prior_mcmc,
    nitt = nitt_val,
    burnin = burnin_val,
    thin = thin_val,
    verbose = FALSE
  )
  
  current_dic <- model_mcmc$DIC
  
  formula_txt <- format(current_formula) 
  
  res_row <- data.frame(
    Model_ID = i,
    Formula = formula_txt,
    DIC = current_dic,
    stringsAsFactors = FALSE
  )
  
  file_name <- paste0(dir_output, "/model_", i, "_dic.rds")
  
  saveRDS(res_row, file = file_name)

  return(NULL)
}

stopCluster(cl)

print("the end")

files <- list.files("partials_MCMCglmm", full.names = T, pattern = "\\.rds$")
files_imp <- list.files("partials_imputed_MCMCglmm", full.names = T, pattern = "\\.rds$")

dic_final <- files %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  arrange(DIC)

run_top_models <- function(dic_df) {
  dic_df$delta_DIC <- dic_df$DIC - min(dic_df$DIC)
  
  top_models_info <- dic_df[dic_df$delta_DIC < 2, ]
  top_ids <- unique(top_models_info$Model_ID)
  
  print(top_models_info)
  
  res <- lapply(top_ids, function(id) {
    current_formula <- all_formulas[[id[1]]]
    print(current_formula)
    
    model_fit <- MCMCglmm(
      fixed = current_formula,
      random = ~ Species,
      family = "categorical",
      data = subdat,
      ginverse = list(Species = Ainv),
      prior = prior_mcmc,
      nitt = nitt_val, 
      burnin = burnin_val,
      thin = thin_val,
      verbose = FALSE
  )
    
    return(list(id = id, model = model_fit, formula = current_formula))
  })
  
  names(res) <- paste0("Mod_", top_ids)
  
  return(res)
}

best_models <- run_top_models(dic_df = dic_final)

coef_comparison <- do.call(rbind, lapply(best_models, function(m) {
  res <- as.data.frame(summary(m$model)$solutions)
  res$Variable <- rownames(res)
  res$Model_ID <- m$id
  return(res)
}))

best_model <- all_formulas[[245]]

dic_imp_final <- files_imp %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  arrange(DIC)

best_model_imp <- all_formulas[[255]]

## repeating for 100 alternative topologies

## non_imputed
dir_output <- "topologies_MCMCglmm"
if (!dir.exists(dir_output)) dir.create(dir_output)

models <- foreach(i = 1:length(trees), 
        .packages = c("MCMCglmm", "phytools", "ape"),
        .verbose = FALSE) %dopar% {
  
  tree_i <- trees[[i]]
  
  tree_i <- keep.tip(tree_i, subdat$Species)

  tree_i <- force.ultrametric(tree_i, message = FALSE)

  inv_phylo <- inverseA(tree_i)$Ainv

  model_run <- try({
    MCMCglmm(fixed = best_model,
             random = ~ Species,
             family = "categorical",
             data = subdat,
             ginverse = list(Species = inv_phylo),
             prior = prior_mcmc,
             nitt = nitt_val,
             burnin = burnin_val,
             thin = thin_val,
             verbose = FALSE)
  }, silent = TRUE)
  
  if(inherits(model_run, "try-error")){
      return(NULL)
    } else {
      saveRDS(model_run, file = paste0(dir_output, "/model_", i, ".rds"))
      return(paste0(dir_output, "/model_", i, ".rds"))
    }
  }

stopCluster(cl)
print("the end")

## imputed
dir_output <- "topologies_imputed_MCMCglmm"
if (!dir.exists(dir_output)) dir.create(dir_output)

models_imp <- foreach(i = 1:length(trees), 
        .packages = c("MCMCglmm", "phytools", "ape"),
        .verbose = TRUE) %dopar% {
  
  tree_i <- trees[[i]]
  subdat_i <- subdat_rep_imp[[i]]
  
  tree_i <- keep.tip(tree_i, subdat_i$Species)

  tree_i <- force.ultrametric(tree_i, message = FALSE)

  inv_phylo <- inverseA(tree_i)$Ainv

  model_run <- try({
    MCMCglmm(fixed = best_model_imp,
             random = ~ Species,
             family = "categorical",
             data = subdat_i,
             ginverse = list(Species = inv_phylo),
             prior = prior_mcmc,
             nitt = nitt_val,
             burnin = burnin_val,
             thin = thin_val,
             verbose = F)
  }, silent = F)
  
  if(inherits(model_run, "try-error")){
      return(NULL)
    } else {
      saveRDS(model_run, file = paste0(dir_output, "/model_", i, ".rds"))
      return(paste0(dir_output, "/model_", i, ".rds"))
    }
  }

stopCluster(cl)
print("the end")

all_models <- list.files("topologies_MCMCglmm", full.names = T, pattern = "\\.rds$")
all_models_imp <- list.files("topologies_imputed_MCMCglmm", full.names = T, pattern = "\\.rds$")

list_Sol <- lapply(all_models, function(arq) {
  mod <- readRDS(arq)
  return(mod$Sol)
})

list_Sol_imp <- lapply(all_models_imp, function(arq) {
  try({mod <- readRDS(arq)
  return(mod$Sol)})
})

list_VCV <- lapply(all_models, function(arq) {
  mod <- readRDS(arq)
  return(mod$VCV)
})

list_VCV_imp <- lapply(all_models_imp, function(arq) {
  try({mod <- readRDS(arq)
  return(mod$VCV)})
})

all_Sol <- do.call(rbind, list_Sol)
all_VCV <- do.call(rbind, list_VCV)

all_Sol_imp <- do.call(rbind, list_Sol_imp)
all_VCV_imp <- do.call(rbind, list_VCV_imp)

summary_fixed <- data.frame(
  Post.mean = colMeans(all_Sol),
  Lower_95_CI = apply(all_Sol, 2, quantile, probs = 0.025),
  Upper_95_CI = apply(all_Sol, 2, quantile, probs = 0.975),
  pMCMC = apply(all_Sol, 2, function(x) {
    if (mean(x) > 0) 2 * (sum(x < 0) / length(x)) else 2 * (sum(x > 0) / length(x))
  })
)

summary_fixed_imp <- data.frame(
  Post.mean = colMeans(all_Sol_imp),
  Lower_95_CI = apply(all_Sol_imp, 2, quantile, probs = 0.025),
  Upper_95_CI = apply(all_Sol_imp, 2, quantile, probs = 0.975),
  pMCMC = apply(all_Sol_imp, 2, function(x) {
    if (mean(x) > 0) 2 * (sum(x < 0) / length(x)) else 2 * (sum(x > 0) / length(x))
  })
)

summary_random <- data.frame(
  Post.mean = colMeans(all_VCV),
  Lower_95_CI = apply(all_VCV, 2, quantile, probs = 0.025),
  Upper_95_CI = apply(all_VCV, 2, quantile, probs = 0.975)
)

summary_random_imp <- data.frame(
  Post.mean = colMeans(all_VCV_imp),
  Lower_95_CI = apply(all_VCV_imp, 2, quantile, probs = 0.025),
  Upper_95_CI = apply(all_VCV_imp, 2, quantile, probs = 0.975)
)

## figure 2

p_range <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = Range_size_log, 
                              fill = Population_trend_binary)) +  
  geom_violin(color = NA) +
  theme_classic() + 
  ylab("Log (Range size (km²))") + 
  guides(fill = FALSE) +
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank()) +
  labs(x=NULL) + 
  geom_boxplot(width=.1,  outliers = FALSE) +
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)

p_SVL <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = SVL_log, 
                            fill = Population_trend_binary)) +
  geom_violin(color = NA) +
  theme_classic() + 
  ylab("Log (BS (mm))") + 
  guides(fill = FALSE) +
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  labs(x=NULL)+geom_boxplot(width=.1,  outliers = FALSE) +
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)

p_prev_bio1 <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, 
                                y = prev_bio1, 
                                fill = Population_trend_binary)) +  
 geom_violin(color = NA) + 
 theme_classic() + 
 ylab("AMT Prevalence") + 
 guides(fill = FALSE)+   
 scale_fill_manual(values = c("#ecc05f", "#8c4733")) +   
 theme(axis.text.y = element_text(angle = 90, vjust = 0.5)) + 
 xlab("Population trend")+
 geom_boxplot(width=.1,  outliers = FALSE)+
 stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)
  
p_prev_moisture <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, 
                                    y = prev_Moisture_scale, 
                                    fill = Population_trend_binary)) +  
  geom_violin(color = NA, scale = "width") +   
  theme_classic() + 
  coord_cartesian(ylim = c(0, 3)) + 
  ylab("CMI Prevalence") + 
  guides(fill = "none")+   
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) +   
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5))+ 
  xlab("Population trend")+
  geom_boxplot(width=.1,  outliers = FALSE)+
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)

p_bio1 <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = bio1_scale, 
                           fill = Population_trend_binary)) +  
  geom_violin(color = NA) +   
  theme_classic() + 
  ylab("AMN (ºC)") + 
  guides(fill = FALSE)+   
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) +   theme(axis.text.y = element_text(angle = 90, vjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ 
  labs(x=NULL)+
  geom_boxplot(width=.1,  outliers = FALSE)+
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)

p_bio7 <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = bio7,
                           fill = Population_trend_binary)) +  
  geom_violin(color = NA) +   
  theme_classic() + 
  ylab("TAR (ºC)") + 
  guides(fill = FALSE)+   
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) +   
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())+ 
  labs(x=NULL)+
  geom_boxplot(width=.1,  outliers = FALSE)+
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)

p_moisture <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = Moisture_scale, 
                               fill = Population_trend_binary)) +  
  geom_violin(color = NA) +   
  theme_classic() + 
  ylab("CMI (kg m−2 month−1)") + 
  guides(fill = FALSE)+   
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) + 
  labs(x=NULL)+
  geom_boxplot(width=.1,  outliers = FALSE)+
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)+  
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())
 
p_lat <- ggplot(subset(subdat_imp, !is.na(Population_trend_binary)), aes(x = Population_trend_binary, y = Latitude_scale, fill = Population_trend_binary)) +  
  geom_violin(color = NA) +   
  theme_classic() + 
  ylab("Absolute latitude midpoint") + 
  guides(fill = FALSE)+   
  scale_fill_manual(values = c("#ecc05f", "#8c4733")) + 
  labs(x=NULL)+
  geom_boxplot(width=.1,  outliers = FALSE)+
  stat_summary(geom = "Point",  fun.y = "mean", col="white", cex=2)+  
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5), axis.text.x=element_blank(), axis.ticks.x=element_blank())


pdf("Fig2_rev1.pdf", width = 8, height = 10)
ggarrange(p_SVL, p_range, 
          p_bio1, p_moisture,
          p_bio7, p_lat,
          p_prev_bio1, p_prev_moisture,
          ncol = 2, nrow = 4)
dev.off()

## figure 3

term_labels <- c(
  "Range_size_log" = "Range Size (log)",
  "bio1_scale" = "AMT",
  "prev_bio1_scale" = "AMT Prevalence",
  "bio7" = "TAR",
  "Moisture_scale" = "CMI",
  "prev_Moisture_scale" = "CMI Prevalence",
  "Latitude_scale" = "Latitude"
)

list_mod <- lapply(all_models, function(arq) {
  mod <- readRDS(arq) 
  stats <- data.frame(
    Term = colnames(mod$Sol),
    Estimate = colMeans(mod$Sol),
    Lower = apply(mod$Sol, 2, quantile, probs = 0.025),
    Upper = apply(mod$Sol, 2, quantile, probs = 0.975),
    row.names = NULL
  )
  stats <- stats[!stats$Term == "(Intercept)", ]
  return(stats)
})

list_mod_imp <- lapply(all_models_imp, function(arq) {
  mod <- readRDS(arq) 
  stats <- data.frame(
    Term = colnames(mod$Sol),
    Estimate = colMeans(mod$Sol),
    Lower = apply(mod$Sol, 2, quantile, probs = 0.025),
    Upper = apply(mod$Sol, 2, quantile, probs = 0.975),
    row.names = NULL
  )
  stats <- stats[!stats$Term == "(Intercept)", ]
  return(stats)
})

p_mod <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

p_mod_imp <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

for (i in 1:100) {
  p_mod <- p_mod + 
    geom_pointrange(data = list_mod[[i]], aes(x = Estimate, y = Term, 
                                              xmin = Lower, xmax = Upper),
                    position = position_dodge(width = 0.5),
                    alpha = 0.1,
                    size = 0.8, fatten = 3) 
  p_mod_imp <- p_mod_imp + 
   geom_pointrange(data = list_mod_imp[[i]], aes(x = Estimate, y = Term, 
                                             xmin = Lower, xmax = Upper),
                   position = position_dodge(width = 0.5),
                   alpha = 0.1,
                   size = 0.8, fatten = 3) 
}

p_mod <- p_mod + 
  scale_y_discrete(labels = term_labels) +
  labs(y = NULL, x = "Estimate", title = "A. Non-imputed") + 
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    plot.title = element_text(face = "bold", size = 14)
  )

p_mod_imp <- p_mod_imp + 
  scale_y_discrete(labels = term_labels) +
  labs(y = NULL, x = "Estimate", title = "B. Imputed") + 
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    plot.title = element_text(face = "bold", size = 14)
  )

pdf("Fig3_rev1.pdf", width = 12, height = 5)
ggarrange(p_mod, p_mod_imp, ncol = 2)
dev.off()
