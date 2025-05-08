libs <- c('dplyr', 'readr', 'ggplot2', 'tidyr', 'brms')
sapply(libs, require, character.only = TRUE)

source("Scripts/utils.R")

# read in the data
starts <- read_csv("Data/starts.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) |> filter(!NetworkFilename %in% omit_videos)
bins <- read_csv("Data/PTZ_bins.csv")
events <- read_csv("Data/PTZ_events.csv")
features <- read_csv("Data/PTZ_features.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
master_df <- read_csv("Data/masterdf.csv") 
master_df <- master_df |> rename(ordAM = seizure)

master_df <- master_df |> mutate(ordMS = case_when(MS == -1 ~ "No", MS <= 2 ~ "Low", MS <= 4 ~ "Medium", TRUE ~ "High")) |> mutate(ordAvg = case_when(AvgScore < 0 ~ "No", AvgScore <= 2 ~ "Low", AvgScore <= 4 ~ "Medium", TRUE ~ "High")) 

master_df$ordAM <- factor(master_df$ordAM, levels = c("No", "Low", "Medium", "High"))
master_df$ordMS <- factor(master_df$ordMS, levels = c("No", "Low", "Medium", "High"))
master_df$ordAvg <- factor(master_df$ordAvg, levels = c("No", "Low", "Medium", "High"))

master_df <- master_df |> mutate(numAM = factor(as.numeric(ordAM), ordered = TRUE), numMS = factor(as.numeric(ordMS), ordered = TRUE), numAvg = factor(as.numeric(ordAvg), ordered = TRUE))

#Use starts to identify the duration of the video 
#This will be used to calculate the proportion of behavior time
starts <- starts |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) |> relocate(MouseID, .after = NetworkFilename)
starts <- starts |> mutate(duration = (endframe - startframe)/(30*60)) #duration in minutes

#Time-normalized features
master_df <- master_df |> inner_join(starts |> select(MouseID, duration), by = "MouseID")
master_df <- master_df |> mutate_at(vars(all_of(feature_set)), list(~ ./duration)) |> select(!c(duration))

ptz <- master_df
ptz <- ptz |> mutate_if(is.numeric, ~ (replace_na(., 0)))
ptz <- ptz |> filter(!Dosereal == 20)

#LDA with PC (Seizure as response variable)
#Figure 3A
X <- ptz |> dplyr::select(any_of(feature_set)) |> as.matrix()
X <- scale(X, center = TRUE, scale = TRUE)
svd_X <- svd(X)
evalues <- svd_X$d^2/(nrow(X) - 1)
X_pc <- svd_X$u%*%diag(svd_X$d)

tmp <- MASS::lda(ordAvg ~ ., data = ptz |> dplyr::select(ordAvg) |> bind_cols(as.data.frame(X_pc[,1:31]))) #34
tmp2 <- predict(tmp)

df_tmp <- ptz |> dplyr::select(ordAvg, Strain, Sex) |> mutate(LD1 = tmp2$x[,1], LD2 = tmp2$x[,2]) 
p1 <- ggplot(df_tmp, aes(x = LD1, y = LD2, color = ordAvg)) + geom_point(size = 4, alpha = 0.80) + theme_classic(base_size = 16) + labs(x = latex2exp::TeX("$LD_1$"), y = latex2exp::TeX("$LD_2$")) + scale_color_brewer(palette = "Set1") + guides(color=guide_legend(title="Seizure Intensity"))

p1
ggsave("lda_ordAvg.pdf", width = 5.2, height = 3.8)

#LDA with PC (Dose as response variable)
#Figure 3B
X <- ptz |> dplyr::select(any_of(feature_set)) |> as.matrix()
X <- scale(X, center = TRUE, scale = TRUE)
svd_X <- svd(X)
evalues <- svd_X$d^2/(nrow(X) - 1)
X_pc <- svd_X$u%*%diag(svd_X$d)

tmp <- MASS::lda(Dosereal ~ ., data = ptz |> dplyr::select(Dosereal) |> bind_cols(as.data.frame(X_pc[,1:31])))
tmp2 <- predict(tmp)

df_tmp <- ptz |> dplyr::select(Dosereal) |> mutate(LD1 = tmp2$x[,1], LD2 = tmp2$x[,2]) 

p2 <- ggplot(df_tmp, aes(x = LD1, y = LD2, color = as.factor(Dosereal))) + geom_point(size = 4, alpha = 0.85) + theme_classic(base_size = 16) + labs(x = latex2exp::TeX("$LD_1$"), y = latex2exp::TeX("$LD_2$")) + scale_color_brewer(palette = "Set1", name = "Dose (mg/kg)")

p2
ggsave("lda_dose.pdf", width = 5.2, height = 3.8)


#Ordinal regression model
ptz <- ptz |> filter(!Dosereal == 20)
X <- ptz |> dplyr::select(c(all_of(feature_set))) |> as.matrix()
X <- scale(X, center = TRUE, scale = TRUE)
y <- ptz |> dplyr::select(ordAvg) |> pull() |> factor(ordered = TRUE)

set.seed(123)
fitcv <- ordinalNet::ordinalNetCV(X, y, family = "cumulative", link = "logit", lambdaMinRatio = 1e-02, maxiterOut = 5000, nFolds = 10, nFoldsCV = 5)

theta_hat <- coef(fitcv$fit)[-c(1:3)] #Plot P1 with -theta_hat for correct direction or better interpretability
tau <- coef(fitcv$fit)[c(1:3)]
df_theta <- data.frame(Features = feature_set, theta = theta_hat[setdiff(names(theta_hat), c("strain", "sex"))])

#Figure 3C
df_theta <- data.frame(Feature = feature_set, theta_hat = theta_hat[setdiff(names(theta_hat), c("strain", "sex"))])
df_theta$Feature <- factor(df_theta$Feature, levels = feature_set)
p2 <- ggplot(df_theta, aes(x = Feature, y = -theta_hat)) + geom_bar(stat = "identity", color = TRUE) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(y = latex2exp::TeX("Coefficients, $\\theta$"))
p2
ggsave("theta_hat_ordAvg.pdf", width = 6.3, height = 4.9)

colMeans(summary(fitcv))

#Figure 3D
df_univariate <- ptz |> dplyr::select(ordAvg) |> mutate(score = as.numeric(X%*%theta_hat))
p1 <- ggplot(df_univariate, aes(x = -score, color = ordAvg)) + geom_density(lwd = 1.1) + scale_color_brewer(palette = "Set1", name = "Seizure") + labs(x = "Univariate Seizure Scale", y = "Density") + scale_y_continuous(expand = c(0,0)) + geom_vline(xintercept = tau[1], linetype = "dashed") + geom_vline(xintercept = tau[2], linetype = "dashed") + geom_vline(xintercept = tau[3], linetype = "dashed") + theme_classic(base_size = 16)
p1 
ggsave("univariate_seizure_scale_enet_ordAvg.pdf", width = 5.2, height = 3.8)

#highscore ('baseline' experiment) (Figure S1F)
ptz <- ptz |> filter(!Dosereal == 20)
X <- ptz |> dplyr::select(c(all_of(feature_set))) |> as.matrix()
X <- scale(X, center = TRUE, scale = TRUE)
y <- ptz |> dplyr::select(highscore) |> pull() |> factor(ordered = TRUE)

set.seed(123)
fitcv <- ordinalNet::ordinalNetCV(X, y, family = "cumulative", link = "logit", lambdaMinRatio = 1e-02, maxiterOut = 5000, nFolds = 10, nFoldsCV = 5)

theta_hat <- coef(fitcv$fit)[-c(1:5)] #Plot P1 with -theta_hat for correct direction or better interpretability
tau <- coef(fitcv$fit)[c(1:5)]
df_theta <- data.frame(Features = feature_set, theta = theta_hat[setdiff(names(theta_hat), c("strain", "sex"))])

#Figure 3C
df_theta <- data.frame(Feature = feature_set, theta_hat = theta_hat[setdiff(names(theta_hat), c("strain", "sex"))])
df_theta$Feature <- factor(df_theta$Feature, levels = feature_set)
p2 <- ggplot(df_theta, aes(x = Feature, y = -theta_hat)) + geom_bar(stat = "identity", color = TRUE) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(y = latex2exp::TeX("Coefficients, $\\theta$"))
p2
ggsave("theta_hat_highscore.pdf", width = 6.3, height = 4.9)

df_univariate <- ptz |> dplyr::select(highscore) |> mutate(score = as.numeric(X%*%theta_hat)) |> mutate(highscore = factor(highscore, ordered = TRUE))
p1 <- ggplot(df_univariate, aes(x = -score, color = highscore)) + geom_density(lwd = 1.1) + ggsci::scale_color_jama(name = "Seizure") + labs(x = "Univariate Seizure Scale", y = "Density") + scale_y_continuous(expand = c(0,0)) + geom_vline(xintercept = tau[1], linetype = "dashed") + geom_vline(xintercept = tau[2], linetype = "dashed") + geom_vline(xintercept = tau[3], linetype = "dashed") + geom_vline(xintercept = tau[4], linetype = "dashed") + geom_vline(xintercept = tau[5], linetype = "dashed") + theme_classic(base_size = 16) + guides(color=guide_legend(title="Racine Scale"))
p1 
ggsave("univariate_seizure_scale_enet_highscore.pdf", width = 14.2, height = 6.6)

p1 | p2
ggsave("univariate_seizure_scale_enet_highscore.pdf", width = 14.2, height = 6.6) #p2 comes from below after changing ordAvg to highscore and knowing a bit of details about the model i.e. coefficients (theta_hat) and tau

#Feature importance/weights across folds

nfolds <- length(fitcv$folds)
set <- 1:nrow(ptz)
df_theta_folds <- data.frame(Features = feature_set)

for (fold in 1:nfolds){

    cat("Fold: ", fold, "\n")
    test <- fitcv$folds[[fold]]
    train <- setdiff(set, test)
    
    Xtrain <- ptz |> dplyr::select(c(all_of(feature_set))) |> filter(row_number() %in% train) |> as.matrix()
    Xtrain <- scale(Xtrain, scale = TRUE, center = TRUE)
    ytrain <- ptz |> dplyr::select(highscore) |> filter(row_number() %in% train) |> pull() |> factor(ordered = TRUE)

    Xtest <- ptz |> dplyr::select(c(all_of(feature_set))) |> filter(row_number() %in% test) |> as.matrix()
    ytest <- ptz |> dplyr::select(highscore) |> filter(row_number() %in% test) |> pull() |> factor(ordered = TRUE)

    fit <- ordinalNet::ordinalNet(Xtrain, ytrain, family = "cumulative", link = "logit", maxiterOut = 5000, standardize = TRUE)

    theta_hat <- coef(fit)[-c(1:5)]
    tau <- coef(fit)[c(1:5)]
    df_theta_folds <- df_theta_folds |> left_join(data.frame(Features = feature_set, fold = theta_hat[setdiff(names(theta_hat), c("strain", "sex"))]), by = "Features")

}

names(df_theta_folds) <- c("Feature", paste0("Fold", 1:nfolds))

df_theta <- df_theta_folds |> tibble() |> rowwise() |> mutate(m = mean(c_across(Fold1:Fold10)), sd = sd(c_across(Fold1:Fold10)), se = sd/sqrt(nfolds)) |> select(Feature, m, sd, se)

df_theta$Feature <- factor(df_theta$Feature, levels = feature_set)
p2 <- ggplot(df_theta, aes(x = Feature, y = -m)) + geom_bar(stat = "identity", color = TRUE) + geom_errorbar(aes(x = Feature, ymin = -m - se, ymax = -m + se), width = 0.2) + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(y = latex2exp::TeX("Coefficients, $\\hat{\\theta}$")) 
p2
ggsave("theta_hat_ordAvg_folds.pdf", width = 6.3, height = 4.9)

ggsave("highscore_results_univariate_feature_importance.pdf", width = 7.4, height = 7.0)


#Interplay between noise and misclassification

X <- ptz |> filter(!Dosereal == 20) |> select(all_of(feature_set)) |> as.matrix()
mu <- colMeans(X) 

Xpca <- prcomp(X)
Xhat <- Xpca$x %*% t(Xpca$rotation) 
Xhat <- scale(Xhat, center = -mu, scale = FALSE)


ff <- function(eps, seed){

    #Gaussian
    Xpc <- Xpca$x + matrix(rnorm(nrow(Xpca$x)*ncol(Xpca$x), 0, eps), nrow = nrow(Xpca$x), ncol = ncol(Xpca$x))
    Xhatn <- Xpc %*% t(Xpca$rotation) 
    Xhatn <- scale(Xhatn, center = -mu, scale = FALSE)
    X <- scale(Xhatn, center = TRUE, scale = TRUE)
    y <- ptz |> filter(!Dosereal == 20) |> dplyr::select(ordAvg) |> pull()
    set.seed(seed)
    fitcv <- ordinalNet::ordinalNetCV(X, y, family = "cumulative", link = "probit", lambdaMinRatio = 1e-02, maxiterOut = 5000, nFolds = 10, nFoldsCV = 5)
    m1 <- colMeans(summary(fitcv))[["misclass"]]
    s1 <- summary(fitcv)[["misclass"]] |> sd()
    coef1 <- coef(fitcv$fit)
    
    
    #Uniform
    Xpc <- Xpca$x + matrix(runif(nrow(Xpca$x)*ncol(Xpca$x), -eps, eps), nrow = nrow(Xpca$x), ncol = ncol(Xpca$x))
    Xhatn <- Xpc %*% t(Xpca$rotation)
    Xhatn <- scale(Xhatn, center = -mu, scale = FALSE)
    X <- scale(Xhatn, center = TRUE, scale = TRUE)
    y <- ptz |> filter(!Dosereal == 20) |> dplyr::select(ordAvg) |> pull()
    set.seed(seed)
    fitcv <- ordinalNet::ordinalNetCV(X, y, family = "cumulative", link = "probit", lambdaMinRatio = 1e-02, maxiterOut = 5000, nFolds = 10, nFoldsCV = 5)

    m2 <- colMeans(summary(fitcv))[["misclass"]]
    s2 <- summary(fitcv)[["misclass"]] |> sd()
    coef2 <- coef(fitcv$fit)

    return(list(err = matrix(c(m1, s1, m2, s2), nrow = 2, ncol = 2, byrow = TRUE), coef1 = coef1, coef2 = coef2))

}

require(future)
err <- as.list(seq(0, 0.5, by = 0.1))
plan(multicore, workers = availableWorkers() |> length())

res_list <- furrr::future_map(err, ~ff(eps = .x, seed = 1121), .options=furrr::furrr_options(seed = TRUE))

dfres <- do.call(rbind, lapply(res_list, `[[`, 1))
dfres <- dfres |> as_tibble() |> mutate(Distribution = rep(c("Gaussian", "Uniform"), length(err)) |> as.factor()) |> rename("cvmean" = V1, "cvsd" = V2) |> mutate(Eps = rep(err, each = 2) |> as.numeric()) |> mutate(cvse = cvsd/sqrt(10))

dfres |> ggplot(aes(x = Eps, y = cvmean, color = Distribution)) + geom_line() + geom_point(size = 2) + theme_classic(base_size = 16) + labs(x = latex2exp::TeX("$Noise (\\sigma_\\epsilon$)"), y = "Misclassification Rate") + scale_color_brewer(palette = "Set1") + guides(color=guide_legend(title="Error Distribution")) + geom_hline(yintercept = 0.75, linetype = 'dashed', color = 'black')
ggsave("noise_validation4.pdf", width = 5.3, height = 3.6)