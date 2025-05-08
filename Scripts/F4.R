libs <- c("dplyr", "readr", "ggplot2", "tidyr", "stringr", "ordinal", "patchwork")
sapply(libs, require, character.only = TRUE)

source("Scripts/utils.R")


# read in the data
starts <- read_csv("Data/starts.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) |> filter(!NetworkFilename %in% omit_videos)
bins <- read_csv("Data/PTZ_bins.csv") |> rename_all(make.names) |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
events <- read_csv("Data/PTZ_events.csv") |> mutate(MouseID = str_extract(NetworkFilename, "[0-9]{6}")) 
features <- read_csv("PTZ_features.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
meta <- read_csv("meta.csv") 
circling <- read_csv("Data/circling_minutebins.csv", col_types = (NetworkFilename = "f")) |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) |> rename_all(make.names) |> relocate(MouseID, .before = "X1") |> mutate(behavior = "Circling")

df_correction <- read_csv("Data/nfilenames_corrected.csv") |> rename("NetworkFilename" = "network_filename", "NetworkFilename_corrected" = "network_filename_corrected") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) #We extract MouseID from Nfilename since starts has Nfilenames
df_correction_20 <- starts |> select(NetworkFilename, MouseID) |> dplyr::slice(86:109) |> mutate(NetworkFilename_corrected = NetworkFilename)
df_correction <- rbind(df_correction, df_correction_20)

circling <- circling |> left_join(df_correction |> select(!NetworkFilename), by = "MouseID") 
circling <- circling |> select(!NetworkFilename) |> rename("NetworkFilename" = "NetworkFilename_corrected") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) 

master_df <- read_csv("Data/masterdf.csv") |> mutate_at("Dosereal", factor)
master_df <- master_df |> drop_na(Score)

meta <- meta |> rename(ordAM = seizure)
meta <- meta |> mutate(ordMS = case_when(MS == -1 ~ "No", MS <= 2 ~ "Low", MS <= 4 ~ "Medium", TRUE ~ "High")) |> mutate(ordAvg = case_when(AvgScore < 0 ~ "No", AvgScore <= 2 ~ "Low", AvgScore <= 4 ~ "Medium", TRUE ~ "High")) 
meta$ordAM <- factor(meta$ordAM, levels = c("No", "Low", "Medium", "High"), ordered = TRUE)
meta$ordMS <- factor(meta$ordMS, levels = c("No", "Low", "Medium", "High"), ordered = TRUE)
meta$ordAvg <- factor(meta$ordAvg, levels = c("No", "Low", "Medium", "High"), ordered = TRUE)
meta <- meta |> mutate(numAM = factor(as.numeric(ordAM), ordered = TRUE), numMS = factor(as.numeric(ordMS), ordered = TRUE), numAvg = factor(as.numeric(ordAvg), ordered = TRUE))

bins <- rbind(bins, circling)
bins_melt <- bins |> filter(!behavior == "Freeze_2cm") |> pivot_longer(names_to = "bin", values_to = "count", cols = contains("X")) |> rename("minute" = "bin") |> mutate(minute = gsub("X", "", minute) |> as.numeric()) |> left_join(meta |> select(!c(NetworkFilename, highscore)), by = "MouseID")

df_ptz <- bins_melt |> pivot_wider(names_from = "behavior", values_from = "count") |> rename("freezing" = "Freeze_001", "leg_splaying" = "Leg_splaying", "side_seizure" = "Side_seizure", "tail_jerk" = "Tail_jerk", "wild_jump" = "Wild_jumping", "circling" = "Circling") 
df_ptz[, c("circling", "freezing", "leg_splaying", "side_seizure", "tail_jerk", "wild_jump")] <- scale(df_ptz[, c("circling", "freezing","leg_splaying", "side_seizure", "tail_jerk", "wild_jump")], center = TRUE, scale = TRUE)
df_ptz$Strain <- as.double(ifelse(df_ptz$Strain == "B6J", 0, 1))
df_ptz$Sex <- as.double(ifelse(df_ptz$Sex == "F", 0, 1))
df_ptz$minute <- factor(df_ptz$minute)

df_ptz <- df_ptz |> filter(!Dosereal == 20)
ptz <- master_df
ptz <- ptz |> filter(!Dosereal == 20)
ptz <- ptz |> filter(!MouseID == "007741")

#Implementing LOOCV in ordinal linear mixed model
mIDs <- ptz |> select(MouseID) |> unique() |> pull()
df_test <- df_ptz |> select(MouseID, ordAvg, AvgScore) |> distinct(MouseID, ordAvg, AvgScore) |> mutate(mean_seizure_intensity = NA, max_seizure_intensity = NA, median_seizure_intensity = NA) 

for (i in 1:length(mIDs)){
    cat("MouseID:", paste0(mIDs[i]), "\n")
    fit0 <- clmm(ordAvg ~ freezing + leg_splaying + side_seizure + tail_jerk + (1|minute), data = df_ptz |> filter(!MouseID %in% mIDs[i]))
    beta_hat <- coefficients(fit0)[-c(1:3)]
    X <- as.matrix(df_ptz |> filter(MouseID == mIDs[i]) |> select(all_of(c("freezing", "leg_splaying", "side_seizure", "tail_jerk"))))
    df_test[df_test$MouseID %in% mIDs[i], "mean_seizure_intensity"] <- mean(X%*%beta_hat + ranef(fit0)$minute[,1], na.rm = TRUE)
    df_test[df_test$MouseID %in% mIDs[i], "max_seizure_intensity"] <- max(X%*%beta_hat + ranef(fit0)$minute[,1], na.rm = TRUE)  
    df_test[df_test$MouseID %in% mIDs[i], "median_seizure_intensity"] <- median(X%*%beta_hat + ranef(fit0)$minute[,1], na.rm = TRUE) 
}

#Figure 4C
p1 <- ggplot(df_test |> filter(mean_seizure_intensity < 20), aes(x = ordAvg, y = mean_seizure_intensity, fill = ordAvg)) + geom_boxplot() + geom_point() + theme_classic(base_size = 16) + labs(x = "Racine Group", y = latex2exp::TeX("$\\hat{y}_{Mean}$")) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none")

p2 <- ggplot(df_test, aes(x = ordAvg, y = max_seizure_intensity, fill = ordAvg)) + geom_boxplot() + geom_point() + theme_classic(base_size = 16) + labs(x = "Racine Group", y = latex2exp::TeX("$\\hat{y}_{Max}$")) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none")

p3 <- ggplot(df_test, aes(x = ordAvg, y = median_seizure_intensity, fill = ordAvg)) + geom_boxplot() + geom_point() + theme_classic(base_size = 16) + labs(x = "Racine Group", y = latex2exp::TeX("$\\hat{y}_{Median}$")) + scale_fill_brewer(palette = "Set1") + theme(legend.position = "none")

(p3 + theme(legend.position = "none") + labs(x = NULL)) / (p1 + theme(legend.position = "none") + labs(x = NULL)) / (p2 + theme(legend.position = "none"))

ggsave("LOOCV_ordAvg.pdf", width = 3.12, height = 7.6)


#Figure 4B
fit <- clmm(ordAvg ~ freezing + leg_splaying + side_seizure + tail_jerk + (1|minute), data = df_ptz)
X <- as.matrix(df_ptz |> select(all_of(c("freezing", "leg_splaying", "side_seizure", "tail_jerk"))))
beta_hat <- coefficients(fit)[-c(1:3)]
df_ptz <- df_ptz |> mutate(seizure_intensity = as.numeric(X%*%beta_hat) + ranef(fit)$minute[,1])

p1 <- ggplot(df_ptz, aes(x = minute, y = seizure_intensity, group = interaction(MouseID, ordAvg))) + geom_step(alpha = 0.2, aes(color = ordAvg)) + theme_classic(base_size = 16) + scale_color_brewer(palette = "Set1") + stat_summary(aes(y = seizure_intensity, group = ordAvg, color = ordAvg), fun = max, geom = "step", size = 1.5, na.rm = TRUE) + labs(y = latex2exp::TeX("$\\hat{y}_{Max}$"), x = "Time Bins (minute)", color = "Racine Group") + theme(legend.position = "top") + facet_grid(cols = vars(ordAvg)) + scale_x_discrete(breaks = seq(0,20,5))

p2 <- ggplot(df_ptz, aes(x = minute, y = seizure_intensity, group = interaction(MouseID, ordAvg))) + geom_step(alpha = 0.2, aes(color = ordAvg)) + theme_classic(base_size = 16) + scale_color_brewer(palette = "Set1") + stat_summary(aes(y = seizure_intensity, group = ordAvg, color = ordAvg), fun = mean, geom = "step", size = 1.5, na.rm = TRUE) + labs(y = latex2exp::TeX("$\\hat{y}_{Mean}$"), x = "Time Bins (minute)", color = "Racine Group") + theme(legend.position = "top") + facet_grid(cols = vars(ordAvg)) + scale_x_discrete(breaks = seq(0,20,5))

p3 <- ggplot(df_ptz, aes(x = minute, y = seizure_intensity, group = interaction(MouseID, ordAvg))) + geom_step(alpha = 0.2, aes(color = ordAvg)) + theme_classic(base_size = 16) + scale_color_brewer(palette = "Set1") + stat_summary(aes(y = seizure_intensity, group = ordAvg, color = ordAvg), fun = median, geom = "step", size = 1.5, na.rm = TRUE) + labs(y = latex2exp::TeX("$\\hat{y}_{Median}$"), x = "Time Bins (minute)", color = "Racine Group") + theme(legend.position = "top") + facet_grid(cols = vars(ordAvg)) + scale_x_discrete(breaks = seq(0,20,5), labels = seq(0,20,5))

(p3 + theme(legend.position = "none") + labs(x = NULL)) / (p2 + theme(legend.position = "none") + labs(x = NULL)) / (p1 + theme(legend.position = "none"))

ggsave("seizure_intensity_time_series_ordAvg.pdf", width = 7, height = 7)

#Qualititave validation 
#Select representative animals from each group and overlay annotated behaviors on top of predicted behavior profiles
df_tmp <- df_test |> arrange(desc(max_seizure_intensity)) |> distinct(ordAvg, .keep_all = TRUE)  
mIDs <- df_tmp$MouseID
ggplot(df_ptz |> filter(MouseID %in% mIDs), aes(x = minute, y = seizure_intensity, group = interaction(MouseID, ordAvg))) + geom_step(alpha = 0.2, aes(color = ordAvg)) + theme_classic(base_size = 16) + scale_color_brewer(palette = "Set1") + stat_summary(aes(y = seizure_intensity, group = ordAvg, color = ordAvg), fun = max, geom = "step", size = 1.5, na.rm = TRUE) + labs(y = latex2exp::TeX("$\\hat{y}_{Max}$"), x = "Time Bins (minute)", color = "Racine Group") + facet_grid(rows = vars(ordAvg)) + scale_x_discrete(breaks = seq(0,20,5)) + theme(legend.position = "none")
dev.print(pdf, "../plots/qualitative_validation1_ordAvg.pdf", width = 7.4, height = 3)



# read in the data
starts <- read_csv("Data/starts.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"), duration = (endframe - startframe)/(30*60)) |> filter(!NetworkFilename %in% omit_videos)
bins <- read_csv("Data/PTZ_bins.csv") 
events <- read_csv("Data/PTZ_events.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
events <- events |> add_row(NetworkFilename = "007666_B6J_F_vehicle_", behavior_key = "Circling", start = 3363, duration = 75, MouseID = "007666") |> add_row(NetworkFilename = "007688_B6J_M_80mg-kg_", behavior_key = "Circling", start = 1003, duration = 75, MouseID = "007668") |> add_row(NetworkFilename = "007710_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 433, duration = 75, MouseID = "007710") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1696, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1799, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1914, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 2685, duration = 75, MouseID = "007707")

events <- events |> inner_join(meta |> select(!NetworkFilename), by = "MouseID") |> mutate(minute = start/(30*60)) |> mutate_at("Dosereal", factor)
events <- events |> mutate(behavior_key = case_when(behavior_key == "Freeze_001" ~ "Freeze", behavior_key == "Tail_jerk" ~ "Straub Tail", behavior_key == "Leg_splaying" ~ "Leg Splaying", behavior_key == "Side_seizure" ~ "Side Seizure", behavior_key == "Wild_jumping" ~ "Wild Jumping", behavior_key == "Circling" ~ "Circling", TRUE ~ behavior_key))

df_master <- events |> filter(MouseID %in% mIDs) |> rename("behavior" = "behavior_key") |> filter(!behavior %in% "Freeze_2cm") 


p1 <- ggplot(df_ptz |> filter(MouseID %in% mIDs), aes(x = minute, y = seizure_intensity, group = interaction(ordAvg, MouseID))) + geom_step(alpha = 0.2, aes(color = ordAvg)) + theme_classic(base_size = 16) + stat_summary(aes(y = seizure_intensity, group = ordAvg, color = ordAvg), fun = max, geom = "step", size = 1.5, na.rm = TRUE) + labs(y = latex2exp::TeX("$\\hat{y}_{Max}$"), x = "Time Bins (minute)", color = "Racine Group") + facet_grid(cols = vars(factor(ordAvg, levels = c("No", "Low", "Medium", "High")))) + scale_x_discrete(breaks = seq(0,20,5)) + theme(legend.position = "none") + geom_jitter(data = df_master |> filter(minute < 20), aes(x = minute, y = 14, color = behavior), height = 1.5, size = 2, alpha = 0.5) + scale_color_manual(values = c("Straub Tail" = "#1B9E77", "Leg Splaying" = "#D95F02", "Side Seizure" = "#7570B3", "Wild Jumping" = "#E7298A", "Freeze" = "#577590", "Circling" = "#E6AB02", "end" = "#000000")) 

ggsave("qualitative_val3_ordAvg.pdf", height = 2.4, width = 7)


p2 <- df_master |> group_by(ordAvg, behavior) |> mutate(Count = n(), behavior = factor(behavior)) |> mutate(behavior = factor(behavior, levels = c("Circling", "Freeze", "Leg Splaying", "Side Seizure", "Straub Tail", "Wild Jumping"))) |> ggplot(aes(y = as.factor(behavior), x = ordAvg, fill = Count)) + geom_tile(col = "white") + scale_fill_viridis_b() + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x = "Seizure Intensity", y = NULL) 
# + geom_text(aes(label = count)) 

ggsave("qualitative_val4_ordAvg.pdf", width = 3.90, height = 3.1)
