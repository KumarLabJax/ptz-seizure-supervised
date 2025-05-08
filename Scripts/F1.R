libs <- c("dplyr", "readr", "ggplot2", "tidyr", "circlize", "ComplexHeatmap", "correlation")
sapply(libs, require, character.only = TRUE)

source("Scripts/utils.R")

# read in the data
starts <- read_csv("Data/starts.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}")) |> filter(!NetworkFilename %in% omit_videos)
bins <- read_csv("Data/PTZ_bins.csv")
events <- read_csv("Data/PTZ_events.csv")
features <- read_csv("Data/PTZ_features.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))

master_df <- read_csv("masterdf.csv") 
master_df <- master_df |> drop_na(Score)
master_df <- master_df |> mutate_at("Dosereal", factor)
master_df <- master_df |> slice(1:85)

#Figure 1
#Figure 1B
ggplot(master_df, aes(x = Strain, fill = Sex, alpha = Dosereal)) + geom_bar() + theme_classic(base_size = 16) + scale_fill_brewer(palette = "Set1", name = "Sex") + labs(alpha = "PTZ Dose", x = "Strain", y = "Count") + scale_y_continuous(expand=c(0,0))
ggsave("counts.pdf", width = 3.2, height = 3.6)

#Figure 1C
ggplot(master_df |> mutate(highscore = as.numeric(as.character(AvgScore))), aes(x = Dosereal, y = AvgScore, fill = Sex)) + geom_violin() + theme_classic(base_size = 16) + scale_fill_brewer(palette = "Set1") + labs(x = "PTZ Dose (mg/kg)", y = "High Score") + facet_wrap(.~Strain, scales = "free_y") + scale_x_discrete(breaks = c("0", "20", "40", "60", "80"), labels = c("0" = "0", "20" = "20", "40" = "40", "60" = "60", "80" = "80")) + expand_limits(y = c(-1, 8)) + guides(fill = guide_legend(title = "Sex"))
ggsave("avgscore_dose.pdf", width = 5.6, height = 2.7)

#Figure 1D
levels(master_df$highscore) <- c(levels(master_df$highscore), "0", "2", "3")
master_df$highscore <- factor(master_df$highscore, levels = as.character(seq(-1,7,1)))
ggplot(master_df, aes(x = highscore, fill = factor(seizure, levels = c("No", "Low", "Medium", "High")))) + geom_bar(stat = "count") + theme_classic(base_size = 16) + labs(x = "High Score", y = "Count") + scale_x_discrete(drop = FALSE) + scale_y_continuous(expand = c(0,0)) + scale_fill_brewer(palette = "Set1") + guides(fill = guide_legend(title = NULL)) + theme(legend.position = "top")
ggsave("highscore.pdf", width = 2.6, height = 2.6)
ggsave("highscore.pdf", width = 4.2, height = 4.2)

#Figure 1E
ptz <- master_df |> janitor::clean_names() |> mutate(seizure_numeric = case_when(seizure == "No" ~ 0, seizure == "Low" ~ 1, seizure == "Medium" ~ 2, seizure == "High" ~ 3)) 

ptz <- ptz |> mutate_if(is.numeric, ~ (replace_na(., 0)))


features <- c(straub_tail_features, leg_splaying_features, side_seizure_features, wild_jumping_features, freeze_features, circle_features, of_features)

ptz$dosereal <- factor(ptz$dosereal, ordered = TRUE)
ptz$dosereal <- as.numeric(as.character(ptz$dosereal))
corr_all <- correlation::correlation(ptz |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) 

ptz_M <- ptz |> filter(sex == "M")
ptz_F <- ptz |> filter(sex == "F")

corr_M <- correlation::correlation(ptz_M |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rM = r, pM = p)
corr_F <- correlation::correlation(ptz_F |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rF = r, pF = p)
df_corr <- inner_join(corr_M, corr_F, by = c("Parameter2", "Type"))

ptz_b6 <- ptz |> filter(strain %in% "B6J")
ptz_b6nj <- ptz |> filter(strain %in% "B6NJ")

corr_b6 <- correlation::correlation(ptz_b6 |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6 = r, pb6 = p)

corr_b6nj <- correlation::correlation(ptz_b6nj |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6nj = r, pb6nj = p)

df_corr <- inner_join(corr_b6, corr_b6nj, by = c("Parameter2", "Type"))

ptz_b6_M <- ptz |> filter(strain %in% "B6J" & sex %in% "M")
ptz_b6_F <- ptz |> filter(strain %in% "B6J" & sex %in% "F")

corr_b6_M <- correlation::correlation(ptz_b6_M |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6 = r, pb6 = p)

corr_b6_F <- correlation::correlation(ptz_b6_F |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6 = r, pb6 = p)

ptz_b6nj_M <- ptz |> filter(strain %in% "B6NJ" & sex %in% "M")
ptz_b6nj_F <- ptz |> filter(strain %in% "B6NJ" & sex %in% "F")

corr_b6nj_M <- correlation::correlation(ptz_b6nj_M |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6nj = r, pb6nj = p)

corr_b6nj_F <- correlation::correlation(ptz_b6nj_F |> select(all_of(c("seizure_numeric", features))), method = "pearson") |> filter(Parameter1 == "seizure_numeric") |> mutate(Type = case_when(Parameter2 %in% straub_tail_features ~ "Straub Tail", Parameter2 %in% leg_splaying_features ~ "Leg Splaying", Parameter2 %in% side_seizure_features ~ "Side Seizure", Parameter2 %in% wild_jumping_features ~ "Wild Jumping", Parameter2 %in% freeze_features ~ "Freeze", Parameter2 %in% circle_features ~ "Circle", Parameter2 %in% of_features ~ "Open Field")) |> select(all_of(c("Parameter2", "r", "p", "Type"))) |> rename(rb6nj = r, pb6nj = p)

#Show all correlations vs only significant ones
show_corr <- "all" # "all" or "significant"
if (show_corr == "all"){
    corr_all <- corr_all |> mutate(tau = ifelse(p > 0.05, 0, r)) |> select(Parameter2, r) |> rename(Features = Parameter2, All = r)
    corr_M <- corr_M |> mutate(tau = ifelse(pM > 0.05, 0, rM)) |> select(rM) |> rename(Males = rM)
    corr_F <- corr_F |> mutate(tau = ifelse(pF > 0.05, 0, rF)) |> select(rF) |> rename(Females = rF)
    corr_b6 <- corr_b6 |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(rb6) |> rename(B6J = rb6)
    corr_b6nj <- corr_b6nj |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(rb6nj) |> rename(B6NJ = rb6nj)
    corr_b6_M <- corr_b6_M |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(rb6) |> rename(B6J_M = rb6)
    corr_b6_F <- corr_b6_F |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(rb6) |> rename(B6J_F = rb6)
    corr_b6nj_M <- corr_b6nj_M |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(rb6nj) |> rename(B6NJ_M = rb6nj)
    corr_b6nj_F <- corr_b6nj_F |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(rb6nj) |> rename(B6NJ_F = rb6nj)
} else {
    corr_all <- corr_all |> mutate(tau = ifelse(p > 0.05, 0, r)) |> select(Parameter2, tau) |> rename(Features = Parameter2, All = tau)
    corr_M <- corr_M |> mutate(tau = ifelse(pM > 0.05, 0, rM)) |> select(tau) |> rename(Males = tau)
    corr_F <- corr_F |> mutate(tau = ifelse(pF > 0.05, 0, rF)) |> select(tau) |> rename(Females = tau)
    corr_b6 <- corr_b6 |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(tau) |> rename(B6J = tau)
    corr_b6nj <- corr_b6nj |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(tau) |> rename(B6NJ = tau)
    corr_b6_M <- corr_b6_M |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(tau) |> rename(B6J_M = tau)
    corr_b6_F <- corr_b6_F |> mutate(tau = ifelse(pb6 > 0.05, 0, rb6)) |> select(tau) |> rename(B6J_F = tau)
    corr_b6nj_M <- corr_b6nj_M |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(tau) |> rename(B6NJ_M = tau)
    corr_b6nj_F <- corr_b6nj_F |> mutate(tau = ifelse(pb6nj > 0.05, 0, rb6nj)) |> select(tau) |> rename(B6NJ_F = tau)
}

df_cor_sex <- cbind(corr_all, corr_b6, corr_b6nj, corr_M, corr_F, corr_b6_M, corr_b6_F, corr_b6nj_M, corr_b6nj_F)
rownames(df_cor_sex) <- df_cor_sex$Features
df_cor_sex <- df_cor_sex[,-1]
col_fun <- colorRamp2(c(min(df_cor_sex, na.rm = TRUE),0,0,max(df_cor_sex, na.rm = TRUE)),c("#045a8d","#f7f7f7","#f7f7f7","#cb181d"))
ht.z <- ComplexHeatmap::Heatmap(t(as.matrix(df_cor_sex)), row_names_gp = gpar(fontsize = 10),row_names_side = "left", column_names_gp = gpar(fontsize = 10), column_names_side = "bottom", heatmap_legend_param = list(at = c(-1,0,1),title = latex2exp::TeX("Pearson correlation"), title_position = "leftcenter-rot", border = "black",legend_height = unit(4, "cm"), labels_gp = gpar(fontsize = 14),just = c("right", "top")), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))})
ht.z

col_fun2 <- colorRamp2(c(1,2,3,4,5,6,7),c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", '#577590', '#ffcb69', '#797d62'), transparency = 0)
my_features <- data.frame(Type = ifelse(features %in% straub_tail_features, 1, ifelse(features %in% leg_splaying_features, 2, ifelse(features %in% side_seizure_features, 3, ifelse(features %in% wild_jumping_features, 4, ifelse(features %in% freeze_features, 5, ifelse(features %in% circle_features, 6, 7)))))))
rownames(my_features) <- features
colnames(my_features) <- ""

ht.cluster <- Heatmap(t(as.matrix(my_features)),cluster_rows = FALSE, cluster_columns = FALSE,col = col_fun2, show_heatmap_legend = FALSE,column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 10),row_names_side = "left")

ht.z %v% ht.cluster
dev.print(pdf, "correlation_heatmap_seizure_significant.pdf", height = 4.45, width = 7.55)


#Extracting legend for Figure 1E
my_features <- data.frame(Type = as.factor(c(rep("Straub Tail", 5), rep("Leg Splaying", 5), rep("Side Seizure", 5), rep("Wild Jumping", 5), rep("Freeze", 8), rep("Circle", 1), rep("Open Field", 9))))
levels(my_features$Type) <- c("Straub Tail", "Leg Splaying", "Side Seizure", "Wild Jumping", "Freeze", "Circle", "Open Field")
p <- ggplot(my_features, aes(x = Type, fill = factor(Type))) + geom_bar() + scale_fill_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", '#577590', '#ffcb69', '#797d62')) + guides(fill = guide_legend(title = "Feature Type")) + theme(legend.position = "top", legend.text = element_text(size = 12))
legend <- cowplot::get_legend(p)
grid.newpage()
grid.draw(legend)
ggsave("feature_types_legend.pdf", height = 1.3, width = 5, plot = legend)


#Figure 1F (Plots 1 and 2 need to be combined in Inkscape)
corrs <- ptz |> select(all_of(features)) |> correlation() |> summary(redundant = TRUE) 
corr_mat <- as.matrix(corrs |> select(!Parameter))
colnames(corr_mat) <- NULL

col_fun <- colorRamp2(c(min(corr_mat, na.rm = TRUE),0,0,max(corr_mat, na.rm = TRUE)),c("#045a8d","#f7f7f7","#f7f7f7","#cb181d"))
ht.z <- ComplexHeatmap::Heatmap(t(as.matrix(corr_mat)), row_names_gp = gpar(fontsize = 10),row_names_side = "left", column_names_gp = gpar(fontsize = 10), column_names_side = "bottom", heatmap_legend_param = list(at = c(-1,0,1),title = latex2exp::TeX("Pearson correlation"), title_position = "leftcenter-rot", border = "black",legend_height = unit(4, "cm"), labels_gp = gpar(fontsize = 14),just = c("right", "top")), col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,cell_fun = function(j, i, x, y, width, height, fill) {grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey"))})
ht.z

col_fun2 <- colorRamp2(c(1,2,3,4,5,6,7),c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", '#577590', '#ffcb69', '#797d62'), transparency = 0)
my_features <- data.frame(Type = ifelse(features %in% straub_tail_features, 1, ifelse(features %in% leg_splaying_features, 2, ifelse(features %in% side_seizure_features, 3, ifelse(features %in% wild_jumping_features, 4, ifelse(features %in% freeze_features, 5, ifelse(features %in% circle_features, 6, 7)))))))
colnames(my_features) <- ""

ht.cluster <- Heatmap(t(as.matrix(my_features)),cluster_rows = FALSE, cluster_columns = FALSE,col = col_fun2, show_heatmap_legend = FALSE,column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 10),row_names_side = "left")
ht.z %v% ht.cluster
dev.print(pdf, "correlation_heatmap.pdf", width = 3.9, height = 3.3) 

ht.cluster_row <- Heatmap((as.matrix(my_features)),cluster_rows = FALSE, cluster_columns = FALSE,col = col_fun2, show_heatmap_legend = FALSE,column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 10),row_names_side = "left")
dev.print(pdf, "correlation_heatmap_rows.pdf", width = 0.25, height = 3.3) #Plot2

dev.print(pdf, "feature_color_code.pdf", width = 6.3, height = 2.5)