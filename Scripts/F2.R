libs <- c("dplyr", "readr", "ggplot2", "tidyr", "circlize", "ComplexHeatmap", "correlation", "patchwork")
sapply(libs, require, character.only = TRUE)

source("Scripts/utils.R")

# read in the data
starts <- read_csv("Data/starts.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"), duration = (endframe - startframe)/(30*60)) |> filter(!NetworkFilename %in% omit_videos)
bins <- read_csv("Data/PTZ_bins.csv") 
events <- read_csv("Data/PTZ_events.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
events <- events |> add_row(NetworkFilename = "007666_B6J_F_vehicle_", behavior_key = "Circling", start = 3363, duration = 75, MouseID = "007666") |> add_row(NetworkFilename = "007688_B6J_M_80mg-kg_", behavior_key = "Circling", start = 1003, duration = 75, MouseID = "007668") |> add_row(NetworkFilename = "007710_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 433, duration = 75, MouseID = "007710") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1696, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1799, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 1914, duration = 75, MouseID = "007707") |> add_row(NetworkFilename = "007707_B6NJ_F_80mg-kg_", behavior_key = "Circling", start = 2685, duration = 75, MouseID = "007707") #adding circling bouts information to the events data frame

features <- read_csv("Data/PTZ_features.csv") |> mutate(MouseID = stringr::str_extract(NetworkFilename, "[0-9]{6}"))
meta <- read_csv("Data/meta.csv")

master_df <- read_csv("Data/masterdf.csv") |> mutate_at("Dosereal", factor)
master_df <- master_df |> drop_na(Score) 

master_df <- master_df |> dplyr::slice(1:85)

starts <- starts |> select(MouseID, startframe, endframe, duration) |> left_join(meta |> select(MouseID, NetworkFilename, highscore, NetworkFilename_corrected, Dosereal, Strain, Sex, Score, seizure, AM, MS, AvgScore), by = "MouseID")
starts <- starts |> rename(start = startframe) |> mutate(minute = duration, behavior_key = "End")

starts <- starts |> select(!endframe)
events <- events |> inner_join(meta |> select(!NetworkFilename), by = "MouseID") |> mutate(minute = start/(30*60)) |> mutate_at("Dosereal", factor)
events <- events |> mutate(behavior_key = case_when(behavior_key == "Freeze_001" ~ "Freeze", behavior_key == "Tail_jerk" ~ "Straub Tail", behavior_key == "Leg_splaying" ~ "Leg Splaying", behavior_key == "Side_seizure" ~ "Side Seizure", behavior_key == "Wild_jumping" ~ "Wild Jumping", behavior_key == "Circling" ~ "Circling", TRUE ~ behavior_key))
events <- rbind(events, starts)
events <- events |> filter(!Dosereal == 20)

events <- events |> group_by(MouseID) |> mutate(dif = minute - minute[behavior_key == "End"]) |> filter(dif <= 0) |> ungroup() |> select(!dif)

p1 <- events |> filter(!behavior_key == "Freeze_2cm") |> filter(minute <= 20) |> ggplot(aes(x = minute, y = MouseID, color = behavior_key)) + geom_point() + theme_classic(base_size = 14) + facet_wrap(.~Dosereal, ncol = 1, scales = "free_y") + labs(x = "Bin", y = "MouseID", color = 'Behavior') + scale_color_brewer(palette = "Set1") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_color_manual(values = c("Straub Tail" = "#1B9E77", "Leg Splaying" = "#D95F02", "Side Seizure" = "#7570B3", "Wild Jumping" = "#E7298A", "Freeze" = "#577590", "Circling" = "#E6AB02", "End" = "#000000"))

p2 <- events |> filter(!behavior_key %in% c("Freeze_2cm", "End")) |> filter(minute <= 20) |> ggplot(aes(x = minute, y = behavior_key, color = behavior_key)) + geom_jitter(size = 1) + theme_classic(base_size = 14) + facet_wrap(.~Dosereal, ncol = 1, scales = "free_y") + labs(x = "Bin", y = "MouseID", color = 'Behavior') + scale_color_brewer(palette = "Set1") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + scale_color_manual(values = c("Straub Tail" = "#1B9E77", "Leg Splaying" = "#D95F02", "Side Seizure" = "#7570B3", "Wild Jumping" = "#E7298A", "Freeze" = "#577590", "Circling" = "#E6AB02", "End" = "#000000"))

p1 + theme(legend.position = "none") + labs(x = "Time (mins)") | p2 + labs(y = NULL) + theme(legend.position = "none") + labs(x = "Time (mins)")

ggsave("fig2.pdf", width = 9.23, height = 11.33)