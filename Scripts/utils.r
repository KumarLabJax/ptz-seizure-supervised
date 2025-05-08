freeze_features <- c("freeze_time_sec", "freeze_behavior_proportion", "freeze_num_bouts", "freeze_avg_bout_lens", "freeze_first_bout") #"freeze_bout_3sec", "freeze_bout_6sec", "freeze_bout_9sec"
circle_features <- c("circle_num_bouts")
of_features <- c("f2_distance_cm", "f2_center_distance_cm", "f2_periphery_distance_cm", "f2_corner_distance_cm", "f2_grooming_number_bouts", "f2_grooming_duration_secs") #"f2_center_time_secs", "f2_periphery_time_secs", "f2_corner_time_secs",
straub_tail_features <- c("straub_tail_time_sec", "straub_tail_behavior_proportion", "straub_tail_num_bouts", "straub_tail_avg_bout_lens", "straub_tail_first_bout")
leg_splaying_features <- c("leg_splaying_time_sec", "leg_splaying_behavior_proportion", "leg_splaying_num_bouts", "leg_splaying_avg_bout_lens", "leg_splaying_first_bout")
side_seizure_features <- c("side_seizure_time_sec", "side_seizure_behavior_proportion", "side_seizure_num_bouts", "side_seizure_avg_bout_lens", "side_seizure_first_bout")
wild_jumping_features <- c("wild_jumping_time_sec", "wild_jumping_behavior_proportion", "wild_jumping_num_bouts", "wild_jumping_avg_bout_lens", "wild_jumping_first_bout")

feature_set <- c(straub_tail_features, leg_splaying_features, side_seizure_features, wild_jumping_features, freeze_features, circle_features, of_features)

omit_videos <- c("007737_B6NJ_M_60mg-kg_MAYBE_WRONG_ID DO NOT USE.avi", "007728-2_B6NJ_M_60mg-kg_maybe part 2 top of hour.avi", "007678_B6J_M_60mg-kg_.avi", "007719 B6NJ M 40 mg-kg NV6_2022-03-22_11-02-26.avi", "007737_B6NJ_M_vehicle_.avi", "007712_B6NJ_F_80mg-kg_.avi") #lacking good qc, racine score