# Agnieszka Goroncy, 2026
# IV. Final stage and summary

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Selection of media for the report: 5 most distant and similar to healthy face
media_rank <- dist_summary %>% filter(MediaID != 1) %>% arrange(desc(mean_dist_total))
top_far_media <- media_rank %>% slice_head(n = 5) %>% pull(MediaID)
top_close_media <- media_rank %>% arrange(mean_dist_total) %>% slice_head(n=5) %>% pull(MediaID)
selected_media <- c(top_far_media, top_close_media)

# Collective ranking of selected media
report_table_main <- dist_summary %>% filter(MediaID %in% selected_media) %>% arrange(desc(mean_dist_total)) %>%
  mutate(media_role=case_when(
      MediaID %in% top_far_media ~ "most different from healthy face",
      MediaID %in% top_close_media ~ "most similar to healthy face",
      TRUE ~ "selected")) %>%
  select(MediaID, media_role, mean_dist_occ, mean_dist_centers, mean_dist_trans, mean_dist_total, median_dist_total, sd_dist_total, n)

# The biggest shifts in state centers
report_table_centers <- center_shift_summary %>% filter(MediaID %in% selected_media) %>% arrange(MediaID, desc(mean_shift)) %>% group_by(MediaID) %>%
  slice_head(n=5) %>% ungroup() %>% select(MediaID, clean_state, n, mean_dx, mean_dy, mean_shift)

# Largest changes in state shares
report_table_occ <- occ_shift_summary %>% filter(MediaID %in% selected_media, n>=10) %>% mutate(abs_delta = abs(mean_delta_occ)) %>% arrange(MediaID, desc(abs_delta)) %>%
  group_by(MediaID) %>% slice_head(n=6) %>% ungroup() %>% select(MediaID, ref_state, type, n, mean_occ_clean, mean_occ_media, mean_delta_occ)
 
# Biggest transition changes
report_table_trans <- transition_shift_summary %>% filter(MediaID %in% selected_media) %>%
  mutate(abs_delta = abs(mean_delta_trans)) %>% arrange(MediaID, desc(abs_delta)) %>% group_by(MediaID) %>%
  slice_head(n=8) %>% ungroup() %>% select(MediaID, clean_from, clean_to, mean_trans_clean, mean_trans_media, mean_delta_trans)

report_table_trans
