# Agnieszka Goroncy
# II. ANALYSIS: What does the gaze pattern in the HMM look like?
# How much does the medium change the gaze pattern of the same person relative to the healthy face?

# For each ParticipantID x MediaID, we count three types of HMM features
# and compare them to the model of the same person for MediaID = 1 (healthy face)
# 1. State share - what percentage of observations fall into State1...State5
# 2. State location (ROI) - State1_mean_x, State1_mean_y, ..., State5_mean_x, State5_mean_y
# 3. Transition structure - the T_Si_Sj matrix

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

max_states <- 5

# State shares for each model: share of observations in each state
state_occ <- df_pelna %>% group_by(ParticipantID, MediaID, state) %>%
  summarise(n=n(), .groups="drop") %>% group_by(ParticipantID, MediaID) %>%
  mutate(prop=n/sum(n)) %>% ungroup() %>%
  select(ParticipantID, MediaID, state, prop) %>%
  tidyr::pivot_wider(names_from=state, values_from=prop, names_prefix="Occ_State",values_fill=0)

# make sure all columns 1 to 5 exist
for (s in 1:max_states) {
  nm <- paste0("Occ_State", s)
  if (!nm %in% names(state_occ)) state_occ[[nm]] <- 0
}

state_occ <- state_occ %>% select(ParticipantID, MediaID, paste0("Occ_State", 1:max_states))
head(state_occ)

# Building a single HMM feature table
# Selecting columns with center and transition features
center_cols <- unlist(lapply(1:max_states, function(s) {
	c(paste0("State", s, "_mean_x"),
	paste0("State", s, "_mean_y"))
}))

sd_cols <- unlist(lapply(1:max_states, function(s) {
  c(paste0("State", s, "_sd_x"),
    paste0("State", s, "_sd_y"))
}))

trans_cols <- as.vector(outer(
  1:max_states, 1:max_states,
  FUN = function(i, j) paste0("T_S", i, "_S", j)
))

hmm_features <- features_df %>%
  left_join(state_occ, by=c("ParticipantID", "MediaID")) %>%
  mutate(across(all_of(c(center_cols, sd_cols, trans_cols, paste0("Occ_State", 1:max_states))),~ replace_na(., 0)))

# Distance functions between model and MediaID = 1
occ_distance <- function(v1, v2) sqrt(sum((v1-v2)^2))

# Distance between states
# Important: we do not penalize missing a state as much as moving an existing state
# Therefore, we calculate the distance between centers only for states that have
# a non-zero contribution in at least one of the two models, and weight it by the average contribution
center_distance_weighted <- function(row1, row2, max_states=5) {
  dsum <- 0
  wsum <- 0
  for (s in 1:max_states) {
    x1 <- row1[[paste0("State", s, "_mean_x")]]
    y1 <- row1[[paste0("State", s, "_mean_y")]]
    x2 <- row2[[paste0("State", s, "_mean_x")]]
    y2 <- row2[[paste0("State", s, "_mean_y")]]
    p1 <- row1[[paste0("Occ_State", s)]]
    p2 <- row2[[paste0("Occ_State", s)]]
    w <- mean(c(p1, p2), na.rm=TRUE)
    if (is.na(w) || w==0) next

    d <- sqrt((x1-x2)^2+(y1-y2)^2)

    dsum <- dsum+w*d
    wsum <- wsum+w
  }
  if (wsum==0) return(NA_real_)
  dsum/wsum
}

# transition distance
transition_distance <- function(row1, row2, trans_cols) {
  v1 <- as.numeric(row1[trans_cols])
  v2 <- as.numeric(row2[trans_cols])
  sqrt(sum((v1-v2)^2))
}
# calculating the distance of each medium from media=1
# base models - MediaID=1
clean_models <- hmm_features %>% filter(MediaID==1) %>%
  select(ParticipantID, all_of(center_cols), all_of(trans_cols), starts_with("Occ_State")) %>%
  rename_with(~ paste0(.x, "_clean"), -ParticipantID)
hmm_compare <- hmm_features %>% left_join(clean_models, by="ParticipantID")
# auxiliary function: extracts a "clean row" for one record
get_clean_row <- function(row) {
  out <- list()
  for (nm in c(center_cols, trans_cols, paste0("Occ_State", 1:max_states))) {
    out[[nm]] <- row[[paste0(nm, "_clean")]]
  }
  as.list(out)
}

# distance calculation
dist_list <- lapply(seq_len(nrow(hmm_compare)), function(i) {
  row_i <- hmm_compare[i,]
  row_current <- as.list(row_i)
  row_clean <- get_clean_row(row_i)
  occ_cols <- paste0("Occ_State", 1:max_states)
  dist_occ <- occ_distance(
    unlist(row_current[occ_cols]),
    unlist(row_clean[occ_cols])
  )
  dist_centers <- center_distance_weighted(
    row1=row_current,
    row2=row_clean,
    max_states=max_states
  )
  dist_trans <- transition_distance(
    row1=row_current,
    row2=row_clean,
    trans_cols=trans_cols
  )
  data.frame(ParticipantID=row_i$ParticipantID, MediaID=row_i$MediaID, n_states=row_i$n_states,
    dist_occ=dist_occ, dist_centers=dist_centers, dist_trans=dist_trans)
})
dist_df <- bind_rows(dist_list)

# building the total distance
# standardization first - because we have different scales
dist_df <- dist_df %>% mutate(z_occ=as.numeric(scale(dist_occ)),
    z_centers=as.numeric(scale(dist_centers)),
    z_trans=as.numeric(scale(dist_trans)),
    dist_total=z_occ+z_centers+z_trans)

# Summary for MediaID:
dist_summary <- dist_df %>% group_by(MediaID) %>%
  summarise(mean_dist_occ=mean(dist_occ, na.rm=TRUE),
    mean_dist_centers=mean(dist_centers, na.rm=TRUE),
    mean_dist_trans=mean(dist_trans, na.rm=TRUE),
    mean_dist_total=mean(dist_total, na.rm=TRUE),
    median_dist_total=median(dist_total, na.rm=TRUE),
    sd_dist_total=sd(dist_total, na.rm=TRUE), n=n(), .groups="drop") %>%
  arrange(desc(mean_dist_total))

dist_summary %>% filter(MediaID!=1) %>% arrange(desc(mean_dist_total))

# Do media differ from healthy face in terms of HMM features?
library(lmerTest)
library(nlme)

dist_df_nonclean <- dist_df %>% filter(MediaID!=1)

mod_dist_total_nlme <- lme(fixed=dist_total~factor(MediaID), random=~1|ParticipantID,
  data=dist_df_nonclean, method="REML")

summary(mod_dist_total_nlme)
anova(mod_dist_total_nlme)

# Simple global test for media
mod_dist_null_nlme <- lme(fixed=dist_total~1, random=~1|ParticipantID, data=dist_df_nonclean, method = "ML")
mod_dist_full_nlme <- lme(fixed=dist_total~factor(MediaID), random=~1|ParticipantID, data=dist_df_nonclean, method="ML")
anova(mod_dist_null_nlme, mod_dist_full_nlme)

#Three analogous models separately:
# dist_occ ~ factor(MediaID)+(1|ParticipantID)
# do media differ in their attention distribution across states?
mod_occ <- lme(fixed=dist_occ~factor(MediaID), random=~1|ParticipantID, data=dist_df_nonclean,
  method="REML")

# dist_centers ~ factor(MediaID) + (1|ParticipantID)
# Do media differ in ROI location?
mod_centers <- lme(fixed=dist_centers~factor(MediaID), random=~1|ParticipantID,
	data=dist_df_nonclean,method="REML")

#dist_trans ~ factor(MediaID) + (1|ParticipantID)
# Do the media differ in their transition sequence?
mod_trans <- lme(fixed=dist_trans~factor(MediaID), random=~1|ParticipantID,
  data=dist_df_nonclean, method="REML")

anova(mod_occ)
anova(mod_centers)
anova(mod_trans)

mod_occ_null <- lme(fixed=dist_occ~1, random=~1|ParticipantID, data=dist_df_nonclean, method="ML")
mod_centers_null <- lme(fixed=dist_centers ~ 1, random=~1|ParticipantID,
  data=dist_df_nonclean, method="ML")
mod_trans_null <- lme(fixed=dist_trans~1, random=~1|ParticipantID, data=dist_df_nonclean, method="ML")

mod_occ_full <- update(mod_occ_null, fixed=dist_occ~factor(MediaID))
mod_centers_full <- update(mod_centers_null, fixed=dist_centers~factor(MediaID))
mod_trans_full <- update(mod_trans_null, fixed=dist_trans~factor(MediaID))

anova(mod_occ_null, mod_occ_full)
anova(mod_centers_null, mod_centers_full)
anova(mod_trans_null, mod_trans_full)
