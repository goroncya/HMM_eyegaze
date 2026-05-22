# Agnieszka Goroncy, 2026
# I. DOES THE DISTRIBUTION OF THE NUMBER OF ROIs DEPEND ON MEDIA ID?

# Individual participant analysis:
# does the number of states for a given medium
# differ from that for "clean face" (MediaID=1) for the same participant?
library(dplyr)
library(tidyr)
library(ordinal)

nstate_wide <- features_df %>% dplyr::select(ParticipantID, MediaID, n_states) %>%
  pivot_wider(names_from=MediaID, values_from=n_states, names_prefix="M")

# for each medium: calculate the difference from the healthy face 
nstate_diff <- features_df %>% dplyr::select(ParticipantID, MediaID, n_states) %>%
  left_join(features_df %>% filter(MediaID==1) %>% dplyr::select(ParticipantID, n_states_clean=n_states), by="ParticipantID") %>%
  mutate(delta_states=n_states-n_states_clean)

nstate_diff %>% group_by(MediaID) %>%
  summarise(mean_delta=mean(delta_states, na.rm=TRUE), median_delta=median(delta_states, na.rm=TRUE), sd_delta=sd(delta_states, na.rm=TRUE),
    n=n(), .groups="drop")

# Mixed model for n_states -tests whether MediaID affects the probability of transitioning to models with a larger number of states, taking into account
# that observations are repeated within a participant 
features_df$n_states_ord <- ordered(features_df$n_states)

# Cumulative Link Mixed Model
# baseline: media=1 (healthy face)
mod_states <- clmm(n_states_ord~factor(MediaID)+(1|ParticipantID), data=features_df)

# Likelihood ratio test
# null model: without MediaID
mod_states_null <- clmm(n_states_ord~1+(1|ParticipantID), data=features_df)
anova(mod_states, mod_states_null)

# Healthy vs lesion:
features_df <- features_df %>% mutate(media_type = ifelse(MediaID == 1, "healthy", "lesion"))
features_df$media_type <- factor(features_df$media_type, levels = c("healthy", "lesion"))

mod_healthy_vs_lesion <- clmm(n_states_ord~media_type+(1|ParticipantID),data=features_df)

# Global effect healthy vs lesion
mod_healthy_vs_lesion_null <- clmm(n_states_ord~1+(1|ParticipantID),data=features_df)
anova(mod_healthy_vs_lesion, mod_healthy_vs_lesion_null)

# Media and n_states
features_df %>% count(MediaID, n_states) %>%
  group_by(MediaID) %>% mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(MediaID), y = prop, fill = factor(n_states))) +
  geom_col(position = "fill") +
  labs(x = "MediaID", y = "Proportion", fill = "n_states") +
  theme_minimal()
