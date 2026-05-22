# Agnieszka Goroncy
# III. In what direction do the HMM features change with respect to MediaID=1?
# We don't assume that State1 in medium = State1 in healthy face, but just first match states between models of the same person,
# and then count shifts and changes
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

max_states <- 5
occ_threshold <- 0.05 # "active" state if it has >= 5% observations
unmatched_penalty <- 0.35 # penalty for not matching state
occ_weight <- 0.20 #weight of the difference in the shares of states in the matching cost

# Extract states from the first line of hmm_features
extract_states_from_row <- function(row_df, suffix="") {
  bind_rows(lapply(1:max_states, function(s) {
    data.frame(state=s,
      mean_x=as.numeric(row_df[[paste0("State", s, "_mean_x", suffix)]][1]),
      mean_y=as.numeric(row_df[[paste0("State", s, "_mean_y", suffix)]][1]),
      occ=as.numeric(row_df[[paste0("Occ_State", s, suffix)]][1])
    )
  }))
}

all_perms <- function(x) {
  if (length(x)==1) return(matrix(x, nrow=1))
  out <- do.call(rbind, lapply(seq_along(x), function(i) {
    rest <- x[-i]
    cbind(x[i], all_perms(rest))
  }))
  out
}

# Optimal fit healthy <-> media
# Returns:
# - matched: matched pairs of states
# - clean_only: healthy (clean) state with no media match
# - media_only: media state with no healthy match

optimal_state_match <- function(clean_states, media_states, occ_threshold=0.05,
				unmatched_penalty=0.35, occ_weight=0.20) {
  clean_act <- clean_states %>%
    filter(!is.na(mean_x), !is.na(mean_y), occ >= occ_threshold)
  media_act <- media_states %>%
    filter(!is.na(mean_x), !is.na(mean_y), occ >= occ_threshold)
  nc <- nrow(clean_act)
  nm <- nrow(media_act)
  if (nc==0 && nm==0) {
    return(list(matched=tibble(),clean_only = tibble(),
      media_only = tibble()))
  }
  if (nc==0) {
    return(list(matched=tibble(), clean_only=tibble(),
      media_only=media_act %>% mutate(type="media_only")
    ))
  }
  if (nm==0) {
    return(list(matched=tibble(), clean_only=clean_act %>% mutate(type = "clean_only"), media_only = tibble()
    ))
  }
  # cost of matching clean_i -> media_j
  cost_mat <- outer(seq_len(nc), seq_len(nm), 
	Vectorize(function(i, j) {
    d_center <- sqrt((clean_act$mean_x[i]-media_act$mean_x[j])^2 +
                     (clean_act$mean_y[i]-media_act$mean_y[j])^2)
    d_occ <- abs(clean_act$occ[i]-media_act$occ[j])
    d_center+occ_weight*d_occ
  }))
  # extension to square matrix via dummy rows/cols  
  n <- max(nc, nm)
  big_cost <- matrix(0, nrow=n, ncol=n)
  big_cost[1:nc, 1:nm] <- cost_mat
  if (nc<n) big_cost[(nc+1):n, 1:nm] <- unmatched_penalty
  if (nm<n) big_cost[1:nc, (nm+1):n] <- unmatched_penalty
  perms <- all_perms(seq_len(n))
  perm_cost <- apply(perms, 1, function(p) sum(big_cost[cbind(seq_len(n), p)]))
  best_perm <- perms[which.min(perm_cost), ]  
  matched <- list()
  clean_only <- list()
  media_used <- integer(0)  
  for (i in seq_len(nc)) {
    j <- best_perm[i]
    if (j<=nm) {
      media_used <- c(media_used, j)
      matched[[length(matched) + 1]] <- tibble(clean_state=clean_act$state[i],
        media_state=media_act$state[j], clean_x=clean_act$mean_x[i],
        clean_y=clean_act$mean_y[i], media_x=media_act$mean_x[j],
        media_y=media_act$mean_y[j], occ_clean=clean_act$occ[i],
        occ_media=media_act$occ[j], match_cost=cost_mat[i, j],
        type="matched"
      )
    } else {
      clean_only[[length(clean_only) + 1]] <- tibble(clean_state=clean_act$state[i],
        media_state=NA_integer_, clean_x=clean_act$mean_x[i],
        clean_y=clean_act$mean_y[i], media_x=NA_real_,
        media_y=NA_real_, occ_clean=clean_act$occ[i],
        occ_media=0, match_cost=unmatched_penalty, type="clean_only"
      )
    }
  }  
  media_only_idx <- setdiff(seq_len(nm), media_used)
  media_only <- lapply(media_only_idx, function(j) {
    tibble(clean_state=NA_integer_, media_state=media_act$state[j],
      clean_x=NA_real_,clean_y = NA_real_,
      media_x=media_act$mean_x[j], media_y=media_act$mean_y[j],
      occ_clean=0, occ_media=media_act$occ[j],
      match_cost=unmatched_penalty, type="media_only"
    )
  })  
  list(matched = bind_rows(matched), clean_only = bind_rows(clean_only),
    media_only = bind_rows(media_only))
}

# Matching states for each participant and media pair
# Healthy (clean) reference for each person
clean_ref <- hmm_features %>% filter(MediaID==1) %>%
  select(ParticipantID, all_of(center_cols), starts_with("Occ_State")) %>%
  rename_with(~paste0(.x, "_clean"), -ParticipantID)

# Matching medium and healthy (clean) face of the same person
paired_hmm <- hmm_features %>% filter(MediaID!=1) %>%
  left_join(clean_ref, by="ParticipantID")

# State maching
state_match_list <- lapply(1:nrow(paired_hmm), function(i) {
  row_i <- paired_hmm[i,]
  clean_states <- extract_states_from_row(row_i, suffix="_clean")
  media_states <- extract_states_from_row(row_i, suffix="")
  mt <- optimal_state_match(clean_states=clean_states,
    media_states=media_states, occ_threshold=occ_threshold,
    unmatched_penalty=unmatched_penalty, occ_weight=occ_weight)
  bind_rows(mt$matched, mt$clean_only, mt$media_only) %>%
    mutate(ParticipantID=row_i$ParticipantID, MediaID=row_i$MediaID,
      n_states=row_i$n_states)
})

state_match_df <- bind_rows(state_match_list)
 
# State shifts after matching
center_shift_summary <- state_match_df %>% filter(type=="matched") %>%
  mutate(dx=media_x-clean_x, dy=media_y-clean_y, shift=sqrt(dx^2+dy^2)) %>%
  group_by(MediaID, clean_state) %>%
  summarise(n=n(), clean_x=mean(clean_x, na.rm=TRUE), clean_y = mean(clean_y, na.rm=TRUE),
    media_x = mean(media_x, na.rm=TRUE), media_y = mean(media_y, na.rm=TRUE),
    mean_dx = mean(dx, na.rm=TRUE), mean_dy = mean(dy, na.rm=TRUE),
    mean_shift = mean(shift, na.rm=TRUE),.groups = "drop")

# Changes in states share after matching
# Here we also take into account disappearing and new states:
# clean_only -> the state was in healthy image, disappeared in lesioned medium
# media_only -> a new state appeared in lesioned medium
occ_shift_df <- state_match_df %>% mutate(ref_state=case_when(type=="matched"~clean_state, 
	type=="clean_only"~clean_state, type=="media_only"~media_state),
    delta_occ = occ_media - occ_clean)

occ_shift_summary <- occ_shift_df %>% group_by(MediaID, ref_state, type) %>%
  summarise(n=n(), mean_occ_clean=mean(occ_clean, na.rm=TRUE),
    mean_occ_media=mean(occ_media, na.rm=TRUE),
    mean_delta_occ=mean(delta_occ, na.rm=TRUE), .groups = "drop")

# state transitions after matching
# Function: transition delta after matching states
compute_transition_delta_matched <- function(row_df, match_sub, max_states=5) {  
  matched <- match_sub %>% filter(type=="matched") %>% select(clean_state, media_state)  
  if (nrow(matched)<2) return(NULL)
  out <- list()
  for (i in seq_len(nrow(matched))) {
    for (j in seq_len(nrow(matched))) {
      c_from <- matched$clean_state[i]
      c_to <- matched$clean_state[j]
      m_from <- matched$media_state[i]
      m_to <- matched$media_state[j]
      clean_col <- paste0("T_S", c_from, "_S", c_to, "_clean")
      media_col <- paste0("T_S", m_from, "_S", m_to)
      t_clean <- as.numeric(row_df[[clean_col]][1])
      t_media <- as.numeric(row_df[[media_col]][1])
      out[[length(out) + 1]] <- tibble(clean_from = c_from, clean_to=c_to,
        media_from=m_from, media_to=m_to, trans_clean=t_clean,
        trans_media=t_media, delta_trans=t_media-t_clean)
    }
  }  
  bind_rows(out)
}

# Append the clean transition matrix to paired_hmm
clean_ref_trans <- hmm_features %>% filter(MediaID==1) %>%
  select(ParticipantID, all_of(trans_cols)) %>%
  rename_with(~ paste0(.x, "_clean"), -ParticipantID)

paired_hmm_trans <- hmm_features %>% filter(MediaID!=1) %>%
  left_join(clean_ref_trans, by="ParticipantID")

# Counting delta in transitions
transition_delta_list <- lapply(1:nrow(paired_hmm_trans), 
	function(i) {
  	row_i <- paired_hmm_trans[i,] 
  	match_sub <- state_match_df %>% filter(ParticipantID==row_i$ParticipantID, MediaID==row_i$MediaID)
  	td <- compute_transition_delta_matched(row_i, match_sub, max_states=max_states)
  	if (is.null(td)) return(NULL)
  	td %>% mutate(ParticipantID=row_i$ParticipantID, MediaID=row_i$MediaID)
	}
)

transition_delta_df <- bind_rows(transition_delta_list)

# Transitions summary
transition_shift_summary <- transition_delta_df %>%
  group_by(MediaID, clean_from, clean_to) %>% summarise(
    mean_trans_clean=mean(trans_clean, na.rm=TRUE),
    mean_trans_media=mean(trans_media, na.rm=TRUE),
    mean_delta_trans=mean(delta_trans, na.rm=TRUE),
    mean_abs_delta=mean(abs(delta_trans), na.rm=TRUE),.groups = "drop")

# biggest shift in centers
top_center_shifts <- center_shift_summary %>% arrange(desc(mean_shift)) %>%
  select(MediaID, clean_state, n, mean_dx, mean_dy, mean_shift)

# Biggest changes in state share
top_occ_shifts <- occ_shift_summary %>% mutate(abs_delta=abs(mean_delta_occ)) %>%
  arrange(desc(abs_delta)) %>% select(MediaID, ref_state, type, n, mean_occ_clean, mean_occ_media, mean_delta_occ)

# Biggest state transition changes
top_transitions <- transition_shift_summary %>% mutate(abs_delta=abs(mean_delta_trans)) %>%
  arrange(desc(abs_delta)) %>% select(MediaID, clean_from, clean_to,
         mean_trans_clean, mean_trans_media, mean_delta_trans)
