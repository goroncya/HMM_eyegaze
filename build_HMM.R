# Agnieszka Goroncy, 2026
library(dplyr)
library(depmixS4)	
str(df)	# data 

## Model settings
max_states <- 5	# max number of states (ROIs)
min_obs_total <- 8	# min number of observations for Participant and Media
min_sd_xy <- 0.01    # min standard deviation
multistart  <- TRUE
nstart_ms <- 10		
initIters_ms <- 20		
store_all_models <- FALSE  # TRUE = zapisuj wszystkie modele dla każdego k

## Auxiliary functions
# build_hmm_multiseq: builds HMM for multi-sequence data (every ModelID is a unique sequence)
# get_state_table: get state parameters
# get_state_mapping: HMM state ordering for comparison, from up to down (mean_y)
				# and left to right (mean_x) 
# extract_state_params: state parameters after ordering
# extract_transition_vector: transition matrix

build_hmm_multiseq <- function(df_pair, n_states, ntimes, use_multistart = multistart,nstart = 10, initIters = 10) {
	data1 <- df_pair %>% dplyr::select(x, y) 

	if (nrow(data1) < 2) stop("Not enough observations to build HMM.")
	if (any(!is.finite(data1$x)) || any(!is.finite(data1$y))) stop("Data contain NA/Inf values in x or y.")
  	if (sd(data1$x) == 0 || sd(data1$y) == 0) stop("Zero variance in x or y.")
 	if (sum(ntimes) != nrow(data1)) stop("sum(ntimes) is not equal to number of observations.")
 
	mod <- depmix(response=list(x~1, y~1), data=data1, nstates=n_states, family=list(gaussian(), gaussian()), ntimes=ntimes)
	if (multistart) {
    		fit_mod <- multistart(mod, nstart = nstart, initIters = initIters)
  	} else {
    	fit_mod <- fit(mod)
  	}	
	return(fit_mod)
}

get_state_table <- function(hmm_model) {
	n <- hmm_model@nstates
	out <- vector("list", n)
	for (s in seq_len(n)) {
		resp_x <- hmm_model@response[[s]][[1]]
		resp_y <- hmm_model@response[[s]][[2]]
		out[[s]] <- data.frame(state_old=s, mean_x=resp_x@parameters$coefficients[1], 
				sd_x=resp_x@parameters$sd[1], mean_y=resp_y@parameters$coefficients[1], sd_y=resp_y@parameters$sd[1])
 	}
	bind_rows(out)
}

get_state_mapping <- function(hmm_model) {
	st <- get_state_table(hmm_model)
	st_ord <- st %>% arrange(mean_y, mean_x) %>% mutate(state_new = row_number())
	map_old_to_new <- st_ord$state_new
	names(map_old_to_new) <- as.character(st_ord$state_old)
	list(state_table=st_ord, map_old_to_new=map_old_to_new)
}

extract_state_params <- function(hmm_model, max_states=5) {
	out_names <- c(paste0("State", rep(1:max_states, each=4), c("_mean_x", "_sd_x", "_mean_y", "_sd_y")))
	res <- setNames(rep(NA_real_, length(out_names)), out_names)
#  if (is.null(hmm_model) || !inherits(hmm_model, "depmix.fitted")) return(res)
	mapping <- get_state_mapping(hmm_model)
#  if (is.null(mapping)) return(res)
	st_ord <- mapping$state_table
	n <- nrow(st_ord)
	for (i in seq_len(min(n, max_states))) {
		base <- (i-1)*4
		res[base+1] <- st_ord$mean_x[i]
    		res[base+2] <- st_ord$sd_x[i]
    		res[base+3] <- st_ord$mean_y[i]
    		res[base+4] <- st_ord$sd_y[i]
  	}
  res
}

extract_transition_vector <- function(hmm_model, df_pair_sorted, max_states=5) {
	nm <- as.vector(outer(1:max_states, 1:max_states, FUN = function(i, j) paste0("T_S", i, "_S", j)))
  	res <- setNames(rep(NA_real_, length(nm)), nm)
	post <- posterior(hmm_model, type = "viterbi")
	mapping <- get_state_mapping(hmm_model)
	state_map <- mapping$map_old_to_new
	n <- hmm_model@nstates
	tmp <- df_pair_sorted %>% mutate(
		state_old=post$state, state=unname(state_map[as.character(state_old)]),
		seq_id=interaction(ParticipantID, MediaID, ModelID, drop=TRUE)) %>%
		group_by(seq_id) %>% mutate(next_state=lead(state)) %>% ungroup() %>% filter(!is.na(next_state))
	mat <- matrix(0, nrow=n, ncol=n)
	if (nrow(tmp) > 0) {
		counts <- tmp %>% count(state, next_state, name="n")
		for (r in seq_len(nrow(counts))) {
      		i <- counts$state[r]
      		j <- counts$next_state[r]
      		mat[i, j] <- counts$n[r]
	    	}
		rs <- rowSums(mat)
		for (i in seq_len(n)) {
      		if (rs[i]>0) mat[i,] <- mat[i,]/rs[i]
		}
  	}
	for (i in seq_len(n)) {
		for (j in seq_len(n)) {
			res[paste0("T_S", i, "_S", j)] <- mat[i,j]
    		}
  	}
  res
}

# prepare results objects:
model_list    <- list()   # all models for all k
best_models   <- list()   # the best models for each pair Participant and Media
features_list <- list()   # the best models features
fit_counter   <- setNames(rep(0, max_states - 1), as.character(2:max_states))

# List of pairs: Participant and Media
pairs <- data %>% distinct(ParticipantID, MediaID) %>% arrange(ParticipantID, MediaID)
dim(pairs)	

## Build one HMM per Participant and Media
                        
for (r in seq_len(nrow(pairs))) {
	p <- pairs$ParticipantID[r]
	m <- pairs$MediaID[r]
	cat(sprintf("\nparticipant=%s media=%s\n", p, m))
	df_sub <- data %>% filter(ParticipantID == p, MediaID == m) %>%
			dplyr::select(ParticipantID, MediaID, ModelID, TrialID, seq_id, x, y)
	# sequence length = number of observations for ModelID
	seq_info <- df_sub %>% count(ModelID, name = "seq_len") %>% arrange(ModelID)
	ntimes <- seq_info$seq_len
  	# check the quality for pairs Participant and Media
	if (nrow(df_sub) < min_obs_total || any(is.na(df_sub$x)) || any(is.na(df_sub$y)) ||
      	sd(df_sub$x, na.rm = TRUE) < min_sd_xy ||
      	sd(df_sub$y, na.rm = TRUE) < min_sd_xy) {
    		cat(" skipped: not enough data OR low variance OR NAs\n")
		next
	}
  	# we need at least 2 states
  	max_k_local <- min(max_states, nrow(df_sub) - 1)
  	if (max_k_local < 2) {
    		cat(" skipped: too few observations for >= 2 states\n")
    		next
  	}
	key <- paste(p, m, sep="_")
	ks <- 2:max_k_local
	aic_candidates <- setNames(rep(Inf, length(ks)), as.character(ks))
	model_candidates <- list()

	for (k in ks) {
    		kname <- as.character(k)
		cat(sprintf("trying k=%d\n", k))
		hmm_k <- tryCatch( {
      	build_hmm_multiseq(df_pair=df_sub, n_states=k, ntimes=ntimes, use_multistart=multistart, 
			nstart = nstart_ms, initIters = initIters_ms)},
    			error = function(e) { 
				cat(sprintf(" ! failed k=%d: %s\n", k, conditionMessage(e)))
 			return(NULL) })
		if (!is.null(hmm_k) && inherits(hmm_k, "depmix.fitted")) {
	      	aic_val <- AIC(hmm_k)
    			aic_candidates[kname] <- aic_val
      		model_candidates[[kname]] <- hmm_k
	      	if (store_all_models) {
      	  		model_list[[paste0(key, "_k", k)]] <- hmm_k
      		}
      		fit_counter[kname] <- fit_counter[kname]+1L
      		cat(sprintf(" * fitted k=%d  AIC=%.2f\n", k, aic_val))
    		} else {
      		aic_candidates[kname] <- Inf
      		cat(sprintf(" ! failed k=%d\n", k))
    		}
  	}
	if (all(!is.finite(aic_candidates))) {
  		cat(" skipped: all candidate models failed\n")
  		next
	}
  	best_k_name <- names(aic_candidates)[which.min(aic_candidates)]
  	best_k <- as.integer(best_k_name)
  	best_model <- model_candidates[[best_k_name]]
  	best_aic <- aic_candidates[best_k_name]
  	best_ll <- logLik(best_model)  
	trans_vec <- extract_transition_vector(best_model, df_sub, max_states = max_states)  
  	state_vec <- extract_state_params(best_model, max_states = max_states)
	feature_named_list <- c(list(ParticipantID = p, MediaID=m, n_obs=nrow(df_sub),
      		n_sequences=length(ntimes), n_states=best_k, logLik=as.numeric(best_ll),
      		AIC=as.numeric(best_aic)), as.list(trans_vec), as.list(state_vec))
  	features_list[[key]] <- as.data.frame(feature_named_list, stringsAsFactors = FALSE)
  	best_models[[key]] <- best_model
}

### MAIN RESULTS
features_df <- bind_rows(features_list)

# Pairs: Participant and Media without fitted model
pairs_df <- df %>% distinct(ParticipantID, MediaID) %>%
  arrange(ParticipantID, MediaID) %>%
  mutate(key=paste(ParticipantID, MediaID, sep="_"), is_model=key %in% names(best_models))

# missing_pairs - pairs for which the models have not been fitted
missing_pairs <- pairs_df %>% filter(!is_model)

# Data with states of best models
# (states are spacially mapped for comparison)

df_pelna_list <- lapply(names(best_models), function(key) {
	model <- best_models[[key]]
	if (is.null(model)) return(NULL)
	ids <- strsplit(key, "_")[[1]]
	p <- as.integer(ids[1])
	m <- as.integer(ids[2])
	df_sub <- df %>% filter(ParticipantID==p, MediaID==m) %>%
		mutate(x=as.numeric(OX), y=as.numeric(OY)) %>%
    		filter(!is.na(x), !is.na(y)) %>% arrange(ModelID, TrialID) %>%
    		mutate(seq_id=interaction(ParticipantID, MediaID, ModelID, drop=TRUE)) %>%
		dplyr::select(ParticipantID, MediaID, ModelID, TrialID, seq_id, x, y)
	post <- depmixS4::posterior(model, type="viterbi")
	if (is.null(post) || !"state" %in% names(post)) return(NULL)
	mapping <- get_state_mapping(model)
	if (is.null(mapping)) return(NULL)
	state_map <- mapping$map_old_to_new
	n_states_val <- features_df %>% filter(ParticipantID==p, MediaID==m) %>% pull(n_states)
		df_sub %>% mutate(n_states=n_states_val, state_raw=post$state, state=unname(state_map[as.character(state_raw)]))
  })

df_pelna <- bind_rows(df_pelna_list)

# FINAL OBJECTS
# best_models - the best HMM for each pair: Participant and Media
length(best_models)

# features_df - HMM model features 
dim(features_df)	

# df_pelna - original observations and HMM states
str(df_pelna)

# fit_counter - number of fits for each number of states state

# HMM DIAGNOSTICS:
summary(features_df$AIC)
summary(features_df$logLik)

apply(features_df[, grep("_sd_", names(features_df))], 2, summary)

# too small SD
sd_cols <- grep("_sd_", names(features_df), value=TRUE)

small_sd_005 <- features_df %>%
  mutate(flag_small_sd=if_any(all_of(sd_cols), ~ !is.na(.) & .<0.005),
         flag_tiny_sd=if_any(all_of(sd_cols), ~ !is.na(.) & .<0.001))

table(small_sd_005$flag_small_sd)
table(small_sd_005$flag_tiny_sd)

small_sd_models <- small_sd_005 %>% filter(flag_small_sd)

nrow(small_sd_models)
head(small_sd_models[, c("ParticipantID", "MediaID", "n_states", "AIC", "logLik")])

# states with small SD: 4 and 5
small_sd_005 %>% count(n_states, flag_small_sd) %>% group_by(n_states) %>%
  mutate(prop=n/sum(n))

small_sd_005 %>% count(n_states, flag_tiny_sd) %>%
  group_by(n_states) %>% mutate(prop = n / sum(n))

used_states_df <- df_pelna %>%
  group_by(ParticipantID, MediaID) %>%
  summarise(n_states_model=first(n_states), n_states_used=n_distinct(state), .groups = "drop")

table(used_states_df$n_states_model, used_states_df$n_states_used)

model_quality <- small_sd_005 %>%
  mutate(quality_flag = case_when(flag_tiny_sd~"suspect", flag_small_sd~"borderline",
      TRUE~"ok"))

table(model_quality$quality_flag)

# "tight" states, little observations 
state_sizes <- df_pelna %>% group_by(ParticipantID, MediaID, n_states, state) %>%
  summarise(n=n(), .groups="drop") %>%
  group_by(ParticipantID, MediaID) %>% mutate(prop=n/sum(n)) %>%
  ungroup()

min_state_prop <- state_sizes %>% group_by(ParticipantID, MediaID) %>%
  summarise(min_prop=min(prop), .groups="drop")

model_diag <- model_quality %>% left_join(min_state_prop, by=c("ParticipantID", "MediaID"))

summary(model_diag$min_prop)
table(model_diag$quality_flag, model_diag$min_prop<0.05)

## Check whether the suspicious states are actually far from each other, or are they "duplicates"
# For each model, we calculate the minimum Euclidean distance between the centers of the states,
# and a separation index: the distance between two states divided by their "typical size"
# We append this to model_quality
library(purrr)

# Extract centers and SD of used states from first row of features_df
get_state_info <- function(row_df) {
  n <- as.integer(row_df$n_states[1])
  out <- lapply(seq_len(n), function(s) {
    data.frame(state=s, mean_x=as.numeric(row_df[[paste0("State", s, "_mean_x")]][1]),
      mean_y=as.numeric(row_df[[paste0("State", s, "_mean_y")]][1]),
      sd_x=as.numeric(row_df[[paste0("State", s, "_sd_x")]][1]),
      sd_y=as.numeric(row_df[[paste0("State", s, "_sd_y")]][1]))
  	})
  bind_rows(out)
}

# Count state pairs and separation
calc_state_pair_metrics <- function(row_df) {
  st <- get_state_info(row_df)
  if (nrow(st)<2) return(NULL)
  pairs <- combn(st$state, 2, simplify = FALSE)
  res <- lapply(pairs, function(p) {
    i <- p[1]
    j <- p[2]
    si <- st%>%filter(state==i)
    sj <- st%>%filter(state==j)
    # distance between centers
    d <- sqrt((si$mean_x-sj$mean_x)^2+(si$mean_y-sj$mean_y)^2)
    # state "size": mean of SD in x and y
    size_i <- mean(c(si$sd_x, si$sd_y), na.rm=TRUE)
    size_j <- mean(c(sj$sd_x, sj$sd_y), na.rm=TRUE)
    # mean size of state pairs
    pooled_size <- mean(c(size_i, size_j), na.rm=TRUE)
    # separation index
    # >1 - centers further than the average "radius" of the state
    # >2 - separation
    sep_index <- d/pooled_size
    # overlap:
    # d<pooled_size - quite strong overlapping
    overlap_flag <- d<pooled_size
    data.frame(state1=i, state2=j, center_dist=d, pooled_size=pooled_size,
      sep_index=sep_index, overlap_flag=overlap_flag)
  })
  bind_rows(res)
}

# Diagnostics for each model
pairwise_diag_list <- lapply(seq_len(nrow(features_df)), function(i) {
  row_i <- features_df[i,]
  pair_metrics <- calc_state_pair_metrics(row_i)
  pair_metrics %>%
    mutate(ParticipantID=row_i$ParticipantID[1], MediaID=row_i$MediaID[1], n_states=row_i$n_states[1])
})

pairwise_diag <- bind_rows(pairwise_diag_list)

# Reduction: one record per model
center_sep_diag <- pairwise_diag %>%
  group_by(ParticipantID, MediaID, n_states) %>%
  summarise(min_center_dist=min(center_dist, na.rm=TRUE),
    median_center_dist=median(center_dist, na.rm=TRUE),
    min_sep_index=min(sep_index, na.rm=TRUE),
    median_sep_index=median(sep_index, na.rm=TRUE),
    any_overlap=any(overlap_flag, na.rm=TRUE),
    n_overlapping_pairs=sum(overlap_flag, na.rm=TRUE),.groups="drop")

# Add model quality
model_diag2 <- model_diag %>%
  left_join(center_sep_diag, by=c("ParticipantID", "MediaID", "n_states"))

# Summary - model quality
model_diag2 %>% group_by(quality_flag) %>%
  summarise(n=n(), mean_min_center_dist=mean(min_center_dist, na.rm=TRUE),
    median_min_center_dist=median(min_center_dist, na.rm=TRUE),
    mean_min_sep_index=mean(min_sep_index, na.rm=TRUE),
    median_min_sep_index=median(min_sep_index, na.rm=TRUE),
    prop_any_overlap=mean(any_overlap, na.rm=TRUE),
    mean_n_overlapping_pairs=mean(n_overlapping_pairs, na.rm=TRUE), .groups="drop"
  )

# separation index
# <1 - very poor separation
# <1.5 - rather poor separation
# >=2 = quite good separation
table(model_diag2$quality_flag, model_diag2$min_sep_index<1, useNA="ifany")
table(model_diag2$quality_flag, model_diag2$min_sep_index<1.5, useNA="ifany")
table(model_diag2$quality_flag, model_diag2$min_sep_index<2, useNA="ifany")

# Most suspicious models
worst_sep_models <- model_diag2 %>% arrange(min_sep_index, min_center_dist) %>%
  dplyr::select(ParticipantID, MediaID, n_states, quality_flag,
         min_center_dist, min_sep_index, any_overlap, n_overlapping_pairs,
         min_prop, flag_small_sd, flag_tiny_sd) %>% head(30)

# The closest state pairs
closest_pairs <- pairwise_diag %>%
  left_join(model_quality %>% dplyr::select(ParticipantID, MediaID, n_states, quality_flag, 
	flag_small_sd, flag_tiny_sd), by = c("ParticipantID", "MediaID", "n_states")) %>%
  arrange(sep_index, center_dist)

