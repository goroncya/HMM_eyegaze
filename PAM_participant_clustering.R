# Agnieszka Goroncy, 2026
library(dplyr)
library(tidyr)
library(purrr)
library(cluster)
library(ggplot2)

max_states <- 5
occ_threshold <- 0.05
unmatched_penalty <- 0.35
occ_weight <- 0.20
sil_threshold <- 0.25

# we do not allow small clusters
min_cluster_size <- 8
max_k_clusters <- 4


# Distance components between two models
pairwise_hmm_distance_components <- function(row_a, row_b, occ_threshold=0.05, unmatched_penalty=0.35, occ_weight=0.20) {
  states_a <- extract_states_from_row(row_a)
  states_b <- extract_states_from_row(row_b)  
  mt <- optimal_state_match(clean_states=states_a, media_states=states_b, occ_threshold=occ_threshold, 
		unmatched_penalty=unmatched_penalty, occ_weight=occ_weight)
  matched <- mt$matched
  a_only  <- mt$clean_only
  b_only  <- mt$media_only
  if (is.null(matched)) matched <- tibble::tibble()
  if (is.null(a_only))  a_only  <- tibble::tibble()
  if (is.null(b_only))  b_only  <- tibble::tibble()
  # occupancy distance
  d_occ_matched <- if (nrow(matched)>0) sum(abs(matched$occ_clean-matched$occ_media)) else 0
  d_occ_unmatched <- sum(a_only$occ_clean, na.rm=TRUE)+sum(b_only$occ_media, na.rm=TRUE)
  d_occ <- d_occ_matched+d_occ_unmatched
  # center distance
  if (nrow(matched) > 0) { 
	center_d <- sqrt((matched$clean_x-matched$media_x)^2+(matched$clean_y-matched$media_y)^2)
    center_w <- (matched$occ_clean+matched$occ_media)/2
    d_centers_matched <- weighted.mean(center_d, w=center_w, na.rm=TRUE)} 
	else {
    d_centers_matched <- 0
  	}  
  # penalty for mismatched states
  unmatched_occ <- sum(a_only$occ_clean, na.rm=TRUE)+sum(b_only$occ_media, na.rm=TRUE)
  d_centers <- d_centers_matched+unmatched_penalty*unmatched_occ
  # transition distance
  trans_cols <- grep("^T_S[0-9]+_S[0-9]+$", names(row_a), value=TRUE)
#  d_trans <- transition_distance(row_a, row_b, matched)
  d_trans <- transition_distance(row_a, row_b, trans_cols) 
  tibble(d_occ = d_occ,d_centers = d_centers,d_trans = d_trans)
}

# distance matrix for 1 medium
build_medium_distance_matrix <- function(hmm_features, media_id, occ_threshold=0.05,
							unmatched_penalty=0.35, occ_weight=0.20) {
  df_m <- hmm_features %>% filter(MediaID == media_id) %>% arrange(ParticipantID)
  n <- nrow(df_m)
  if (n<2) return(NULL)
  combs <- combn(seq_len(n), 2, simplify=FALSE)
  pair_df <- bind_rows(lapply(combs, function(idx) {
    i <- idx[1]
    j <- idx[2] 
    comps <- pairwise_hmm_distance_components(row_a=df_m[i,, drop=FALSE], row_b=df_m[j,, drop=FALSE],
      occ_threshold=occ_threshold, unmatched_penalty=unmatched_penalty, occ_weight=occ_weight)
    tibble(i=i, j=j, ParticipantA=df_m$ParticipantID[i], ParticipantB=df_m$ParticipantID[j], d_occ=comps$d_occ,
      d_centers=comps$d_centers, d_trans=comps$d_trans)
  }))
 # scaling components
  sd_occ <- sd(pair_df$d_occ, na.rm=TRUE)
  sd_centers <- sd(pair_df$d_centers, na.rm=TRUE)
  sd_trans <- sd(pair_df$d_trans, na.rm=TRUE)  
  if (is.na(sd_occ)||sd_occ==0) sd_occ <- 1
  if (is.na(sd_centers)||sd_centers==0) sd_centers <- 1
  if (is.na(sd_trans)||sd_trans==0) sd_trans <- 1
  pair_df <- pair_df %>% mutate(d_occ_s=d_occ/sd_occ, d_centers_s=d_centers/sd_centers,
		d_trans_s=d_trans/sd_trans, d_total=d_occ_s+d_centers_s+d_trans_s)
  # dinstance matrix
  D <- matrix(0, nrow=n, ncol=n)
  for (r in seq_len(nrow(pair_df))) {
    i <- pair_df$i[r]
    j <- pair_df$j[r]
    D[i, j] <- pair_df$d_total[r]
    D[j, i] <- pair_df$d_total[r]
  }
  rownames(D) <- df_m$ParticipantID
  colnames(D) <- df_m$ParticipantID
  list(media_id=media_id, df_models=df_m, pair_df=pair_df, D=D)
}

# Clustering media
# Partitioning Around Medoids (PAM), choice of k using silhouette
cluster_medium_patterns <- function(dist_obj,max_k_clusters=4, sil_threshold=0.25, min_cluster_size=8) {
  D <- dist_obj$D
  n <- nrow(D)
  if (n<3) {
    return(list(
      summary=tibble(MediaID=dist_obj$media_id, n_models=n, best_k=1, best_sil=NA_real_, final_k=1, cluster_sizes="all in one cluster"),
      assignments = tibble(ParticipantID=dist_obj$df_models$ParticipantID, MediaID=dist_obj$media_id, cluster=1), pam_fits = NULL))
  }
  k_candidates <- 2:min(max_k_clusters, n-1)
  pam_results <- lapply(k_candidates, function(k) {
    fit <- pam(as.dist(D), k=k, diss=TRUE)
    sizes <- table(fit$clustering)
    avg_sil <- fit$silinfo$avg.width
    tibble(k=k, avg_sil=avg_sil, min_size=min(sizes), valid=min(sizes)>= min_cluster_size) %>%
      mutate(fit = list(fit))
  })
  pam_tbl <- bind_rows(pam_results)
# selecting the best k among solutions with a reasonable cluster size
  pam_valid <- pam_tbl %>% filter(valid)  
  if (nrow(pam_valid)==0 || max(pam_valid$avg_sil, na.rm=TRUE)<sil_threshold) {
    final_k <- 1
    assignments <- tibble(ParticipantID=dist_obj$df_models$ParticipantID,
      MediaID=dist_obj$media_id, cluster=1)
    cluster_sizes <- "all in one cluster"
    best_k <- if (nrow(pam_tbl)>0) pam_tbl$k[which.max(pam_tbl$avg_sil)] else 1
    best_sil <- if (nrow(pam_tbl)>0) max(pam_tbl$avg_sil, na.rm = TRUE) else NA_real_
  } else {
    best_row <- pam_valid %>% slice_max(order_by=avg_sil, n=1, with_ties=FALSE)
    best_fit <- best_row$fit[[1]]
    final_k <- best_row$k
    assignments <- tibble(ParticipantID=dist_obj$df_models$ParticipantID, MediaID=dist_obj$media_id, cluster=best_fit$clustering)
    cluster_sizes <- paste(as.integer(table(best_fit$clustering)), collapse="/")
    best_k <- best_row$k
    best_sil <- best_row$avg_sil
  }
  list(summary=tibble(MediaID=dist_obj$media_id, n_models=n, best_k=best_k, best_sil=best_sil, final_k=final_k,
      cluster_sizes=cluster_sizes), assignments=assignments, pam_fits=pam_tbl)
}

media_ids <- sort(unique(hmm_features$MediaID))

dist_objects <- lapply(media_ids, function(m) {
  cat("Building distance matrix for Media", m, "\n")
  build_medium_distance_matrix(hmm_features=hmm_features, media_id=m, occ_threshold=occ_threshold, unmatched_penalty=unmatched_penalty, occ_weight=occ_weight)})

names(dist_objects) <- paste0("M", media_ids)

cluster_results <- lapply(dist_objects, function(obj) {
  cat("Clustering Media", obj$media_id, "\n")
  cluster_medium_patterns(dist_obj=obj, max_k_clusters=max_k_clusters, sil_threshold=sil_threshold, min_cluster_size=min_cluster_size)})

cluster_summary <- bind_rows(lapply(cluster_results, `[[`, "summary")) %>% arrange(MediaID)

cluster_assignments <- bind_rows(lapply(cluster_results, `[[`,"assignments"))


# Medoids: cluster representatives
get_cluster_medoids <- function(dist_obj, assignments_df, media_id) {
  D <- dist_obj$D
  ass <- assignments_df %>% filter(MediaID==media_id)
  out <- lapply(sort(unique(ass$cluster)), function(cl) {
    ids <- ass$ParticipantID[ass$cluster==cl]
    idx <- which(rownames(D) %in% ids)
    if (length(idx)==1) {
      medoid_id <- rownames(D)[idx]
    } else {
      subD <- D[idx, idx, drop=FALSE]
      medoid_id <- rownames(subD)[which.min(rowSums(subD))]
    }
    tibble(MediaID=media_id, cluster=cl, n=length(idx), medoid_participant=as.integer(medoid_id))
  })
  bind_rows(out)
}

cluster_medoids <- bind_rows(lapply(media_ids, function(m) {
  get_cluster_medoids(dist_obj=dist_objects[[paste0("M", m)]],
    assignments_df=cluster_assignments, media_id=m)}))
