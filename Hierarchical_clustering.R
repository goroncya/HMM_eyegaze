# Agnieszka Goroncy, 2026

library(dplyr)
library(cluster)

sil_threshold <- 0.25
min_cluster_size <- 8
max_k_clusters <- 4

cluster_medium_hierarchical <- function(dist_obj, method="average", max_k_clusters=4,
	sil_threshold=0.25, min_cluster_size=8) {
  D <- dist_obj$D
  n <- nrow(D)
  if (n<3) {
    return(list(summary=tibble(MediaID=dist_obj$media_id, linkage=method,
        n_models=n, best_k=1, best_sil=NA_real_, final_k=1, cluster_sizes="all in one cluster"),
      assignments = tibble(ParticipantID=dist_obj$df_models$ParticipantID,
		MediaID=dist_obj$media_id, linkage=method,cluster=1),
      fit = NULL,eval_table = NULL))
  }
  # Agglomerative Nesting (Hierarchical Clustering)
  fit <- agnes(as.dist(D), diss=TRUE, method=method)
  hc <- as.hclust(fit)
  k_candidates <- 2:min(max_k_clusters, n-1)
  eval_table <- bind_rows(lapply(k_candidates, function(k) {
    cl <- cutree(hc, k=k)
    sizes <- table(cl)
    sil <- silhouette(cl, dist=as.dist(D))
    avg_sil <- mean(sil[, "sil_width"])
    tibble(k=k, avg_sil=avg_sil, min_size=min(sizes),valid=min(sizes)>=min_cluster_size, cluster_sizes=paste(as.integer(sizes), collapse="/"))
  }))
  eval_valid <- eval_table %>% filter(valid)
  if (nrow(eval_valid)==0 || max(eval_valid$avg_sil, na.rm=TRUE) < sil_threshold) {
    final_k <-1
    best_k <- if (nrow(eval_table)>0) eval_table$k[which.max(eval_table$avg_sil)] else 1
    best_sil <- if (nrow(eval_table)>0) max(eval_table$avg_sil, na.rm = TRUE) else NA_real_
    cluster_sizes <- "all in one cluster"
    assignments <- tibble( ParticipantID=dist_obj$df_models$ParticipantID,
      MediaID = dist_obj$media_id, linkage = method, cluster = 1)
  } else {
    best_row <- eval_valid %>% slice_max(order_by=avg_sil, n=1, with_ties=FALSE)
    final_k <- best_row$k
    best_k <- best_row$k
    best_sil <- best_row$avg_sil
    cluster_sizes <- best_row$cluster_sizes
    cl <- cutree(hc, k=final_k)
    assignments <- tibble(ParticipantID=dist_obj$df_models$ParticipantID,
      MediaID=dist_obj$media_id, linkage=method, cluster=cl)
  }
  list(summary=tibble(MediaID=dist_obj$media_id, linkage=method, n_models=n,
      best_k=best_k, best_sil=best_sil, final_k=final_k, cluster_sizes=cluster_sizes),
    assignments=assignments, fit=hc, eval_table=eval_table)
}

# average method (UPGMA)
hc_results_avg <- lapply(dist_objects, function(obj) {
  cat("Hierarchical clustering (average) for Media", obj$media_id, "\n")
  cluster_medium_hierarchical(dist_obj=obj, method="average", max_k_clusters=max_k_clusters, sil_threshold=sil_threshold,
    min_cluster_size=min_cluster_size)})

hc_summary_avg <- bind_rows(lapply(hc_results_avg, `[[`, "summary")) %>% arrange(MediaID)

hc_assignments_avg <- bind_rows(lapply(hc_results_avg, `[[`, "assignments"))

# complete linkage
hc_results_complete <- lapply(dist_objects, function(obj) {
  cat("Hierarchical clustering (complete) for Media", obj$media_id, "\n")
  cluster_medium_hierarchical(dist_obj=obj, method="complete", max_k_clusters=max_k_clusters, sil_threshold=sil_threshold,
    min_cluster_size=min_cluster_size)})
