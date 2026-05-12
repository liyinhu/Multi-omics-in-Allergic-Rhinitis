library(cluster)
library(factoextra)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(umap)
library(vegan)

# Import target multi-omics data
data <- read.csv(file.choose(), header = TRUE, row.names = 1)

# === Step 1: Determine the optimal number of clusters ===
max_k <- 10
wss <- numeric(max_k)
wss[1] <- sum(scale(data, center = TRUE, scale = FALSE)^2)

set.seed(666)
for (k in 2:max_k) {
  km <- kmeans(data, centers = k, nstart = 5)
  wss[k] <- km$tot.withinss
}

plot(1:max_k, wss, type = "b", pch = 19, col = "#58ABDB",
     xlab = "Number of Clusters (K)", ylab = "Within-cluster sum of squares", main = "Elbow Curve")
abline(v = which.min(diff(diff(wss))), col = "red", lty = 2) 

# === Step 2: Clustering (Optimal k can be adjusted according to the above figure)===
best_k <- 5
set.seed(123)
km_res <- kmeans(data, centers = best_k, nstart = 10)

data_clustered <- data.frame(Sample = rownames(data),
                              Cluster = factor(km_res$cluster),
                              data)

# === Step 3: PCoA visualization ===

dist_mat <- dist(data)

pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(data),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  Cluster = factor(km_res$cluster)
)

system_colors <- structure(
    RColorBrewer::brewer.pal(length(unique(km_res$cluster)), "Set2"),
    names = unique(km_res$cluster)
)

eig <- pcoa_res$eig

ggplot(pcoa_df)+ geom_point(aes(x = PCoA1, y = PCoA2, color = Cluster))+stat_ellipse(aes(x = PCoA1, y = PCoA2, group=Cluster, color=Cluster), level = 0.95)+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) + labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+scale_color_manual(values = system_colors)


# === Step 4: Remove clusters with ≤ 5 samples ===
min_cluster_size <- 5
cluster_size <- table(data_clustered$Cluster)
valid_clusters <- names(cluster_size[cluster_size >= min_cluster_size])

filtered_df <- data_clustered %>% filter(Cluster %in% valid_clusters)

long_df <- filtered_df %>%
  pivot_longer(cols = -c(Sample, Cluster), names_to = "Feature", values_to = "Value")

# === Step 5: Kruskal-Wallis + DunnTest Analysis ===
diff_res <- long_df %>%
  group_by(Feature) %>%
  summarise(p = kruskal.test(Value ~ Cluster)$p.value) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

sig_feats <- diff_res %>% filter(p_adj < 0.05) %>% pull(Feature)

dunn_results <- list()
for (f in sig_feats) {
  subset_data <- long_df %>% filter(Feature == f)
  dunn <- dunnTest(Value ~ Cluster, data = subset_data, method = "bh")
  dunn_df <- dunn$res
  dunn_df$Feature <- f
  dunn_results[[f]] <- dunn_df
}
dunn_all <- do.call(rbind, dunn_results)

# === Step 6: Isolate representative features for the clusters ===
cluster_means <- long_df %>%
  filter(Feature %in% sig_feats) %>%
  group_by(Feature, Cluster) %>%
  summarise(mean_val = mean(Value), .groups = "drop")

top_cluster_feats <- cluster_means %>%
  group_by(Feature) %>%
  filter(mean_val == max(mean_val)) %>%
  ungroup()

rep_features <- top_cluster_feats %>%
  group_by(Cluster) %>%
  summarise(Representative_Features = paste(Feature, collapse = ", "))

print(rep_features)

# === Step 7: Visualize the representative features  ===

for (f in unique(top_cluster_feats$Feature)) {
    p <- ggplot(long_df %>% filter(Feature == f),
                aes(x = Cluster, y = Value, fill = Cluster)) +
        geom_violin(trim = FALSE, alpha = 0.6) +
        geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +theme_bw()+
        scale_fill_manual(values = system_colors) +
        labs(title = paste("Violin Plot of", f)) +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
    print(p)
}

# === Step 8: Export results ===
write.csv(diff_res, "Kruskal_Wallis_pvalues.csv", row.names = FALSE)
write.csv(dunn_all, "DunnTest_pairwise_results.csv", row.names = FALSE)
write.csv(rep_features, "Cluster_Representative_Features.csv", row.names = FALSE)
write.csv(data_clustered, "Sample_clustered.csv", row.names = FALSE)