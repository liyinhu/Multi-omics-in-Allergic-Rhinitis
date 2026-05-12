library(qgraph)
library(bootnet)
library(dplyr)
library(tidyr)
library(tidyverse)
library(corrplot)


## 1. Import organ matrix
organ <- read.csv(file.choose(), row.names = 1, check.names = FALSE)

# Import group information
info <- read.delim(file.choose())  # Sample | Group (AR / HC)

# Select common samples
common_samples <- intersect(rownames(organ), info$ID)
organ <- organ[common_samples, ]
info  <- info[match(common_samples, info$ID), ]

group <- info$Group

## 2. Z-score normalization
organ_z <- scale(organ)
organ_z <- as.data.frame(organ_z)


## 3. Grouping
organ_AR <- organ_z[group == "2AR", ]
organ_HC <- organ_z[group == "1HC", ]

## 4. Bootstrap
## 4.1 Network estimation
estimate_glasso <- function(data) {
  EBICglasso(cor(data), n = nrow(data), gamma = 0.5)
}

## 4.2 Bootstrap
set.seed(123)

boot_AR <- bootnet(
  organ_AR,
  default = "EBICglasso",
  nBoots = 1000,
  type = "nonparametric",
  tuning = 0.75
)

boot_HC <- bootnet(
  organ_HC,
  default = "EBICglasso",
  nBoots = 1000,
  type = "nonparametric",
  tuning = 0.75
)

## 5. Isolate bootstrap results
## 5.1 Isolate bootTable for AR and HC  
boot_AR_edges <- boot_AR$bootTable %>% 
  filter(type == "edge") %>%
  rename(value_AR = value)

boot_HC_edges <- boot_HC$bootTable %>% 
  filter(type == "edge") %>%
  rename(value_HC = value)

## 5.2 Summarize the bootstrap statistics for each edge
edge_AR <- boot_AR_edges %>%
  group_by(node1, node2) %>%
  summarise(
    mean = mean(value_AR, na.rm = TRUE),
    sd   = sd(value_AR, na.rm = TRUE),
    CI_lower = quantile(value_AR, 0.025, na.rm = TRUE),
    CI_upper = quantile(value_AR, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

edge_HC <- boot_HC_edges %>%
  group_by(node1, node2) %>%
  summarise(
    mean = mean(value_HC, na.rm = TRUE),
    sd   = sd(value_HC, na.rm = TRUE),
    CI_lower = quantile(value_HC, 0.025, na.rm = TRUE),
    CI_upper = quantile(value_HC, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## 6. Difference network between AR and HC
## 6.1 Merge the bootstrap edge tables of AR and HC
boot_AR_edges <- boot_AR$bootTable %>% 
  filter(type == "edge") %>%
  rename(value_AR = value)

boot_HC_edges <- boot_HC$bootTable %>% 
  filter(type == "edge") %>%
  rename(value_HC = value)

## 6.2 Paired by edge and bootstrap
boot_diff <- inner_join(
  boot_AR_edges,
  boot_HC_edges,
  by = c("node1", "node2", "name")
)

## 6.3 Calculate bootstrap differences
boot_diff <- boot_diff %>%
  mutate(diff = value_AR - value_HC)

## 6.4 Calculate the bootstrap p-value for each edge
edge_diff <- boot_diff %>%
  group_by(node1, node2) %>%
  summarise(
    mean_diff = mean(diff),
    p = {
      prop_neg <- mean(diff <= 0)
      prop_pos <- mean(diff >= 0)
      p_raw <- 2 * min(prop_neg, prop_pos)
      min(p_raw, 1)
    },
    CI_lower = quantile(diff, 0.025),
    CI_upper = quantile(diff, 0.975),
    .groups = "drop"
  )

## 6.5 Adjustement for multiple test
edge_diff$padj <- p.adjust(edge_diff$p, method = "BH")


## 7. Export results
write.csv(edge_AR, "glasso_edges_AR_bootstrap.csv", row.names = FALSE)
write.csv(edge_HC, "glasso_edges_HC_bootstrap.csv", row.names = FALSE)
write.csv(edge_diff, "glasso_edges_AR_vs_HC_diff.csv", row.names = FALSE)

## 8. Convert edge_diff into a matrix

edge_diff_rev <- edge_diff %>%
    rename(
        node1_tmp = node1,
        node2_tmp = node2
    ) %>%
    transmute(
        node1 = node2_tmp,
        node2 = node1_tmp,
        mean_diff  = mean_diff,
        p     = p,
        CI_lower = CI_lower,
        CI_upper = CI_upper,
        padj  = padj
)

edge_diff_bidir <- bind_rows(edge_diff, edge_diff_rev)

organs <- sort(unique(c(edge_diff_bidir$node1, edge_diff_bidir$node2)))

full_pairs <- expand.grid(
  node1 = organs,
  node2 = organs,
  stringsAsFactors = FALSE
) %>%
  filter(node1 != node2)

edge_full <- full_pairs %>%
  left_join(
    edge_diff_bidir %>%
      select(node1, node2, mean_diff, p, padj),
    by = c("node1", "node2")
  ) %>%
  mutate(
    diff = ifelse(is.na(mean_diff), 0, mean_diff),
    p    = ifelse(is.na(p), 1, p),
    padj = ifelse(is.na(padj), 1, padj)
  )


make_matrix <- function(edge_full, value_col, organs) {
    edge_full %>%
        mutate(
            node1 = factor(node1, levels = organs),
            node2 = factor(node2, levels = organs)
        ) %>%
        select(node1, node2, {{ value_col }}) %>%
        pivot_wider(
            names_from  = node2,
            values_from = {{ value_col }}
        ) %>%
        arrange(node1) %>%
        column_to_rownames("node1") %>%
        as.matrix()
}

diff_mat  <- make_matrix(edge_full, diff,  organs)
p_mat     <- make_matrix(edge_full, p,     organs)
padj_mat  <- make_matrix(edge_full, padj,  organs)


diff_mat  <- (diff_mat + t(diff_mat))/2 
p_mat     <- (p_mat + t(p_mat))/2
padj_mat  <- (padj_mat + t(padj_mat))/2


diag(diff_mat) <- 0
diag(p_mat)    <- 1
diag(padj_mat) <- 1

## 9. Export results
write.csv(diff_mat, "diff_edges_mat.csv")
write.csv(p_mat, " diff_edges_p.csv")
write.csv(padj_mat, "diff_edges_padj.csv")


##10.Plotting
diff_mat_scaled <- diff_mat / 0.2

corrplot(diff_mat_scaled,
         p.mat= p_mat,
         method = 'circle',
         order = 'AOE', 
         type = 'lower',
         diag = FALSE,
         insig = 'label_sig', sig.level = c(0.01, 0.05, 0.10), pch.cex = 0.5,
         col = COL2('RdBu', 8) 
)


col2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(6)

corrplot(as.matrix(diff_mat_scaled),
         p.mat = as.matrix(p_mat),
         order = 'FPC', 
         type = 'lower',
         diag = FALSE,
         insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5,
         col = col2
)