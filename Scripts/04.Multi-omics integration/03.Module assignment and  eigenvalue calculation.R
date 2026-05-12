library(MetaNet)
library(dplyr)
library(tibble)
library(WGCNA)
library(tidyr)

#Import edges from DIABOLO and mmVec models
edge_list1 <- read.csv(edges_from_DIABOLO.csv)
edge_list2 <- read.csv(edges_from_mmVec.csv)

edge_list <- rbind(edge_list1, edge_list2)

#Construct multi-omics network
dnet <- c_net_from_edgelist(edge_list, direct = F)

#Module detection
co_net_modu <- module_detect(dnet, method = "cluster_fast_greedy")
node_modu <- get_v(co_net_modu)[, c("name", "module")]
write.csv(node_modu, "Node_module.csv")

#Plot
plot_module_tree(co_net_modu, label.size = 0.6)

## Covert to long format: rows and columns represent features and samples, respectively
expr_mat <- read.csv(multi_omics_signatures.csv)
expr_long <- as.data.frame(t(expr_mat))
expr_long$Feature <- rownames(expr_long)

#Omics integration
expr_annot <- expr_long %>%
  left_join(feature_type, by = "Feature")

#Z-score normalization by features
expr_scaled_long <- expr_annot %>%
  group_by(Type) %>%
  mutate(across(where(is.numeric), scale)) %>%
  ungroup()

#Convert to wide format
expr_scaled <- expr_scaled_long %>%
  select(-Type) %>%
  column_to_rownames("Feature") %>%
  as.data.frame()

expr_scaled <- t(expr_scaled) 

node_modu$module <- as.character(node_modu$module)

#Ensure feature consistency
expr_scaled <- expr_scaled[intersect(rownames(expr_scaled), node_modu$name), ]
node_modu <- node_modu %>% filter(name %in% rownames(expr_scaled))

#Calculate the eigengene by modules
modules_list <- split(node_modu$name, node_modu$module)

datExpr <- t(expr_scaled)

#Calcualte module eigengenes
ME_list <- lapply(names(modules_list), function(mod) {
  subset_features <- modules_list[[mod]]
  if (length(subset_features) >= 2) {
    dat_mod <- datExpr[, subset_features, drop = FALSE]

    # Process NA
    dat_mod <- as.matrix(dat_mod)
    storage.mode(dat_mod) <- "numeric"
    dat_mod <- dat_mod[, colSums(is.na(dat_mod)) == 0, drop = FALSE]
    
    # moduleEigengenes
    if (ncol(dat_mod) >= 2) {
      ME <- moduleEigengenes(dat_mod, colors = rep(mod, ncol(dat_mod)))$eigengenes
      colnames(ME) <- mod
      return(ME)
    }
  }
  return(NULL)
})

#Export the results
ME_matrix <- do.call(cbind, ME_list)
write.csv(ME_matrix, "ME_matrix.csv")