library(igraph)
library(dplyr)
library(tidyr)
library(purrr)

#Import data
edges <- read.csv(multi-omics network.edge, header = T)
nodes <- read.csv(multi-omics network.nodes, header = T)

#Network construction
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

#Module isolation
modules <- unique(nodes$module)

#Initialization
all_paths <- list()

for (mod in modules) {
    cat("Processing module:", mod, "\n")
    
    #Subpath isolation
    mod_nodes <- nodes %>% filter(module == mod)
    mod_edges <- edges %>% filter(Source %in% mod_nodes$ID & Target %in% mod_nodes$ID)
    
    g_mod <- graph_from_data_frame(mod_edges, vertices = mod_nodes, directed = FALSE)
       
    micro_nodes <- V(g_mod)[Group == "Microbiome"]$name
    
    for (start_node in micro_nodes) {
        paths_mod <- all_simple_paths(g_mod, from = start_node, cutoff = 5)
    }
    
        paths_named <- lapply(paths_mod, function(p) names(p))
    
        keep_paths <- list()
    
    for (p in paths_named) {
        p_group <- mod_nodes$Group[match(p, mod_nodes$ID)]
        
        if (any(p_group == "Microbiome") &&
            any(p_group %in% c("Metabolome", "Lipidome")) &&
            any(p_group == "Proteome")) {
                        
            if (length(p) >= 3) {
                if (length(p) == 3) {
                    triplets <- matrix(p, nrow = 1)
                } else {
                    triplets <- embed(rev(p), 3)[, 3:1, drop = FALSE]  
                }
                
                for (row in seq_len(nrow(triplets))) {
                    trip <- triplets[row, ]
                    g_trip <- mod_nodes$Group[match(trip, mod_nodes$ID)]
                    
                    valid <- any(
                        all(g_trip == c("Microbiome", "Metabolome", "Proteome")),
                        all(g_trip == c("Proteome", "Metabolome", "Microbiome")),
                        all(g_trip == c("Microbiome", "Lipidome", "Proteome")),
                        all(g_trip == c("Proteome", "Lipidome", "Microbiome"))
                    )
                    
                    if (valid) {
                        mid <- trip[2]
                        ends <- sort(c(trip[1], trip[3]))
                        path_id <- paste(mid, paste(ends, collapse = "_"), sep = "|")
                        keep_paths[[path_id]] <- trip
                    }
                }
            }
        }
    }
    all_paths[[as.character(mod)]] <- keep_paths 
}



# Export data.frame
final_path1 <- unique(all_paths[[as.character(1)]])
df_paths1 <- bind_rows(lapply(final_path1, function(x) {
    data.frame(Source = x[1], Middle = x[2], Target = x[3])
}))

write.csv(df_paths1, "Subpath.csv", row.names = FALSE)