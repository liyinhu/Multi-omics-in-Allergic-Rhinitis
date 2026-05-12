library(lme4)
library(mediation)
library(dplyr)

results_list <- list()

#Improt isolated subpath
triplets_df <- read.csv(Subpath.csv, header = T)
all_data <- read.csv(multi_omics_signatures.csv, header = T, row.names = 1)

#双向中介效应分析
for (i in 1:nrow(triplets_df)) {
    microbe <- triplets_df$Source[i]
    metabolite <- triplets_df$Middle[i]
    protein <- triplets_df$Target[i]
    
    dat <- data.frame(
        microbe = all_data[[microbe]],
        metabolite = all_data[[metabolite]],
        protein = all_data[[protein]],
        Group = all_data$Group
    )
    
    # Model A：microbe → metabolite → protein

    model.m1 <- lm(metabolite ~ microbe, dat)
    model.y1 <- lm(protein ~ microbe + metabolite, data = dat)
    
    med1 <- mediate(model.m1, model.y1, treat = "microbe", mediator = "metabolite", sims = 500)
    
    # Model B：protein → metabolite → microbe
    model.m2 <- lm(metabolite ~ protein, dat)
    model.y2 <- lm(microbe ~ protein + metabolite, dat)

    med2 <- mediate(model.m2, model.y2, treat = "protein", mediator = "metabolite", sims = 500)
    
    # Results
    extract_summary <- function(med, direction) {
        data.frame(
            Direction = direction,
            ACME = med$d0,
            ACME_CI_Lower = med$d0.ci[1],
            ACME_CI_Upper = med$d0.ci[2],
            ACME_p = med$d0.p,
            
            ADE = med$z0,
            ADE_CI_Lower = med$z0.ci[1],
            ADE_CI_Upper = med$z0.ci[2],
            ADE_p = med$z0.p,
            
            Total_Effect = med$tau.coef,
            Total_CI_Lower = med$tau.ci[1],
            Total_CI_Upper = med$tau.ci[2],
            Total_p = med$tau.p,
            
            Prop_Mediated = med$n0,
            Prop_p = med$n0.p
        )
    }
    
    res_A <- extract_summary(med1, "Microbe→Metabolite→Protein")
    res_B <- extract_summary(med2, "Protein→Metabolite→Microbe")
    
    result_row <- bind_rows(res_A, res_B) %>%
        mutate(
            Microbe = microbe,
            Metabolite = metabolite,
            Protein = protein,
            Triplet = paste(microbe, metabolite, protein, sep = " | ")
        )
    
    results_list[[i]] <- result_row
}

final_results <- bind_rows(results_list)

write.csv(final_results, "Mediation_results.csv", row.names = FALSE)