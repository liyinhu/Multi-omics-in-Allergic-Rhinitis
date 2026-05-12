library(mixOmics) 
set.seed(123)

#Import proteome, metabolome and lipidome
proteome <- read.csv(file.choose, header = T, row.names = 1) 
metabolome <- read.csv(file.choose, header = T, row.names = 1)
lipidome <- read.csv(file.choose, header = T, row.names = 1)

data = list(Proteome = proteome, 
            Metabolome = metabolome,
            Lipidome = lipidome)

#Import group information for samples
Y = read.csv(file.choose, header = T, row.names = 1) 

#Pairwise PLS Comparisons
#Generate three pairwise PLS models
list.keepX = c(25, 25) 
list.keepY = c(25, 25)

pls1 <- spls(data[["Metabolome"]], data[["Lipidome"]], 
             keepX = list.keepX, keepY = list.keepY) 
pls2 <- spls(data[["Proteome"]], data[["Metabolome"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["Proteome"]], data[["Lipidome"]], 
             keepX = list.keepX, keepY = list.keepY)


# Plot features of PLS
plotVar(pls1, cutoff = 0.5, title = "(a) Metabolome vs Lipidome", 
        legend = c("Metabolome", "Lipidome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('#D86C76', '#14976F'))


plotVar(pls2, cutoff = 0.5, title = "(b) Proteome vs Metabolome", 
        legend = c("Proteome", "Metabolome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('#58ABDB', '#D86C76'))


plotVar(pls3, cutoff = 0.5, title = "(c) Proteome vs Lipidome", 
        legend = c("Proteome", "Lipidome"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('#58ABDB', '#14976F'))

#Initial DIABLO Model
#For square matrix filled with 0.1s
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

#Form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) 


#Tuning the number of components
#Run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning

#Set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
#Show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

#Tuning the number of features
#Set grid of values for each component to test
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))


#Run the feature selection tuning
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")

#Final DIABLO Model
#Set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

final.diablo.model$design #Design matrix for the final model


#Plots
plotDiablo(final.diablo.model, ncomp = 1)

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

#Export multi-omics relations
net <- network(final.diablo.model, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
edges <- net$gR$edges
write.csv(edges, "DIABLO_network_edges.csv")



