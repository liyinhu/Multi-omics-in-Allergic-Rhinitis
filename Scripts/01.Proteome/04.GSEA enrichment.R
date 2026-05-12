library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(forcats)
library(dplyr)
library(ggplot2)

#Import differential protein analysis results
data <- read.csv(file.choose(), header = T)

#Isolate protein ID, Estimate size and P-value
df <- data[,c(1,2,5)]
df$SYMBOL <- df$ProteinID

#Protein ID transformation
entrez <- bitr(df$ProteinID,
fromType = 'SYMBOL',
toType = c('ENTREZID'),
OrgDb = 'org.Hs.eg.db')
df <- merge(df, entrez, by.y="SYMBOL")

#导出gene list并排序
geneList <- df$AR_Estimate
names(geneList) <- df$ENTREZID
geneList <- sort(geneList,decreasing = T)

#GSEA KEGG enrichment analysis
kegg <- gseKEGG(geneList,
organism = 'hsa',
pvalueCutoff = 0.2,
pAdjustMethod = 'BH',
minGSSize = 5,
maxGSSize = 200
)
kegg <- append_kegg_category(kegg)

#Convert to readable mode
kegg <- setReadable(kegg,
OrgDb= org.Hs.eg.db,
keyType= "ENTREZID"
)

#Exprot enrichment results
kegg_result<- kegg@result
write.csv(kegg_result, "AR_gsea_KEGG.csv")


#Plotting 
kegg_result2 <- kegg_result %>% mutate(Description= fct_reorder(Description, NES))

ggplot(data = kegg_result2, aes(x = NES, y = Description))+
geom_point(aes(alpha=0.7, size= -log10(pvalue), color = NES))+
theme_bw() +
theme(panel.grid.major= element_line(color = "white"),panel.grid.minor =element_line(color= "white"),legend.title = element_blank())+
scale_color_gradient2(low='#205D9B',high='#D4423A', limits=c(-1.8,1.8))+
labs(x = "Normalized Enrichment Score", y = "", title = "AR vs HC") +
scale_size("-log10(P-value)", range=c(3,7)
)


