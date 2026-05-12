#Under qiime2-2020 environment

#Import microbiome, proteome, metabolome and lipidome
biom convert -i ./Microbiome.sample_aligned.txt -o ./Microbiome.biom --table-type="Taxon table" --to-json
biom convert -i ./Proteome.sample_aligned.txt -o ./Proteome.biom --table-type="Metabolite table" --to-json
biom convert -i ./Metabolome.sample_aligned.txt -o ./Metabolome.biom --table-type="Metabolite table" --to-json
biom convert -i ./Lipidome.sample_aligned.txt -o ./Lipidome.biom --table-type="Metabolite table" --to-json

#File transformation
qiime tools import \
        --input-path ./Microbiome.biom \
        --input-format BIOMV100Format \
        --output-path ./Microbiome.qza \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path ./Proteome.biom \
        --input-format BIOMV100Format \
        --output-path ./Proteome.qza \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path ./Metabolome.biom \
        --input-format BIOMV100Format \
        --output-path ./Metabolome.qza \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path ./Lipidome.biom \
        --input-format BIOMV100Format \
        --output-path ./Lipidome.qza \
        --type FeatureTable[Frequency]

##mmVec model training
#Microbiome-Proteome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Proteome.qza \
	    --p-epochs 100 \
	    --p-latent-dim 3 \
	    --p-learning-rate 0.001 \
        --p-summary-interval 1 \
	    --p-no-arm-the-gpu \
        --output-dir Microbiome_Proteome_Summary

#Microbiome-Metabolome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Metabolome.qza \
        --p-epochs 100 \
        --p-latent-dim 3 \
	    --p-learning-rate 0.001 \
        --p-summary-interval 1 \
        --p-no-arm-the-gpu \
        --output-dir Microbiome_Metabolome_Summary

#Microbiome-Lipidome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Lipidome.qza \
        --p-epochs 100 \
        --p-latent-dim 3 \
	    --p-learning-rate 0.001 \
        --p-summary-interval 1 \
        --p-no-arm-the-gpu \
        --output-dir Microbiome_Lipidome_Summary

##Export feature relations
#Microbiome-Proteome
qiime metadata tabulate \
        --m-input-file ./Microbiome_Proteome_Summary/conditionals.qza \
        --o-visualization ./Microbiome_Proteome_conditionals_viz.qzv

qiime emperor biplot \
        --i-biplot ./Microbiome_Proteome_Summary/conditional_biplot.qza \
        --m-sample-metadata-file ./metadata/Proteome_metadata.txt \
	    --m-feature-metadata-file ./metadata/Microbiome_metadata.txt \
	    --p-number-of-features 15 \
	    --p-ignore-missing-samples \
        --o-visualization ./Microbiome_Proteome_emperor.qzv

qiime mmvec summarize-single \
        --i-model-stats ./Microbiome_Proteome_Summary/model_stats.qza \
	    --o-visualization ./Microbiome_Proteome_model_summary.qzv

#Microbiome-Metabolome
qiime metadata tabulate \
        --m-input-file ./Microbiome_Metabolome_Summary/conditionals.qza \
        --o-visualization ./Microbiome_Metabolome_conditionals_viz.qzv

qiime emperor biplot \
        --i-biplot ./Microbiome_Metabolome_Summary/conditional_biplot.qza \
        --m-sample-metadata-file ./metadata/Metabolome_metadata.txt \
        --m-feature-metadata-file ./metadata/Microbiome_metadata.txt \
        --p-number-of-features 15 \
        --p-ignore-missing-samples \
        --o-visualization ./Microbiome_Metabolome_emperor.qzv

qiime mmvec summarize-single \
        --i-model-stats ./Microbiome_Metabolome_Summary/model_stats.qza \
        --o-visualization ./Microbiome_Metabolome_model_summary.qzv


#Microbiome-Lipidome
qiime metadata tabulate \
        --m-input-file ./Microbiome_Lipidome_Summary/conditionals.qza \
        --o-visualization ./Microbiome_Lipidome_conditionals_viz.qzv

qiime emperor biplot \
        --i-biplot ./Microbiome_Lipidome_Summary/conditional_biplot.qza \
        --m-sample-metadata-file ./metadata/Lipidome_metadata.txt \
        --m-feature-metadata-file ./metadata/Microbiome_metadata.txt \
        --p-number-of-features 15 \
        --p-ignore-missing-samples \
        --o-visualization ./Microbiome_Lipidome_emperor.qzv

qiime mmvec summarize-single \
        --i-model-stats ./Microbiome_Lipidome_Summary/model_stats.qza \
        --o-visualization ./Microbiome_Lipidome_model_summary.qzv


##Q2 calculation
#Microbiome-Proteome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Proteome.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir Proteome_Summary

qiime mmvec summarize-paired \
        --i-model-stats ./Microbiome_Proteome_Summary/model_stats.qza \
        --i-baseline-stats ./Proteome_Summary/model_stats.qza \
        --o-visualization Proteome_paired-summary.qzv

#Microbiome-Metabolome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Metabolome.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir Metabolome_Summary

qiime mmvec summarize-paired \
        --i-model-stats ./Microbiome_Metabolome_Summary/model_stats.qza \
        --i-baseline-stats ./Metabolome_Summary/model_stats.qza \
        --o-visualization Metabolome_paired-summary.qzv

#Microbiome-Lipidome
qiime mmvec paired-omics \
        --i-microbes ./Microbiome.qza \
        --i-metabolites ./Lipidome.qza \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir Lipidome_Summary

qiime mmvec summarize-paired \
        --i-model-stats ./Microbiome_Lipidome_Summary/model_stats.qza \
        --i-baseline-stats ./Lipidome_Summary/model_stats.qza \
        --o-visualization Lipidome_paired-summary.qzv