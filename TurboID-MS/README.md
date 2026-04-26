This workflow included a comprehensive pipeline for TurboID-MS analysis with human culture cell samples after finishing the database searching phase.

Please check with the required packages for the R script.

The pipeline dealt with the data as below,

1. Collect all the numeric quantity calculate and draw the distribution
2. For all N/A value fill with 1
3. For each condition perform normalization within group, normalized by median to the largest one
4. Perform log2 transformation to the matrix generate in step 3
5. Perform PCA based on log_matrix
6. With limma package, compare each experimental condition to blank control; for each comparation, pre-filter the cases that 6 values are all less than the last five percent of pre-defined distribution in step 1; generate a Differential table, a volcano plot, and perform GO enrichment analysis based on the upregulated proteins in experimental conditions as output.
7. With limma package, compare each variants to WT condition; for each comparation, pre-filter the cases that 6 values are all less than the last five percent of pre-defined distribution in step 1; generate a differential table and a volcano plot as output.
8. Gather upregulated ones together, and generate a Venn plot to identify common upregulated proteins; Gather downregulated ones together, and generate a Venn plot to identify common downregulated proteins.
9. For common upregulated ones, try GO/KEGG enrichment and PPI, for common downregulated ones, try GO/KEGG enrichment and PPI.