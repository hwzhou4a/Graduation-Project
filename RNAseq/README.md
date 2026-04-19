This workflow included a comprehensive pipeline for RNA-seq analysis with human culture cell samples, starting from raw Fastq data to the downstream analysis including DEA, PCA and enrichment analysis.

Please pre-build all the required environment with demanded packages for the Linux step, and check with the required packages for the python and R scripts.

First with run_fastp.sh to preprocess and aggregate the raw Fastq data.

Second with run_rna_seq_pipeline.sh to perform genome alignment, multiple post-alignment QC and quantification.

(The first 2 step are under Linux environment, the current code is for HPC Platform@PKU, do adjust accordingly.)

Third with generate_gene_matrix.py to generate the Counts and TPM expression matrix on gene-level. (Utilize merge_isoforms.py instead if aiming for isoform-level study.)

Finally with Graduation_downstream_PCA_DEG_enrichment.r to further tidy up the expression matrix and perform DEA, PCA and enrichment analysis.
