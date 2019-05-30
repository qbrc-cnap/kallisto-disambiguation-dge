library("sleuth")
library("ggplot2")

args <- commandArgs(TRUE)
DESIGN_FILE <- args[1] # path to a design matrix file- see below for format
TRANSCRIPT_TO_GENE_MAPPING_FILE <- args[2] # file mapping the transcript to gene names.  Format below
OUTPUT_FILE <- args[3] # filepath to write the sleuth differential results to
NC_FILE <- args[4] # filepath to write the table of normalized counts
PLOTS_DIR <- args[5] # path to a directory for all the plots
PLOT_N_TOP_HITS <- args[6] # integer giving the max number of transcripts to plot
QVAL_THRESHOLD <- args[7] # hits are considered significant if less than this value

# cast to appropriate types:
QVAL_THRESHOLD <- as.double(QVAL_THRESHOLD)
PLOT_N_TOP_HITS <- as.integer(PLOT_N_TOP_HITS)

# 3 column file-- condition, path to kallisto results, and sample ID.  need to be named 'path', 'condition', and 'sample'
# 'path', 'sample' names are required by sleuth, and 'condition' is required for our model equation below
# The file has headers and the name of the condition column should be 'condition' for the model formula below
s2c <- read.table(DESIGN_FILE, header=T, stringsAsFactors=F)

# map the ENST IDs to common gene names.  The files used have (in order):
# - ensembl gene ID
# - ensembl transcript ID
# - common gene name
# sleuth requires that this has a column named 'target_id' which is the identifier used in kallisto (the ensembl transcript ID in our case)
transcript_to_gene_map <- read.table(TRANSCRIPT_TO_GENE_MAPPING_FILE, header=T, sep='\t')
colnames(transcript_to_gene_map) <- c('ensemble_gene_id', 'target_id','gene_name')

so <- sleuth_prep(s2c, target_mapping = transcript_to_gene_map, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, formula= ~ condition, fit_name='full')
so <- sleuth_fit(so, formula = ~1, fit_name='reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
results_table <-sleuth_results(so, 'reduced:full', test_type = 'lrt')
write.table(results_table, OUTPUT_FILE, sep='\t', quote=F, row.names=F)

# get the normalized counts:
nc <- kallisto_table(so, normalized = TRUE, include_covariates = TRUE)
write.table(nc, NC_FILE, sep='\t', quote=F, row.names=F)

# create a directory for for the figures
dir.create(PLOTS_DIR)

# plot PCA:
pca_filepath <- sprintf('%s/pca.png', PLOTS_DIR)
png(pca_filepath)
plot_pca(so, color_by='condition')
dev.off()

# for the top N transcripts (if significant), produce a boxplot for each
significant_hits_table <- na.omit(results_table[results_table[,'qval'] <= QVAL_THRESHOLD,])
N <- if (dim(significant_hits_table)[1] < PLOT_N_TOP_HITS) dim(significant_hits_table)[1] else PLOT_N_TOP_HITS
if (N == 0){
	print(sprintf('No significant hits at %s threshold.', QVAL_THRESHOLD))
}else{
	print(sprintf('Found %s significant hits', N))
	for(i in 1:N){
		transcript_id = significant_hits_table[i, 'target_id']
		outfile = sprintf('%s/%s.png', PLOTS_DIR, transcript_id)
		gg <- plot_bootstrap(so, transcript_id, units = "est_counts", color_by = "condition")
		ggsave(filename=outfile, gg)
	}
}

