Report for RNA-Seq transcript-level differential expression analysis
---

This document discusses the steps that were performed in the analysis pipeline.  It also describes the format of the output files and some brief interpretation.  For more detailed questions about interpretation of results, consult the documentation of the various tools.

The process used for this pipeline is a combination of the Kallisto and Sleuth tools created by the Pachter group (See references).  Kallisto is one of a generation of new "alignment-free" quantification methods that skip the computationally expensive alignment against the reference transcriptome.  Rather, these methods focus on comparing the sequence content of the raw sequence reads, which allows much faster quantification, allegedly without loss of result quality.

An important note is that the Sleuth process quantifies and tests differential expression at the *transcript* level, *not* gene level.  Where applicable, we have provided columns which map the transcript ID (in Ensembl's notation) to the common gene name.  

For further details, please see the references below or visit the [Kallisto](https://pachterlab.github.io/kallisto/about), and [Sleuth](https://pachterlab.github.io/sleuth/about) sites which provide additional information.

## Results:

We summarize some brief results in this section.  Full results can be found in the files, as described in the *Outputs* section.  


The following contrasts were performed, yielding the differentially expressed transcript counts shown below.  The threshold for significance was set such that the multiple-test corrected p-value is less than {{qval}}.  When referenced, the "top" genes refers to the {{num_hits}} genes with lowest multiple-test corrected p-value (q-value).   

|Experimental condition | Base condition| Upregulated | Downregulated | Result table | Heatmap of top DE transcripts |
|---|---|---|---|---|---|
{% for item in contrast_display %}
|{{item.exp_condition}} | {{item.base_condition}} | {{item.up_counts}}|{{item.down_counts}}| [Table](differential_expression/{{item.contrast_name}}.{{sleuth_output_suffix}}) | [Figure](differential_expression/{{item.contrast_name}}.{{top_heatmap_suffix}})|
{% endfor %}

## Outputs:

This section describes the contents of the delivered results.

#### Main results

The main results are contained in a zip-archive and should be downloaded an "unzipped" on your local computer.  It contains several sub-directories which contain files produced in each step of the pipeline.

- **QC**
    - This directory contains an interactive HTML-based QC report which summarizes read quality, pseudo-alignment information, and other metrics.  It was produced by MultiQC, and information can be found at <https://multiqc.info/>.
- **Quantifications**
    - Quantification tables, which give estimates of the mRNA quantity, as expressed in units of *transcripts per million* (TPM).  Files are tab-delimited.  These may be opened with your software of choice, including spreadsheet software such as Excel (note: <https://doi.org/10.1186/s13059-016-1044-7>).
- **Logs**
    - This contains logs and summaries produced by the various tools.  These can be used for troubleshooting, if necessary.
- **Differential expression results**

    Differential expression results performed by Sleuth are organized by the contrasts requested.  Thus, each folder, named by the contrast, contains the following files:

    - The differential expression results are summarized in a table in the file ending with `{{sleuth_output_suffix}}`.  It is saved in a tab-delimited text format.  Each row contains information about a particular transcript and the results of the statistical test performed..

    Although the details of the model and differential expression testing by Sleuth are different, it can helpful to recall the Student's t-test when thinking about interpretation of the columns.  The parameter estimates and tests are different, but the ideas are very similar.  In particular, Sleuth compares two models of transcript expression-- one model accounts for differences in treatment/phenotype in accordance with the contrast you are performing; the other model ignores that distinction.  The resulting test determines if the model including the effect of the treatment/phenotype "explains" the data better.  

    The columns and brief interpretations are:
      - **ensemble_gene_id**:  The gene identifier using the {{genome}} annotation.	
      - **gene_name**: The common gene name.	
      - **target_id**:  The transcript ID, using the {{genome}} annotation	
      - **pval**: The raw p-value of the statistical test.  Lower values indicate more evidence to reject the null hypothesis.	
      - **qval**: The "adjusted" p-value, which adjusts for the large number of statistical tests performed.  This addresses issues encountered in the "multiple-testing problem".  Corrections are based on the Benjamini-Hochberg procedure. 
      - **tpm_***: The mean TPM values from each group of samples.  Useful for assessing the direction of differential expression.	
      - **remaining columns**:  Various statistics and measures related to the statistical testing.  Left in the file for transparency, but requires thorough understanding of Sleuth modeling for interpretation.

  - Normalized expression table

    A table of normalized expressions, in *transcripts per million* (TPM) units is provided in `{{normalized_counts_suffix}}`.  It is saved in a tab-delimited text file.

- **Figures**

    Principal component analysis (PCA) of the estimated abundance matrix was performed for each contrast. The first two components, PC1 and PC2, are shown.  Each sample group is represented with a different color

    We provide boxplots of the top differentially expressed transcipts (if any) for quick reference.  Note that these figures are used for quick inspection, and it is not expected that such figures will be "publication-ready". 

    Finally, there is a heatmap showing the top differentially expressed transcripts.  It displays either the number of significant transcripts or the top {{num_hits}} by multiple-test corrected p-value.  


## Methods:

Input fastq-format files are processed using Kallisto software {{kallisto_version}}, which performs its "pseudo-alignment" against the {{genome}} reference genome.  Inferential variability in transcript quantification was assessed by performing N={{num_bootstraps}} bootstrap resamplings.  

Quantification estimates from Kallisto were subsequently used for differential transcript expression testing using Sleuth software.  For the basic contrast, we performed a likelihood ratio test, comparing the full model (which accounts for the experimental grouping) versus a reduced model which does not model the group membership.

Quality-control software included FastQC ({{fastqc_version}}) and MultiQC ({{multiqc_version}}).  Please see the respective references for interpretation of output information and figures.

The R `sessionInfo()` produced the following output.  We print here so that the same combination of packages/software may be recreated, if necessary.

```
{{session_info}}
```

## Inputs:
The inputs to the workflow were given as:

The inputs to the workflow were given as:

Samples and sequencing fastq-format files:

{% for obj in file_display %}
  - {{obj.sample_name}}
    - R1 fastq: {{obj.r1}}
    - R2 fastq: {{obj.r2}}
{% endfor %}

Sample annotations file: `{{annotations_file}}`

Parsed sample and condition table:

|Sample|Condition|
|---|---|
{% for item in annotation_objs %}
|{{item.name}} | {{item.condition}} |
{% endfor %}


## Version control:
To facilitate reproducible analyses, the analysis pipeline used to process the data is kept under git-based version control.  The repository for this workflow is at 

<{{git_repo}}>

and the commit version was {{git_commit}}.

This allows us to run the *exact* same pipeline at any later time, discarding any updates or changes in the process that may have been added. 


#### References:

[1] Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525–527 (2016), doi:10.1038/nbt.3519

[2] Harold J. Pimentel, Nicolas Bray, Suzette Puente, Páll Melsted and Lior Pachter, Differential analysis of RNA-Seq incorporating quantification uncertainty, Nature Methods (2017), advanced access http://dx.doi.org/10.1038/nmeth.4324
