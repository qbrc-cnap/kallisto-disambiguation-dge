import "single_sample_kallisto.wdl" as single_sample_kallisto
import "multiqc.wdl" as multiqc
import "fastqc.wdl" as fastqc
import "sleuth.wdl" as sleuth
#import "figures.wdl" as figures
#import "report.wdl" as reporting


workflow KallistoAndSleuthWorkflow{
    # This workflow is a 'super' workflow that parallelizes
    # RNA-seq analysis over multiple samples

    Array[File] r1_files
    Array[File] r2_files
    File sample_annotations
    Array[String] base_conditions
    Array[String] experimental_conditions
    String genome
    File kallisto_index_path
    File transcript_to_gene_mapping
    Int kallisto_bootstraps = 500
    String output_zip_name
    String git_repo_url
    String git_commit_hash

    Float padj_threshold = 0.01
    Float lfc_threshold = 1.5

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)
    Array[Pair[String, String]] contrast_pairs = zip(base_conditions, experimental_conditions)


    scatter(item in fastq_pairs){

        call fastqc.run_fastqc as fastqc_for_read1 {
            input:
                fastq = item.left
        }

        call fastqc.run_fastqc as fastqc_for_read2 {
            input:
                fastq = item.right
        }

        call single_sample_kallisto.SingleSampleKallistoWorkflow as single_sample_process{
            input:
                r1_fastq = item.left,
                r2_fastq = item.right,
                kallisto_index_path = kallisto_index_path,
                kallisto_bootstraps = kallisto_bootstraps
        }
    }

    call sleuth.sleuth_dge as sleuth_dge {
        input:
            abundance_h5_files = single_sample_process.abundance_h5,
            abundance_tsv_files = single_sample_process.abundance_tsv,
            run_info_files = single_sample_process.run_info,
            annotations = sample_annotations,
            transcript_to_gene_mapping = transcript_to_gene_mapping
    }

    call multiqc.create_qc as experimental_qc {
        input:
            kallisto_stdout = single_sample_process.kallisto_stdout,
            r1_fastqc_zips = fastqc_for_read1.fastqc_zip,
            r2_fastqc_zips = fastqc_for_read2.fastqc_zip
    }

    output {
        File sleuth_results = sleuth_dge.sleuth_results
        File norm_counts = sleuth_dge.norm_counts
        Array[File] plots = sleuth_dge.sleuth_plots
    }

    meta {
        workflow_title : "Kallisto + Sleuth RNA-Seq differential expression"
        workflow_short_description : "For determining differential expression of transcript abundances using Kallisto and Sleuth"
        workflow_long_description : "Use this workflow for performing pseudo-alignments with Kallisto, which produces estimated transcript abundances.  Differential expression is performed using the companion Sleuth tool."
    }
}