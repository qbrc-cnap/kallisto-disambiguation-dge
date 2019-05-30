task sleuth_dge {
    Array[File] abundance_h5_files
    Array[File] abundance_tsv_files
    Array[File] run_info_files
    File annotations
    File transcript_to_gene_mapping

    Int disk_size = 100
    String sleuth_annotations = "sleuth_annotations.tsv"
    String sleuth_output = "sleuth_differential_expression_results.tsv"
    String sleuth_norm_counts = "normalized_counts.tsv"
    String plot_dir = "sleuth_plots"
    Int max_transcripts = 30 # the maximum number of transcripts to plot
    Float qval_threshold = 0.05

    command {
        python3 /opt/software/scripts/move_files.py ${sep=" " abundance_h5_files}
        python3 /opt/software/scripts/move_files.py ${sep=" " abundance_tsv_files}
        python3 /opt/software/scripts/move_files.py ${sep=" " run_info_files}

        python3 /opt/software/scripts/create_sleuth_annotation_file.py -o ${sleuth_annotations}

        Rscript /opt/software/scripts/sleuth.R \
            ${sleuth_annotations} \
            ${transcript_to_gene_mapping} \
            ${sleuth_output} \
            ${sleuth_norm_counts} \
            ${plot_dir} \
            ${max_transcripts} \
            ${qval_threshold}
    }

    output {
        File norm_counts = "${sleuth_norm_counts}"
        File sleuth_results = "${sleuth_output}"
        Array[File] sleuth_plots = glob("${plot_dir}/*")
    }

    runtime {
        zones: "us-east4-c"
        docker: "docker.io/blawney/kallisto:v0.0.1"
        cpu: 4
        memory: "3 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}