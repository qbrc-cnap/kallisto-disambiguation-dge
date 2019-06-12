task sleuth_dge {
    Array[File] abundance_h5_files
    File annotations
    File transcript_to_gene_mapping
    String base_group
    String experimental_group
    String versus_sep
    String normalized_counts_suffix
    String sleuth_output_suffix
    String pca_suffix
    String top_heatmap_suffix
    Float qval_threshold
    Int max_transcripts


    Int disk_size = 100
        
    String contrast_name = experimental_group + versus_sep + base_group
    String sleuth_output = contrast_name + "." + sleuth_output_suffix
    String sleuth_normalized_counts = contrast_name + "." + normalized_counts_suffix
    String output_figures_dir = contrast_name + "_figures"
    String pca_output = contrast_name + "." + pca_suffix
    String heatmap_output = contrast_name + "." + top_heatmap_suffix
    String normalized_counts = "nc.tsv"
    String sleuth_annotations = 'sleuth_annotationts.tsv'

    command {
        echo "Moving abundance.h5 files..."

        /usr/bin/python3 /opt/software/move_files.py ${sep=" " abundance_h5_files}

        echo "Completed moving files..."
        /usr/bin/python3 /opt/software/create_sleuth_annotation_file.py -i ${annotations} -o ${sleuth_annotations}
        echo "Created annotations for sleuth process..."
        
        Rscript /opt/software/sleuth.R \
            ${sleuth_annotations} \
            ${transcript_to_gene_mapping} \
            ${sleuth_output} \
            ${normalized_counts} \
            ${output_figures_dir} \
            ${max_transcripts} \
            ${qval_threshold} \
            ${base_group} \
            ${experimental_group} \
            ${pca_output} \
            ${heatmap_output}

        echo "Completed sleuth"

        /usr/bin/python3 /opt/software/pivot_normalized_counts.py \
            -i ${normalized_counts} \
            -o ${sleuth_normalized_counts}
    }

    output {
        File norm_counts = "${sleuth_normalized_counts}"
        File sleuth_results = "${sleuth_output}"
        Array[File] sleuth_plots = glob("${output_figures_dir}/*.png")
    }

    runtime {
        docker: "docker.io/blawney/kallisto:v0.0.1"
        cpu: 4
        memory: "5 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
