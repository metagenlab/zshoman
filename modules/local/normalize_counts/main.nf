process NORMALIZE_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.22.1--py312hcfdcdd7_2':
        'biocontainers/pysam:0.22.1--py312hcfdcdd7_2' }"

    input:
    tuple val(meta), path(aligned_reads), path(motus_profile)

    output:
    tuple val(meta), path("*_genes_per_cell.csv"), emit: gene_counts

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'count_reads.py'
}
