process FILTER_SCAFFOLDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(scaffolds), path(classification)

    output:
    tuple val(meta), path("${meta.id}.scaffolds.min500.fasta"), emit: all_scaffolds
    tuple val(meta), path("${meta.id}.scaffolds.min500_eukaryotes.fasta"), emit: euk_scaffolds
    tuple val(meta), path("${meta.id}.scaffolds.min500_prokaryotes.fasta"), emit: prok_scaffolds

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'scaffold_filter.py'
}
