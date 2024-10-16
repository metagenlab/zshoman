process FILTERSAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://metagenlab/filtersam:1.0"

    input:
    tuple val(meta), path(aligned_reads)

    output:
    tuple val(meta), path('*_filtered.bam'), emit: reads
    tuple val(meta), path("*.filtersam.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filtersam \\
    $args \\
    -p $task.cpus \\
    -o ${prefix}_filtered.bam \\
    $aligned_reads \\
    &> ${prefix}.filtersam.log
    """
}
