process ASSEMBLY_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/assembly-stats:1.0.1--hc9558a2_3':
        'biocontainers/assembly-stats:1.0.1--hc9558a2_3' }"

    input:
    tuple val(meta), path(scaffolds)

    output:
    tuple val(meta), path('*.stats'), emit: stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assembly-stats \\
        -l 500 \\
        $args \\
        -t \\
        <(cat $scaffolds) \\
        > ${prefix}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assembly-stats: \$(assembly-stats -v)
    END_VERSIONS
    """
}
