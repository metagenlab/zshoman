process GET_HEADERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path('*.headers'), emit: headers

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    grep "^>" $alignment | \
    cut -f 2 -d ">" | \
    cut -f 1 -d " " > ${prefix}.headers
    """
}
