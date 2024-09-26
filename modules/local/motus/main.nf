process MOTUS_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/motus:3.1.0--pyhdfd78af_0 ':
        'biocontainers/motus:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.motus'), emit: motus
    tuple val(meta), path("*.motus.log"), emit: log
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-s ${reads[0]}" : "-f ${reads[0]} -r ${reads[1]} -s ${reads[2]},${reads[3]}"
    """
    motus downloadDB  # does nothing if DB is already downloaded
    motus profile \\
        $input \\
        -n $meta.id \\
        $args \\
        -o ${prefix}.motus \\
        -t $task.cpus\\
        &> ${prefix}.motus.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        motus: \$(motus --version)
    END_VERSIONS
    """
}
