process BBMAP_BBDUK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.06--h92535d8_1':
        'biocontainers/bbmap:39.06--h92535d8_1' }"

    input:
    tuple val(meta), path(reads)
    path contaminants
    val keep_singletons

    output:
    tuple val(meta), path("${prefix}{_[12],}.fastq.gz")                   , emit: reads
    tuple val(meta), path('*.log')                                        , emit: log
    path "versions.yml"                                                   , emit: versions
    tuple val(meta), path("${prefix}_singletons.fastq.gz"), optional: true, emit: singletons

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    outs = keep_singletons ? "outs=${prefix}_singletons.fastq.gz" : ""
    trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.fastq.gz out2=${prefix}_2.fastq.gz ${outs}"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        &> ${prefix}.bbduk.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_command  = meta.single_end ? "echo '' | gzip > ${prefix}.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastq.gz ; echo '' | gzip > ${prefix}_2.fastq.gz"
    """
    touch ${prefix}.bbduk.log
    $output_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
