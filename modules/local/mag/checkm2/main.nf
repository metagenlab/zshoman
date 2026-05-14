process CHECKM2 {

    tag "$meta.id"
    label 'process_high'

    conda "bioconda::checkm2"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("checkm2"), emit: checkm2

    script:
    """
    checkm2 predict \
        --input ${bins} \
        --output-directory checkm2 \
        --threads ${task.cpus}
    """
}
