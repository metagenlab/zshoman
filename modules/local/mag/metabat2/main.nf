process METABAT2 {

    tag "$meta.id"
    label 'process_high'

    conda "bioconda::metabat2=2.15"

    input:
    tuple val(meta), path(scaffolds), path(depth)

    output:
    tuple val(meta), path("bins"), emit: bins

    script:
    """
    mkdir -p bins

    metabat2 \
        -i ${scaffolds} \
        -a ${depth} \
        -o bins/${meta.id}

    """
}
