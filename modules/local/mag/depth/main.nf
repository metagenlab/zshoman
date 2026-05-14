process MAG_DEPTH {

    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::metabat2=2.15"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths \
        --outputDepth ${meta.id}.depth.txt \
        ${bam}
    """
}
