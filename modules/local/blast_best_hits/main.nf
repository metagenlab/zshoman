process BLAST_BEST_HITS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d7126100b0eb7cb53dfb50291707ea8dda3b9738b76551ab73605d0acbe114b/data':
        'community.wave.seqera.io/library/pandas:2.3.3--5a902bf824a79745' }"

    input:
    tuple val(meta), path(hits)

    output:
    tuple val(meta), path("*_best_hits.csv"), emit: best_hits
    tuple val(meta), path("*_best_hits.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'blast_best_hits.py'
}
