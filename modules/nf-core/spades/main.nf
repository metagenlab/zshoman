process SPADES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:4.0.0--h5fb382e_1' :
        'biocontainers/spades:4.0.0--h5fb382e_1' }"

    input:
    tuple val(meta), path(illumina), path(pacbio), path(nanopore)
    path yml
    path hmm

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')      , optional:true, emit: contigs
    tuple val(meta), path('*.transcripts.fa.gz')  , optional:true, emit: transcripts
    tuple val(meta), path('*.gene_clusters.fa.gz'), optional:true, emit: gene_clusters
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.scaffolds.paths.gz') , optional:true, emit: assembly_paths
    tuple val(meta), path('*.warnings.log')       , optional:true, emit: warnings
    tuple val(meta), path('*.spades.log')         , emit: log
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "--pe-1 1 ${illumina[0]} --pe-2 1 ${illumina[1]} --pe-m 1 ${illumina[2]} --pe-s 1 ${illumina[3]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$illumina_reads $pacbio_reads $nanopore_reads"
    """
    spades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $custom_hmms \\
        $reads \\
        -o ./
    mv spades.log ${prefix}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${prefix}.scaffolds.fa
        gzip -n ${prefix}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${prefix}.contigs.fa
        gzip -n ${prefix}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${prefix}.transcripts.fa
        gzip -n ${prefix}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
        gzip -n ${prefix}.assembly.gfa
    fi

    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
        gzip -n ${prefix}.gene_clusters.fa
    fi

    if [ -f scaffolds.paths ]; then
        mv scaffolds.paths ${prefix}.scaffolds.paths
        gzip -n ${prefix}.scaffolds.paths
    fi

    if [ -f warnings.log ]; then
        mv warnings.log ${prefix}.warnings.log
    fi

    rm -rf tmp K21 K33 K55 K77 misc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$illumina_reads $pacbio_reads $nanopore_reads"
    """
    echo "" | gzip > ${prefix}.scaffolds.fa.gz
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > ${prefix}.transcripts.fa.gz
    echo "" | gzip > ${prefix}.gene_clusters.fa.gz
    echo "" | gzip > ${prefix}.assembly.gfa.gz
    touch ${prefix}.spades.log
    touch ${prefix}.warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """
}
