process CLASSIFY_4CAC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://metagenlab/4cac:1.0"

    input:
    tuple val(meta), path(contigs), path(graph), path(paths)

    output:
    tuple val(meta), path('4CAC_classification.fasta'), emit: classification
    tuple val(meta), path("*.4cac.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c $contigs > scaffolds.fasta

    python ${params.fourcac_dir}/classify_xgb.py \\
    -f scaffolds.fasta \\
    -p $task.cpus \\
    $args \\
    &> ${prefix}.4cac.log

    gunzip -c $graph > assembly_graph_with_scaffolds.gfa
    gunzip -c $paths > scaffolds.paths

    python ${params.fourcac_dir}/classify_4CAC.py \\
    --assembler metaSPAdes \\
    --asmdir . \\
    $args2 \\
    &>> ${prefix}.4cac.log
    """
}
