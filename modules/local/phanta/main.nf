process PHANTA_PROFILE {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "docker://metagenlab/phanta:1.1"

    input:
    tuple val(meta), path(reads)
    path(phanta_db)

    output:
    tuple val(meta), path('final_merged_outputs/*.txt'), emit: classification
    tuple val(meta), path("*.phanta.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffixes = meta.single_end ? "--fwd _filtered" : "--fwd host_filtered_1 --rev host_filtered_2"
    def input = meta.single_end ? "-s ${reads[0]}" : "-f ${reads[0]} -r ${reads[1]} -s ${reads[2]},${reads[3]}"
    """
    python ${params.phanta_dir}/run_phanta.py \\
    -i ./ \\
    -p $params.phanta_dir \\
    -d $phanta_db \\
    -o ./ \\
    -t $task.cpus \\
    --run \\
    $suffixes \\
    $args \\
    &> ${prefix}.phanta.log
    """
}
