#!/usr/bin/env nextflow

include { BBMAP_BBDUK as BBDUK_TRIM_ADAPTERS } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_FILTER_PHIX } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_QUALITY_FILTERING } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_SINGLE_END } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_PAIRED_END_PAIRS } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_PAIRED_END_SINGLETONS } from './modules/nf-core/bbmap/align/main'
include { BBMAP_INDEX as BBMAP_INDEX_HOST } from './modules/nf-core/bbmap/index/main'
include { BBMAP_BBMERGE as BBMAP_MERGE_PAIRS } from './modules/nf-core/bbmap/bbmerge/main'
include { MOTUS_PROFILE } from './modules/local/motus/main'


process assemble_paired_end {
    cpus = 20
    memory '250 GB'
    input:
        tuple val(sample_name), path(merged), path(paired_1), path(paired_2), path(singleton_reads)

    output:
        tuple val(sample_name), path("assembly")

    script:
        """
        mkdir assembly
        metaspades.py -t ${task.cpus} -m ${task.memory.toGiga()} --only-assembler --pe-1 1 $paired_1 --pe-2 1 $paired_2 --pe-m 1 $merged --pe-s 1 $singleton_reads -o assembly
        """
    }

process filter_short_contigs {
    cpus = 1
    input:
        tuple val(sample_name), path(assembly)

    output:
        tuple val(sample_name), path("${sample_name}.scaffolds.min500.fasta")

    script:
        """
        python $params.scaffold_filtering_script $sample_name scaffolds $assembly/scaffolds.fasta .
        """
    }

process get_assembly_stats {
    cpus = 1
    input:
        tuple val(sample_name), path(filtered_assembly)

    output:
        tuple val(sample_name), path("${sample_name}.stats")

    script:
        """
        assembly-stats -l 500 -t <(cat $filtered_assembly) > ${sample_name}.stats
        """
    }

process call_genes {
    cpus = 1
    input:
        tuple val(sample_name), path(filtered_assembly)

    output:
        tuple val(sample_name), path("${sample_name}.faa"), path("${sample_name}.fna"), path("${sample_name}.gff")

    script:
        """
        prodigal -a ${sample_name}.faa -d ${sample_name}.fna -f gff \
        -o ${sample_name}.gff -c -q -p meta -i $filtered_assembly
        """
    }

process make_gene_catalog {
    cpus = 20
    input:
        path(amino_acids)
        path(nucleotides)

    output:
        path("cdhit9590")

    script:
        """
        mkdir cdhit9590
        cd-hit-est -i $nucleotides -o cdhit9590/gene_catalog_cdhit9590.fasta \
        -c 0.95 -T 20 -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0

        grep "^>" cdhit9590/gene_catalog_cdhit9590.fasta | \
        cut -f 2 -d ">" | \
        cut -f 1 -d " " > cdhit9590/cdhit9590.headers
        seqtk subseq $amino_acids cdhit9590/cdhit9590.headers > cdhit9590/gene_catalog_cdhit9590.faa
        """
    }

process align_reads {
    cpus = 20
    input:
        path(gene_catalog)
        tuple val(sample_name), path(merged), path(paired_1), path(paired_2), path(singleton_reads)

    output:
        tuple val(sample_name), path("alignments")

    script:
        """
        mkdir alignments
        bwa index $gene_catalog/gene_catalog_cdhit9590.fasta

        bwa mem -a -t 20 $gene_catalog/gene_catalog_cdhit9590.fasta $paired_1 \
        | samtools view -F 4 -bh - > alignments/${sample_name}_r1.bam

        bwa mem -a -t 20 $gene_catalog/gene_catalog_cdhit9590.fasta $paired_2 \
        | samtools view -F 4 -bh - > alignments/${sample_name}_r2.bam

        bwa mem -a -t 20 $gene_catalog/gene_catalog_cdhit9590.fasta $merged \
        | samtools view -F 4 -bh - > alignments/${sample_name}_merged.bam

        bwa mem -a -t 20 $gene_catalog/gene_catalog_cdhit9590.fasta $singleton_reads \
        | samtools view -F 4 -bh - > alignments/${sample_name}_singleton.bam
        """
    }

process filter_reads {
    cpus = 5
    input:
        tuple val(sample_name), path(alignments)

    output:
        tuple val(sample_name), path("${sample_name}_filtered.bam")

    script:
        """
        samtools merge ${sample_name}.bam ${alignments}/${sample_name}_merged.bam ${alignments}/${sample_name}_r1.bam ${alignments}/${sample_name}_r2.bam ${alignments}/${sample_name}_singleton.bam
        samtools view -F 256 -e '(qlen-sclen)>45' -O BAM -o filtered_primary.bam ${sample_name}.bam
        filtersam -i 95 -p 5 -o ${sample_name}_filtered.bam filtered_primary.bam
        """
    }

process count_reads {
    cpus = 1
    input:
        path(read_counter)
        tuple val(sample_name), path(filtered_reads), path(motus)

    output:
        tuple val(sample_name), path("${sample_name}_read_counts.csv")

    script:
        """
        python $read_counter $filtered_reads $motus ${sample_name}_read_counts.csv
        """
    }

process rpsblast_COG {
    cpus = 4
    input:
        tuple (file(cog_db), file(seq))

    output:
        path result_file

    script:
        n = seq.name
        result_file = "${n}.tab"
        """
        rpsblast -db cog/cog_db -query $seq -outfmt 6 -evalue 0.001 \
                -num_threads 4 > ${result_file}
        """
}

process merge_cogs {
    cpus = 1
    input:
        path(cog_db)
        path(cog_hit_files)

    output:
        path("cogs.csv")

    script:
        cdd_to_cog = "${cog_db}/cdd_to_cog"
        """
        python $params.merge_cogs $cdd_to_cog $cog_hit_files
        """
}

process download_cog_definitions {
    cpus = 1
    publishDir params.publish_dir

    output:
        tuple (path("cog-20.def.tab"), path("fun-20.tab"))

    script:
        """
        wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
        wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab
        """
}

process execute_kofamscan {

    cpus = 4

    input:
        tuple (path(ko_db), path(seq))

    output:
        path '*tab'

    script:
        n = seq.name

        """
        exec_annotation ${n} -p ${ko_db}/profiles/prokaryote.hal \
            -k ${ko_db}/ko_list --cpu 4 -o ${n}.tab
        """
}

process get_ko_table {
    publishDir params.publish_dir
    cpus = 1
    input:
       path(ko_hit_files)

    output:
        path("kos.csv")

    script:
        """
        python $params.prepare_ko_table $ko_hit_files
        """
}

process get_ko_mappings {
    publishDir params.publish_dir
    cpus = 1

    output:
        tuple( path("ko_module_classes.csv"),
               path("ko_modules.csv"),
               path("ko_pathways.csv"),
               path("ko_definitions.csv"),
               path("ko_to_pathway.csv"),
               path("ko_to_module.csv")
            )

    script:
        """
        python $params.get_ko_mappings
        """
}

process collect_motus_outputs {
    cpus = 1
    publishDir params.publish_dir

    input:
        tuple val(sample_name), path(motus_file)

    output:
        path(motus_file)

    script:
        """
        """
}

process collect_assemblies {
    cpus = 1
    publishDir params.publish_dir

    input:
        tuple val(sample_name), path(assembly)

    output:
        path(assembly)

    script:
        """
        """
}

process collect_counts {
    cpus = 1
    publishDir params.publish_dir

    input:
        tuple val(sample_name), path(counts)

    output:
        path(counts)

    script:
        """
        """
}

process collect_cogs {
    cpus = 1
    publishDir params.publish_dir

    input:
        path(cogs)

    output:
        path(cogs)

    script:
        """
        """
}

process collect_gene_catalog {
    cpus = 1
    publishDir params.publish_dir

    input:
        path(gene_catalog)

    output:
        path(gene_catalog)

    script:
        """
        """
}

workflow {
    input_file = Channel.fromPath(params.input)
    samples = input_file
        .splitCsv(header: false, strip: false, limit: 1662)

    // prepare metadata for samples to distinguish single from paired-end reads
    samples = samples.map( {
        row ->
            if (!row[2].strip()) {
                return new Tuple ([ id:row[0], single_end:true ], [ row[1] ])
            } else {
                return new Tuple ([ id:row[0], single_end:false ], [ row[1], row[2] ])
            }
        })

    trimmed_reads = BBDUK_TRIM_ADAPTERS(samples, params.references.adapters).reads
    phix_filtered_reads = BBDUK_FILTER_PHIX(trimmed_reads, params.references.phix).reads
    quality_filtered_reads = BBDUK_QUALITY_FILTERING(phix_filtered_reads, []).reads

    // here we split up paired-end and single-end samples, as we will need to
    // treat them separately because we kept the singletons for the paired-end
    // reads.
    quality_filtered_reads = quality_filtered_reads.branch(
        {
            single_end: it[0].single_end
            paired_end: !it[0].single_end
        })

    // Sort the files for paired-end reads -> (*_1, *_2, singletons)
    def compare_files = { a, b ->
        if (a.toString().endsWith("_1.fastq.gz")) {
            return -1
        } else if (a.toString().endsWith("singleton_reads_qf.fastq.gz")) {
            return 1
        } else if (b.toString().endsWith("_1.fastq.gz")) {
            return 1
        } else {
            return -1
        }
    }

    paired_end = quality_filtered_reads.paired_end.map( { new Tuple (it[0], it[1].sort(compare_files)) } )
    paired_end = paired_end.multiMap(
        {
            pairs: new Tuple (it[0], [it[1][0], it[1][1]])
            singletons: new Tuple (it[0] + [single_end:true], [it[1][2]])
        })

    host_index = BBMAP_INDEX_HOST(params.references.masked_human)
    single_end_reads = BBMAP_FILTER_HOST_SINGLE_END(quality_filtered_reads.single_end, host_index.index).reads
    paired_end_reads_pairs = BBMAP_FILTER_HOST_PAIRED_END_PAIRS(paired_end.pairs, host_index.index).reads
    paired_end_reads_singletons = BBMAP_FILTER_HOST_PAIRED_END_SINGLETONS(paired_end.singletons, host_index.index).reads

    paired_end_reads_merged = BBMAP_MERGE_PAIRS(paired_end_reads_pairs, false)

    // prepare a single channel with elements of the form
    // (meta, [*_1_unmerged.fastq.gz, *_2_unmerged.fastq.gz, *_merged.fastq.gz, "*_singletons.fastq.gz"])
    // for paired-end and (meta, *.fastq.gz) for single-end samples
    paired_end_reads = paired_end_reads_singletons.map({ new Tuple (it[0] + [single_end:false], it[1]) })
    paired_end_reads = paired_end_reads.join(paired_end_reads_merged.merged.join(paired_end_reads_merged.unmerged))
    paired_end_reads = paired_end_reads.map({ new Tuple (it[0], it[3] + [it[2]] + [it[1]]) })
    preprocessed_samples = single_end_reads.mix(paired_end_reads)

    MOTUS_PROFILE(preprocessed_samples)
    /*
    assembly = assemble_paired_end(preprocessed_paired)

    filtered_assembly = filter_short_contigs(assembly)

    assembly_stats = get_assembly_stats(filtered_assembly)

    genes = call_genes(filtered_assembly)

    amino_acids = genes.collectFile( {row ->  [ "genes.faa", row[1] ]} )
    nucleotides = genes.collectFile( {row ->  [ "genes.fna", row[2] ]} )
    gene_catalog = make_gene_catalog(amino_acids, nucleotides)

    aligned_reads = align_reads(gene_catalog, preprocessed_paired)

    filtered_reads = filter_reads(aligned_reads)

    reads_and_motus = filtered_reads.join(motus_paired_end.out)
    counts = count_reads(params.count_reads, reads_and_motus)

    gene_catalog_aa = gene_catalog.first() + "/gene_catalog_cdhit9590.faa"
    split_aa_seqs = gene_catalog_aa.splitFasta( by: 300, file: "chunk_" )

    // COGS
    cog_db = Channel.fromPath("$params.cog_db", type: "dir")
    cogs = rpsblast_COG(cog_db.combine(split_aa_seqs))
    merged_cogs = merge_cogs(cog_db, cogs.collect())
    cog_def_files = download_cog_definitions()

    // KOs
    ko_db = Channel.fromPath("$params.ko_db", type: "dir")
    ko_hits = execute_kofamscan(ko_db.combine(split_aa_seqs))
    get_ko_table(ko_hits.collect())
    get_ko_mappings()

    // Collection should be replaced by publishing directly from
    // the processes
    collect_motus_outputs(motus_paired_end.out)
    collect_assemblies(filtered_assembly)
    collect_counts(counts)
    collect_cogs(merged_cogs)
    collect_gene_catalog(gene_catalog_aa)
    */
}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
