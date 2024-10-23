#!/usr/bin/env nextflow

include { BBMAP_BBDUK as BBDUK_TRIM_ADAPTERS } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_FILTER_PHIX } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_QUALITY_FILTERING } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_SINGLETONS } from './modules/nf-core/bbmap/align/main'
include { BBMAP_INDEX as BBMAP_INDEX_HOST } from './modules/nf-core/bbmap/index/main'
include { BBMAP_BBMERGE as BBMAP_MERGE_PAIRS } from './modules/nf-core/bbmap/bbmerge/main'
include { MOTUS_PROFILE } from './modules/local/motus/main'
include { SPADES } from './modules/nf-core/spades/main'
include { FILTER_SCAFFOLDS } from './modules/local/filter_scaffolds/main'
include { ASSEMBLY_STATS } from './modules/local/assembly_stats/main'
include { PHANTA_PROFILE } from './modules/local/phanta/main'
include { CLASSIFY_4CAC } from './modules/local/4CAC/main'
include { PRODIGAL } from './modules/nf-core/prodigal/main'
include { METAEUK_EASYPREDICT } from './modules/nf-core/metaeuk/easypredict/main'
include { CDHIT_CDHITEST } from './modules/nf-core/cdhit/cdhitest/main'
include { GET_HEADERS } from './modules/local/seq_headers/main'
include { SEQTK_SUBSEQ } from './modules/nf-core/seqtk/subseq/main'
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { BWA_INDEX as BWA_INDEX_GC } from './modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_SAMPLES } from './modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_GC } from './modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_SAMPLES } from './modules/nf-core/bwa/mem/main'
include { FILTERSAM as FILTERSAM_GC } from './modules/local/filtersam/main'
include { FILTERSAM as FILTERSAM_SAMPLES } from './modules/local/filtersam/main'
include { NORMALIZE_COUNTS as NORMALIZE_COUNTS_GC } from './modules/local/normalize_counts/main'
include { NORMALIZE_COUNTS as NORMALIZE_COUNTS_SAMPLES } from './modules/local/normalize_counts/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_GC } from './modules/nf-core/eggnogmapper/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_SAMPLES } from './modules/nf-core/eggnogmapper/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_1 } from './modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_2 } from './modules/nf-core/pigz/compress/main'
include { CAT_CAT as CAT_AA } from './modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_NT } from './modules/nf-core/cat/cat/main'


workflow {
    // Read input file
    input_file = Channel.fromPath(params.input)
    samples = input_file
        .splitCsv(header: false, strip: false, limit: 1662)

    // Prepare metadata for samples to distinguish single from paired-end reads
    samples = samples.map( {
        row ->
            if (!row[2].strip()) {
                return new Tuple ([ id:row[0], single_end:true ], [ row[1] ])
            } else {
                return new Tuple ([ id:row[0], single_end:false ], [ row[1], row[2] ])
            }
        })


    ///////////////////
    // PREPROCESSING //
    ///////////////////

    trimmed_reads = BBDUK_TRIM_ADAPTERS(samples, params.references.adapters, false).reads
    phix_filtered_reads = BBDUK_FILTER_PHIX(trimmed_reads, params.references.phix, false).reads
    qf_reads = BBDUK_QUALITY_FILTERING(phix_filtered_reads, [], true).reads
    qf_singletons = BBDUK_QUALITY_FILTERING.out.singletons.map({ new Tuple (it[0] + [single_end:true], it[1]) })

    host_index = BBMAP_INDEX_HOST(params.references.masked_human)

    hf_reads = BBMAP_FILTER_HOST(qf_reads, host_index.index).reads
    hf_singletons = BBMAP_FILTER_HOST_SINGLETONS(qf_singletons, host_index.index).reads

    // We only merge for paired-end reads
    hf_reads_split = hf_reads.branch({
        single: it[0].single_end
        paired: !it[0].single_end
        })
    paired_end_reads_merged = BBMAP_MERGE_PAIRS(hf_reads_split.paired, false)

    // prepare a single channel preprocessed_samples with elements of the form
    // (meta, [*_1_unmerged.fastq.gz, *_2_unmerged.fastq.gz, *_merged.fastq.gz, "*.fastq.gz"])
    // for paired-end and (meta, *.fastq.gz) for single-end samples
    paired_end_reads = hf_singletons.map({ new Tuple (it[0] + [single_end:false], it[1]) })
    paired_end_reads = paired_end_reads.join(paired_end_reads_merged.merged.join(paired_end_reads_merged.unmerged))
    paired_end_reads = paired_end_reads.map({ new Tuple (it[0], it[3] + [it[2]] + [it[1]]) })
    preprocessed_samples = hf_reads_split.single.mix(paired_end_reads)
    // Make sure we have a single fastq file for all reads per sample
    reads = CAT_FASTQ(preprocessed_samples, true).reads

    /////////////////////////
    // Taxonomic Profiling //
    /////////////////////////

    motus_profiles = MOTUS_PROFILE(preprocessed_samples, params.motus_db).motus

    if (!params.skip_phanta) {
        // we cannot use the singletons nor the merged reads here so he use hf_reads instead.
        phanta = PHANTA_PROFILE(hf_reads, params.phanta_db)
    }


    ///////////////////////////////
    // Assembly and gene calling //
    ///////////////////////////////

    scaffolds = SPADES(preprocessed_samples.map({ new Tuple (it[0], it[1], [], []) }), [], []).scaffolds

    assembly_graph_and_paths = scaffolds.join(SPADES.out.gfa).join(SPADES.out.assembly_paths)
    contig_classification = CLASSIFY_4CAC(assembly_graph_and_paths).classification

    filtered_assembly = FILTER_SCAFFOLDS(scaffolds.join(contig_classification)).all_scaffolds
    assembly_stats = ASSEMBLY_STATS(filtered_assembly)

    prokaryotic_genes = PRODIGAL(FILTER_SCAFFOLDS.out.prok_scaffolds, "gff")
    eukaryotic_genes = METAEUK_EASYPREDICT(FILTER_SCAFFOLDS.out.euk_scaffolds, params.metaeuk_db)


    // we need to compress the output from MetaEuk so that both eukaryotic
    // and prokaryotic genes are compressed
    eukaryotic_genes_aa = PIGZ_COMPRESS_1(eukaryotic_genes.faa).archive
    eukaryotic_genes_nt = PIGZ_COMPRESS_2(eukaryotic_genes.codon).archive

    //////////////////
    // Gene catalog //
    //////////////////

    // gather all amino acids and nucleotides from prokaryotic and eukaryotic genes
    all_amino_acids = prokaryotic_genes.amino_acid_fasta.mix(eukaryotic_genes_aa)
                    .collectFile( {row ->  [ "genes.faa.gz", row[1] ]} )
                    .map( { new Tuple([id: 'all'], it )} )
    all_nucleotides = prokaryotic_genes.nucleotide_fasta.mix(eukaryotic_genes_nt)
                    .collectFile( {row ->  [ "genes.fna.gz", row[1] ]} )
                    .map( { new Tuple([id: 'all'], it )} )

    gene_catalog_nt = CDHIT_CDHITEST(all_nucleotides).fasta

    headers = GET_HEADERS(gene_catalog_nt).headers
    gene_catalog_aa = SEQTK_SUBSEQ(all_amino_acids, headers.map( { it[1] } ).first()).sequences

    catalog_index = BWA_INDEX_GC(gene_catalog_nt).index
    aligned_reads = BWA_MEM_GC(reads.combine(catalog_index).map( { new Tuple(it[0], it[1], it[3]) } ), false).bam
    filtered_reads = FILTERSAM_GC(aligned_reads).reads
    NORMALIZE_COUNTS_GC(filtered_reads.join(motus_profiles))


    ///////////////////////////
    // Functional annotation //
    ///////////////////////////

    EGGNOGMAPPER_GC(gene_catalog_aa, params.eggnog_db, params.eggnog_dbdir, new Tuple([:], params.eggnog_dmnd))

    /////////////////////////////////////
    // Genes counts individual samples //
    /////////////////////////////////////

    // Let's gather the prokaryotic genes and eukaryotic genes together
    amino_acids = CAT_AA(prokaryotic_genes.amino_acid_fasta.mix(eukaryotic_genes_aa).groupTuple()).file_out
    nucleotides = CAT_NT(prokaryotic_genes.nucleotide_fasta.mix(eukaryotic_genes_nt).groupTuple()).file_out

    catalog_index = BWA_INDEX_SAMPLES(nucleotides).index
    reads_index = reads.join(catalog_index)
    aligned_reads = BWA_MEM_SAMPLES(reads_index, false).bam
    filtered_reads = FILTERSAM_SAMPLES(aligned_reads).reads
    NORMALIZE_COUNTS_SAMPLES(filtered_reads.join(motus_profiles))


    ///////////////////////////////////
    // Functional individual samples //
    ///////////////////////////////////

    EGGNOGMAPPER_SAMPLES(amino_acids, params.eggnog_db, params.eggnog_dbdir, new Tuple([:], params.eggnog_dmnd))

}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
