#!/usr/bin/env nextflow

include { BBMAP_BBDUK as BBDUK_TRIM_ADAPTERS } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_FILTER_PHIX } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_QUALITY_FILTERING_SINGLE_END } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_QUALITY_FILTERING_PAIRED_END } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_SINGLE_END } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_PAIRED_END_PAIRS } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_PAIRED_END_SINGLETONS } from './modules/nf-core/bbmap/align/main'
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
include { BWA_INDEX } from './modules/nf-core/bwa/index/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { FILTERSAM } from './modules/local/filtersam/main'
include { NORMALIZE_COUNTS } from './modules/local/normalize_counts/main'
include { EGGNOGMAPPER } from './modules/nf-core/eggnogmapper/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_1} from './modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_2 } from './modules/nf-core/pigz/compress/main'


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

    trimmed_reads = BBDUK_TRIM_ADAPTERS(samples, params.references.adapters).reads
    phix_filtered_reads = BBDUK_FILTER_PHIX(trimmed_reads, params.references.phix).reads

    // here we split up paired-end and single-end samples, as we will need to
    // treat them separately because we keep the singletons for the paired-end
    // reads.
    phix_filtered_reads = phix_filtered_reads.branch(
        {
            single_end: it[0].single_end
            paired_end: !it[0].single_end
        })

    qf_reads_single = BBDUK_QUALITY_FILTERING_SINGLE_END(phix_filtered_reads.single_end, []).reads
    qf_reads_paired = BBDUK_QUALITY_FILTERING_PAIRED_END(phix_filtered_reads.paired_end, []).reads

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

    qf_reads_paired = qf_reads_paired.map( { new Tuple (it[0], it[1].sort(compare_files)) } )
    qf_reads_paired = qf_reads_paired.multiMap(
        {
            pairs: new Tuple (it[0], [it[1][0], it[1][1]])
            singletons: new Tuple (it[0] + [single_end:true], [it[1][2]])
        })

    host_index = BBMAP_INDEX_HOST(params.references.masked_human)
    single_end_reads = BBMAP_FILTER_HOST_SINGLE_END(qf_reads_single, host_index.index).reads
    paired_end_reads_pairs = BBMAP_FILTER_HOST_PAIRED_END_PAIRS(qf_reads_paired.pairs, host_index.index).reads
    paired_end_reads_singletons = BBMAP_FILTER_HOST_PAIRED_END_SINGLETONS(qf_reads_paired.singletons, host_index.index).reads

    paired_end_reads_merged = BBMAP_MERGE_PAIRS(paired_end_reads_pairs, false)

    // prepare a single channel with elements of the form
    // (meta, [*_1_unmerged.fastq.gz, *_2_unmerged.fastq.gz, *_merged.fastq.gz, "*_singletons.fastq.gz"])
    // for paired-end and (meta, *.fastq.gz) for single-end samples
    paired_end_reads = paired_end_reads_singletons.map({ new Tuple (it[0] + [single_end:false], it[1]) })
    paired_end_reads = paired_end_reads.join(paired_end_reads_merged.merged.join(paired_end_reads_merged.unmerged))
    paired_end_reads = paired_end_reads.map({ new Tuple (it[0], it[3] + [it[2]] + [it[1]]) })
    preprocessed_samples = single_end_reads.mix(paired_end_reads)


    /////////////////////////
    // Taxonomic Profiling //
    /////////////////////////

    motus_profiles = MOTUS_PROFILE(preprocessed_samples, params.motus_db).motus

    if (!params.skip_phanta) {
        phanta = PHANTA_PROFILE(paired_end_reads_pairs, params.phanta_db)
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


    //////////////////
    // Gene catalog //
    //////////////////

    // gather all amino acids and nucleotides from prokaryotic and eukaryotic genes
    // we first need to compress the output from MetaEuk
    eukaryotic_genes_aa = PIGZ_COMPRESS_1(eukaryotic_genes.faa).archive
    eukaryotic_genes_nt = PIGZ_COMPRESS_2(eukaryotic_genes.codon).archive
    amino_acids = prokaryotic_genes.amino_acid_fasta.mix(eukaryotic_genes_aa)
                    .collectFile( {row ->  [ "genes.faa.gz", row[1] ]} )
                    .map( { new Tuple([id: 'all'], it )} )
    nucleotides = prokaryotic_genes.nucleotide_fasta.mix(eukaryotic_genes_nt)
                    .collectFile( {row ->  [ "genes.fna.gz", row[1] ]} )
                    .map( { new Tuple([id: 'all'], it )} )

    gene_catalog_nt = CDHIT_CDHITEST(nucleotides).fasta

    headers = GET_HEADERS(gene_catalog_nt).headers
    gene_catalog_aa = SEQTK_SUBSEQ(amino_acids, headers.map( { it[1] } ).first()).sequences

    // Make sure we have a single fastq file for all reads per sample
    reads = CAT_FASTQ(preprocessed_samples, true).reads

    catalog_index = BWA_INDEX(gene_catalog_nt).index
    aligned_reads = BWA_MEM(reads, catalog_index.first(), gene_catalog_nt.first(), false).bam
    filtered_reads = FILTERSAM(aligned_reads).reads
    NORMALIZE_COUNTS(filtered_reads.join(motus_profiles))


    ///////////////////////////
    // Functional annotation //
    ///////////////////////////

    EGGNOGMAPPER(gene_catalog_aa, params.eggnog_db, params.eggnog_dbdir, new Tuple([:], params.eggnog_dmnd))
}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
