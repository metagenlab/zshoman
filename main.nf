#!/usr/bin/env nextflow

include { ASSEMBLY_STATS } from './modules/local/assembly_stats/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST } from './modules/nf-core/bbmap/align/main'
include { BBMAP_ALIGN as BBMAP_FILTER_HOST_SINGLETONS } from './modules/nf-core/bbmap/align/main'
include { BBMAP_BBDUK as BBDUK_FILTER_PHIX } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_QUALITY_FILTERING } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBDUK as BBDUK_TRIM_ADAPTERS } from './modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBMERGE as BBMAP_MERGE_PAIRS } from './modules/nf-core/bbmap/bbmerge/main'
include { BBMAP_INDEX as BBMAP_INDEX_HOST } from './modules/nf-core/bbmap/index/main'
include { BLAST_BLASTP } from './modules/nf-core/blast/blastp/main'
include { BLAST_MAKEBLASTDB } from './modules/nf-core/blast/makeblastdb/main'
include { BWA_INDEX as BWA_INDEX_GC } from './modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_SAMPLES } from './modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_MEM_GC } from './modules/nf-core/bwa/mem/main'
include { BWA_MEM as BWA_MEM_SAMPLES } from './modules/nf-core/bwa/mem/main'
include { CAT_CAT as CAT_AA } from './modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_NT } from './modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_R1 } from './modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_R2 } from './modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_SINGLE_END } from './modules/nf-core/cat/cat/main'
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { CDHIT_CDHITEST } from './modules/nf-core/cdhit/cdhitest/main'
include { CLASSIFY_4CAC } from './modules/local/4CAC/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_GC } from './modules/nf-core/eggnogmapper/main'
include { EGGNOGMAPPER as EGGNOGMAPPER_SAMPLES } from './modules/nf-core/eggnogmapper/main'
include { FILTER_SCAFFOLDS } from './modules/local/filter_scaffolds/main'
include { FILTERSAM as FILTERSAM_GC } from './modules/local/filtersam/main'
include { FILTERSAM as FILTERSAM_SAMPLES } from './modules/local/filtersam/main'
include { GET_HEADERS } from './modules/local/seq_headers/main'
include { METAEUK_EASYPREDICT } from './modules/nf-core/metaeuk/easypredict/main'
include { MOTUS_PROFILE } from './modules/local/motus/main'
include { NORMALIZE_COUNTS as NORMALIZE_COUNTS_GC } from './modules/local/normalize_counts/main'
include { NORMALIZE_COUNTS as NORMALIZE_COUNTS_SAMPLES } from './modules/local/normalize_counts/main'
include { PHANTA_PROFILE } from './modules/local/phanta/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_1 } from './modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_COMPRESS_2 } from './modules/nf-core/pigz/compress/main'
include { PRODIGAL } from './modules/nf-core/prodigal/main'
include { SEQTK_SUBSEQ } from './modules/nf-core/seqtk/subseq/main'
include { SPADES } from './modules/nf-core/spades/main'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

import java.nio.file.Paths
import java.nio.file.Files

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


workflow {
    outdir_abs = Paths.get(params.outdir).toAbsolutePath().toString()
    // Create a new channel of metadata from the sample sheet passed to the pipeline through the --input parameter
    samples = samplesheetToList(params.input, "assets/schema_input.json")

    // We first need to deal with multilane samples.
    // We collect the lanes into lists of files, concatenate them into a single
    // file if there are several files and then recombine the channels.
    samples = samples.collect({[it[0], it[1..5].findAll(), it[6..10].findAll()]})

    samples = Channel.fromList(samples).branch({
        single_lane: it[1].size()==1 && it[2].size()<=1
        multi_lane_paired: it[1].size()>1 && it[2].size()>1
        multi_lane_single: it[1].size()>1 && it[2].size()==0
        other: true
        })

    samples.other.count().map({
        if (it > 0) {
            error "We do not support samples with a single forward but multiple backward read files."
        }
    })

    r1_multi = samples.multi_lane_paired.map({
            [it[0], it[1]]
        })
    r2_multi = samples.multi_lane_paired.map({
            [it[0], it[2]]
        })

    r1_multi = CAT_R1(r1_multi).file_out
    r2_multi = CAT_R2(r2_multi).file_out

    single_end_multi = samples.multi_lane_single.map({
            [it[0], it[1]]
        })
    single_end_multi = CAT_SINGLE_END(single_end_multi).file_out

    samples = samples.single_lane.map({
            [it[0], it[1][0], it[2][0]]
        }).mix(r1_multi.join(r2_multi)).mix(single_end_multi.map({
            [it[0], it[1], it[2]]
        }))

    // Update the metadata with the single_end parameter and put reads files in a list
    samples = samples.map( {
        row ->
            if (!row[2]) {
                return new Tuple (row[0] + [ single_end:true ], [ row[1] ])
            } else {
                return new Tuple (row[0] + [ single_end:false ], [ row[1], row[2] ])
            }
        })

    ///////////////////
    // PREPROCESSING //
    ///////////////////

    // We will skip preprocessing for samples which already have the preprocessed
    // folder in the ouput
    samples = samples.branch({
        already_preprocessed: params.resume_from_output && Files.isDirectory(Paths.get(outdir_abs, it[0].id, "preprocessed_reads"))
        to_preprocess: true
        })

    trimmed_reads = BBDUK_TRIM_ADAPTERS(samples.to_preprocess, params.adapters, false).reads
    phix_filtered_reads = BBDUK_FILTER_PHIX(trimmed_reads, params.phix, false).reads
    qf_reads = BBDUK_QUALITY_FILTERING(phix_filtered_reads, [], true).reads
    qf_singletons = BBDUK_QUALITY_FILTERING.out.singletons.map({ new Tuple (it[0] + [single_end:true], it[1]) })

    host_index = BBMAP_INDEX_HOST(params.masked_human)

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

    // Add the samples that had already been preprocessed
    preprocessed_samples = preprocessed_samples.mix(
        samples.already_preprocessed.map({
            it[0].single_end ?
            new Tuple (
                it[0],
                [Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_host_filtered.fastq.gz")]
                ):
            new Tuple (
                it[0],
                [Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_1_unmerged.fastq.gz"),
                 Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_2_unmerged.fastq.gz"),
                 Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_merged.fastq.gz"),
                 Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_host_filtered_singletons.fastq.gz")]
                )
    }))

    hf_reads = hf_reads.mix(
        samples.already_preprocessed.map({
            it[0].single_end ?
            new Tuple (
                it[0],
                [Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_host_filtered.fastq.gz")]
                ):
            new Tuple (
                it[0],
                [Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_host_filtered_1.fastq.gz"),
                 Paths.get(outdir_abs, it[0].id, "preprocessed_reads", "${it[0].id}_host_filtered_2.fastq.gz")]
                )
    }))

    // Make sure we have a single fastq file for all reads per sample
    reads = CAT_FASTQ(preprocessed_samples, true).reads

    /////////////////////////
    // Taxonomic Profiling //
    /////////////////////////

    if (!params.skip_motus) {
        // Skip samples for which motus has already been run and readd afterwards
        samples_motus = preprocessed_samples.branch({
            done: params.resume_from_output && Files.isDirectory(Paths.get(outdir_abs, it[0].id, "motus"))
            to_do: true
            })

        motus_profiles = MOTUS_PROFILE(samples_motus.to_do, params.motus_db).motus

        motus_profiles = motus_profiles.mix(
            samples_motus.done.map({
                new Tuple (
                    it[0],
                    Paths.get(outdir_abs, it[0].id, "motus", "${it[0].id}.motus")
                    )
                })
            )
    }

    if (!params.skip_phanta) {
        // we cannot use the singletons nor the merged reads here so he use hf_reads instead.
        PHANTA_PROFILE(
            hf_reads.filter({(!params.resume_from_output) || Files.notExists(Paths.get(outdir_abs, it[0].id, "phanta"))}),
            params.phanta_db)
    }

    if (!params.skip_assembly) {
        ///////////////////////////////
        // Assembly and gene calling //
        ///////////////////////////////

        // If we are not making a gene catalog we can skip samples for which we
        // already have the annotations and gene counts.
        if (params.skip_gene_catalog) {
            preprocessed_samples = preprocessed_samples.filter({
                (!params.resume_from_output) ||
                Files.notExists(Paths.get(outdir_abs, it[0].id, "annotations")) ||
                Files.notExists(Paths.get(outdir_abs, it[0].id, "gene_counts"))
                })
        }

        // Skip samples for which spades has already been run and readd afterwards
        samples_spades = preprocessed_samples.branch({
            done: params.resume_from_output && Files.isDirectory(Paths.get(outdir_abs, it[0].id, "assembly"))
            to_do: true
            })

        SPADES(samples_spades.to_do.map({ new Tuple (it[0], it[1], [], []) }), [], [])

        scaffolds = SPADES.out.scaffolds.mix(
            samples_spades.done.map({new Tuple (
                it[0],
                Paths.get(outdir_abs, it[0].id, "assembly", "${it[0].id}.scaffolds.fa.gz"))
            }))
        graphs = SPADES.out.gfa.mix(
            samples_spades.done.map({new Tuple (
                it[0],
                Paths.get(outdir_abs, it[0].id, "assembly", "${it[0].id}.assembly.gfa.gz"))
            }))
        paths = SPADES.out.assembly_paths.mix(
            samples_spades.done.map({new Tuple (
                it[0],
                Paths.get(outdir_abs, it[0].id, "assembly", "${it[0].id}.scaffolds.paths.gz"))
            }))

        assembly_graph_and_paths = scaffolds.join(graphs).join(paths)
        contig_classification = CLASSIFY_4CAC(assembly_graph_and_paths).classification

        filtered_assembly = FILTER_SCAFFOLDS(scaffolds.join(contig_classification)).all_scaffolds
        assembly_stats = ASSEMBLY_STATS(filtered_assembly)

        prokaryotic_genes = PRODIGAL(FILTER_SCAFFOLDS.out.prok_scaffolds, "gff")
        // We only predict eukaryotic genes if there are any (i.e. file is not empty)
        // this channel will therefore not contain all the samples, but it's fine
        // as we mix it back with the prokaryotic_genes channel afterwards.
        eukaryotic_genes = METAEUK_EASYPREDICT(
            FILTER_SCAFFOLDS.out.euk_scaffolds.filter({ it[1].size() > 0}),
            params.metaeuk_db)

        // we need to compress the output from MetaEuk so that both eukaryotic
        // and prokaryotic genes are compressed
        eukaryotic_genes_aa = PIGZ_COMPRESS_1(eukaryotic_genes.faa).archive
        eukaryotic_genes_nt = PIGZ_COMPRESS_2(eukaryotic_genes.codon).archive

        if (!params.skip_gene_catalog) {
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
        }

        if (!params.skip_per_sample) {
            /////////////////////////////////////
            // Genes counts individual samples //
            /////////////////////////////////////

            // We first gather the prokaryotic genes and eukaryotic genes together
            nt_tuples = prokaryotic_genes.nucleotide_fasta.mix(eukaryotic_genes_nt).groupTuple(size: 2, remainder: true)
            // Avoid redoing the mapping and count calculation if it was already done
            nt_tuples = nt_tuples.filter({
                (!params.resume_from_output) || Files.notExists(Paths.get(outdir_abs, it[0].id, "gene_counts"))
                })
            nucleotides = CAT_NT(nt_tuples).file_out

            catalog_index = BWA_INDEX_SAMPLES(nucleotides).index
            reads_index = reads.join(catalog_index)
            aligned_reads = BWA_MEM_SAMPLES(reads_index, false).bam
            filtered_reads = FILTERSAM_SAMPLES(aligned_reads).reads
            NORMALIZE_COUNTS_SAMPLES(filtered_reads.join(motus_profiles))


            ///////////////////////////////////
            // Functional individual samples //
            ///////////////////////////////////

            // We first gather the prokaryotic genes and eukaryotic genes together
            aa_tuples = prokaryotic_genes.amino_acid_fasta.mix(eukaryotic_genes_aa).groupTuple(size: 2, remainder: true)
            amino_acids = CAT_AA(aa_tuples).file_out

            if (params.custom_annotation_db) {
                // Avoid redoing the annotations if it was already done
                amino_acids_blast = amino_acids.filter({
                    (!params.resume_from_output) || Files.notExists(Paths.get(outdir_abs, it[0].id, "custom_annotations"))
                    })
                annot_db = BLAST_MAKEBLASTDB(new Tuple( [id: 'custom_annotations'], params.custom_annotation_db )).db
                hits = BLAST_BLASTP(amino_acids_blast, annot_db).tsv
            }
            // Avoid redoing the annotations if it was already done
            amino_acids_eggnog = amino_acids.filter({
                (!params.resume_from_output) || Files.notExists(Paths.get(outdir_abs, it[0].id, "annotations"))
                })
            EGGNOGMAPPER_SAMPLES(amino_acids_eggnog, params.eggnog_db, params.eggnog_dbdir, new Tuple([:], params.eggnog_dmnd))
        }
    }
}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
