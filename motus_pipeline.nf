#!/usr/bin/env nextflow

process index_host_genome {

    output:
        path("host_bbmap_ref")

    script:
        """
        # Filtering out human host reads
        bbmap.sh -Xmx23g ref=$params.references.masked_human path="host_bbmap_ref" 2>> prepare_index.log
        """
    }


process preprocess_paired_end {

    input:
        tuple val(sample_name), path(file1), path(file2)
        path(host_bbmap_ref)

    output:
        tuple val(sample_name), path("merged.fq.gz"), path("R1.fq.gz"), path("R2.fq.gz"), path("singleton_reads.fq.gz")

    script:
        """
        # adapter trimming
        bbduk.sh -Xmx1G usejni=t in=$file1 in2=$file2 out=stdout.fq \
        refstats=adapter_trim.stats statscolumns=5 overwrite=t ref=$params.references.adapters \
        ktrim=r k=23 mink=11 hdist=1 2>> preprocessing.log | \

        # contaminant filtering
        bbduk.sh -Xmx1G usejni=t interleaved=true overwrite=t \
        in=stdin.fq out=stdout.fq ref=$params.references.phix \
        k=31 hdist=1 refstats=phix.stats statscolumns=5 2>> preprocessing.log | \

        # quality filtering. We also retain singleton reads that pass the filter
        # i.e. when only one of a paired end read passes,
        # we retain just that one in singleton_reads_qf.fq.gz
        bbduk.sh -Xmx1G usejni=t overwrite=t interleaved=true \
        in=stdin.fq out=qf.fasta.gz minlength=45 qtrim=rl maq=20 maxns=1 \
        stats=qc.stats statscolumns=5 trimq=14 outs=singleton_reads_qf.fq.gz 2>> preprocessing.log

        # Filtering out human host reads
        bbmap.sh -Xmx24g usejni=t interleaved=true overwrite=t \
        qin=33 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast \
        minhits=2 path=$host_bbmap_ref qtrim=rl trimq=15 untrim in=qf.fasta.gz \
        out=stdout.fq t=8 2>> removeHost.log |

        # Merging overlapping paired-end reads
        bbmerge.sh -Xmx24G interleaved=true in=stdin.fq out=merged.fq.gz \
        outu1=R1.fq.gz outu2=R2.fq.gz minoverlap=16 usejni=t \
        ihist=Sample1.merge.hist &> merge.log

        # Filtering out human host reads for singleton reads
        bbmap.sh -Xmx24g usejni=t threads=24 overwrite=t qin=33 minid=0.95 maxindel=3 \
        bwr=0.16 bw=12 quickmatch fast minhits=2 path=$host_bbmap_ref qtrim=rl trimq=15 \
        untrim=t in=singleton_reads_qf.fq.gz outu=singleton_reads.fq.gz 2>> out.rmHost.log
        """
    }

process preprocess_single_end {

    input:
        tuple val(sample_name), path(file1)
        path(host_bbmap_ref)

    output:
        tuple val(sample_name), path("singleton_reads.fq.gz")

    script:
        """
        # adapter trimming
        bbduk.sh -Xmx1G usejni=t in=$file1 out=stdout.fq \
        refstats=adapter_trim.stats statscolumns=5 overwrite=t ref=$params.references.adapters \
        ktrim=r k=23 mink=11 hdist=1 2>> preprocessing.log | \

        # contaminant filtering
        bbduk.sh -Xmx1G usejni=t interleaved=false overwrite=t \
        in=stdin.fq out=stdout.fq ref=$params.references.phix \
        k=31 hdist=1 refstats=phix.stats statscolumns=5 2>> preprocessing.log | \

        # quality filtering.
        # we retain just that one in singleton_reads_qf.fq.gz
        bbduk.sh -Xmx1G usejni=t overwrite=t interleaved=false \
        in=stdin.fq out=qf.fasta.gz minlength=45 qtrim=rl maq=20 maxns=1 \
        stats=qc.stats statscolumns=5 trimq=14 2>> preprocessing.log

        # Filtering out human host reads
        bbmap.sh -Xmx24g usejni=t interleaved=false overwrite=t \
        qin=33 minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast \
        minhits=2 path=$host_bbmap_ref qtrim=rl trimq=15 untrim in=qf.fasta.gz \
        outu=singleton_reads.fq.gz t=8 2>> removeHost.log
        """
    }

process motus_paired_end {

    input:
        tuple val(sample_name), path(merged), path(paired_1), path(paired_2), path(singleton_reads)

    output:
        tuple val(sample_name), path("${sample_name}.motus")

    script:
        """
        motus profile -f $paired_1 -r $paired_2 -s $merged,$singleton_reads\
        -n $sample_name -c -k mOTU -q -p -o ${sample_name}.motus
        """
    }

process motus_single_end {

    input:
        tuple val(sample_name), path(singleton_reads)

    output:
        tuple val(sample_name), path("${sample_name}.motus")

    script:
        """
        motus profile -s $singleton_reads \
        -n $sample_name -c -k mOTU -q -p -o ${sample_name}.motus
        """
    }

process assemble_paired_end {

    input:
        tuple val(sample_name), path(merged), path(paired_1), path(paired_2), path(singleton_reads)

    output:
        tuple val(sample_name), path("assembly")

    script:
        """
        mkdir assembly
        metaspades.py -t 4 -m 10 --only-assembler --pe-1 1 $paired_1 --pe-2 1 $paired_2 --pe-m 1 $merged --pe-s 1 $singleton_reads -o assembly
        """
    }

process filter_short_contigs {

    input:
        tuple val(sample_name), path(assembly)

    output:
        tuple val(sample_name), path("filtered_${sample_name}.scaffoldss.min1000.fasta")

    script:
        """
        python scaffold_filter.py $sample_name scaffolds $assembly filtered_$sample_name
        """
    }

process get_assembly_stats {

    input:
        tuple val(sample_name), path(filtered_assembly)

    output:
        tuple val(sample_name), path("${sample_name}.stats")

    script:
        """
        assembly-stats -l 500 -t <(cat $filtered_assembly) > $sample_name.stats
        """
    }

process call_genes {

    input:
        tuple val(sample_name), path(filtered_assembly)

    output:
        tuple val(sample_name), path("${sample_name}.faa"), path("${sample_name}.fna"), path("${sample_name}.gff")

    script:
        """
        prodigal -a $sample_name.faa -d $sample_name.fna -f gff \
        -o $sample_name.gff -c -q -p meta -i $filtered_assembly
        """
    }

process make_gene_catalog {

    input:
        tuple path(amino_acids), path(nucleotides)

    output:
        tuple path("cdhit9590")

    script:
        """
        mkdir cdhit9590
        cd-hit-est -i $nucleotides -o cdhit9590/gene_catalog_cdhit9590.fasta \
        -c 0.95 -T 20 -M 0 -G 0 -aS 0.9 -g 1 -r 1 -d 0

        grep "^>" cdhit9590/gene_catalog_cdhit9590.fasta | \
        cut -f 2 -d ">" | \
        cut -f 1 -d " " > cdhit9590/cdhit9590.headers
        seqtk subseq gene_catalog_all.faa cdhit9590/cdhit9590.headers > cdhit9590/gene_catalog_cdhit9590.faa
        """
    }


workflow {
    input_file = Channel.fromPath(params.input)
    samples = input_file
        .splitCsv(header: false, strip: false, limit: 1662)

    paired_end_samples = samples.filter( { it[2].strip() } )
    single_end_samples = samples.filter( { !it[2].strip() } ).map( {row -> new Tuple (row[0], row[1])})

    host_bbmap_ref = index_host_genome()

    preprocessed_paired = preprocess_paired_end(paired_end_samples, host_bbmap_ref)
    motus_paired_end(preprocessed_paired)

    preprocessed_single = preprocess_single_end(single_end_samples, host_bbmap_ref)

    motus_single_end(preprocessed_single)

    assembly = assemble_paired_end(preprocessed_paired)

    filtered_assembly = filter_short_contigs(assembly)

    assembly_stats = get_assembly_stats(filtered_assembly)

    genes = call_genes(filtered_assembly)

    amino_acids = genes.collectFile( {row ->  [ "genes.faa", row[1] ]} )
    nucleotides = genes.collectFile( {row ->  [ "genes.fna", row[2] ]} )
    gene_catalog = make_gene_catalog(amino_acids, nucleotides)
}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
