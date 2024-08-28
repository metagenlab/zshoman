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
        tuple val(sample_name), path("merged.fq.gz"), path("R1.fq.gz"), path("R2.fq.gz"), path("singleton_reads.fq.gz")

    output:
        tuple val(sample_name), path("${sample_name}.motus")

    script:
        """
        motus profile -f R1.fq.gz -r R2.fq.gz -s merged.fq.gz,singleton_reads.fq.gz \
        -n $sample_name -c -k mOTU -q -p -o ${sample_name}.motus
        """
    }

process motus_single_end {

    input:
        tuple val(sample_name), path("singleton_reads.fq.gz")

    output:
        tuple val(sample_name), path("${sample_name}.motus")

    script:
        """
        motus profile -s singleton_reads.fq.gz \
        -n $sample_name -c -k mOTU -q -p -o ${sample_name}.motus
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
}

workflow.onComplete {
    println "Done!!!!!!!!!!!!"
}
