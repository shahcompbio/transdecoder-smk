# extract cDNA from transcripts
rule extract_cDNA:
    input:
        gtf = transcript_gtf
    output:
        transcripts = os.path.join(out_dir, "gffread", "{sample}", "transcripts.fa")
    params:
        guide_fa = ref_genome
    container:
        "docker://quay.io/biocontainers/gffread:0.9.8--0",
    shell:
        """
        gffread -w {output.transcripts} -g {params.guide_fa} {input.gtf}
        """

# extract longest ORFs with transdecoder
rule extract_orfs:
    input:
        fasta = os.path.join(out_dir, "gffread", "{sample}", "transcripts.fa"),
    params:
        out_dir = os.path.join(out_dir, "transdecoder", "{sample}"),
    output:
        dat = os.path.join(out_dir, "transdecoder", "{sample}","transcripts.fa.transdecoder_dir", "base_freqs.dat"),
        cds = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.cds"),
        gff3 = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.gff3"),
        peps = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.pep")

    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "TransDecoder.LongOrfs -t {input.fasta} -O {params.out_dir}"

# run diamond blast against swissprot
# prep database if necessary
rule prep_diamond_db:
    input:
        sp_fasta = sp_fasta
    params:
        diamond_dir = diamond_dir,
        db_name = diamond_db_prefix
    output:
        sp_db = os.path.join(diamond_dir, diamond_db_prefix+".dmnd")
    resources:
        mem_mb = 50000,
        time = 20,
        partition = "cpushort"
    threads: 4
    container:
        "docker://quay.io/biocontainers/diamond:2.1.11--h5ca1c30_0"
    shell:
        """
        cd {diamond_dir} &&
        diamond makedb --in {input.sp_fasta} -d {params.db_name}
        """

rule run_diamond_blast:
    input:
        peps = os.path.join(out_dir,"transdecoder","{sample}","transcripts.fa.transdecoder_dir","longest_orfs.pep"),
        sp_db = os.path.join(diamond_dir, diamond_db_prefix+".dmnd")
    output:
        diamond_out = os.path.join(out_dir,"diamond_blast","{sample}","blastp.tsv"),
    params:
        diamond_dir = diamond_dir,
        db_name = diamond_db_prefix
    container:
        "docker://quay.io/biocontainers/diamond:2.1.11--h5ca1c30_0"
    resources:
        mem_mb = 50000,
        time = 20,
        partition = "cpushort"
    threads: 10
    shell:
        """
        diamond blastp --db {input.sp_db} --query {input.peps} \
        --max-target-seqs 1 --outfmt 6 \
        --evalue 1e-5 --threads 20 --out {output.diamond_out}
        """

# predict CDS regions with trandecoder
rule predict_cds:
    input:
        peps = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.pep"),
        fasta = os.path.join(out_dir, "gffread", "{sample}", "transcripts.fa"),
        diamond_blast= os.path.join(out_dir,"diamond_blast","{sample}","blastp.tsv")
    output:
        bed = os.path.join(out_dir, "transdecoder", "{sample}","transcripts.fa.transdecoder.bed"),
        cds = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.cds"),
        gff = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.gff3"),
        pep = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.pep"),
    params:
        out_dir = os.path.join(out_dir, "transdecoder", "{sample}"),
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        """
        TransDecoder.Predict -t {input.fasta} -O {params.out_dir} \
        --retain_blastp_hits {input.diamond_blast}
        """

# write protein fasta
rule generate_proteome:
    input:
        fasta = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.pep"),
    params:
        gtf = transcript_gtf,
        sp_fasta = sp_fasta,
    output:
        protein_fasta = os.path.join(out_dir, "proteome", "{sample}", "proteins.fa"),
        protein_tsv = os.path.join(out_dir, "proteome", "{sample}", "protein_transcript_info.tsv")
    container:
        "docker://quay.io/preskaa/biopython:v241011a"
    script:
        "../scripts/generate_proteome.py"