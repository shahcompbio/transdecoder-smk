# convert to gtf to bed file
rule gtf2bed:
    input:
        gtf = transcript_gtf,
    output:
        bed = os.path.join(out_dir, "visualization", "{sample}", "transcripts.bed"),
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "gtf_to_bed.pl {input.gtf} > {output.bed}"
# convert the predictions to a genome based bed as well
rule gtf2gff3:
    input:
        gtf = transcript_gtf,
    output:
        gff3 = os.path.join(out_dir, "visualization", "{sample}", "transcripts.gff3")
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "gtf_to_alignment_gff3.pl {input.gtf} > {output.gff3}"
# generate genome-based coding region annotation
rule genome_aligned_gff3:
    input:
        gff3 = os.path.join(out_dir,"visualization","{sample}","transcripts.gff3"),
        cds_regions = os.path.join(out_dir,"transdecoder","{sample}","transcripts.fa.transdecoder.gff3"),
        transcripts_fa = os.path.join(out_dir, "gffread", "{sample}", "transcripts.fa")
    output:
        genome_gff3 = os.path.join(out_dir,"visualization","{sample}","transcripts.transdecoder.genome.gff3")
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "cdna_alignment_orf_to_genome_orf.pl {input.cds_regions} {input.gff3} "
        "{input.transcripts_fa} > {output.genome_gff3}"
# generate bed file
rule gff3tobed:
    input:
        gff3 = os.path.join(out_dir,"visualization","{sample}","transcripts.transdecoder.genome.gff3")
    output:
        bed = os.path.join(out_dir,"visualization","{sample}","transcripts.fasta.transdecoder.genome.bed")
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "gff3_file_to_bed.pl {input.gff3} > {output.bed}"
