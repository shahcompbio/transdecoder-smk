# convert to gtf to bed file
rule gtf2bed:
    input:
        gtf = transcript_gtf,
    output:
        bed = os.path.join(out_dir, "visualization", "{sample}", "transcripts.bed"),
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "gtf_to_bed.pl {input.gtf} {output.bed}"
# convert the predictions to a bed file as well
rule gff3tobed:
    input:
        gff3 = os.path.join(out_dir,"transdecoder","{sample}","transcripts.fa.transdecoder.gff3")
    output:
        bed = os.path.join(out_dir,"visualization","{sample}","transcripts.fasta.transdecoder.genome.bed")
    container:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "gff3_file_to_bed.pl {input.gff3} > {output.bed}"
