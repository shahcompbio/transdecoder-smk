## extract longest ORFs with transdecoder
rule extract_orfs:
    input:
        fasta = _fetch_cDNA,
    output:
        dat = os.path.join(out_dir, "transdecoder", "{sample}","transcripts.fa.transdecoder_dir", "base_freqs.dat"),
        cds = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.cds"),
        gff3 = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.gff3"),
        peps = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.pep")
    params:
        out_dir = os.path.join(out_dir, "transdecoder", "{sample}"),
    singularity:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "TransDecoder.LongOrfs -t {input.fasta} -O {params.out_dir}"

## predict CDS regions with trandecoder
rule predict_cds:
    input:
        peps = os.path.join(out_dir, "transdecoder", "{sample}", "transcripts.fa.transdecoder_dir", "longest_orfs.pep"),
        fasta = os.path.join(out_dir, "merged_cDNA", "{sample}", "transcripts.fa")
    output:
        bed = os.path.join(out_dir, "transdecoder", "{sample}","transcripts.fa.transdecoder.bed"),
        cds = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.cds"),
        gff = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.gff3"),
        pep = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.pep"),
    params:
        out_dir = os.path.join(out_dir, "transdecoder", "{sample}"),
    singularity:
        "docker://quay.io/biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0"
    shell:
        "TransDecoder.Predict -t {input.fasta} -O {params.out_dir} --single_best_only"

# ## output final fasta file for msfragger ...
# rule write_protein_fasta:
#     input:
#         pep = os.path.join(out_dir,"transdecoder", "{sample}","transcripts.fa.transdecoder.pep"),
#         jaffa_csv = _get_jaffal_fusion_csv(samples[0]),
#     output:
#         protein_fasta=os.path.join(out_dir, "proteome", "{sample}", "proteins.fa")
#     singularity:
#         "docker://quay.io/preskaa/biopython:v240911"
#     script:
#         "../scripts/write_fasta.py"
