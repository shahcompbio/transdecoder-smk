import os
import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
"""
write protein fasta from transdecoder results
"""

def extract_complete_orfs(input_file):
    """
    filter fasta file based on fusion contig name ....
    :param input_file: Transcoder peptide sequences
    :return: dictionary of sequences with complete ORFs
    """
    seq_dict = {}
    with open(input_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            ## just take complete ORFs
            if "ORF type:complete" in record.description:
                ## remove asterisk at end of sequence if its there
                if record.seq.endswith("*"):
                    record.seq = record.seq[:-1]
                    seq_dict[record.description] = str(record.seq)
    return seq_dict

def write_fasta(seq_dict, output_file, novel_prefix="Bambu"):
    """
    :param seq_dict: complete ORFs predicted by transdecoder
    :param output_file: output fasta in msfragger format
    :return: writes fasta file in msfragger format
    """
    ## for fusions; may not be necessary
    fusion_id = 1
    fusion_ids = []
    contigs = []
    fusion_genes = []
    with open(output_file, "w") as outfile:
        for header, seq in seq_dict.items():
            ## write out non-canonical protein sequences ...
            ## write out canonical protein sequences ....
            if header.startswith("ENST"):
                protein = header.split(" ")[0]
                protein_id = protein.split(".p")[0]
                outfile.write(f'>%s PE=1\n' % protein_id)
                outfile.write(f'%s\n' % seq)
            else:
                protein = header.split(" ")[0]
                protein_id = protein.split(".p")[0]
                outfile.write(f'>%s PE=2\n' % protein_id)
                outfile.write(f'%s\n' %seq)
    return


## format fasta for msfragger; if fusions are included collect and rename

# ## inputs
# transdecoder_peps = sys.argv[1]
# jaffacsv_path = sys.argv[2]
# ## params
# fusion_calling = sys.argv[3]
# fusion_stats_dir = sys.argv[4]
# ## outputs
# filtered_fasta = sys.argv[5]

####
transdecoder_peps = snakemake.input["pep"]
jaffacsv_path = snakemake.input["jaffa_csv"]
## params
fusion_calling = snakemake.params["fusions"]
fusion_stats_dir = snakemake.params["fusion_stats_dir"]
## outputs
filtered_fasta = snakemake.output["protein_fasta"]


##filter for complete ORFs returned from transdecoder
seq_dict = extract_complete_orfs(transdecoder_peps)
## write the fasta and return a fusion dataframe
write_fasta(seq_dict, filtered_fasta)
