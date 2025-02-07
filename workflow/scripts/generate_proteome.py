import pandas as pd
from gtfparse import read_gtf
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')
###
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

"""
generate proteome from transdecoder results in UniProt format
"""
def extract_features(line):
    """
    extract proteoform features
    """
    terms = line.split(" ")
    temp = terms[0]
    id, start, end = temp.split(":")
    _, ORFid = id.split("|")
    _, txid = ORFid.split("_")
    if start < end:
        strand = "+"
    else:
        strand = "-"
    return start, end, strand, ORFid, txid

def make_seqdict(input_file):
    """
    make a dict of sequences from fasta file
    :param input_file: fasta file
    :return: dictionary of sequences with ORFs
    """
    seq_dict = {}
    with open(input_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            record.seq = record.seq[:-1]
            seq_dict[record.description] = str(record.seq)
    return seq_dict

### inputs
input_fasta = snakemake.input["fasta"]
### params
gtfpath = snakemake.params["gtf"]
sp_fasta = snakemake.params["sp_fasta"]
### outputs
output_file = snakemake.output["protein_fasta"]
protein_tsv = snakemake.output["protein_tsv"]


## make sequence dictionary
seqdict = make_seqdict(input_fasta)

### load in de novo transcript assembly
gtf = read_gtf(gtfpath)
txdb = gtf[gtf["feature"] == "transcript"]
## extract gene info
tx2gene = dict(zip(txdb["transcript_id"], txdb["gene_id"]))
# load in swissprot fasta headers to give these proteoforms the exact swissprot headers
sp_headers = []
sp_ids = []
with open(sp_fasta, "r") as f:
    for line in f:
        if line.startswith(">"):
            sp_id = line.split("|")[1]
            sp_ids.append(sp_id)
            sp_headers.append(line)

sp_dict = dict(zip(sp_ids, sp_headers))

## make sequence dictionary
seqdict = make_seqdict(input_fasta)
# load in swissprot fasta headers to give these proteoforms the exact swissprot headers
tx_ids = []
protein_ids = []
gene_ids = []
og_headers = []
### proteoform count for novel proteoforms (with PG accession numbers)
count = 1

with open(output_file, "w+") as outfile:
    for header, seq in seqdict.items():
        temp = header.split(" ")
        tx_id = temp[0].split(".p")[0]
        # fetch gene
        gene = tx2gene[tx_id]
        if "sp|" in header:
            protein_id = temp[5].split("|")[1]
            sp_header = sp_dict[protein_id]
            outfile.write(sp_header)
        else:
            new_header = f">tr|PG{count}|PG{count}_HUMAN Predicted ORF from transcript ({tx_id}) OS=Homo sapiens OX=9606 GN={gene} PE=2\n"
            count = count + 1
            # give unique protein id
            protein_id = f"PG{count}"
            outfile.write(new_header)
        # write sequences
        outfile.write(f"{seq}\n")
        # save to create index table
        tx_ids.append(tx_id)
        protein_ids.append(protein_id)
        gene_ids.append(gene)
        og_headers.append(header)

### write protein csv file
descriptdat = pd.DataFrame(zip(protein_ids, tx_ids, gene_ids, og_headers),
                           columns=["Protein", "Transcript", "Gene", "transdecoder_header"])

descriptdat.to_csv(protein_tsv, sep="\t", index=None)