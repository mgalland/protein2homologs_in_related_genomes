import sys
import os
import subprocess
from snakemake.utils import R
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from collections import Counter

#########################
## Pipeline configuration
##########################
configfile:"config.json"


# Genomes fasta files
GENOMES_DIR = config["genomes"]["dir"]

# Blast configuration
TBLASTN_PARAMS = " ".join(list(config["blast"]["tblastn"].values()))

# THREADS
THREADS = 8


# Desired outputs
SCAFF_CHR = "blast/protein2genome.outfmt6"

rule all:
	input:
		SCAFF_CHR 



###################
## Download genomes
###################

#rule download_genomes_from_ncbi:


#################
## Compare 
##################



###########################################################################################################
## Find scaffold(s) or chromosome for protein of interest (in the species you are working with)
###########################################################################################################
rule find_scaffold_or_chromosome:
    input:
        ref = config["genomes"]["dir"] + config["genomes"]["species_of_interest"],
        protein = config["query_protein"]
    output:"blast/protein2genome.outfmt6"
    message:"finding scaffold or chromosome where the protein of interest lies"
    shell:
        "makeblastdb -in {input.ref} -parse_seqids -dbtype nucl -out {input.ref} "
        "tblastn -num_threads {THREADS} {TBLASTN_PARAMS} -in {input.protein} -out {output}"
