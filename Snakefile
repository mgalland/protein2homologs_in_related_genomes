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
LIST_OF_SPECIES = list(config["genomes"]["target_species"].keys())
#LIST_OF_SPECIES_FASTAS = list(str(config["genomes"]["dir"]+ config["genomes"]["target_species"].values()))

# Blast configuration
TBLASTN_PARAMS = " ".join(list(config["blast"]["tblastn"].values()))
BLAST_HEADER = "\t".join(config["blast"]["header"].split(" ")) # to add a tab separated header to blast results

# THREADS
THREADS = 8


# Desired outputs
SCAFF_CHR = "blast/protein2genome.outfmt6"

rule all:
	input:
		SCAFF_CHR,	
		"blast/query_for_all_genomes.fasta"



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
rule extract_fasta_sequences_of_genomic_regions:
    input:
        ids = "blast/list_of_regions.txt",
        ref = config["genomes"]["dir"] + config["genomes"]["species_of_interest"]
    output:
        "blast/query_for_all_genomes.fasta"
    message:"retrieving scaffold(s)/chromosome(s) fasta file that match protein of interest"
    shell:
        "blastdbcmd -db {input.ref} -entry_batch {input.ids} > {output}"

rule extract_genomic_ids_of_matches_sequences:
    input:"blast/protein2genome.outfmt6"
    output:"blast/list_of_regions.txt"
    message:"extracting scaffold(s) / chromosome(s) matched by protein of interest"
    shell:
        """
        tail -n +2 {input} |awk '{{print $3}}' |sort|uniq > {output} 
        """ 
      
rule find_scaffold_or_chromosome:
    input:
        ref = config["genomes"]["dir"] + config["genomes"]["species_of_interest"],
        protein = config["query_protein"]
    output:"blast/protein2genome.outfmt6"
    message:"In species of interest, finding scaffold or chromosome where the protein of interest lies"
    shell:
        """
        makeblastdb -in {input.ref} -parse_seqids -dbtype nucl -out {input.ref};
        tblastn -num_threads {THREADS} {TBLASTN_PARAMS} -query {input.protein} -db {input.ref} -out {output};
        sed -i '1s/^/{BLAST_HEADER}\\n/' {output} 
        """    


