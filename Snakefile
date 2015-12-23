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
SPECIES2GENOMES_DICT = config["genomes"]["target_species"]
SPECIES = list(config["genomes"]["target_species"].keys())
SPECIES_GENOMES = [str(GENOMES_DIR + f) for f in list(config["genomes"]["target_species"].values())]

# Blast configuration
TBLASTN_PARAMS = " ".join(list(config["blast"]["tblastn"].values()))
BLAST_HEADER = "\t".join(config["blast"]["header"].split(" ")) # to add a tab separated header to blast results

# THREADS
THREADS = 8


# Desired outputs
SCAFF_CHR = "blast/protein2genome.outfmt6"
QUERY = "blast/query.fasta"
ALNS = expand("aln/{species}.delta",species = SPECIES)

rule all:
	input:
		SCAFF_CHR,	
		QUERY,
		ALNS



###################
## Download genomes
###################

#rule download_genomes_from_ncbi:


#######################
## Compare with genomes
#######################
rule align_genomes_with_nucmer:
    input:"blast/query.fasta"
    output:
        delta = "aln/{species}.delta"      
    message:"aligning {input} against {wildcards.species} genome"
    run:
        for species in SPECIES:
            shell("nucmer -p aln/" + species + " " + SPECIES2GENOMES_DICT[species] + " " + "{input)")
            
###########################################################################################################
## Find scaffold(s) or chromosome harbouring the protein of interest (in the species you are working with)
###########################################################################################################
rule extract_fasta_sequences_of_genomic_regions:
    input:
        ids = "blast/list_of_regions.txt",
        ref = config["genomes"]["dir"] + config["genomes"]["species_of_interest"]
    output:
        "blast/query.fasta"
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


