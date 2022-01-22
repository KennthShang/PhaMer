#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import pickle as pkl
import subprocess
import argparse
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#############################################################
########################  Parameters  #######################
#############################################################

parser = argparse.ArgumentParser(description="""PhaMer is a python library for identifying bacteriophages from metagenomic data. 
                                 PhaMer is based on a Transorfer model and rely on protein-based vocabulary to convert DNA sequences into sentences.""")
parser.add_argument('--contigs', help='FASTA file of contigs',  default = 'test_contigs.fa')
parser.add_argument('--len', help='minimun length of contigs', type=int, default=3000)
parser.add_argument('--midfolder', help='folder to store the intermediate files', type=str, default='phamer')
inputs = parser.parse_args()


#############################################################
######################  Check folders  ######################
#############################################################

out_fn = inputs.midfolder
transformer_fn = inputs.midfolder

if not os.path.isdir(out_dir):
    os.makedirs(out_fn)

#############################################################
##################  Filter short contigs  ###################
#############################################################
rec = []
for record in SeqIO.parse(inputs.contigs, 'fasta'):
    if len(record.seq) > inputs.len:
        rec.append(record)
SeqIO.write(rec, out_fn+'filtered_contigs.fa', 'fasta')

#############################################################
####################  Prodigal translation  #################
#############################################################


prodigal_cmd = f'prodigal -i {out_fn}/filtered_contigs.fa -a {out_fn}/test_protein.fa -f gff -p meta'
print("Running prodigal...")
_ = subprocess.check_call(prodigal_cmd, shell=True)



#############################################################
####################  DIAMOND BLASTP  #######################
#############################################################

print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
print("Creating Diamond database and running Diamond...")

try:
    # create database
    make_diamond_cmd = f'diamond makedb --threads 8 --in database/database.fa -d {out_fn}/database.dmnd'
    print("Creating Diamond database...")
    _ = subprocess.check_call(make_diamond_cmd, shell=True)
    # running alignment
    diamond_cmd = f'diamond blastp --threads 8 --sensitive -d {out_fn}/database.dmnd -q {out_fn}/test_protein.fa -o {out_fn}/results.tab -k 1'
    print("Running Diamond...")
    _ = subprocess.check_call(diamond_cmd, shell=True)
    diamond_out_fp = f"{out_fn}/results.tab"
    database_abc_fp = f"{out_fn}/results.abc"
    _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
except:
    print("create database failed")
    exit(1)




#############################################################
####################  Contig2Sentence  ######################
#############################################################


# Load dictonary and BLAST results
proteins_df = pd.read_csv('database/proteins.csv')
proteins_df.dropna(axis=0, how='any', inplace=True)
pc2wordsid = {pc: idx for idx, pc in enumerate(sorted(set(proteins_df['cluster'].values)))}
protein2pc = {protein: pc for protein, pc in zip(proteins_df['protein_id'].values, proteins_df['cluster'].values)}
blast_df = pd.read_csv(f"{out_fn}/results.abc", sep=' ', names=['query', 'ref', 'evalue'])

# Parse the DIAMOND results
contig2pcs = {}
for query, ref, evalue in zip(blast_df['query'].values, blast_df['ref'].values, blast_df['evalue'].values):
    conitg = query.rsplit('_', 1)[0]
    idx    = query.rsplit('_', 1)[1]
    pc     = pc2wordsid[protein2pc[ref]]
    try:
        contig2pcs[conitg].append((idx, pc, evalue))
    except:
        contig2pcs[conitg] = [(idx, pc, evalue)]

# Sorted by position
for contig in contig2pcs:
    contig2pcs[contig] = sorted(contig2pcs[contig], key=lambda tup: tup[0])



# Contigs2sentence
contig2id = {contig:idx for idx, contig in enumerate(contig2pcs.keys())}
id2contig = {idx:contig for idx, contig in enumerate(contig2pcs.keys())}
sentence = np.zeros((len(contig2id.keys()), 300))
sentence_weight = np.ones((len(contig2id.keys()), 300))
for row in range(sentence.shape[0]):
    contig = id2contig[row]
    pcs = contig2pcs[contig]
    for col in range(len(pcs)):
        try:
            _, sentence[row][col], sentence_weight[row][col] = pcs[col]
            sentence[row][col] += 1
        except:
            break

# Corresponding Evalue weight
#sentence_weight[sentence_weight<1e-200] = 1e-200
#sentence_weight = -np.log(sentence_weight)/200

# propostion
rec = []
for key in blast_df['query'].values:
    name = key.rsplit('_', 1)[0]
    rec.append(name)
counter = Counter(rec)
mapped_num = np.array([counter[item] for item in id2contig.values()])

rec = []
for record in SeqIO.parse(f'{out_fn}/test_protein.fa', 'fasta'):
    name = record.id
    name = name.rsplit('_', 1)[0]
    rec.append(name)
counter = Counter(rec)
total_num = np.array([counter[item] for item in id2contig.values()])
proportion = mapped_num/total_num


# Store the parameters
pkl.dump(sentence,        open(transformer_fn + 'sentence.feat', 'wb'))
pkl.dump(id2contig,       open(transformer_fn + 'sentence_id2contig.dict', 'wb'))
pkl.dump(proportion,      open(transformer_fn + 'sentence_proportion.feat', 'wb'))
pkl.dump(pc2wordsid,      open(transformer_fn + 'pc2wordsid.dict', 'wb'))
