#!/bin/env/python

# < IMPORTS >

import pandas as pd
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# < GLOBAL VARIABLES (default values) >
# -- global --

SUB_MATRICES_PATH='../utils/substitution_matrices.json'
data=json.load(open(SUB_MATRICES_PATH))
SUB_MATRICES = {k:data[k]['type'] for k in list(data.keys())}

# -- user defined variables: --

# mandatory
FASTA_1 = '../../example/protein_seq1.fasta'
FASTA_2 = '../../example/protein_seq2.fasta'

# optional
SEQ_TYPE='protein'
WINDOW_SIZE=10
THRESHOLD=23
SCORE_MATRIX='blosum62'