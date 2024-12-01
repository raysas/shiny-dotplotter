# from default_params import *

# < IMPORTS >

import pandas as pd
import json
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from plotly import express as px

# set path to this file parent directory
import sys,os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# < GLOBAL VARIABLES (default values) >
# -- global --

SUB_MATRICES_PATH='utils/substitution_matrices.json'
data=json.load(open(SUB_MATRICES_PATH))
SUB_MATRICES = {k:data[k]['type'] for k in list(data.keys())}

# -- user defined variables: --

# mandatory 
FASTA_1 = '../example/protein_seq1.fasta'
FASTA_2 = '../example/protein_seq2.fasta'

# optional (set to default  if not given)
SEQ_TYPE='protein'
WINDOW_SIZE=10
THRESHOLD=23
SCORE_MATRIX='blosum62'

def validate_path(path_str):
    path=Path(path_str)
    if not path.exists():
        print(f'<!> Path {path} does not exist.')
        return False
    return True

def is_protein(fasta_file):
    try:
        seq = SeqIO.read(fasta_file, 'fasta')
        seq = str(seq.seq)
        if all([aa in 'ACDEFGHIKLMNPQRSTVWY' for aa in seq]):
            return True
        else:
            return False
    except:
        print(f'<!> File {fasta_file} is not a valid fasta file.')
        return False
    
def is_dna(fasta_file):
    try:
        seq = SeqIO.read(fasta_file, 'fasta')
        seq = str(seq.seq)
        if all([nt in 'ACGT' for nt in seq]):
            return True
        else:
            return False
    except:
        print(f'<!> File {fasta_file} is not a valid fasta file.')
        return False
    
def get_type(fasta_file):
    if is_protein(fasta_file):
        return 'protein'
    elif is_dna(fasta_file):
        return 'dna'
    else:
        return None
    
def get_seq_from_fasta(fasta_file):
    t= get_type(fasta_file)
    seq = SeqIO.read(fasta_file, 'fasta')
    seq.id = seq.id + f' {t}'
    return seq

def validate_substitution_matrix(matrix_name='blosum62', sequence_type='protein'):
    types=['protein','dna']
    matrix_name = matrix_name.lower()
    sequence_type = sequence_type.lower()


    if sequence_type not in types:
        print('<!> Invalid sequence type. Must be one of: {}\n\t-- Taking protein as default--'.format(types))
        sequence_type='protein'

    keys = [k for k,v in SUB_MATRICES.items() if v==sequence_type]
    # print(keys)
    
    if matrix_name not in keys:
        if sequence_type=='protein':
            print('<!> Invalid substitution matrix. Must be one of: {}\n\t'.format(keys))
            matrix_name='blosum62'
        else:
            print('<!> Invalid substitution matrix. Must be one of: {}'.format(keys))
            matrix_name='DNAFull'
        
    print(f'\n\t -- Taking {matrix_name} as substitution matrix for {sequence_type} sequences --')
    return matrix_name, sequence_type

def read_submat_from_json(matrix_name, json_file=SUB_MATRICES_PATH):

    json_file = open(json_file)
    data = json.load(json_file)
    df = pd.DataFrame(data[matrix_name]['matrix'])
    df.index = df.columns
    return df

def validate_positive_integers(w=WINDOW_SIZE,t=THRESHOLD):
    if not all([isinstance(i,int) for i in [w,t]]):
        print(f'<!> Window size and threshold must be integers. \n\t-- Taking default values')
        w=10
        t=23
    return w,t

def get_parameters(w=WINDOW_SIZE, t=THRESHOLD, sm=SCORE_MATRIX, stype=SEQ_TYPE, sub_mat_path=SUB_MATRICES_PATH):
    SCORE_MATRIX, SEQ_TYPE=validate_substitution_matrix(sm, stype)
    # THRESHOLD=validate_threshold(w, t)
    WINDOW_SIZE, THRESHOLD=validate_positive_integers(w, t)
    SUBSTITUTION_MATRIX=read_submat_from_json(sm, sub_mat_path)

    return WINDOW_SIZE, THRESHOLD, SCORE_MATRIX, SEQ_TYPE, SUBSTITUTION_MATRIX



def get_parameters_description(seq_type=SEQ_TYPE, window_size=WINDOW_SIZE, threshold=THRESHOLD, score_matrix=SCORE_MATRIX):
    return f'''
    -- Parameters --
    Sequence type: {seq_type}
    Window size: {window_size}
    Threshold: {threshold}
    Score matrix: {score_matrix}
    '''

def get_parameters_description_df(seq_type=SEQ_TYPE, window_size=WINDOW_SIZE, threshold=THRESHOLD, score_matrix=SCORE_MATRIX, seq1=None, seq2=None):
    data = {
        'Parameter': ['Sequence type', 'Window size', 'Threshold', 'Score matrix'],
        'Value': [seq_type, window_size, threshold, score_matrix]
    }
    df = pd.DataFrame(data, columns=['Parameter', 'Value'])
    df=df.set_index('Parameter')
    return df



def validation(w=WINDOW_SIZE, t=THRESHOLD, sm=SCORE_MATRIX, st=SEQ_TYPE, fasta1=FASTA_1, fasta2=FASTA_2):
    '''takes input paramters, validates them and returns the parameters, in this order:  
    WINDOW_SIZE, THRESHOLD, SCORE_MATRIX, SEQ_TYPE, SUBSTITUTION_MATRIX, SEQ_1, SEQ_2'''
    window, thresh, scor_mat, seqtype, sub_mat = get_parameters(w, t, sm, st)
    seq1 = get_seq_from_fasta(fasta1)
    seq2 = get_seq_from_fasta(fasta2)

    # print(get_parameters_description())
    # get_parameters_description_df()
    return window, thresh, scor_mat, seqtype, sub_mat, seq1, seq2

# < ALGORITHM >

def dotmatrix_wo_threshold(S1,S2, score_matrix,w):

    seq1=str(S1.seq)
    seq2=str(S2.seq)

    n=len(seq1); m=len(seq2) #seq class shpuld have len and iter

    rownames=[char for char in seq1]  
    colnames=[char for char in seq2]

    matrix=pd.DataFrame(
        [[0 for _ in range(m)] for _ in range(n)],
        index=rownames,
        columns=colnames
    )
    # matrix
    matrix=matrix.iloc[::-1] #reverse the rows to make the seq start from bottom up, left to right

    for res1_i in range(len(rownames)):
        for res2_i in range(len(colnames)):
            # print(matrix.iloc[0][0]) #success
            # print(matrix.iloc[res1_i][res2_i])

            w_score=0
            w_i=0; w_score=0
            while w_i in range(w) and res1_i+w_i<len(rownames) and res2_i+w_i<len(colnames):
                tmp_res1=rownames[res1_i+w_i]
                tmp_res2=colnames[res2_i+w_i]
                # print(identity_matrix[tmp_res1][tmp_res2]) #success
                w_score+=score_matrix[tmp_res1][tmp_res2]
                w_i+=1
            matrix.iloc[res1_i][res2_i]=w_score
    return matrix
        
def dotmatrix_w_threshold(matrix, threshold):
    '''takes output of dotmatrix_wo_threshold and applies threshold'''
    rownames=list(matrix.index)
    colnames=list(matrix.columns)
    new_matrix=matrix
    for res_i in range(len(rownames)):
        for res_j in range(len(colnames)):
            if new_matrix.iloc[res_i][res_j] <=threshold:
                new_matrix.iloc[res_i][res_j]=0
    return new_matrix

def get_dotmatrix(S1, S2, score_matrix, window, threshold):
    matrix=dotmatrix_wo_threshold(S1, S2, score_matrix, window)
    matrix=dotmatrix_w_threshold(matrix, threshold)
    return matrix

def plot_dotplot(matrix, s1, s2, seq_type=None):

    if not seq_type:
        seq_type = is_protein(s1) #needs fixing to take pbject of type Seq
        if not seq_type:
            seq_type = is_dna(s1)

    df=matrix
    df.columns=[i for i in range(len(df.columns))]
    df.index=[i for i in range(len(df.index))]

    df=df.stack().reset_index()
    df.columns=['x','y','value']
    df=df[df['value']!=0]

    # -- testing each seq on which axis --
    # x_axis=str(len(s1.seq))
    # y_axis=str(len(s2.seq))
    x_axis=f'Sequence 1:{s1.id}'
    y_axis=f'Sequence 2:{s2.id}'
    title=f'Dotplot of {seq_type} sequences'
    # return x_axis, y_axis, title


    fig = px.scatter(df, x='x', y='y', color='value', color_continuous_scale='cividis')
    fig.update_layout(
        title=title,
        xaxis_title=x_axis,
        yaxis_title=y_axis,
        paper_bgcolor='rgba(0,0,0,0)',  #background
        plot_bgcolor='rgba(0,0,0,0)',   #plot area
        template='plotly_dark', 
        # color_continuous_scale='Viridis' 
    )    
    # fig.show()
    return fig

def workflow():
    '''workflow to be called from main, returns fig object'''
    # -- get parameters from user..

    # -- validation --
    window, thresh, scor_mat, seqtype, sub_mat, seq1, seq2 = validation()
    metadata=get_parameters_description_df(seq_type=seqtype, window_size=window, threshold=thresh, score_matrix=scor_mat, seq1=seq1, seq2=seq2)
    
    # -- get dotmatrix and ploting --
    matrix = get_dotmatrix(seq1, seq2, sub_mat, window, thresh)
    return plot_dotplot(matrix, seq1, seq2, seqtype)


if __name__ == '__main__':
    # print(os.getcwd())
    fig=workflow()
    fig.write_html('../output/v3.html')
    