import plotly.express as px
from shiny.express import ui, input, render
from shinywidgets import render_plotly

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from dotplot import validation, plot_dotplot, get_parameters_description_df, workflow
from dotplot import WINDOW_SIZE, THRESHOLD, SCORE_MATRIX, SEQ_TYPE, FASTA_1, FASTA_2

with ui.sidebar():
    # input file
    # ui.input_file("fasta_file1", "Fasta File 1")
    # ui.input_file("fasta_file2", "Fasta File 2")
    # slider
    ui.input_slider("window_size", "Window Size", 1, 50, 1)
    ui.input_slider("threshold", "Threshold", 1, 50, 1)
    # dropdown
    ui.input_select("score_matrix", "Score Matrix", ["blosum62", "pam120", "dnafull"])
    ui.input_select("seq_type", "Sequence Type", ["protein", "dna"])

    # -- set default values
    # input.fasta_file1(FASTA_1)
    # input.fasta_file2(FASTA_2)
    # input.window_size(WINDOW_SIZE)
    # input.threshold(THRESHOLD)
    # input.score_matrix(SCORE_MATRIX)
    # input.seq_type(SEQ_TYPE)





# -- parameters (still testing on default ones)
window_size, threshold, score_matrix, seq_type, sub_mat,seq1, seq2 = validation(
    w=WINDOW_SIZE,
    t=THRESHOLD,
    sm=SCORE_MATRIX,
    st=SEQ_TYPE,
    fasta1=FASTA_1,
    fasta2=FASTA_2
)

# @render_plotly
# def parameters():
#     # -- description
#     df = get_parameters_description_df(
#         seq_type=seq_type,
#         window_size=window_size,
#         threshold=threshold,
#         score_matrix=score_matrix,
#         seq1=input.fasta_file1(),
#         seq2=input.fasta_file2()
#     )
#     fig = px.table(df)
#     return fig

@render_plotly
def dotplot():

    # -- plot the dotplot
    fig = workflow(
        w= window_size,
        t= threshold,
        sm= score_matrix,
        st= seq_type,
        fasta1= FASTA_1,
        fasta2= FASTA_2
    )
    return fig