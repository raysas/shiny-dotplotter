import plotly.express as px
import plotly.graph_objects as go
from shiny.express import ui, input, render
from shinywidgets import render_plotly
from shiny import reactive
import pandas as pd

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from dotplot import validation, plot_dotplot, get_parameters_description_df, workflow
from dotplot import WINDOW_SIZE, THRESHOLD, SCORE_MATRIX, SEQ_TYPE, FASTA_1, FASTA_2

# ui.tags.head(
#         ui.tags.link(rel="stylesheet", href="styles.css")
# )
ui.input_dark_mode()
ui.tags.style(
            """
            body {
                background-color: black;
                color: white; /* Optional: Set text color to white for better contrast */
            }
            """
    )

# with ui.sidebar():
#     # input file
#     # ui.input_file("fasta_file1", "Fasta File 1")
#     # ui.input_file("fasta_file2", "Fasta File 2")
#     # slider
#     ui.input_slider("window_size", "Window Size", 1, 50, 1)
#     ui.input_slider("threshold", "Threshold", 1, 50, 1)
#     # dropdown
#     ui.input_select("score_matrix", "Score Matrix", ["blosum62", "pam120", "dnafull"])
#     ui.input_select("seq_type", "Sequence Type", ["protein", "dna"])

#     # -- set default values
#     # input.fasta_file1(FASTA_1)
#     # input.fasta_file2(FASTA_2)
#     # input.window_size(WINDOW_SIZE)
#     # input.threshold(THRESHOLD)
#     # input.score_matrix(SCORE_MATRIX)
#     # input.seq_type(SEQ_TYPE)


with ui.sidebar():
    ui.input_file("file1", "Input first fasta:", multiple=False)
    ui.input_file("file2", "Input second fasta:", multiple=False)
    ui.input_slider("window_size", "Window Size", min=0, max=50, value=10) 
    ui.input_slider("threshold", "Threshold", min=0, max=50, value=23)
    ui.input_select("score_matrix", "Score Matrix", ["blosum62", "pam120", "dnafull"])
    ui.input_select("seq_type", "Sequence Type", ["protein", "dna"])
        
    ui.input_action_link("example_sequences", "Run Example")
    @reactive.effect
    @reactive.event(input.example_sequences)
    def run_example():
        # make window_size and threshold 10 and 23
        global window_size, threshold, score_matrix, seq_type, sub_mat, file1, file2
        window_size = WINDOW_SIZE
        threshold = THRESHOLD
        score_matrix = SCORE_MATRIX
        seq_type = SEQ_TYPE
        file1 = FASTA_1
        file2 = FASTA_2

    def input_taken():
        return input.file1() and input.file2()


    @render.text
    def log():
        global window_size, threshold, score_matrix, seq_type, sub_mat, file1, file2
        if input.example_sequences():
            window_size, threshold, score_matrix, seq_type, sub_mat,file1, file2 = validation(
            w=WINDOW_SIZE,
            t=THRESHOLD,
            sm=SCORE_MATRIX,
            st=SEQ_TYPE,
            fasta1=FASTA_1,
            fasta2=FASTA_2
        )
            return f"Running example on orthologous hemoglobin sequences"
        elif input_taken():
            window_size, threshold, score_matrix, seq_type, sub_mat,file1, file2 = validation(
                w=input.window_size(),
                t=input.threshold(),
                sm=input.score_matrix(),
                st=input.seq_type(),
                fasta1=input.file1(),
                fasta2=input.file2()
            )
            return "Input taken"
        return ""
        
    # useless?
    # @render.text
    # def value():
    #     return f"{input.slider()}"


    # << Test on parameters output: success >>

    # -- parameters (still testing on default ones)
    # window_size, threshold, score_matrix, seq_type, sub_mat,seq1, seq2 = validation(
    #     w=WINDOW_SIZE,
    #     t=THRESHOLD,
    #     sm=SCORE_MATRIX,
    #     st=SEQ_TYPE,
    #     fasta1=FASTA_1,
    #     fasta2=FASTA_2
    # )

    # if input.example_sequences():
    #     window_size, threshold, score_matrix, seq_type, sub_mat,file1, file2 = validation(
    #         w=WINDOW_SIZE,
    #         t=THRESHOLD,
    #         sm=SCORE_MATRIX,
    #         st=SEQ_TYPE,
    #         fasta1=FASTA_1,
    #         fasta2=FASTA_2
    #     )
    # elif input_taken():
    #     window_size, threshold, score_matrix, seq_type, sub_mat,file1, file2 = validation(
    #         w=input.window_size(),
    #         t=input.threshold(),
    #         sm=input.score_matrix(),
    #         st=input.seq_type(),
    #         fasta1=input.file1(),
    #         fasta2=input.file2()
    #     )

    @render.table
    def table():
        # -- default param for example
        if input.example_sequences():
            window_size, threshold, score_matrix, seq_type, sub_mat,file1, file2 = validation(
                w=WINDOW_SIZE,
                t=THRESHOLD,
                sm=SCORE_MATRIX,
                st=SEQ_TYPE,
                fasta1=FASTA_1,
                fasta2=FASTA_2
            )
            df= get_parameters_description_df(
                seq_type=seq_type,
                window_size=window_size,
                threshold=threshold,
                score_matrix=score_matrix,
                seq1=file1,
                seq2=file2
            )
            
        elif input_taken():
            df= pd.DataFrame({
                "Parameter": ["Sequence Type","Window Size", "Threshold", "Score Matrix", "Sequence 1 length", "Sequence 2 length"],
                "Value": [ input.seq_type(),input.window_size(), input.threshold(), input.score_matrix(), input.file1(), input.file2()]
            })
        else:
            data = {
                'Parameter': ['Sequence type', 'Window size', 'Threshold', 'Score matrix', 'Sequence 1 length', 'Sequence 2 length'],
                'Value': ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']}
            df = pd.DataFrame(data, columns=['Parameter', 'Value'])
            
        
        # fig = go.Figure(data=[go.Table(
        # header=dict(values=list(df.columns),
        #             fill_color='#08051D',
        #             align='left'),
        # cells=dict(values=[df[col] for col in df.columns],
        #         fill_color='#08051D',
        #         align='left'))
        # ])
        # fig.update_layout(
        #     title="Parameters",
        #     autosize=False,
        #     width=500,
        #     height=300,
        #     paper_bgcolor='rgba(0,0,0,0)',
        #     plot_bgcolor='rgba(0,0,0,0)',
        #     font=dict(
        #         family="Courier New, monospace",
        #         size=12,
        #         color="white"
        #     ),
        #     # theme='plotly_dark'
        # )

        return df




# -- parameters (still testing on default ones)
# window_size, threshold, score_matrix, seq_type, sub_mat,seq1, seq2 = validation(
#     w=WINDOW_SIZE,
#     t=THRESHOLD,
#     sm=SCORE_MATRIX,
#     st=SEQ_TYPE,
#     fasta1=FASTA_1,
#     fasta2=FASTA_2
# )

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
    if input.example_sequences() or input_taken():

        fig = workflow(
            w= window_size,
            t= threshold,
            sm= score_matrix,
            st= seq_type,
            fasta1= FASTA_1,
            fasta2= FASTA_2
        )
        return fig
    return None