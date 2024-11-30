# Memo


## set up

First, input:

Second, validating parameters; we have the following: 
- fasta file 1 (mandatory)  
- fasta file 2 (mandatory)   
- substitution matrix (optional, default is blosum62)
    - sequence type (derived)  
    - score matrix (the matrix object itself)  
- window size (optional, default is 10)  
- threshold (optional, default is 27)  

From fasta files we need to validate that  
1. the path of the file exists  
2. the files are actually fasta files  
3. only one sequence is present in each file (didnt do it yet)  

Any error in validation will not proceed to the following steps, so user has to make sure these are met in order to get a dotplot.

After that we would transform it into a Bio.Seq.Seq object, and add to the id its sequence type (whether protein or dna), which is something we check for from the alphabet used in the sequence (this step is useful later on when we want to access information about the sequence type). Now the sequences are ready to be used in the dotplot.  

Next we want to validate the actual algorithm parameters:  
1. the window size   
    - is a positive integer (didnt do yet)  
2. the threshold is a positive integer (didnt do yet)  
3. the substitution matrix:  
    - Got default matrices (only BLOSUM62 and PAM120) for protein sequences from [ncbi ftp site](ftp://ftp.ncbi.nih.gov/blast/matrices), and the DNAFull matrix from [rosalind](https://rosalind.info/glossary/dnafull/). These were converted to json and saved in one file in utils, adding to them, the sequence type (protein or dna) and the matrix object itself (to ease out incoming checks).  
    - The user when choosing it, 2 things to validate, 1. that the matrix actually is from those in the list (in json file) and 2. that the sequence type of the matrix is the same as the sequence type of the sequences (cant allow doing blosum for DNA)   
    - ensure parsing from json to a pandas dataframe   

Any error within substitution matrix will set seqtype to protein and matrix to blosum62, this might be problematic if the user inputs a dna seq, might think about adjusting the default value based on the actual fasta type. So no error will actually be raised to allow running the dotplot.


## algorithm  

First thing would be to initialize a matrix of rownames = seq1 and colnames = seq2, and fill it with zeros. Normally, we invert the rows, to make it in a way the 1st seq starts from bottom up and the 2nd from left to right.  
Then after than would need to fill it follwoing the dotplot algorithm:  
- iterate over the sequences, and for each window of size w, we would calculate the score of the window, no thresholding will be applied here.
- regoing over the matrix, we would apply the thresholding, and if the score is above the threshold, we would mark with the score (initially we would only mark it with a dot, but here we care about producing a continuous color map in the final visualization, so we would mark it with the score to allow for that)  
- after that generating the plot out of the matrix. For interactive plots, plotly will be used - particularly to get a heatmap with hover information

## shiny app  
