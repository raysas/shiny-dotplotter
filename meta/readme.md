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

First thing would be to initialize a matrix of rownames = seq1 and colnames = seq2, and fill it with zeros. Normally, we invert the rows, to make it in a way the 1st seq starts from bottom up and the 2nd from left to right (consider putting seq 1 horizontally instead).  
Then after than would need to fill it follwoing the dotplot algorithm:  
- iterate over the sequences, and for each window of size w, we would calculate the score of the window, no thresholding will be applied here.
- regoing over the matrix, we would apply the thresholding, and if the score is above the threshold, we would mark with the score (initially we would only mark it with a dot, but here we care about producing a continuous color map in the final visualization, so we would mark it with the score to allow for that)  
- after that generating the plot out of the matrix. For interactive plots, plotly will be used - particularly to get a heatmap with hover information

## visualization

Using plotly for interactive visualization, figured that the best way is to do actually a scatter plot, because the heatmap is going to miss a lot of datapoints due to some errors in plotly rendering when theres a lot of data points (it worked for seaborn and matplotlib but not plotly, couldn't figure out why or ho to fix it)  
So to do the scatter plot, we would need to sequeeze the matrix into a list of tuples, where each tuple is (x, y, score) - so coordinates - and filter out null values. Then we would plot these points, coloring tehm using continuous color map to give it a heatmap look and focus on score matching regions.  
The plot would be interactive, with hover information showing the score of the point.  

Now comes a step to polish the plot, background, theme, colors, labels. Here the importance of saving up teh seq type in a variable, we would use the seq id of each fasta as an axis labels and then add the seq type to the title. For scaling the width and length, this will be relative to the length/width ratio because we dont want to distort the plot ( we want to keep it orthonormal): so starting from length = len(seq1) and width = len(seq2) (not sure about order).  
To fix sizes, define a ratio, based on a specified width or length (better user input), lets say of size s pixels, we wanna have lets say new_width=s, so we do ratio: new_width/width.  And now compute new_width=width * ratio and new_length =length * ratio (didnt do, still set to auto)

## kniting

need to make a script that knits all the steps together as one function that takes paramaters as input and returns the plot.  
Then we would use this function in the shiny app to generate the plot.

## shiny app  
