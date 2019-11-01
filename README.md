<h1>Analysis of Protein Complexes in Cancer</h1>
This directory includes code and files for the analysis of the role of disrupted protein complexes in cancer.

<h2>Jupyter Notebooks</h2>

<h3>Analysis_Functions.ipynb</h3> Analysis and plotting functions used in other notebooks. Includes data and package
imports, functions involved in performing statistical tests on ratios within protein complexes, and visualization functions.

<h3>Plotting_Functions.ipynb</h3> Functions used to create various histograms and boxplots of the data as well as perform basic t-tests. Contains most of the functions used in pancancer analysis.

<h2>Data Files</h2>

<h3>proteinGroups_cleaned.txt</h3> Output of Clean_Dataset.ipynb. This is the file we used in our analyses, and contains
proteomics data for both tumor and normal samples, some of which are matched.

<h3>proteinGroups_simplified.txt</h3> The original data file, before it was cleaned up.

<h3>allComplexes.txt</h3> Data about protein complexes downloaded from CORUM (http://mips.helmholtz-muenchen.de/corum/#download)

<h3>trrust_rawdata.human.tsv.txt</h3> Data on transcription factors and their target genes downloaded from TRRUST v.2 
(https://www.grnpedia.org/trrust/downloadnetwork.php)

<h3>spliceosome_proteins.csv</h3> Information on which version of the spliceosome (A complex, B complex, etc) each protein is in
