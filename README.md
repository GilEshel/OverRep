# OverRep
A simple and flexable Python program to find over-represented annotations (GO, KEGG, etc.) in a subset of interesting features (e.g. genes).

## Running OverRep
\# You can invoke the help option to see all required and optional arguments:

**python overrep.py -h**

\# Running over-repressentated analysis for the provided example (use $'\\t' for tab-delimiters):

**python overrep.py -i ./Example_files/Aratha_interesting.list -b ./Example_files/Aratha_background.list -g ./Example_files/Aratha_gene2go_with_parents.txt -gd $'\t' -ad ',' -s ./Example_files/GO_descriptions.txt -sd $'\t' -a GO_term > analysis.log** # Analyses usually finish within one minute...

### Dependencies:
- Python (https://www.python.org/downloads/)
- SciPy (https://scipy.org/install.html). 'pip install scipy' worked for me.
- I only tested this program on Mac and Linux (Not sure if works on Windows - let me know...)

### Some notes before you start
1. Examples for input and output files can be found in the Example_files folder.
2. The features **do not** have to be genes: They can be anything that you like, the script just need a unique id (i.e. a string) for each feature. At least some of these features need to have at least one assigned annotation. Features can be genes, ortholog groups, etc.
3. Annotations **do not** have to be GO terms: They can be anything that you like, as long as you provide a list of annotation ids per feature (can be 'NA'). Avoid having spaces within annotation ids, probably will cause problems. Examples for common annotations are GO terms, KEGG ids, InterPro ids, etc. You can also create your own annotations and test for overrepresentation, It is very flexable. 
4. Currently, you have to provide an annotation description file. First column will be the annotation id (e.g. GO:0006950) and second column will be a short description (e.g. response to stress). This is added to the output table. For GO terms, using the utility add_parent_GOs.py will generate the GO description file. If you don't have one for your annotation, you can modify the code to avoid it (or you can write an issue and I will add that option...).

### Adding parental GO terms:
For GO term overrepresentation analysis, you should make sure your gene2GO annotation file includes also the parental GO terms for the assigned GO terms. Since some genes will be annotated with more deeper level GO terms, and some with higher levels, we want to make sure we are counting all the genes for a given GO term. Use the utility add_parent_GOs.py to achive that:

\# Use the help option to see all input arguments:

**python ./utils/add_parent_GOs.py -h**

\# Run the script on the example files:

**python ./utils/add_parent_GOs.py -i ./Example_files/Aratha_gene2go.txt -d $'\t' -gd ',' -o Aratha_gene2go_with_parents.txt -r**


## Running OverCompare
\# I provide another util program that can compare multipe overrepresentation result tables, for two groups of species. It count the number of species with significant overrepresentation for a given annotation id, from each group, and run Fisher's Exact test to find overrepresented annotations that are more frequent in one group over the other.

\# Use the help option to see all input arguments:

**python ./utils/overcompare.py -h**

\# Run the script on the example files:

**python ./utils/overcompare.py -p1 Example_files/group_comparison/Atacama -p2 Example_files/group_comparison/Sister -fe _overrepresentation_table.txt -c 8 -t 0.05 -pt padj -s1 Example_files/group_comparison/Atacama/Atacama_species.txt -s2 Example_files/group_comparison/Sister/Sister_species.txt -o test_example -a GO_term**

- It needs a path to a folder containing all the overrepresentation tables for group 1 (-p1), and group 2 (-p2).
- It needs a common end of file for all overrepresentation tables, to be able to find them (-fe). Default = '_overrepresentation_table.txt'.
- It needs to know the column in the overrepresentation table that represents significant value (-c). Default = '8' (the p-adjusted column).
- It needs to know the significance threshold to define overrepresentation (-t). Default = '0.05'.
- It needs to know which type of significance value is used (-pt). Default = 'padj'. This will be appended to the output file names for your records.
- It needs a file that contain the species short label names (e.g. Aratha) and full names (e.g. Arabidopsis thaliana) for group 1 (-s1). This file should be tab-delimited, and contain no header. The species short label should be at the begining of the overrepresentation table file name (e.g. Aratha_something_something_overrepresentation_table.txt).
- It needs a file that contain the species short label names (e.g. Zeamay) and full names (e.g. Zea mays) for group 2 (-s2).
- It needs a short string to add to output file names (-o), for your records, and so that you can run multiple comarisons without overwritting the output files. Default = '10_top_expressed' (I was looking at 10% top expressed genes. It should reflect the way the genes of interest were selected for the original annotation overrepresentation analysis...).
- Lastly it needs the type of annotation (-a), to add to output file names and columns. Default = 'GO_Term'.

#### I was using these scripts and many others for my own research, and thought it might be useful for others. I haven't tested these scripts under different environments, so I can't guarantee that they will work. If you are getting truble, and cann't figure the code, please post an issue, and hopefully I can help you.


**Cheers,**

**Gil Eshel**
