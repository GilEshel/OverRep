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
3. Annotations **do not** have to be GO terms: They can be anything that you like, as long as you provide a list of annotation ids per feature (can be 'NA'). Avoid having spaces within annotation ids, probably will cause problems. Examples for common annotations are GO terms, KEGG ids, InterPro ids, etc. You can also create your own annotations and test for overrepressentation, It is very flexable. 
4. Currently, you have to provide an annotation description file. First column will be the annotation id (e.g. GO:0006950) and second column will be a short description (e.g. response to stress). This is added to the output table. For GO terms, using the utility add_parent_GOs.py will generate the GO description file. If you don't have one for your annotation, you can modify the code to avoid it (or you can write an issue and I will add that option...).

### Adding parental GO terms:
For GO term overrepressentation analysis, you should make sure your gene2GO annotation file includes also the parental GO terms for the assigned GO terms. Since some genes will be annotated with more deeper level GO terms, and some with higher levels, we want to make sure we are counting all the genes for a given GO term. Use the utility add_parent_GOs.py to achive that:

\# Use the help option to see all input arguments:

**python add_parent_GOs.py -h**

\# Run the script on the example files

**python ./utils/add_parent_GOs.py -i ./Example_files/Aratha_gene2go.txt -d $'\t' -gd ',' -o Aratha_gene2go_with_parents.txt -r**
