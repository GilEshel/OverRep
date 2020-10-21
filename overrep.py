#!/usr/bin/env python

#########################################################################################
### The MIT License (MIT)
###
### Copyright (c) 2019 Gil Eshel
###
### Permission is hereby granted, free of charge, to any person obtaining a copy
### of this software and associated documentation files (the "Software"), to deal
### in the Software without restriction, including without limitation the rights
### to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
### copies of the Software, and to permit persons to whom the Software is
### furnished to do so, subject to the following conditions:
###
### The above copyright notice and this permission notice shall be included in all
### copies or substantial portions of the Software.
###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
#########################################################################################
### DATE CREATED: April 12, 2019
### AUTHOR: Gil Eshel
### CONTACT1: ge30@nyu.edu
### CONTACT2: giltu1@gmail.com
###
### CITATION: To be added
###
#########################################################################################
### DESCRIPTION:
### 
### To pass a tab-delimitor character to the relevant options, type: $'\t'
### 
### Required input: 
### A file containing a list of interesting gene identifiers (one per line) - should also be included in the background gene list file
### A file containing a list of background gene identifiers (one per line) - should correspond to the identifiers in the gene2Annotation file
### A gene2Annotation file: Two columns: (I) Gene_ID, (II) list of annotation ids (delimiters can be defined using the -gd and -ad options)
### A Annotation2Description file: Two columns: (I) Annotation_ID, (II) Annotation description (delimiter can be defined using the -sd option)
### 
### usage: python overrep.py [-h] [-i INTERESTING_GENE_FILE]
###                          [-b BACKGROUND_GENE_FILE] [-g GENE_TO_ANNOT_FILE]
###                          [-gd GENE_TO_ANNOT_DELIMITER]
###                          [-ad GENE_TO_ANNOT_ANNOT_DELIMITER]
###                          [-s ANNOT_DESCRIPTION_FILE]
###                          [-sd ANNOT_DESCRIPTION_DELIMITER]
###                          [-a ANNOTATION_TYPE]
### 
### A python script to conduct an over-repressentation analysis (for any given
### functional annotation [GO, KEGG, etc.]), for a list of 'interesting' genes
### given a list of background genes. It needs a predefined gene2Annotation
### mapping file (see ./Example_files/Aratha_gene2go_with_parents.txt as an
### example), and an Annotation2Decription files (see
### ./Example_files/GO_descriptions.txt as an example). It performs a 2-sides
### Fisher's Exact test using the stats.fisher_exact function (https://docs.scipy.
### org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html) in SciPy (a
### dependancy), and conduct a multiple testing FDR correction based on the
### Benjamini-Hochberg Method (1995) method. See https://scipy.org/install.html
### for how to install SciPy. Typing 'pip install scipy' in your command line
### terminal should work.
### 
### optional arguments:
###   -h, --help            show this help message and exit
###   -i INTERESTING_GENE_FILE, --interesting-gene-file INTERESTING_GENE_FILE
###                         [required] Specify a file containing a gene list of
###                         interest (one gene id per line). Expecting no column
###                         header.
###   -b BACKGROUND_GENE_FILE, --background-gene-file BACKGROUND_GENE_FILE
###                         [required] Specify a file containing the background
###                         gene list (one gene id per line). Expecting no column
###                         header.
###   -g GENE_TO_ANNOT_FILE, --gene-to-annot-file GENE_TO_ANNOT_FILE
###                         [required] Specify a file containing the gene to
###                         annotation (e.g. GO ids) mapping. Expecting two
###                         columns: (I) Gene_ID, (II) list of annotation ids
###                         (delimiters can be defined using the -gd and -ad
###                         options)
###   -gd GENE_TO_ANNOT_DELIMITER, --gene-to-annot-delimiter GENE_TO_ANNOT_DELIMITER
###                         Indicate the delimiter of the gene-to-annot file [e.g
###                         ' ', ',', '\t'. Also 't', 'w', 'c' will work]. Default
###                         is '\t')
###   -ad GENE_TO_ANNOT_ANNOT_DELIMITER, --gene-to-annot-annot-delimiter GENE_TO_ANNOT_ANNOT_DELIMITER
###                         Indicate the delimiter of the annotation list for each
###                         gene in the gene-to-annot file [e.g ' ', ',', '\t',
###                         ';']. Default is ',')
###   -s ANNOT_DESCRIPTION_FILE, --annot-description-file ANNOT_DESCRIPTION_FILE
###                         [required] Specify a file containing the annotation
###                         (e.g. GO ids) ids and their respective annotation
###                         (e.g. GO description/name), for all the included
###                         annotations (Two columns, not expecting headers)
###   -sd ANNOT_DESCRIPTION_DELIMITER, --annot-description-delimiter ANNOT_DESCRIPTION_DELIMITER
###                         Indicate the delimiter of the annot-description file
###                         [e.g ' ', ',', '\t', ' | ', ';'.]. Default is '\t')
###   -a ANNOTATION_TYPE, --annotation-type ANNOTATION_TYPE
###                         [optional] Specify the type of annotation used.
###                         Default = 'GO_Term'
### 
#########################################################################################

import scipy.stats as stats
import sys, argparse
import operator as op
from decimal import Decimal
from collections import Counter
def main():
	###	Define the arguments ###
	parser = argparse.ArgumentParser(description="A python script to conduct an over-repressentation analysis (for any given functional annotation [GO, KEGG, etc.]), for a list of 'interesting' genes given a list of background genes. It needs a predefined gene2Annotation mapping file (see ./Example_files/Aratha_gene2go_with_parents.txt as an example), and an Annotation2Decription files (see ./Example_files/GO_descriptions.txt as an example).\nIt performs a 2-sides Fisher's Exact test using the stats.fisher_exact function (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html) in SciPy (a dependancy), and conduct a multiple testing FDR correction based on the Benjamini-Hochberg Method (1995) method. See https://scipy.org/install.html for how to install SciPy. Typing 'pip install scipy' in your command line terminal should work.", prog="python overrep.py")
	parser.add_argument("-i","--interesting-gene-file", help="[required] Specify a file containing a gene list of interest (one gene id per line). Expecting no column header.")
	parser.add_argument("-b","--background-gene-file", help="[required] Specify a file containing the background gene list (one gene id per line). Expecting no column header.")
	parser.add_argument("-g","--gene-to-annot-file", help="[required] Specify a file containing the gene to annotation (e.g. GO ids) mapping.  Expecting two columns: (I) Gene_ID, (II) list of annotation ids (delimiters can be defined using the -gd and -ad options)")
	parser.add_argument("-gd","--gene-to-annot-delimiter",help="Indicate the delimiter of the gene-to-annot file [e.g ' ', ',', '\\t']. For tab delimiter, type $'\\t'. Default is '\\t')", default="\t")
	parser.add_argument("-ad","--gene-to-annot-annot-delimiter",help="Indicate the delimiter of the annotation list for each gene in the gene-to-annot file [e.g ' ', ',', '\\t', ';']. Should be different from the gene-to-annot-delimiter. For tab delimiter, type $'\\t'. Default is ',')", default=",")
	parser.add_argument("-s","--annot-description-file", help="[required] Specify a file containing the annotation (e.g. GO ids) ids and their respective annotation (e.g. GO description/name), for all the included annotations (Two columns, not expecting headers)")
	parser.add_argument("-sd","--annot-description-delimiter",help="Indicate the delimiter of the annot-description file [e.g ' ', ',', '\\t', ' | ', ';'.]. Default is '\\t')", default="\t")
	parser.add_argument("-a","--annotation-type", action="store", default="GO_Term", help="[optional] Specify the type of annotation used. Default = 'GO_Term'")
	#
	if len(sys.argv[1:])==0:
		parser.print_usage()
		parser.exit()
	args = parser.parse_args()
	#
	print('Loaded input arguments')
	#
	### Functions ###
	# A function to conduct a p-value correction based on FDR using the Benjamini-Hochberg Method (1995):
	# Brief explanation:
	# (1) It take all the p-values and order them from the smallest (most significant) to the largest. (and their position (index + 1 ) is their rank
	# (2) It assign corrected p-values as follows:
	# 		- the largest p-value gets the same original p-value
	#		- it loops from the second largest to the smallest p-value. each ajust-p-val will be the smallest of two options: (1) the previous adjust p-value, or (2) the current p-value *(total_num_of_p-values / the current p-value rank) 
	def fdr_corr(pval_list):
		# Make sure that the returned adjust p-values will be in the same original order:
		original_p_order = range(len(pval_list))	# Record the original order of the p-values
		p_zip = zip(pval_list,original_p_order)	# Zip the original p-values with their original positions
		p_reverse_sorted = sorted(p_zip, key=lambda tup: tup[0], reverse=True)	# Reverse sort the p-values together with original positions
		new_p_order = [x[1] for x in p_reverse_sorted]	# Extract back the p-value positions, with the new order (after, zip this list with the new adjust p-values and sort by the postions (back to original positions))
		original_p_reverse_sorted = [x[0] for x in p_reverse_sorted]	# Extract back the p-value
		# Now start the actual calculations:
		adj_p_reverse_order = [original_p_reverse_sorted[0]]	# Store the adjusted p-values, start with the largest p-value = adjust p-value
		rank = len(original_p_reverse_sorted) - 1	# Rank of each p-value, starting with the second largest rank (for the second largest p-value)
		previous_adj = original_p_reverse_sorted[0]	# Need to compare the the adjust p-values with the previous adjust-pvalues in the loop - start with the largest p-val and keep update in the loop
		for p in original_p_reverse_sorted[1:]: # Loop from the second largest
			p_adj = float(p) * len(original_p_reverse_sorted) / rank
			rank = rank - 1	# Update for the next iteration...
			adj_p_reverse_order.append(min(p_adj,previous_adj))	# Take the min of the two values
			previous_adj = min(p_adj,previous_adj)	# Update for the next iteration...
		p_zip_adj = zip(adj_p_reverse_order,new_p_order)	# Zip the adjust p-values and their reversed ordered original positions
		p_adj_original_order = sorted(p_zip_adj, key=lambda tup: tup[1])	# Sort the adjust p-values to their original positions	
		return [x[0] for x in p_adj_original_order]	# Return the list of adjust p-values in their original positions	
	#
	print('Loaded script functions')
	#
	### Input ###
	#
	interesting_genes = []
	if args.interesting_gene_file is not None:
		with open(str(args.interesting_gene_file), 'r') as intFile:
			for line in intFile:
				interesting_genes.append(line.strip())
	else:
		sys.exit("Error: Please provide a file containing the genes of interest (one gene id per line), using the '-i' option. Type: python overrep.py -h for more instructions of usage")
	#
	print('Parsed ' + str(args.interesting_gene_file) + ' file')
	#
	background_genes = []
	if args.background_gene_file is not None:
		with open(str(args.background_gene_file), 'r') as backFile:
			for line in backFile:
				background_genes.append(line.strip())
	else:
		sys.exit("Error: Please provide a file containing the background genes (one gene id per line), using the '-b' option. Type: python overrep.py -h for more instructions of usage")
	#
	print('Parsed ' + str(args.background_gene_file) + ' file')
	#
	gene_to_annot_dict = {}
	gene_to_annot_dict_for_phylobrowse = {}	# Store the gene id and the whole line as the value (will use this at the end of the script...)
	if args.gene_to_annot_file is not None:
		with open(str(args.gene_to_annot_file), 'r') as annFile:
			for annot in annFile:
				annotation = annot.strip().split(args.gene_to_annot_delimiter)
				if len(annotation) != 2:
					sys.exit("Error: Expecting 2-column gene-to-annotation file. The file " + str(args.gene_to_annot_file) + " has " + str(len(annotation)) + " columns! Separator used: " + str(args.gene_to_annot_delimiter))
				else:
					if annotation[1] != 'NA':
						gene_to_annot_dict[annotation[0]] = annotation[1].split(args.gene_to_annot_annot_delimiter)
						gene_to_annot_dict_for_phylobrowse[annotation[0]] = annot
	else:
		sys.exit("Error: Please provide a gene to annotation (e.g. GO ids) mapping, for all the included model organisms. Using the '-g' option, and the delimiter options '-gd' and 'ad'. No header line is expected. Type: python overrep.py -h for more instructions of usage")
	#
	print('Parsed ' + str(args.gene_to_annot_file) + ' file')
	#
	#
	annot_desc_dict = {}
	annot_desc_dict_for_phylobrowse = {}	# Store the annotation id and the whole line as the value (will use this at the end of the script...)
	if args.annot_description_file is not None:
		with open(str(args.annot_description_file), 'r') as desFile:
			for annot in desFile:
				annotation = annot.strip().split(str(args.annot_description_delimiter))
				if len(annotation) != 2:
					sys.exit("Error: Expecting 2-column annotation description file. The file " + str(args.annot_description_file) + " has " + str(len(annotation)) + " columns! Separator used: " + str(args.annot_description_delimiter))
				else:
					annot_desc_dict[annotation[0]] = annotation[1]
					annot_desc_dict_for_phylobrowse[annotation[0]] = annot
	else:
		sys.exit("Error: Please provide a annotation (e.g. GO terms) description file for all included annotation ids in the gene2annot file. Using the '-s' option, and the delimiter using '-sd' option. No header line is expected. Type: python overrep.py -h for more instructions of usage")
	#
	print('Parsed ' + str(args.annot_description_file) + ' file')
	#
	#
	print('Starting overrepresentation analysis')
	#
	### Start the overrepresentation analysis ###
	#
	## Get background and interesting genes annotation counts
	#
	# Background:
	gene_to_annot_dict_background_genes = { your_key: gene_to_annot_dict[your_key] for your_key in background_genes if your_key in gene_to_annot_dict }
	num_genes_with_annot_universe = len(gene_to_annot_dict_background_genes)
	universe_annot_ids = sum(list(gene_to_annot_dict_background_genes.values()),[])	# Get all annotation ids for the background genes
	universe_annot_id_frequencies = dict(Counter(universe_annot_ids))	# For every unique annotation id, record the frequency in the background
	print('There are ' + str(num_genes_with_annot_universe) + ' genes with a ' + str(args.annotation_type) + ' annotation in the background list')
	print('There are ' + str(len(universe_annot_id_frequencies)) + ' ' + str(args.annotation_type) + 's with at least one gene count in the background')
	#
	# Genes of interest:
	gene_to_annot_dict_interesting_genes = { your_key: gene_to_annot_dict[your_key] for your_key in interesting_genes if your_key in gene_to_annot_dict }
	num_of_interesting_genes_with_annot = len(gene_to_annot_dict_interesting_genes)
	interesting_annot_ids = sum(list(gene_to_annot_dict_interesting_genes.values()),[])	# Get all annotation ids for the interesting genes
	interesting_annot_id_frequencies = dict(Counter(interesting_annot_ids))	# For every unique annotation id, record the frequency in the genes of interest
	print('There are ' + str(num_of_interesting_genes_with_annot) + ' genes with a ' + str(args.annotation_type) + ' annotation in the list of interesting genes')
	print('There are ' + str(len(interesting_annot_id_frequencies)) + ' ' + str(args.annotation_type) + 's with at least one gene count in the list of interesting genes')
	#
	annot_id_interesting_gene_ids_dict = {}	# Record the ids of interesting genes having a particular annotation id
	for i in list(gene_to_annot_dict_interesting_genes.keys()):
		for a in gene_to_annot_dict_interesting_genes[i]:
			if a in list(annot_id_interesting_gene_ids_dict.keys()):
				annot_id_interesting_gene_ids_dict[a].append(i)
			else:
				annot_id_interesting_gene_ids_dict[a] = [i]
	#
	## Run Fisher's Exact test:
	results = {}
	pvalues = []
	for a in list(interesting_annot_id_frequencies.keys()):	# Run analysis for each annotation id that have been assigned to a gene of interest
		interesting_genes_with_this_annotid = interesting_annot_id_frequencies[a]	# The number of interesting genes that contain this annotation id (needed for contingency table)
		interesting_genes_without_this_annotid = num_of_interesting_genes_with_annot - interesting_genes_with_this_annotid	# The number of interesting genes that have an annotation, but not this annotation id - need for contingency table
		total_with_this_annotid = universe_annot_id_frequencies[a]	# Total because it includes the interesting genes as well...
		total_without_this_annotid = num_genes_with_annot_universe - total_with_this_annotid	
		universe_with_this_annotid = total_with_this_annotid - interesting_genes_with_this_annotid # Need it for the contingency table
		universe_without_this_annotid = total_without_this_annotid - interesting_genes_without_this_annotid	# Need it for the contingency table
		p_val = stats.fisher_exact([[interesting_genes_with_this_annotid,universe_with_this_annotid],[interesting_genes_without_this_annotid,universe_without_this_annotid]])[1]
		pvalues.append(p_val)	# Collect the p-value
		results[a] = '\t'.join([a,annot_desc_dict[a],str(interesting_genes_with_this_annotid),str(num_of_interesting_genes_with_annot),str(total_with_this_annotid),str(num_genes_with_annot_universe),str(p_val)]) # Store results
	#
	print('Running BH FDR correction')
	adj_pval = fdr_corr(pvalues)
	pos = 0	# Track the position of the adj_p-val
	annot_id_p_adj = {}
	annot_id_p_val = {}
	for a in list(interesting_annot_id_frequencies.keys()):	# Add the adj_p_val to the results list in the same order
		annot_id_p_val[a] = float(results[a].split('\t')[6])
		annot_id_p_adj[a] = float(adj_pval[pos])
		results[a] = results[a] + '\t' + str(adj_pval[pos]) + '\t' + ';'.join(annot_id_interesting_gene_ids_dict[a])
		pos = pos + 1	
	num_overrepresented_terms_pos = sum(i <= 0.05 for i in adj_pval)  
	print('Generating overrepresentation table: ' + str(num_overrepresented_terms_pos) + ' overrepresented (p-adj <= 0.05) ' + str(args.annotation_type) + 's were identified')
	#
	### Sort GOs by the p-val (lowest to highest) and filter by p-adj <= 0.05, p-val <= 0.05 and p-val <= 0.01
	annot_id_p_val_sorted_top = sorted(annot_id_p_val.items(), key = lambda kv:(kv[1], kv[0]))
	annot_id_p_adj_sorted_top = sorted(annot_id_p_adj.items(), key = lambda kv:(kv[1], kv[0]))
	annot_ids_sorted_by_p_val = [i[0] for i in annot_id_p_val_sorted_top]
	annot_ids_sorted_by_p_adj = [i[0] for i in annot_id_p_adj_sorted_top]
	annot_ids_sorted_by_p_val_below_0_5 = [i[0] for i in annot_id_p_val_sorted_top if i[1] <= 0.05]
	annot_ids_sorted_by_p_val_below_0_1 = [i[0] for i in annot_id_p_val_sorted_top if i[1] <= 0.01]
	annot_ids_sorted_by_p_adj_below_0_5 = [i[0] for i in annot_id_p_val_sorted_top if annot_id_p_adj[i[0]] <= 0.05] # keeping the order by p-val for consistency...
	#
	### Save overrepresentation tables:
	# Save the full sorted overrepresentation table (without adj_pval filtering)
	with open(args.interesting_gene_file.rsplit('.',1)[0] + '_' + str(args.annotation_type) + '_overrepresentation_table.txt' ,'w') as enrich_table:
		enrich_table.write('annotid' + '\t' + 'AnnotDesc' + '\t' + 'Interesting_genes_with_this_annotid' + '\t' + 'Interesting_genes_with_any_annotation' + '\t' + 'BG_with_annotid' + '\t' + 'BG_with_any_annotation' + '\t' + 'pvalue' + '\t' + 'p_adj_BH' + '\t' + 'Interesting_gene_ids' + '\n')
		for an in annot_ids_sorted_by_p_val:
			enrich_table.write(results[an] + '\n')
	#
	# Save the filtered overrepresentation table (adj_pval <= 0.05):
	with open(args.interesting_gene_file.rsplit('.',1)[0] + '_' + str(args.annotation_type) + '_overrepresentation_table_padj_0_05.txt' ,'w') as filt_enrich_table1:
		filt_enrich_table1.write('annotid' + '\t' + 'AnnotDesc' + '\t' + 'Interesting_genes_with_this_annotid' + '\t' + 'Interesting_genes_with_any_annotation' + '\t' + 'BG_with_annotid' + '\t' + 'BG_with_any_annotation' + '\t' + 'pvalue' + '\t' + 'p_adj_BH' + '\t' + 'Interesting_gene_ids' + '\n')
		for f1 in annot_ids_sorted_by_p_adj_below_0_5:
			filt_enrich_table1.write(results[f1] + '\n')
	#
	# Save the filtered overrepresentation table (pval <= 0.05):
	with open(args.interesting_gene_file.rsplit('.',1)[0] + '_' + str(args.annotation_type) + '_overrepresentation_table_pval_0_05.txt' ,'w') as filt_enrich_table2:
		filt_enrich_table2.write('annotid' + '\t' + 'AnnotDesc' + '\t' + 'Interesting_genes_with_this_annotid' + '\t' + 'Interesting_genes_with_any_annotation' + '\t' + 'BG_with_annotid' + '\t' + 'BG_with_any_annotation' + '\t' + 'pvalue' + '\t' + 'p_adj_BH' + '\t' + 'Interesting_gene_ids' + '\n')
		for f2 in annot_ids_sorted_by_p_val_below_0_5:
			filt_enrich_table2.write(results[f2] + '\n')
	#
	# Save the filtered overrepresentation table (pval <= 0.01):
	with open(args.interesting_gene_file.rsplit('.',1)[0] + '_' + str(args.annotation_type) + '_overrepresentation_table_pval_0_01.txt' ,'w') as filt_enrich_table3:
		filt_enrich_table3.write('annotid' + '\t' + 'AnnotDesc' + '\t' + 'Interesting_genes_with_this_annotid' + '\t' + 'Interesting_genes_with_any_annotation' + '\t' + 'BG_with_annotid' + '\t' + 'BG_with_any_annotation' + '\t' + 'pvalue' + '\t' + 'p_adj_BH' + '\t' + 'Interesting_gene_ids' + '\n')
		for f3 in annot_ids_sorted_by_p_val_below_0_1:
			filt_enrich_table3.write(results[f3] + '\n')
	#
	### Filter the gene2Annot and Annot2Decription files for PhyloBrowse (to contain only the relevant genes/annotations - i.e. the background gene list (instead of all the genes in the gene2go file) ###
	# gene2Annot:
	#with open(args.gene_to_annot_file.rsplit('.',1)[0] + '_filtered.' + args.gene_to_annot_file.rsplit('.',1)[1], 'w') as filtGene2Ann:
	#	for b in background_genes:
	#		if b in list(gene_to_annot_dict_for_phylobrowse.keys()):
	#			filtGene2Ann.write(gene_to_annot_dict_for_phylobrowse[b])
	#
	# Annot2Decription:
	#with open(args.annot_description_file.rsplit('.',1)[0] + '_filtered.' + args.annot_description_file.rsplit('.',1)[1], 'w') as filtAnn2Desc:
	#	for a in list(universe_annot_id_frequencies.keys()):
	#		filtAnn2Desc.write(annot_desc_dict_for_phylobrowse[a])

if __name__ == '__main__':
	main()

