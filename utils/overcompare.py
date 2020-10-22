# -*- coding: utf-8 -*-
#!/usr/bin/env python
#########################################################################################
### The MIT License (MIT)
###
### Copyright (c) 2020 Gil Eshel
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
### DATE CREATED: October 12, 2020
### AUTHOR: Gil Eshel
### CONTACT1: ge30@nyu.edu
### CONTACT2: giltu1@gmail.com
###
### CITATION: To be added
###
#########################################################################################
### DESCRIPTION:
### 

# Gil Eshel, Oct. 12, 2020
# Compare the overrepresented (padj≤0.05) GO terms between Atacama and Sister groups
# I first got a list of all overrepresented GO terms from all the species:
# cat */*_GO_Term_overrepresentation_table_padj_0_05.txt | cut -f1 | grep -v "annotid" | sort | uniq > All_overrepresentationed_GO_terms_padj_0_05.txt
# Will loop over the *GO_Term_overrepresentation_table_padj_0_05.txt (this is the default, can be changed with option -fe) for each group, and count the number of species with overrepresentation per GO term
# Will compare the proportion of species with overrepresentation between the two groups using Fisher Exact, and will run FDR correction
# Desired Outputs:
# 1. Matrix with padj values per species (columns - ordered like figure 4), per overrepresented GO term (rows)
# 2. A table with the following columns: "GO_term\tGO_desc\tNumber_overrepresented_species_Atacama\tNumber_overrepresented_species_Sister\tpval\tpadj\tAtacama_overrepresented_species\tSister_overrepresented_species\n"
# Usage: python overcompare.py -p1 Atacama -p2 Sister -s1 /Users/gileshel/Documents/NYU_postdoc/EvoNet_project/EvoNet_v1_paper/Dec_2_version_from_Rodrigo/Atacama_species.txt -s2 /Users/gileshel/Documents/NYU_postdoc/EvoNet_project/EvoNet_v1_paper/Dec_2_version_from_Rodrigo/New_expression_from_Tomas_10072020/RSEM_expression_to_phylogenomic_contig_sets/Sister/Sister_18_species.txt

import scipy.stats as stats
import argparse, glob, sys
import operator as op
from decimal import Decimal
from collections import Counter

def main():
	###	Define the arguments ###
	parser = argparse.ArgumentParser(description="OverCompare is a Python program that compares overrepresentation files among species from two defined groups, and performs a 2-sides Fisher's Exact test to find overrepresented annotations that are more frequent in one group over the other. It uses the stats.fisher_exact function from SciPy (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html), so SciPy needs to be installed ('pip install scipy'). The script uses the Benjamini-Hochberg method (1995) to correct for multiple testing.\nOutput1: a significance matrix (rows = overrepressented annotations, columns = species (group 1 first and then group 2)). This can be easily plotted as heatmap in R.\nOutput2: Fisher's Exact test results (p-val and p-adj values, as well as the number and name of speceis overrepresented for each annotation in each group.", prog="python overcompare.py")
	parser.add_argument("-p1","--path-to-group1-overrepresentation-tables", help="[required] Specify the path to a folder containing the overrepresentation tables for all the species in group 1. overrepresentation table files are expected to start with species_label + '_', and to have a file header.")
	parser.add_argument("-p2","--path-to-group2-overrepresentation-tables", help="[required] Specify the path to a folder containing the overrepresentation tables for all the species in group 2. overrepresentation table files are expected to start with species_label + '_', and to have a file header.")
	parser.add_argument("-fe","--overrepresentation-tables-file-end", help=" Specify a common string to look for in the overrepresentation table file names. Default='_overrepresentation_table.txt'", default='_overrepresentation_table.txt')
	parser.add_argument("-c","--column-number-for-significance", help="Specify which column (number) should be checked for significance (the p-adj column). Default='8'", default='8')
	parser.add_argument("-t","--significance-threshold", help="Specify threshold to use to call significant GO overrepresentation (when checking in individual GO overrepresentation files, at the column that was specified using the '-c/--column-number-for-significance option). Default='0.05'", default='0.05')
	parser.add_argument("-pt","--significance-p-type", help="Specify a string for the type of significance (p-value, p-adjust, etc.) used for calling significant GO overrepresentation. It will be added to the output file name. Default='padj'", default='padj')
	parser.add_argument("-s1","--group1-species-file", help="[required] Specify the /path/file/name that contains group 1 species labels (e.g. Aratha) and species full name (e.g. Arabidopsis thaliana). One species per line, tab-delimited file (let me know if you need other formats...). The order of the species will be used for building the overrepresentation matrix, so order the file as you wish (e.g. phylogenetically ordered).")
	parser.add_argument("-s2","--group2-species-file", help="[required] Specify the /path/file/name that contains group 2 species labels (e.g. Aratha) and species full name (e.g. Arabidopsis thaliana). One species per line, tab-delimited file (let me know if you need other formats...). The order of the species will be used for building the overrepresentation matrix, so order the file as you wish (e.g. phylogenetically ordered).")
	parser.add_argument("-o","--output-string", help="Specify a string that will be added to the output files. This should reflect the way the genes of interest were selected for the original annotation overrepressentation analysis (for me it was the 10% top expressed genes from each species). Default = '10_top_expressed'", default='10_top_expressed')
	parser.add_argument("-a","--annotation-type", action="store", default="GO_Term", help="[optional] Specify the type of annotation used. Default = 'GO_Term'")
	#
	if len(sys.argv[1:])==0:
		#parser.print_help()
		parser.print_usage() # for just the usage line
		parser.exit()
	args = parser.parse_args()
	#
	#
	### Function for FDR correction ###
	def fdr_corr(pval_list):
		# First we need to make sure that the returned adjust p-values will be in the same original order:
		original_p_order = range(len(pval_list))	# need to record the original order of the p-values, and to return the adjust p-values in the same original order....
		p_zip = zip(pval_list,original_p_order)	# zip the original p-values and their original positions
		p_reverse_sorted = sorted(p_zip, key=lambda tup: tup[0], reverse=True)	# reverse sort the p-values together with original positions
		new_p_order = [x[1] for x in p_reverse_sorted]	# extract back the p-value positions, with the new order (after, zip this list with the new adjust p-values and sort by the postions (back to original positions))
		original_p_reverse_sorted = [x[0] for x in p_reverse_sorted]	# extract back the p-value
		# Now start the actual calculations:
		adj_p_reverse_order = [original_p_reverse_sorted[0]]	# store the adjusted p-values, start with the largest p-value = adjust p-value
		rank = len(original_p_reverse_sorted) - 1	# need the rank of each p-value, starting with the second largest rank (for the second largest p-value)
		previous_adj = original_p_reverse_sorted[0]	# Need to compare the the adjust p-values with the previous adjust-pvalues in the loop - start with the largest p-val and keep update in the loop
		for p in original_p_reverse_sorted[1:]: # loop from the second largest
			p_adj = float(p) * len(original_p_reverse_sorted) / rank
			rank = rank - 1	# update for the next iteration...
			adj_p_reverse_order.append(min(p_adj,previous_adj))	# take the min of either the original p-value or the 
			previous_adj = min(p_adj,previous_adj)	# update for the next iteration...
		p_zip_adj = zip(adj_p_reverse_order,new_p_order)	# zip the adjust p-values and their reversed ordered original positions
		p_adj_original_order = sorted(p_zip_adj, key=lambda tup: tup[1])	# sort the adjust p-values to their original positions	
		return [x[0] for x in p_adj_original_order]	# return the list of adjust p-values in their original positions
	#
	### Input ###
	# Group 1 path to overrepresentation files:
	if args.path_to_group1_overrepresentation_tables is not None:
		if str(args.path_to_group1_overrepresentation_tables).endswith('/'):
			group1_padj_GO_files = glob.glob(str(args.path_to_group1_overrepresentation_tables) + '*' + str(args.overrepresentation_tables_file_end))
		else:
			group1_padj_GO_files = glob.glob(str(args.path_to_group1_overrepresentation_tables) + '/*' + str(args.overrepresentation_tables_file_end))
		if len(group1_padj_GO_files) == 0:
			sys.exit("Error: Couldn't find any overrepresentation tables for group 1. Path provided: " + str(args.path_to_group1_overrepresentation_tables) + " and file pattern provided: *" + str(args.overrepresentation_tables_file_end))
	else:
		sys.exit("Error: Please provide a path to a folder containing the " + "*" + str(args.overrepresentation_tables_file_end) + " overrepresentation tables, for all the species in group 1, using the '-p1' option. Type: python overcompare.py -h for more instructions of usage.")
	#
	# Group 2 path to overrepresentation files:
	if args.path_to_group2_overrepresentation_tables is not None:
		if str(args.path_to_group2_overrepresentation_tables).endswith('/'):
			group2_padj_GO_files = glob.glob(str(args.path_to_group2_overrepresentation_tables) + str(args.overrepresentation_tables_file_end))
		else:
			group2_padj_GO_files = glob.glob(str(args.path_to_group2_overrepresentation_tables) + str(args.overrepresentation_tables_file_end))
		if len(group2_padj_GO_files) == 0:
			sys.exit("Error: Couldn't find any overrepresentation tables for group 2. Path provided: " + str(args.path_to_group2_overrepresentation_tables) + " and file pattern provided: *" + str(args.overrepresentation_tables_file_end))
	else:
		sys.exit("Error: Please provide a path to a folder containing the " + "*" + str(args.overrepresentation_tables_file_end) + " overrepresentation tables, for all the species in group 2, using the '-p2' option. Type: python overcompare.py -h for more instructions of usage")
	#
	# Get group 1 species label to name conversion file (species_short_label\tspecies_full_name\n). Mine is ordered according to phylogeny for to get the matrix already ordered accordingly"
	if args.group1_species_file is not None:
		group1_species_label_names_file = str(args.group1_species_file)
	else:
		sys.exit("Error: Please provide a path+file to group 1 species_label to species_name file (ordered as desired for the overrepresentation matrix), using the '-s1' option. Type: python overcompare.py -h for more instructions of usage")
	#
	# Get group 2 species label to name conversion file (species_short_label\tspecies_full_name\n). Mine is ordered according to phylogeny for to get the matrix already ordered accordingly"
	if args.group2_species_file is not None:
		group2_species_label_names_file = str(args.group2_species_file)
	else:
		sys.exit("Error: Please provide a path+file to group 2 species_label to species_name file (ordered as desired for the overrepresentation matrix), using the '-s2' option. Type: python overcompare.py -h for more instructions of usage")
	#
	# Get column number in overrepresentation tables to look for significance (hard coded as ≤ 0.05 -> can add argument to define that (as well as argument for overrepresentation file name pattern, instead of '*GO_Term_overrepresentation_table_padj_0_05.txt')...
	sign_column = int(str(args.column_number_for_significance)) - 1	# reduce by 1 to get the correct column in python
	significance_thresh = float(str(args.significance_threshold))	# get the significance cutoff, to call if a species was overrepresented for particular GO term
	significance_type = str(args.significance_p_type).replace('-','_').replace('.','_') # get the significance type (p-val, p-adj...), to add to the output.
	#
	#
	### Parse overrepresentation tables to get a significant GO terms matrix (p-adj per species, per GO term) ###
	go_term_species_overrepresented_dict = {}	# for each GO term (key), store the species label and p-adj values. If not overrepresented, put 'NA' (later will exclude 'NA' only GO terms before Fisher's Exact test)
	go_term_id_desc_dict = {}	# store the GO ids and description for the output table
	for e in group1_padj_GO_files + group2_padj_GO_files:
		species = e.rsplit('/',1)[1].split('_')[0]
		with open(e,'r') as enrichFile:
			next(enrichFile)	# expecting header file
			for line in enrichFile:
				go_term = line.strip().split('\t')[0]	# expecting first column to be the GO term ID
				go_desc = line.strip().split('\t')[1]	# expecting second column to be the GO term description
				if go_term not in go_term_id_desc_dict:
					go_term_id_desc_dict[go_term] = go_desc
				if float(line.strip().split('\t')[sign_column]) <= significance_thresh:
					if go_term in go_term_species_overrepresented_dict:
						go_term_species_overrepresented_dict[go_term][species] = float(line.strip().split('\t')[sign_column])
					else:
						go_term_species_overrepresented_dict[go_term] = {species:float(line.strip().split('\t')[sign_column])}
				else:
					if go_term in go_term_species_overrepresented_dict:
						go_term_species_overrepresented_dict[go_term][species] = 'NA'
					else:
						go_term_species_overrepresented_dict[go_term] = {species:'NA'}
	#
	### Loop over GO terms and build a GO_term species matrix (for all GO terms with at least 1 overrepresented species), and count the number of overrepresented species per GO term, per group ###
	# Get the species label to full name dictionaries:
	#
	group1_sp_label_name_dict = {}
	group1_sp_label_ordered = []
	with open(group1_species_label_names_file, 'r') as group1_spp_namesFile:
		for sp in group1_spp_namesFile:
			group1_sp_label_name_dict[sp.strip().split('\t')[0]] = sp.strip().split('\t')[1]	# expecting two column tab-delim file
			group1_sp_label_ordered.append(sp.strip().split('\t')[0])	# store the species label in their order
	#
	group2_sp_label_name_dict = {}
	group2_sp_label_ordered = []
	with open(group2_species_label_names_file, 'r') as group2_spp_namesFile:
		for sp in group2_spp_namesFile:
			group2_sp_label_name_dict[sp.strip().split('\t')[0]] = sp.strip().split('\t')[1]	# expecting two column tab-delim file
			group2_sp_label_ordered.append(sp.strip().split('\t')[0])	# store the species label in their order
	# str(args.annotation_type)
	out_mat = str(args.annotation_type) + '\tDescription\t' + '\t'.join([group1_sp_label_name_dict[g1] for g1 in group1_sp_label_ordered]) + '\t' + '\t'.join([group2_sp_label_name_dict[g2] for g2 in group2_sp_label_ordered]) + '\n'
	group_counts_per_enrich_GO_dict = {}	# Store the data for Fisher's Exact test and output table
	tested_go_terms = []	# will append only the Fisher's Exact tested (between group1 and group2) GO terms (to keep the order)
	pvalues = []	# will append the Fisher's Exact p-values (between group1 and group2) for the tested GO terms (to keep the order)
	for g in go_term_species_overrepresented_dict:
		group1_overrepresented_spp = []
		group2_overrepresented_spp = []
		for s in go_term_species_overrepresented_dict[g]:
			if go_term_species_overrepresented_dict[g][s] != 'NA':	# if significant...
				if s in group1_sp_label_ordered:
					group1_overrepresented_spp.append(s)
				elif s in group2_sp_label_ordered:
					group2_overrepresented_spp.append(s)
		#
		num_group1_overrepresented = len(group1_overrepresented_spp)
		num_group2_overrepresented = len(group2_overrepresented_spp)
		prop_group1_overrepresented = float(float(num_group1_overrepresented)/float(len(group1_sp_label_ordered)))
		prop_group2_overrepresented = float(float(num_group2_overrepresented)/float(len(group2_sp_label_ordered)))
		group1_overrepresented_spp_full_names = [group1_sp_label_name_dict[g1] for g1 in group1_overrepresented_spp]
		if len(group1_overrepresented_spp_full_names) == 0:
			group1_overrepresented_spp_full_names.append('NA')
		group2_overrepresented_spp_full_names = [group2_sp_label_name_dict[g2] for g2 in group2_overrepresented_spp]
		if len(group2_overrepresented_spp_full_names) == 0:
			group2_overrepresented_spp_full_names.append('NA')
		if num_group1_overrepresented + num_group2_overrepresented > 0:	# if at least one species is overrepresented for a given GO term, include in analysis (can change that to more stringent criterion....)
			# Run Fisher's Exact test of number of overrepresented species, between the two groups (if I will want to test the proportions, I need to use a Test for Two Proportions (https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_Two_Proportions.pdf)):
			#group_1_prop_with_overrepresented_GO = prop_group1_overrepresented
			#group_1_prop_without_overrepresented_GO = 1 - prop_group1_overrepresented
			#group_2_prop_with_overrepresented_GO = prop_group2_overrepresented
			#group_2_prop_without_overrepresented_GO = 1 - prop_group2_overrepresented
			group_1_number_spp_with_overrepresented_GO = num_group1_overrepresented
			group_1_number_spp_without_overrepresented_GO = len(group1_sp_label_ordered) - num_group1_overrepresented
			group_2_number_spp_with_overrepresented_GO = num_group2_overrepresented
			group_2_number_spp_without_overrepresented_GO = len(group2_sp_label_ordered) - num_group2_overrepresented
			p_val = stats.fisher_exact([[group_1_number_spp_with_overrepresented_GO,group_2_number_spp_with_overrepresented_GO],[group_1_number_spp_without_overrepresented_GO,group_2_number_spp_without_overrepresented_GO]])[1]
			pvalues.append(p_val)
			tested_go_terms.append(g)
			group_counts_per_enrich_GO_dict[g] = {'GO_term':g,'GO_desc':go_term_id_desc_dict[g],'Number_overrepresented_species_Group1':num_group1_overrepresented,'Number_overrepresented_species_Group2':num_group2_overrepresented,'Group1_overrepresented_species':';'.join(group1_overrepresented_spp_full_names),'Group2_overrepresented_species':';'.join(group2_overrepresented_spp_full_names),'Prop_overrepresented_species_Group1':prop_group1_overrepresented,'Prop_overrepresented_species_Group2':prop_group2_overrepresented,'p_val':p_val}
			# Store in Species_GO_overrepresentation_matrix:
			out_mat = out_mat + '\t'.join([g,go_term_id_desc_dict[g]])
			for g1 in group1_sp_label_ordered:
				if g1 in group1_overrepresented_spp:
					out_mat = out_mat + '\t' + str(go_term_species_overrepresented_dict[g][g1])
				else:
					out_mat = out_mat + '\tNA'
			for g2 in group2_sp_label_ordered:
				if g2 in group2_overrepresented_spp:
					out_mat = out_mat + '\t' + str(go_term_species_overrepresented_dict[g][g2])
				else:
					out_mat = out_mat + '\tNA'
			out_mat = out_mat + '\n'
	#
	# Save the Species GO overrepresentation matrix:
	mat_out_name = 'Species_' + str(args.output_string) + '_'  + str(args.annotation_type) + '_overrepresentation_matrix_for_' + significance_type + '_' + str(significance_thresh).replace('.','_') + '.txt'
	with open(mat_out_name,'w') as MatFile:
		MatFile.write(out_mat)
	#
	### Run BH FDR correction ###
	adj_pval = fdr_corr(pvalues)	# output is a list in the same order as pvalues and tested_go_terms
	#
	### Build the output table ###
	#
	out_dict = {}	# store the lines, p-val and p-adj as dict, to sort and output by order of significance:
	for i in range(len(tested_go_terms)):	# add the p-adj
		out_dict[tested_go_terms[i]] = {'line':'\t'.join([tested_go_terms[i],go_term_id_desc_dict[tested_go_terms[i]],str(group_counts_per_enrich_GO_dict[tested_go_terms[i]]['Number_overrepresented_species_Group1']),str(len(group1_sp_label_ordered)),str(group_counts_per_enrich_GO_dict[tested_go_terms[i]]['Number_overrepresented_species_Group2']),str(len(group2_sp_label_ordered)),str(group_counts_per_enrich_GO_dict[tested_go_terms[i]]['p_val']),str(adj_pval[i]),group_counts_per_enrich_GO_dict[tested_go_terms[i]]['Group1_overrepresented_species'],group_counts_per_enrich_GO_dict[tested_go_terms[i]]['Group2_overrepresented_species']]) + '\n','p_val':group_counts_per_enrich_GO_dict[tested_go_terms[i]]['p_val'],'adj_pval':adj_pval[i]}
	#
	# Sort by p_val (because usually many GO terms will have the same adj_pval):
	out_p_val_sorted_top = sorted(out_dict.items(), key = lambda x: x[1]['p_val'])
	#
	# Save Fisher's Exact table:
	out = str(args.annotation_type) + '\tDescription\tNumber_overrepresented_species_group1\tNumber_of_species_group1\tNumber_overrepresented_species_group2\tNumber_of_species_group2\tpval\tpadj\tGroup1_overrepresented_species\tGroup2_overrepresented_species\n'
	for j in range(len(list(out_p_val_sorted_top))):
		out = out + out_p_val_sorted_top[j][1]['line']
	#
	fisher_table_out_name = 'Between_groups_overrepresented_annotations_comparison_' + str(args.output_string) + '_Fisher_Exact_for_' + significance_type + '_' + str(significance_thresh).replace('.','_') + '.txt'
	with open(fisher_table_out_name,'w') as FisherTable:
		FisherTable.write(out)

if __name__ == '__main__': # if this script is run directly (and not imported by another script...) # than run the main
	main()

