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
### DATE CREATED: June 22, 2010
### AUTHOR: Gil Eshel
### CONTACT1: ge30@nyu.edu
### CONTACT2: giltu1@gmail.com
###
### CITATION: To be added
### 
### Code for parsing the obo file was adopted from Damian Kao (http://blog.nextgenetics.net/?e=6)
### I modified it to include both is_a and part_of as parent GO_ids. Thank you Damian Kao for a great post!
###
#########################################################################################
### DESCRIPTION:
### 
### usage: python add_parent_GOs.py [-h] [-i INPUT_GENE2GO_FILE]
###                                 [-d GENE2GO_DELIMITER]
###                                 [-gd GENE2GO_GOIDS_DELIMITER]
###                                 [-o OUTPUT_GENE2GO_FILE] [-r]
### 
### A python script to include parental GO terms (based on is_a and part_of
### relationships) to a provided gene2go annotation file containing the genes and
### their GO term ids. This is a recommended step before running a GO
### overrepressentation analysis using OverRep, to ensure a proper count of genes
### per GO term). It will also produce a GO_description file (GO id and its
### description), and GO_category file (GO id and the GO spece/lineage, e.g.
### 'Biological Process').
### 
### optional arguments:
###   -h, --help            show this help message and exit
###   -i INPUT_GENE2GO_FILE, --input-gene2go-file INPUT_GENE2GO_FILE
###                         [required] Specify a file containing a gene2go mapping
###                         file. Expecting two columns: (I) Gene_ID, (II) list of
###                         GO ids (delimiters can be defined using the -gd and
###                         -ad options). Not expecting file header (see
###                         OverRep/Example_files/Aratha_gene2go.txt as an
###                         example).
###   -d GENE2GO_DELIMITER, --gene2go-delimiter GENE2GO_DELIMITER
###                         Indicate the delimiter seperating the two columns of
###                         the input-gene2go-file [e.g ' ', ',', '\t']. For tab delimiter
###                         , type $'\t'. Default is '\t')
###   -gd GENE2GO_GOIDS_DELIMITER, --gene2go-goids-delimiter GENE2GO_GOIDS_DELIMITER
###                         Indicate the delimiter seperating the GO ids (second
###                         column) in the input-gene2go-file [e.g ' ', ',', '\t',
###                         ';'], should be different than the --gene2go-delimiter
###                         delimiter. For tab delimiter, type $'\t'. Default is ',')
###   -o OUTPUT_GENE2GO_FILE, --output-gene2go-file OUTPUT_GENE2GO_FILE
###                         Specify an output file name for the gene2go mapping
###                         that will include parental terms. If not specified, it
###                         will add '_with_parents' to the input file name.
###   -r, --remove-obo      A flag to delete the 'go-basic.obo' file upon
###                         complition. The go-basic.obo file is quite large, so
###                         deleting it will be a good idea.
### 
#########################################################################################
#
#
import sys, subprocess, itertools, os, argparse
#
def main():
	###	Define the arguments ###
	parser = argparse.ArgumentParser(description="A python script to include parental GO terms (based on is_a and part_of relationships) to a provided gene2go annotation file containing the genes and their GO term ids. This is a recommended step before running a GO overrepressentation analysis using OverRep, to ensure a proper count of genes per GO term). It will also produce a GO_description file (GO id and its description), and GO_category file (GO id and the GO spece/lineage, e.g. 'Biological Process').", prog="python add_parent_GOs.py")
	parser.add_argument("-i","--input-gene2go-file", help="[required] Specify a file containing a gene2go mapping file. Expecting two columns: (I) Gene_ID, (II) list of GO ids (delimiters can be defined using the -gd and -ad options). Not expecting file header (see OverRep/Example_files/Aratha_gene2go.txt as an example).")	
	parser.add_argument("-d","--gene2go-delimiter",help="Indicate the delimiter seperating the two columns of the input-gene2go-file [e.g ' ', ',', '\\t']. For tab delimiter, type $'\\t'. Default is '\\t')", default="\t")
	parser.add_argument("-gd","--gene2go-goids-delimiter",help="Indicate the delimiter seperating the GO ids (second column) in the input-gene2go-file [e.g ' ', ',', '\\t', ';'], should be different than the --gene2go-delimiter delimiter. For tab delimiter, type $'\\t'. Default is ',')", default=",")
	parser.add_argument("-o","--output-gene2go-file", help="Specify an output file name for the gene2go mapping that will include parental terms. If not specified, it will add '_with_parents' to the input file name.")
	parser.add_argument("-r","--remove-obo", action='store_true', help=" A flag to delete the 'go-basic.obo' file upon complition. The go-basic.obo file is quite large, so deleting it will be a good idea.")
	#
	if len(sys.argv[1:])==0:
		parser.print_usage()
		parser.exit()
	args = parser.parse_args()
	#
	# Check that input file was provided:
	if args.input_gene2go_file is not None:
		print('Provided input file: ' + str(args.input_gene2go_file))
	else:
		sys.exit("Error: Please provide a gene to GO annotation mapping file (see OverRep/Example_files/Aratha_gene2go.txt as an example), using the '-i' option, and the delimiter options '-d' and 'gd'. No header line is expected. Type: python add_parent_GOs.py -h for more instructions of usage")
	#
	### get most recent go-basic.obo file
	if not os.path.exists('go-basic.obo'):
		try:
			subprocess.call(" ".join(["wget", "http://purl.obolibrary.org/obo/go/go-basic.obo"]), shell=True)
		except:
			print("go-basic.obo could not be download")
	
		print('downloaded go-basic.obo file')
	#
	### parse the go-basic.obo to retrieve parent GO terms:
	# few parsing functions:
	def getTerm(stream):
		block = []
		for line in stream:
			if line.strip() == "[Term]" or line.strip() == "[Typedef]":
				break
			else:
				if line.strip() != "":
					block.append(line.strip())
		return(block)
	#
	def parseTagValue(term):
		data = {}
		for line in term:
			if line.startswith('relationship: part_of'):
				tag = line.split(' ')[1] # tag = 'part_of'
				value = line.split(' ')[2] # collect the GO_id
			elif line.startswith('is_a'):
				tag = line.split(': ',1)[0] # tag = 'is_a'
				value = line.split(': ',1)[1].split()[0] # collect the GO_id
			else:
				tag = line.split(': ',1)[0] # other tags like 'id:', 'name:', 'namespace:' etc. are collected
				value = line.split(': ',1)[1]
			if tag not in data:
				data[tag] = []
			data[tag].append(value)
		return(data)
	#
	# A recursive function that take a GO term id and look for all parents (full path until the root) in the terms dictionary (created when we parsed the go-basic.obo file).
	# Return a list of parent ids for the input GO id:
	def find_parents(go_id, terms_dict, go_term_set=[], ret=False):
		for term2 in terms_dict[go_id]['p']:
			#print term2
			go_term_set.append(term2)
			# Recurse on term to find all parents
			if term2 == 'GO:0003674' or term2 == 'GO:0008150' or term2 == 'GO:0005575':
				ret=True
			else:
				find_parents(term2, terms_dict, go_term_set)          
		return(go_term_set)
	#
	#
	oboFile = open('go-basic.obo','r')
	print('Parsing the go-basic.obo file')
	#
	terms = {} # declare a blank dictionary - keys are the GOids and the values will be their parent and Child GO_ids
	getTerm(oboFile) #skip the file header lines
	go_names_dict = {}	# to store the 'GO_id<\t>GO_term' info
	go_cat_dict = {}	# to store the GO category 'GO_id<\t>GO_category' info (GO_categories: biological_process (BP), molecular_function (MF), cellular_component (CC))
	while 1: #infinite loop to go through the obo file. Breaks when the term returned is empty, indicating end of file
		term = parseTagValue(getTerm(oboFile)) #get the term using the two parsing functions
		if len(term) != 0:
			termID = term['id'][0]
			multi_termID = [termID]
			if 'alt_id' in term:
				multi_termID.extend(term['alt_id'])
			termParents = []
			for i in list(term):	# loop through the term fields
				if i == 'is_a' or i == 'part_of':
					termParents.extend(term[i])
			for t in multi_termID:
				if t not in terms:
					terms[t] = {'p':[],'c':[]} # each GO_id will have two arrays of parents and children
				terms[t]['p'] = termParents # append parents of the current term
				for j in termParents: #for every parent term, add this current term as children
					if j not in terms:
						terms[j] = {'p':[],'c':[]}
					terms[j]['c'].append(t)
			#elif term.has_key('name'):
			for g in multi_termID:
				if g.startswith('GO:'):
					if g not in go_names_dict:
						go_names_dict[g] = ''.join(g) + '\t' + ''.join(term['name'])	# store the 'GO_id<\t>GO_term' info, e.g. 'GO:0000008<\t>obsolete thioredoxin' 
						go_cat_dict[g] = ''.join(g) + '\t' + ''.join(term['namespace'])	# store the 'GO_id<\t>GO_category' info, e.g. 'GO:0000008<\t>molecular_function' 
		else:
			break
	#
	oboFile.close()
	#
	# For each GO term, get a list of all parent GO terms (until the root):
	all_parent_terms = {}
	for g in terms:
		if len(terms[g]['p']) == 0:
			all_parent_terms[g] = ''
		elif g == 'negatively_regulates' or g == 'positively_regulates':
			continue
		else:
			all_parent_terms[g] = list(set(find_parents(g,terms,go_term_set=[], ret=False)))
	#
	print('Parsing input gene2go file')
	#
	### Parse the user provided gene2go.txt file to add the parent GO_ids:
	gene_to_go = [] # a list to store the gene<\t>GO_all<\n> strings
	check_file_sep = 'fine'	# Check that the delimiter make sense, if not, raise a warning.
	check_GO_sep = ''	# Check that the delimiter make sense, if not, raise a warning.
	if args.input_gene2go_file is not None:
		with open(str(args.input_gene2go_file), 'r') as goFile:
			for line in goFile:
				li=line.strip()
				if len(li.split(str(args.gene2go_delimiter)))!= 2:	# if the separator didn't split in exactly 2, in at least one line, something is probably wrong...
					check_file_sep = 'not_fine'
				if li.count(str(args.gene2go_goids_delimiter)) > 0:	# if we couldn't find the seperation in any of the lines, it is probably wrong...
					check_GO_sep = 'fine'
				gene = li.split(str(args.gene2go_delimiter))[0]
				GO = li.split(str(args.gene2go_delimiter))[1].split(str(args.gene2go_goids_delimiter))
				GO_parents = []
				for go in GO:
					if go in all_parent_terms:
						GO_parents.extend(all_parent_terms[go])
				GO_all = str(args.gene2go_goids_delimiter).join(list(set(','.join(GO).split(',') + GO_parents))) # merge the GO list (original) with the GO_parents list to a unique list of GO_ids - join them to a comma-separated string
				gene_to_go.append(gene + str(args.gene2go_delimiter) + GO_all)	
	#
	if check_file_sep == 'not_fine':
		print('WARNING: Check your input file and/or your command. At least one line was not seperated using the provided seperator --gene2go-delimiter (' + str(args.gene2go_delimiter)  + ').')
	#
	if check_GO_sep != 'fine':
		print('WARNING: Check your input file and/or your command. The provided --gene2go-goids-delimiter seperator (' + str(args.gene2go_goids_delimiter)  + ') was not found in any of the lines.')
	#	
	#
	if args.output_gene2go_file is not None:
		outname = str(args.output_gene2go_file)
	else:
		outname = str(args.input_gene2go_file).rsplit('.',1)[0] + '_with_parents.' + str(args.input_gene2go_file).rsplit('.',1)[1]
	#
	print('Saving the gene2go file with parental GO terms as ' + outname)
	### Save the annotation file with parental GO_ids:
	with open(outname, 'w') as output:
		output.write('\n'.join(gene_to_go) + '\n')
	#
	print('Saving GO description and GO catergory files')
	### open an output file to print into the GO_id to GO description file
	og_desc = ''
	for d in go_names_dict:
		og_desc = og_desc + go_names_dict[d] + '\n'
	#
	with open(outname.rsplit('.',1)[0] + '_GODesc.txt', 'w') as output:
		output.write(og_desc)
	#
	### open an output file to print into the GO_id to GO category file
	og_cat = ''
	for c in go_cat_dict:
		og_cat = og_cat + go_cat_dict[c] + '\n'
	#
	with open(outname.rsplit('.',1)[0] + '_GOCategory.txt', 'w') as output:
		output.write(og_cat)
	#
	### Delete the go-basic.obo file (very large...):
	if str(args.remove_obo) == 'True':
		print('Deleting go-basic.obo from working directory')
		os.remove('go-basic.obo')
	#
	print('Done')

if __name__ == '__main__':
	main()
