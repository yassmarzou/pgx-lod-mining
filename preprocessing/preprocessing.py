import sys
import csv
import re
from operator import itemgetter
import errno

import json
import requests

import os, os.path
import errno
# Taken from http://stackoverflow.com/a/600612/119527
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def safe_open_w(path):
    ''' Open "path" for writing, creating any parent directories as needed.
    '''
    mkdir_p(os.path.dirname(path))
    return open(path, 'w')


csv.field_size_limit(sys.maxsize)

# filenames
original_annotations = 'inputs/clinical_ann_metadata.tsv'
original_relationships = 'inputs/relationships.tsv'
simplified_annotations = 'temp/clinical_ann_simplified.tsv'
simplified_relationships = 'temp/relationships_simplified.tsv'
annotations_with_associations = 'temp/annotations_with_associations.tsv'
training_set = 'outputs/train_test_sets/train_set.tsv'
test_set = 'outputs/train_test_sets/test_set.tsv'
negatives = 'outputs/negative_pairs.tsv'
genes_tsv = "inputs/genes.tsv"
drugs_tsv = "inputs/drugs.tsv"
dgidb_output = "temp/dgidb_output.tsv"
interactions_file = "inputs/interactions.tsv"
api_path = "http://dgidb.genome.wustl.edu/api/v1/interactions.json"



# remove non necessary fields. Output : clinical_ann_simplified.tsv
def simplify_annotaions(input_file, output_file):
	# regex \(([^\)]+)\) to match parenthesis.
	with open(input_file) as original, safe_open_w(output_file) as output:
		reader = csv.DictReader(original, dialect="excel-tab")
		fieldnames = ['location', 'gene', 'level', 'drug', 'disease']
		writer = csv.DictWriter(output, dialect='excel-tab', fieldnames=fieldnames)
		writer.writeheader()
		for row in reader:
			for gene in re.findall(r'\(([^\)]+)\)', row['Gene']) or ['']:
				for drug in re.findall(r'\(([^\)]+)\)', row['Related Drugs']) or ['']:
					for disease in re.findall(r'\(([^\)]+)\)', row['Related Diseases']) or ['']:
						writer.writerow({'location': row['Location'], 'gene': gene, 'level': row['Level of Evidence'], 'drug': drug, 'disease': disease })

simplify_annotaions(original_annotations, simplified_annotations)


# remove gene-disease and drug-disease relationships. Output : relationships_simplified.tsv
def simplify_relationships(input_file, output_file):
	with open(input_file) as original, safe_open_w(output_file) as output:
		reader = csv.DictReader(original, dialect="excel-tab")
		fieldnames = ['gene OR variantlocation OR haplotype', 'type', 'drug', 'association']
		writer = csv.DictWriter(output, dialect='excel-tab', fieldnames=fieldnames)
		writer.writeheader()
		for row in reader:
			if row['Entity1_type'] == 'Drug' and row['Entity2_type'] != 'Disease':
				writer.writerow({'gene OR variantlocation OR haplotype': row['Entity2_id'], 'type': row['Entity2_type'], 'drug': row['Entity1_id'], 'association': row['Association'] })
			elif row['Entity2_type'] == 'Drug' and row['Entity1_type'] != 'Disease':
				writer.writerow({'gene OR variantlocation OR haplotype': row['Entity1_id'], 'type': row['Entity1_type'], 'drug': row['Entity2_id'], 'association': row['Association'] })

simplify_relationships(original_relationships, simplified_relationships)


# check relationships_simplified integrity : every entity-drug pair in relationships should have only one association label. should output nothing.
def check_relationships(input_file):
	with open(input_file) as relationships:
		relationships_reader = csv.DictReader(relationships, dialect="excel-tab")
		pairs = {}
		for relationship in relationships_reader:
			if relationship['gene OR variantlocation OR haplotype']+'_'+relationship['drug'] in pairs:
				pairs[relationship['gene OR variantlocation OR haplotype']+'_'+relationship['drug']].append(relationship['association'])
			else:
				pairs[relationship['gene OR variantlocation OR haplotype']+'_'+relationship['drug']] = [relationship['association']]
		for pair in pairs:
			if len(set(pairs[pair]))!=1:
				print pair, pairs[pair]

check_relationships(simplified_relationships)


# add associations to annotations
def cross_reference(annotaions_input, relationships_input, output_file):
	with open(annotaions_input) as annotations, open(relationships_input) as relationships, safe_open_w(output_file) as output:
		annotaions_reader = csv.DictReader(annotations, dialect="excel-tab")
		relationships_reader = csv.DictReader(relationships, dialect="excel-tab")
		fieldnames = ['location', 'gene', 'level', 'drug', 'disease', 'location_association', 'gene_association']
		output_writer = csv.DictWriter(output, dialect='excel-tab', fieldnames=fieldnames)
		output_writer.writeheader()
		# do not forget relationships.seek(0) after looping !
		# previous integrity check is useful for building relationships_dict
		relationships_dict = {relationship['gene OR variantlocation OR haplotype']+'_'+relationship['drug'] : relationship['association'] for relationship in relationships_reader}
		for annotation in annotaions_reader:
			row = {'location': annotation['location'], 'gene': annotation['gene'], 'level': annotation['level'], 'drug': annotation['drug'], 'disease': annotation['disease'], 'location_association': '', 'gene_association': '' }
			if relationships_dict.has_key(annotation['location']+'_'+annotation['drug']):
				row['location_association'] = relationships_dict[annotation['location']+'_'+annotation['drug']]
			if relationships_dict.has_key(annotation['gene']+'_'+annotation['drug']):
				row['gene_association'] = relationships_dict[annotation['gene']+'_'+annotation['drug']]
			output_writer.writerow(row)

cross_reference(simplified_annotations, simplified_relationships, annotations_with_associations)


# check if all gene_drug pairs have a gene_asociation label.
def check_annotations_with_associations(input_file):
	with open(input_file) as annotations:
		annotaions_reader = csv.DictReader(annotations, dialect="excel-tab")
		pair_association = {}
		for annotation in annotaions_reader:
			if annotation['gene'] != '':
				if annotation['gene']+'_'+annotation['drug'] in pair_association:
					pair_association[annotation['gene']+'_'+annotation['drug']].append(annotation['gene_association'])
				else:
					pair_association[annotation['gene']+'_'+annotation['drug']] = [annotation['gene_association']]
		print 'number of pairs with no association label:', len(set(pair for pair in pair_association if '' in pair_association[pair]))

check_annotations_with_associations(annotations_with_associations)


def generate_dgidb_interactions(genes_file, drugs_file, interactions_file, output_file):
	with open(genes_file) as genes, open(drugs_file) as drugs, open(interactions_file) as interactions, safe_open_w(output_file) as output:
		genes_reader = csv.DictReader(genes, dialect="excel-tab")
		drugs_reader = csv.DictReader(drugs, dialect="excel-tab")
		interactions_reader = csv.DictReader(interactions, dialect="excel-tab")
		fieldnames = ['PharmGKB gene Id', 'Symbol', 'PharmGKB drugs Ids', 'Drugs Names']
		output_writer = csv.DictWriter(output, dialect='excel-tab', fieldnames=fieldnames)
		output_writer.writeheader()
		

		drugs_dict = {}
		for drug in drugs_reader:
			drugs_dict[drug['Name'].lower()] = drug['PharmGKB Accession Id']
		interactions_dict = {}
		for interaction in interactions_reader:
			if interaction['entrez_gene_symbol'] in interactions_dict:
				interactions_dict[interaction['entrez_gene_symbol']].add(interaction['drug_primary_name'].lower())
			else:
				interactions_dict[interaction['entrez_gene_symbol']] = {interaction['drug_primary_name'].lower()}
		for gene in genes_reader:
			if gene['Symbol'] in interactions_dict:
				drugs_names = ', '.join(interactions_dict[gene['Symbol']])
				drugs_ids = ', '.join([drugs_dict[drug] for drug in interactions_dict[gene['Symbol']] if drug in drugs_dict])
			else:
				drugs_names = ''
				drugs_ids = ''
			gene_interactions = {'PharmGKB gene Id': gene['PharmGKB Accession Id'], 'Symbol': gene['Symbol'], 'PharmGKB drugs Ids': drugs_ids, 'Drugs Names': drugs_names}
			output_writer.writerow(gene_interactions)
	# with open(genes_file) as genes, open(drugs_file) as drugs, safe_open_w(interactions_file) as interactions:
	# 	genes_reader = csv.DictReader(genes, dialect="excel-tab")
	# 	drugs_reader = csv.DictReader(drugs, dialect="excel-tab")
	# 	fieldnames = ['PharmGKB gene Id', 'Symbol', 'Status', 'PharmGKB drugs Ids', 'Drugs Names']
	# 	interactions_writer = csv.DictWriter(interactions, dialect='excel-tab', fieldnames=fieldnames)
	# 	interactions_writer.writeheader()
		
	# 	drugs_dict = {}
	# 	for drug in drugs_reader:
	# 		drugs_dict[drug['Name'].lower()] = drug['PharmGKB Accession Id']
	# 	for gene in genes_reader:
	# 		response = json.loads(requests.get(api_path+"?genes="+gene['Symbol']).content)
	# 		drugs_names = ''
	# 		drugs_ids = ''
	# 		status = 'found_no_interactions'
	# 		# if (not response['matchedTerms']) and (not response['unmatchedTerms']):
	# 		# 	status = 'found_no_interactions'
	# 		if (response['unmatchedTerms']):
	# 			status = 'not_found'
	# 		elif (response['matchedTerms']):
	# 			drugs_names_set = set()
	# 			drugs_ids_set = set()
	# 			status = 'found_interactions'
	# 			for interaction in response['matchedTerms'][0]['interactions']:
	# 				drugs_names_set.add(interaction['drugName'])
	# 				if interaction['drugName'].lower() in drugs_dict:
	# 					drugs_ids_set.add(drugs_dict[interaction['drugName'].lower()])	
	# 			drugs_names = ', '.join(drugs_names_set)
	# 			drugs_ids = ', '.join(drugs_ids_set)	
	# 		gene_interactions = {'PharmGKB gene Id': gene['PharmGKB Accession Id'], 'Symbol': gene['Symbol'], 'Status': status, 'PharmGKB drugs Ids': drugs_ids, 'Drugs Names': drugs_names}
	# 		print gene_interactions
	# 		interactions_writer.writerow(gene_interactions)


generate_dgidb_interactions(genes_tsv, drugs_tsv, interactions_file, dgidb_output)


# description
def add_instances_to_set(input_file,instances,append):
	l = []
	if append:
		with open(input_file) as pairs:
			pairs_reader = csv.DictReader(pairs, dialect="excel-tab")
			l = list(pairs_reader)
	for instance in instances:
		l.append({'gene': instance[0], 'drug':instance[1], 'association': instance[2]})
	sorted_list = sorted(l, key=itemgetter('gene','drug','association'))
	with safe_open_w(input_file) as pairs:
		fieldnames = ['gene', 'drug', 'association']
		pairs_writer = csv.DictWriter(pairs, dialect='excel-tab', fieldnames=fieldnames)
		pairs_writer.writeheader()
		for element in sorted_list:
			pairs_writer.writerow(element)

def create_datasets(input_file, training_ouput, test_output):
	with open(input_file) as annotations:
		annotaions_reader = csv.DictReader(annotations, dialect="excel-tab")
		training_pairs = set()
		test_pairs = set()
		for annotation in annotaions_reader:
			if annotation['gene'] != '' and annotation['gene_association'] in ['associated','not associated'] and (annotation['gene'], annotation['drug'], annotation['gene_association']) not in training_pairs:
				if annotation['level'] in ['1A', '1B', '2A', '2B']:
					training_pairs.add((annotation['gene'], annotation['drug'], annotation['gene_association']))
				else:
					test_pairs.add((annotation['gene'], annotation['drug'], annotation['gene_association']))
		add_instances_to_set(training_ouput, training_pairs, False)
		add_instances_to_set(test_output, test_pairs, False)


for x in range(10):
	create_datasets(annotations_with_associations, 'outputs/train_test_sets/train_subset_{0}.tsv'.format(x), test_set)
create_datasets(annotations_with_associations, 'outputs/train_test_sets/train_subset_final.tsv', test_set)
#create_datasets(annotations_with_associations, training_set, test_set)



import random
random.seed(0)
def generate_random_negatives(genes_file, drugs_file, interactions_file, training_set, number):
	with open(genes_file) as genes, open(drugs_file) as drugs, open(interactions_file) as interactions:
		
		genes_reader = csv.DictReader(genes, dialect="excel-tab")
		drugs_reader = csv.DictReader(drugs, dialect="excel-tab")
		interactions_reader = csv.DictReader(interactions, dialect="excel-tab")
		genes_set = [gene['PharmGKB Accession Id'] for gene in genes_reader]
		drugs_set = [drug['PharmGKB Accession Id'] for drug in drugs_reader]
		interactions_dict = {interaction['PharmGKB gene Id']:interaction['PharmGKB drugs Ids'] for interaction in interactions_reader}
		res = []
		while len(res)<number:
			gene1 = random.choice(genes_set)
			drug1 = random.choice(drugs_set)
			if drug1 not in interactions_dict[gene1] and gene1+"_"+drug1 not in res:
				res.append(gene1+"_"+drug1)
		instances = [(element.split('_')[0], element.split('_')[1], 'not associated') for element in res]
		add_instances_to_set(training_set, instances, True)

		# fieldnames = ['gene', 'drug', 'association']
		# pairs_writer = csv.DictWriter(negatives, dialect='excel-tab', fieldnames=fieldnames)
		# pairs_writer.writeheader()
		# for element in res:
		# 	pairs_writer.writerow({'gene':element.split('_')[0] , 'drug':element.split('_')[1], 'association':'not associated'})

# with open(training_set) as foo:
#     size = len(foo.readlines())
for x in range(10):
	generate_random_negatives(genes_tsv, drugs_tsv, dgidb_output, 'outputs/train_test_sets/train_subset_{0}.tsv'.format(x), 91)
generate_random_negatives(genes_tsv, drugs_tsv, dgidb_output, 'outputs/train_test_sets/train_subset_final.tsv', 91)




# def generate_all_negatives(genes_file, drugs_file, interactions_file, training_set):
# 	with open(genes_file) as genes, open(drugs_file) as drugs, open(interactions_file) as interactions, open(training_set, 'a') as pairs:
		
# 		fieldnames = ['gene', 'drug', 'association']
# 		pairs_writer = csv.DictWriter(pairs, dialect='excel-tab', fieldnames=fieldnames)
# 		genes_reader = csv.DictReader(genes, dialect="excel-tab")
# 		drugs_reader = csv.DictReader(drugs, dialect="excel-tab")
# 		interactions_reader = csv.DictReader(interactions, dialect="excel-tab")
# 		genes_set = [gene['PharmGKB Accession Id'] for gene in genes_reader]
# 		drugs_set = [drug['PharmGKB Accession Id'] for drug in drugs_reader]
# 		interactions_dict = {interaction['PharmGKB gene Id']:interaction['PharmGKB drugs Ids'] for interaction in interactions_reader}
# 		print(len(genes_set))
# 		print(len(drugs_set))
# 		res = []
# 		i = 0
# 		for gene1 in genes_set:
# 			i += 1
# 			if i%10==0:
# 				print(i)
# 				for drug1 in drugs_set:
# 					if drug1 not in interactions_dict[gene1]:
# 						res.append({'gene': gene1, 'drug':drug1, 'association': 'not associated'})
# 		pairs_writer.writerows(res)

# generate_all_negatives(genes_tsv, drugs_tsv, dgidb_output, training_set)
import shutil
shutil.rmtree('temp')
