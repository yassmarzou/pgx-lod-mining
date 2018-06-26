from SPARQLWrapper import SPARQLWrapper, JSON
import rdflib
import csv
from pymantic import sparql
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

def getIncludedFromExcluded(to_exclude,endpoint):
	server = sparql.SPARQLServer(endpoint)
	result = server.query('select distinct ?p where {?s ?p ?o}')
	to_include = [b['p']['value'] for b in result['results']['bindings']]
	for element in to_exclude:
		to_include.remove(element)
	to_include = [(p,'') for p in to_include]
	return to_include


def joinPairs(folder):
	# with open('outputs/train_test_sets/train_set.tsv') as training:
	# 		training_reader = csv.DictReader(training, dialect="excel-tab")
	# 		rdf= '''<?xml version="1.0" encoding="utf-8" ?>
	# <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
	#          xmlns:ns0="http://bio2rdf.org/pharmgkb_vocabulary:">
	# '''
	# 		for pair in training_reader:
	# 			label = 'positive' if pair['association']=='associated' else 'negative'
	# 			rdf +='''
	#   <rdf:Description rdf:about="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''">
	#     <ns0:train_set>'''+label+'''</ns0:train_set>
	#     <ns0:gene rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'''"/>
	#     <ns0:drug rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''"/>
	#   </rdf:Description>
	# '''
	# 		with safe_open_w('outputs/'+folder+'/pairs.rdf') as out:
	# 			out.write(rdf)
	with open('outputs/train_test_sets/test_set.tsv') as training:
			training_reader = csv.DictReader(training, dialect="excel-tab")
			rdf= '''<?xml version="1.0" encoding="utf-8" ?>
	<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
	         xmlns:ns0="http://bio2rdf.org/pharmgkb_vocabulary:">
	'''
			for pair in training_reader:
				rdf +='''
	  <rdf:Description rdf:about="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''">
	    <ns0:test_set>unknown</ns0:test_set>
	    <ns0:gene rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'''"/>
	    <ns0:drug rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''"/>
	  </rdf:Description>
	'''
			with safe_open_w('outputs/'+folder+'/pairs.rdf') as out:
				out.write(rdf)
	for x in range(10):		
		with open('outputs/train_test_sets/train_subset_{0}.tsv'.format(x)) as training:
			training_reader = csv.DictReader(training, dialect="excel-tab")
			rdf = ''
			for pair in training_reader:
				label = 'positive' if pair['association']=='associated' else 'negative'
				rdf +='''
	  <rdf:Description rdf:about="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''">
	    <ns0:train_subset_{0}>'''.format(x)+label+'''</ns0:train_subset_{0}>
	    <ns0:gene rdf:resource="http://bio2rdf.org/pharmgkb:'''.format(x)+pair['gene']+'''"/>
	    <ns0:drug rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''"/>
	  </rdf:Description>
	'''
			with open('outputs/'+folder+'/pairs.rdf','a') as out:
				out.write(rdf)
	with open('outputs/train_test_sets/train_subset_final.tsv') as training:
			training_reader = csv.DictReader(training, dialect="excel-tab")
			rdf = ''
			for pair in training_reader:
				label = 'positive' if pair['association']=='associated' else 'negative'
				rdf +='''
	  <rdf:Description rdf:about="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''">
	    <ns0:train_subset_final>'''+label+'''</ns0:train_subset_final>
	    <ns0:gene rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['gene']+'''"/>
	    <ns0:drug rdf:resource="http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''"/>
	  </rdf:Description>
	'''
	  		rdf +='''
	</rdf:RDF>'''
			with open('outputs/'+folder+'/pairs.rdf','a') as out:
				out.write(rdf)


def createUndirectedSubgraph(to_include,endpoint):
	print 'number of included predicates: ', len(to_include)
	exceptions = {}
	wrapper = SPARQLWrapper(endpoint)

	for i,p in enumerate(to_include):
		print i, 'querying : '+p[0].encode('utf-8')
		inverse = p[0]+'_inverse'
		try:
			wrapper.setQuery('construct {?s <'+p[0]+'> ?o . ?o <'+inverse+'> ?s} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
			result = wrapper.query()._convertN3()
		except Exception as ex:
			print ex.msg
			exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
			wrapper.setQuery('construct {?s <http://bio2rdf.org/predicate_number_'+str(i)+'> ?o . ?o <http://bio2rdf.org/predicate_number_'+str(i)+'_inverse> ?s} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
			result = wrapper.query()._convertN3()
		
		with safe_open_w('outputs/undirected_subgraph/{0}.rdf'.format(i)) as out:
			out.write(result)
	with safe_open_w('outputs/undirected_subgraph/undirected_subgraph.log') as out:
		out.write('**** Exceptions ****\n\n')
		for key in exceptions:
			out.write(p+' --> '+exceptions[key]+'\n')
		out.write('**** Included predicates ****\n\n')
		for i,p in enumerate(to_include):
			out.write(str(i)+'\t'+p[0].encode('utf-8')+'\t'+p[1].encode('utf-8')+'\n')


def createDirectedSubgraph(to_include,endpoint):
	print 'number of included predicates: ', len(to_include)
	exceptions = {}
	wrapper = SPARQLWrapper(endpoint)

	for i,p in enumerate(to_include):
		print i, 'querying : '+p[0].encode('utf-8')
		try:
			wrapper.setQuery('construct {?s <'+p[0]+'> ?o .} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
			result = wrapper.query()._convertN3()
		except Exception as ex:
			print ex.msg
			exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
			wrapper.setQuery('construct {?s <http://bio2rdf.org/predicate_number_'+str(i)+'> ?o .} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
			result = wrapper.query()._convertN3()
		
		with safe_open_w('outputs/directed_subgraph/{0}.rdf'.format(i)) as out:
			out.write(result)
	with safe_open_w('outputs/directed_subgraph/directed_subgraph.log') as out:
		out.write('**** Exceptions ****\n\n')
		for key in exceptions:
			out.write(p+' --> '+exceptions[key]+'\n')
		out.write('**** Included predicates ****\n\n')
		for i,p in enumerate(to_include):
			out.write(str(i)+'\t'+p[0].encode('utf-8')+'\t'+p[1].encode('utf-8')+'\n')

def createUndirectedFullgraph(to_include,endpoint):
	print 'number of included predicates: ', len(to_include)
	exceptions = {}
	wrapper = SPARQLWrapper(endpoint)

	for i,p in enumerate(to_include):
		print i, 'querying : '+p[0].encode('utf-8')
		inverse = p[0]+'_inverse'
		try:
			wrapper.setQuery('construct {?s <'+p[0]+'> ?o . ?o <'+inverse+'> ?s} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
			result = wrapper.query()._convertN3()
		except Exception as ex:
			print ex.msg
			exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
			wrapper.setQuery('construct {?s <http://bio2rdf.org/predicate_number_'+str(i)+'> ?o . ?o <http://bio2rdf.org/predicate_number_'+str(i)+'_inverse> ?s} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
			result = wrapper.query()._convertN3()
		
		with safe_open_w('outputs/undirected_fullgraph/{0}.rdf'.format(i)) as out:
			out.write(result)
	with safe_open_w('outputs/undirected_fullgraph/undirected_fullgraph.log') as out:
		out.write('**** Exceptions ****\n\n')
		for key in exceptions:
			out.write(p[0]+' --> '+exceptions[key]+'\n')
		out.write('**** Included predicates ****\n\n')
		for i,p in enumerate(to_include):
			out.write(str(i)+'\t'+p[0].encode('utf-8')+'\t'+p[1].encode('utf-8')+'\n')


def createDirectedFullgraph(to_include,endpoint):
	print 'number of included predicates: ', len(to_include)
	exceptions = {}
	wrapper = SPARQLWrapper(endpoint)

	for i,p in enumerate(to_include):
		print i, 'querying : '+p[0].encode('utf-8')
		try:
			wrapper.setQuery('construct {?s <'+p[0]+'> ?o .} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
			result = wrapper.query()._convertN3()
		except Exception as ex:
			print ex.msg
			exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
			wrapper.setQuery('construct {?s <http://bio2rdf.org/predicate_number_'+str(i)+'> ?o .} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
			result = wrapper.query()._convertN3()
		
		with safe_open_w('outputs/directed_fullgraph/{0}.rdf'.format(i)) as out:
			out.write(result)
	with safe_open_w('outputs/directed_fullgraph/directed_fullgraph.log') as out:
		out.write('**** Exceptions ****\n\n')
		for key in exceptions:
			out.write(p[0]+' --> '+exceptions[key]+'\n')
		out.write('**** Included predicates ****\n\n')
		for i,p in enumerate(to_include):
			out.write(str(i)+'\t'+p[0].encode('utf-8')+'\t'+p[1].encode('utf-8')+'\n')

def createCustomSubgraph(to_include, to_inverse, endpoint):
	print 'number of included predicates: ', len(to_include)
	exceptions = {}
	wrapper = SPARQLWrapper(endpoint)

	for i,p in enumerate(to_include):
		print i, 'querying : '+p[0].encode('utf-8')
		inverse = p[0]+'_inverse'
		if p in to_inverse:
			try:
				wrapper.setQuery('construct {?o <'+inverse+'> ?s} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
				result = wrapper.query()._convertN3()
			except Exception as ex:
				print ex.msg
				exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
				wrapper.setQuery('construct {?o <http://bio2rdf.org/predicate_number_'+str(i)+'_inverse> ?s} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
				result = wrapper.query()._convertN3()
		else:
			try:
				wrapper.setQuery('construct {?s <'+p[0]+'> ?o .} WHERE {?s <'+p[0]+'> ?o .'+p[1]+'}')
				result = wrapper.query()._convertN3()
			except Exception as ex:
				print ex.msg
				exceptions[p]='http://bio2rdf.org/predicate_number_'+str(i)
				wrapper.setQuery('construct {?s <http://bio2rdf.org/predicate_number_'+str(i)+'> ?o .} WHERE {?s <'+p[0]+'> ?o . '+p[1]+'}')
				result = wrapper.query()._convertN3()

		with safe_open_w('outputs/custom_subgraph/{0}.rdf'.format(i)) as out:
			out.write(result)
	with safe_open_w('outputs/custom_subgraph/custom_subgraph.log') as out:
		out.write('**** Exceptions ****\n\n')
		for key in exceptions:
			out.write(p+' --> '+exceptions[key]+'\n')
		out.write('**** Included predicates ****\n\n')
		for i,p in enumerate(to_include):
			out.write(str(i)+'\t'+p[0].encode('utf-8')+'\t'+p[1].encode('utf-8')+'\n')

# def createDirectedFullgraph():
# 	fullgraph = rdflib.Graph()
# 	with open('outputs/training_pairs.tsv') as training:
# 		training_reader = csv.DictReader(training, dialect="excel-tab")
# 		index = 0
# 		for pair in training_reader:
# 			label = 'positive' if pair['association']=='associated' else 'negative'
# 			sparql = SPARQLWrapper("http://localhost:9999/blazegraph/sparql")
# 			sparql.setQuery('''
# PREFIX gas: <http://www.bigdata.com/rdf/gas#>
# construct { <http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''> <http://bio2rdf.org/pharmgkb_vocabulary:association> "'''+label+'''" .
# 			<http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''> <http://bio2rdf.org/pharmgkb_vocabulary:gene> <http://bio2rdf.org/pharmgkb:'''+pair['gene']+'''>.
# 			<http://bio2rdf.org/pharmgkb:'''+pair['gene']+'_'+pair['drug']+'''> <http://bio2rdf.org/pharmgkb_vocabulary:drug> <http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''>.
# 			?out ?p ?o }
# {
#   {SERVICE gas:service {
#      gas:program gas:gasClass "com.bigdata.rdf.graph.analytics.BFS" .
#      gas:program gas:in <http://bio2rdf.org/pharmgkb:'''+pair['gene']+'''> . # one or more times, specifies the initial frontier.
#      gas:program gas:out ?out . # exactly once - will be bound to the visited vertices.
#      gas:program gas:out1 ?depth . # exactly once - will be bound to the depth of the visited vertices.
#      gas:program gas:out2 ?predecessor . # exactly once - will be bound to the predecessor.
#      gas:program gas:maxIterations '''+'{0}'.format(17)+''' . # optional limit on breadth first expansion.
#   }
#   ?out ?p ?o} union {SERVICE gas:service {
#      gas:program gas:gasClass "com.bigdata.rdf.graph.analytics.BFS" .
#      gas:program gas:in <http://bio2rdf.org/pharmgkb:'''+pair['drug']+'''> . # one or more times, specifies the initial frontier.
#      gas:program gas:out ?out . # exactly once - will be bound to the visited vertices.
#      gas:program gas:out1 ?depth . # exactly once - will be bound to the depth of the visited vertices.
#      gas:program gas:out2 ?predecessor . # exactly once - will be bound to the predecessor.
#      gas:program gas:maxIterations '''+'{0}'.format(17)+''' . # optional limit on breadth first expansion.
#   }
#   ?out ?p ?o}
# }
# ''')
# 			index +=1
# 			print 'Querying '+pair['gene']+'_'+pair['drug']+', index :'+str(index)
# 			pair_data = sparql.query()._convertN3()
# 			fullgraph.parse(data=pair_data)
# 	print 'Serializing...'
# 	graphfilename = "outputs/directed_fullgraph.rdf"
# 	fullgraph.serialize(destination=graphfilename)


to_include = [
('http://bio2rdf.org/pharmgkb_vocabulary:x-ncbigene','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-uniprot','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-gene','') ,
('http://bio2rdf.org/clinvar_vocabulary:Variant_Gene','') ,
('http://bio2rdf.org/clinvar_vocabulary:assertion','') ,
('http://bio2rdf.org/clinvar_vocabulary:Variant_Phenotype','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-medgen','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-sequence_ontology','') ,
('http://semanticscience.org/resource/SIO_000062','') ,
('http://semanticscience.org/resource/SIO_000095','') ,
('http://semanticscience.org/resource/SIO_000008','') ,
('http://biodb.jp/mappings/medispan_to_sider','') ,
('http://biodb.jp/mappings/clinvar_to_medispan','') ,
('http://biodb.jp/mappings/clinvar_to_sider','') ,
('http://bio2rdf.org/sider_vocabulary:side-effect','') ,
('http://bio2rdf.org/sider_vocabulary:indication','') ,
('http://orpailleur.fr/medispan/side_effect','') ,
('http://orpailleur.fr/medispan/indication','') ,
('http://bio2rdf.org/sider_vocabulary:pubchem-flat-compound-id','') ,
('http://bio2rdf.org/sider_vocabulary:pubchem-stereo-compound-id','') ,
('http://biodb.jp/mappings/medispan_to_umls','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-umls','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-pubchemcompound','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-pubchemcompound','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-uniprot','') ,
('http://bio2rdf.org/drugbank_vocabulary:target','') ,
('http://bio2rdf.org/drugbank_vocabulary:carrier','') ,
('http://bio2rdf.org/drugbank_vocabulary:enzyme','') ,
('http://bio2rdf.org/drugbank_vocabulary:transporter','') ,
('http://bio2rdf.org/drugbank_vocabulary:action','') ,
('http://bio2rdf.org/drugbank_vocabulary:drug','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-atc','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-atc','') ,
('http://bio2rdf.org/drugbank_vocabulary:category','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-pharmgkb','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-drugbank','') ,
('http://bio2rdf.org/bio2rdf_vocabulary:x-identifiers.org','?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C16612> .') ,
('http://semanticscience.org/resource/SIO_000628','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> .') ,
('http://semanticscience.org/resource/SIO_000001','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> .') ,
('http://www.w3.org/2004/02/skos/core#exactMatch','?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C7057> . filter ( regex(str(?o),"http://bio2rdf.org/medgen:") || regex(str(?o),"http://orpailleur.fr/medispan/") || regex(str(?o),"http://bio2rdf.org/umls:") )') ,
]

to_keep = [
('http://bio2rdf.org/pharmgkb_vocabulary:x-ncbigene','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-uniprot','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-gene','') ,
('http://bio2rdf.org/clinvar_vocabulary:Variant_Gene','') ,
('http://bio2rdf.org/clinvar_vocabulary:assertion','') ,
('http://bio2rdf.org/clinvar_vocabulary:Variant_Phenotype','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-medgen','') ,
('http://bio2rdf.org/clinvar_vocabulary:x-sequence_ontology','') ,
('http://semanticscience.org/resource/SIO_000062','') ,
('http://semanticscience.org/resource/SIO_000095','') ,
('http://semanticscience.org/resource/SIO_000008','') ,
('http://biodb.jp/mappings/medispan_to_sider','') ,
('http://biodb.jp/mappings/clinvar_to_medispan','') ,
('http://biodb.jp/mappings/clinvar_to_sider','') ,
('http://bio2rdf.org/sider_vocabulary:side-effect','') ,
('http://bio2rdf.org/sider_vocabulary:indication','') ,
('http://orpailleur.fr/medispan/side_effect','') ,
('http://orpailleur.fr/medispan/indication','') ,
('http://bio2rdf.org/sider_vocabulary:pubchem-flat-compound-id','') ,
('http://bio2rdf.org/sider_vocabulary:pubchem-stereo-compound-id','') ,
('http://biodb.jp/mappings/medispan_to_umls','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-umls','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-pubchemcompound','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-pubchemcompound','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-uniprot','') ,
('http://bio2rdf.org/drugbank_vocabulary:target','') ,
('http://bio2rdf.org/drugbank_vocabulary:carrier','') ,
('http://bio2rdf.org/drugbank_vocabulary:enzyme','') ,
('http://bio2rdf.org/drugbank_vocabulary:transporter','') ,
('http://bio2rdf.org/drugbank_vocabulary:action','') ,
('http://bio2rdf.org/drugbank_vocabulary:drug','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-atc','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-atc','') ,
('http://bio2rdf.org/drugbank_vocabulary:category','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-pharmgkb','') ,
('http://bio2rdf.org/pharmgkb_vocabulary:x-drugbank','') ,
('http://bio2rdf.org/bio2rdf_vocabulary:x-identifiers.org','?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C16612> .') ,
('http://semanticscience.org/resource/SIO_000628','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> . ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C16612> .') ,
('http://semanticscience.org/resource/SIO_000628','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> . ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C7057> .') ,
('http://semanticscience.org/resource/SIO_000001','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> .') ,
('http://www.w3.org/2004/02/skos/core#exactMatch','?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C7057> . filter ( regex(str(?o),"http://bio2rdf.org/medgen:") || regex(str(?o),"http://orpailleur.fr/medispan/") || regex(str(?o),"http://bio2rdf.org/umls:") )') ,
]

to_inverse = [
('http://bio2rdf.org/drugbank_vocabulary:drug','') ,
('http://bio2rdf.org/drugbank_vocabulary:x-pharmgkb','') ,
('http://biodb.jp/mappings/medispan_to_umls','') ,
('http://semanticscience.org/resource/SIO_000628','?s <http://semanticscience.org/resource/SIO_000253>/<http://purl.org/ontology/wi/core#evidence> <http://rdf.disgenet.org/v3.0.0/void/source_evidence_curated> . ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C16612> .') ,
('http://bio2rdf.org/clinvar_vocabulary:x-gene','') ,
('http://bio2rdf.org/clinvar_vocabulary:Variant_Gene','') ,
]


#createUndirectedSubgraph(to_include,'http://localhost:9999/blazegraph/sparql')
#joinPairs('undirected_subgraph')
#createDirectedSubgraph(to_include,'http://localhost:9999/blazegraph/sparql')
#joinPairs('directed_subgraph')
createCustomSubgraph(to_keep,to_inverse,'http://localhost:9999/blazegraph/sparql')
joinPairs('custom_subgraph')

#all_uris = getIncludedFromExcluded([],'http://localhost:9999/blazegraph/sparql')
# createUndirectedFullgraph(all_uris,'http://localhost:9999/blazegraph/sparql')
# joinPairs('undirected_fullgraph')
# createDirectedFullgraph(all_uris,'http://localhost:9999/blazegraph/sparql')
# joinPairs('directed_fullgraph')


# NOTE : to compress output folder, you can use lrztar, example : lrztar -z directed_subgraph/