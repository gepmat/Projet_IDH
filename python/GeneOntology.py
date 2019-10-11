#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import Graph as gr # Graph library from part 1 of the project
import time

parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-o', '--obo', required=True, help='Obo file')#
parser.add_argument('-g', '--goa', required=True, help='Goa file')
parser.add_argument('-s', '--species', required=True, help='Species name')
param = parser.parse_args()

def load_OBO(filename):
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""
	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'): return
		# otherwise create node
		id = lines.pop(0)[4:].rstrip()
		term = gr.add_node(g,id)
		term['id'] = id
		term['type'] = 'GOTerm'
		for line in lines:
			# attributes (name, namespace, def)
			if line.startswith('name: '): term['name'] = line[6:]
			elif line.startswith('namespace: '): term['namespace'] = line[11:]
			elif line.startswith('def: '): term['def'] = line[5:]
			elif line.startswith('alt_id: '): g['alt_id'][ line[8:] ] = id # alternate ids
			# relationships
			elif line.startswith('is_a:'): # is_a
				parent = line[6:line.index('!')].rstrip()
				e = gr.add_edge(g,id, parent)
				e['type'] = 'is_a'
			elif line.startswith('relationship: part_of '): # part_of
				line = line[line.index('GO:'):]
				dest = line[:line.index(' ')]
				e = gr.add_edge(g,id, dest)
				e['type'] = 'part_of'
	#
	g=gr.create_graph(directed=True, weighted=False)
	g['alt_id'] = {} # alternate GO ids
	with open(filename) as f: 
		line = f.readline().rstrip()
		# skip header to reach 1st Term
		while not line.startswith('[Term]'): 
			line = f.readline().rstrip()
		buff = []
		line = f.readline()
		stop = False
		while line and not stop:
			# buffer lines until the next Term is found
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				# next Term found: create corresponding node and edges in parseTerm and empty buffer
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	return g


def load_GOA(go, filename):
	"""
	parse GOA file and add annotated gene products to previsouly loaded graph go

	Extract of a file to be parsed:
	!gaf-version: 2.1
	!GO-version: http://purl.obolibrary.org/obo/go/releases/2016-10-29/go.owl
	UniProtKB  A5A605  ykfM      GO:0006974  PMID:20128927   IMP              P  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20100901  EcoCyc
	UniProtKB  A5A605  ykfM      GO:0016020  GO_REF:0000037  IEA              C  Uncharacterized protein YkfM    YKFM_ECOLI|ykfM|b4586         protein taxon:83333  20161029  UniProt
	UniProtKB  P00448  sodA      GO:0004784  GO_REF:0000003  IEA  EC:1.15.1.1 F  Superoxide dismutase [Mn]       SODM_ECOLI|sodA|JW3879|b3908  protein taxon:83333  20161029  UniProt
	UniProtKB  P00393  ndh  NOT  GO:0005737  PMID:6784762    IDA              C  NADH dehydrogenase              DHNA_ECOLI|ndh|JW1095|b1109   protein taxon:83333  20100621  EcoliWiki
		0        1       2   3       4             5          6        7      8             9                              10
				 id    name        go_id               evidence-codes                     desc                           aliases
	"""
	names = {}
	go['names'] = names # gene names or gene product names (column 3)
	with open(filename) as f: 
		line = f.readline()
		while line:
			if not line.startswith('!'):
				cols = line.rstrip().split('\t')
				id = cols[1]
				go_id = cols[4]
				if go_id not in go['nodes']: # GOTerm not found search alternate ids
					if go_id in go['alt_id']: # success
						go_id = go['alt_id'][go_id] # replace term
					else: # warn user
						print('Warning: could not attach a gene product (%s) to a non existing GO Term (%s)' % (id, go_id))
				if go_id in go['nodes']:
					# create node for gene product if not already present
					if id not in go['nodes']:
						g = gr.add_node(go,id)
						g['id'] = id
						g['type'] = 'GeneProduct'
						names[cols[2]] = id
					# create or update gene product attributes
					gp = go['nodes'][id]            
					gp['name'] = cols[2]
					gp['desc'] = cols[9]
					gp['aliases'] = cols[10]
					# attach gene product to GOTerm
					go_term = go['nodes'][go_id]
					e = gr.add_edge(go, id, go_id)
					e['type'] = 'annotation'
					if 'evidence-codes' not in e: e['evidence-codes'] = []
					e['evidence-codes'].append( cols[6] )
				else: # go_id or alt_id not found in GOTerms
					print('Error: could not attach a gene product (%s) to non existing GO Term (%s)' % (id, go_id))
			line = f.readline()


def max_depth(Go, node):
	"""
	go graph traversal to find the longest path length from node (GOTerm) to a leaf (node with no successor)
	Returns the length of the longest path
	"""
	go_T = gr.transpose_g(Go) # on le fait sur la transposée car on part des trois ontologies
	dfs = gr.DFS_node_p(go_T, node) #utilisation du DFS modifié avec calcul de la distance
	depth = max(dfs['distance'].values()) - 1 #on prend le max de la distance -1 car on veut le nombre d'edges pas de sommets
	return(depth)

def max_depth2(Go, node):
	"""
	go graph traversal to find the longest path length from node (GOTerm) to a leaf (node with no successor)
	Returns the length of the longest path
	"""
	go_T = gr.transpose_g(Go)
	bfs = gr.BFS2(go_T, node) # utilisation du BFS modifié avec calcul de la distance (comparaison)
	depth = max(bfs['distance'].values()) - 1
	return(depth)
	
	
def recherche_id(Go,gp):
	key  =[]
	for id_gp in Go['nodes']:
		if (Go['nodes'][id_gp]['name'] == gp): #on récupère tous les id du nom de gène en entrée
			key.append(id_gp)
	return(key)
	

def GOTerms(Go, key, all=True, evidence_code=None):
	"""
	return the GOTerms associated to the provided gene product (gp)
	ggo_T['names']['scaf4b']
	ggo['edges']['R4GE22'].keys()
	go: Gene Ontology graph
	gp: gene product
	all: if True, all the GOTerms and their ancestors will be returned, otherwise only the GOTerms directly associated to the gene product will be returned.
	evidence_code: ignored for the moment
	
	Returns a list of GOTerms identifiers, e.g. ['GO:0005215','GO:0005515','GO:0006810','GO:0006974','GO:0008643']
	"""
	list_GO_Terms = []
	liste_id = recherche_id(Go,key)
	for gene in liste_id:
		GO = list(Go['edges'][gene].keys()) # pour tous les id on extend les GOterms liés directement aux id
		list_GO_Terms.extend(GO) # on stack les GOterms
	list_GO_Terms = list(set(list_GO_Terms)) # on vérifie les doublons
	if all == True :
		for GOs in list_GO_Terms:
			dfs = gr.DFS_node(Go, GOs) # si all=True on  fait un DFS à partir de chaque GO directs et on stack ces GO
			for debuts in dfs['debut']:
				if debuts not in list_GO_Terms:
					list_GO_Terms.append(debuts)
		return(list_GO_Terms)
	else :
		return(list_GO_Terms)


def GeneProducts(go_T, term, all=True, evidence_code=None):
	"""
	return the gene products anotated by the provided GOTerm
	
	go: Gene Ontology graph
	term: GOTerm id
	all: if True, all the gene products directly and undirectly annotated (linked to a descendant of GOTerm) will be return, otherwise only the gene products directly associated to the GOTerm will be returned.
	evidence_code: ignored for the moment

	Returns a list of gene products identifiers, e.g. ['P0AAG5', 'P0AFY6', 'P10907', 'P16676', 'P23886']
	"""
	
	 # on fait la transposée pour que les edges pointent vers les Gene products.
	
	identifiers = []
	GO_terms = [term]
	for u in list(go_T['edges'][term]):
		if u.startswith('GO:'):
			GO_terms.append(u) # on stack d'un coté les GO terms et de l'autres les id des Gene products directs
		else:
			identifiers.append(u)
	if all == False:
		return(identifiers) # si on veut les directs on renvoie alors directement les GO trouvés
		
	while GO_terms:
		for u in GO_terms:
			liste2 = list(go_T['edges'][u]) # si all = True on fait une liste d'attente tant qu'il y a des GO liés à d'autres GO et on stack les id 
			GO_terms.pop(0)
			for g in liste2:
				if g.startswith('GO:') and g not in GO_terms:
					GO_terms.append(g)
				elif not(g.startswith('GO:')) and g not in identifiers : # on vérifie de ne pas mettre deux fois le même id de gène 
					identifiers.append(g)

	return(identifiers)
	
 #   (not(GO_terms))==False
def annotation(Go,Term): # cette fonction permet le tri d'une liste de GO term afion de les répartir dans les trois domaines de la GO
	biological_process = []
	molecular_function = []
	cellular_component = []
	for GO in Term:
		if Go['nodes'][GO]['namespace'] == 'biological_process':
			biological_process.append(GO)
		elif Go['nodes'][GO]['namespace'] == 'molecular_function':
			molecular_function.append(GO)
		elif Go['nodes'][GO]['namespace'] == 'cellular_component':
			cellular_component.append(GO)
	return(biological_process,len(biological_process),molecular_function,len(molecular_function),cellular_component,len(cellular_component))
	
	
def liste_trois_ontologies(liste1,liste2,liste3): # on trouve ici les Gene products présents dans les 3 ontologies
	seen = liste1 + liste2 + liste3
	seen = sorted(seen)
	present = []
	for i in range(len(seen)-2):
		   if seen[i] == seen[i+1] and seen[i+1] == seen[i+2]:
					 present.append(seen[i])
	return(list(present))
	
	
def liste_deux_ontologies(liste1, liste2): # on trouve ici les Gene products présents dans 2 ontologies
	seen = liste1 + liste2
	seen = sorted(seen)
	present = []
	for i in range(len(seen)-1):
		if seen[i] == seen[i+1] and seen[i+1] != seen[i+2]:
			present.append(seen[i])
	g = list(set(present))
	return(g)
	
	
def origine_gp(Go,gp): # trouve les fonctions d'un nom de gène ( les différentes implications dans la GO)
	iden = recherche_id(Go,gp)
	fonction=[]
	for i in iden :
		fonction.append(Go['nodes'][i]['desc'])

	nombre =len(fonction)
	return(list(set(fonction)),nombre)

def origine_GO(Go,GO): # trouve le processus et la fonction d'un GO term
   processus = Go['nodes'][GO]['namespace']
   fonction = Go['nodes'][GO]['name']
   return(processus,fonction)



	##### TEST SUR LE GRAPH ENTIER  #######
def projet():
	
	print('__________Loading Fichier OBO...__________')
	obo_file = param.obo
	start_time = time.time()
	GO = load_OBO(obo_file)
	print('Temps d exécution de load_OBO:',time.time()-start_time,'s')
	go = GO
	
	
	goa_file=param.goa
	print('__________Ajout des Genes Products du Zebrafish...__________')
	start_time = time.time()
	load_GOA(go,goa_file)

	go_T = gr.transpose_g(go)
	print('Temps d exécution de load_GOA:',time.time()-start_time,'s')
	
	print('__________Etude Gene Products...__________')
	start_time = time.time()
	Nombre_Products_Biological_Process = GeneProducts(go_T,'GO:0008150')
	print('Temps d exécution de GeneProducts:',time.time()-start_time,'s')
	Nombre_Products_molecular_function = GeneProducts(go_T,'GO:0003674')
	Nombre_Products_cellular_component = GeneProducts(go_T,'GO:0005575')
	total = Nombre_Products_Biological_Process + Nombre_Products_molecular_function + Nombre_Products_cellular_component
	total = list(set(total))
	Gene_Products_total_GOA = 31379 
	print('Le nombre de Gene products annotés sans redondance est',len(total),'et le nombre de Gene Products total est', Gene_Products_total_GOA)
	print('La Proportion de Gene products annotés dans au moins une des trois ontologie est :',len(total)/Gene_Products_total_GOA)
	trois_ontologies = liste_trois_ontologies(Nombre_Products_Biological_Process,Nombre_Products_molecular_function,Nombre_Products_cellular_component)
	print('Il y a', len(trois_ontologies),'Gene Products présents et donc annotés dans les trois ontologies en même temps et la proportion correspondante est :',len(trois_ontologies)/(Gene_Products_total_GOA))
	proportion_process_function = liste_deux_ontologies(Nombre_Products_Biological_Process,Nombre_Products_molecular_function)
	proportion_process_component = liste_deux_ontologies(Nombre_Products_Biological_Process,Nombre_Products_cellular_component)
	proportion_function_component = liste_deux_ontologies(Nombre_Products_cellular_component,Nombre_Products_molecular_function)

	print('Il y a',len(Nombre_Products_Biological_Process) ,'Genes products présents dans Biological process soit une proportion de :', len(Nombre_Products_Biological_Process)/(Gene_Products_total_GOA))
	print('I y a',len(Nombre_Products_cellular_component) ,'Genes products présents dans Cellular Component soit une proportion de :', len(Nombre_Products_cellular_component)/(Gene_Products_total_GOA))
	print('Il y a',len(Nombre_Products_molecular_function) ,'Genes products présents dans Molecular Function soit une proportion de :', len(Nombre_Products_molecular_function)/(Gene_Products_total_GOA))

	print('Il y a',len(proportion_process_function) ,'Genes products présents dans Biological process et Molecular Function soit une proportion de :', len(proportion_process_function)/(Gene_Products_total_GOA))
	print('I y a',len(proportion_process_component) ,'Genes products présents dans Biological process et Cellular Component soit une proportion de :', len(proportion_process_component)/(Gene_Products_total_GOA))
	print('Il y a',len(proportion_function_component) ,'Genes products présents dans Cellular Component et Molecular Function soit une proportion de :', len(proportion_function_component)/(Gene_Products_total_GOA))

		
	outd = goa_file,"_direct.sets"
	outi = goa_file,"_indirect.sets"
	outputnamed = ''.join(outd)
	outputnamei = ''.join(outi)
	outfildistrib = ''.join('distribution.txt')
	date = "1/10/19"
	species = param.species

	with open(outputnamed, 'w') as outfiledirect, open(outputnamei, 'w') as outfileindirect,open(outfildistrib, 'w') as outfiledistrib:
		outfiledirect.write("# format: sets\n")
		outfiledirect.write("# version: 1.2\n")
		outfiledirect.write("# strain: ")
		outfileindirect.write(species)
		outfileindirect.write("\n")
		outfiledirect.write("# date: %s \n" % date )
		outfiledirect.write("# comment: Gene Ontology terms and direct genes \n")

		outfileindirect.write("# format: sets\n")
		outfileindirect.write("# version: 1.2\n")
		outfileindirect.write("# strain: ")
		outfileindirect.write(species)
		outfileindirect.write("\n")
		outfileindirect.write("# date: %s \n" % date )
		outfileindirect.write("# comment: Gene Ontology terms and indirect genes \n")
		distribution=[]
		for Goterm in go['nodes']:
			if go['nodes'][Goterm]['type']=='GOTerm':
				print(Goterm)
				processus,fonction = origine_GO(go_T,Goterm)

				direct_protein = GeneProducts(go_T, Goterm, all=False, evidence_code=None)
				distribution.append(len(direct_protein))
				indirect_protein = GeneProducts(go_T, Goterm, all=True, evidence_code=None)

				gene_product_id =[]
				gene_product_indirect = []

				for identifiers in direct_protein:
					gene_product_id.append(go_T['nodes'][identifiers]['name'])
				if direct_protein:
					outfiledirect.write("%s\t%s%s%s\t%s\n" % (Goterm,processus,': ',fonction,"\t".join(gene_product_id)))
					

				for identifiers in indirect_protein:
					gene_product_indirect.append(go_T['nodes'][identifiers]['name'])
				if indirect_protein:
					outfileindirect.write("%s\t%s%s%s\t%s\n" % (Goterm,processus,': ',fonction,"\t".join(gene_product_indirect)))
		# print(distribution)
		for dis in distribution:			
			outfiledistrib.write("%s\n" % dis)

##### lib tests #####
if __name__ == "__main__":
	print('GeneOntology lib tests')
	projet()


	
	
	
		# outputname = goa_file,".sets"
	# with open(outputname, 'w') as outfile:
	#     outfile.write("# format: sets")
	#     outfile.write("# version: 1.2")
	#     outfile.write("# strain: Gallus gallus")
	#     outfile.write("%s %s\n" % "# date: ",date.today())
	#     outfile.write("# comment: Gene Ontology terms")
	#     for GOtrm in ggo['nodes']:
	#         if GOtrm.startswith('GO'):
	#             print(GOtrm)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	