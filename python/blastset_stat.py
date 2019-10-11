#!/usr/bin/env python
# Copyright 2018 BARRIOT Roland
# This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
import numpy as np
from os.path import isfile
from scipy.stats import binom, hypergeom, chi2_contingency
import matplotlib.pyplot as plt
import seaborn as sn
import pandas as pd

# SCRIPT PARAMETERS
# e.g. ./blastset.py --sets EcolA.biocyc.sets --query 'ALAS ARGS ASNS ASPS CYSS GLTX GLYQ GLYS HISS ILES'
parser = argparse.ArgumentParser(description='Search enriched terms/categories in the provided (gene) set')
parser.add_argument('-q', '--query', required=True, help='Query set.')#
parser.add_argument('-t', '--sets', required=True, help='Target sets (categories).')
parser.add_argument('-a', '--alpha', required=False, type=float, default=0.05, help='Significance threshold.')
param = parser.parse_args()

class ComparedSet(object):
	def __init__(self, id, name = '', common = 0, size = 0, pvalue = 1, elements = [], common_elements = []):
		self.id = id
		self.name = name
		self.common = common
		self.size = size
		self.pvalue = pvalue
		self.elements = elements
		self.common_elements = common_elements
#


# LOAD QUERY il regarde si c'est un nom de fichier et pour chaque ligne il découpe la ligne avec tab en mot et tout ces mots finissent dans la query
text = param.query
query = set() #c'est un type d'objet et il fourni des opérateurs en plus (des unions etc)
if isfile(text):
	with open(text) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			if l!='':
				query |= set(l.split())
else: # parse string
	query |= set(text.split()) #s'il n'y a pas de fichier à ce nom alors l'entrée est deirect les gènes

# LOAD REFERENCE SETS, charge en mémoire le set  cible (un pathway cible vers un dictionnaire de nom:liste d'éléments du pathway
def load_sets(filename):
	sets = {}
	ids = set()
	with open( filename ) as f:
		content = f.read()
		lines = content.split('\n')
		for l in lines:
			words = l.split('\t')
			if len(words) > 2 and not words[0].startswith('#'):
				id = words.pop(0)
				name = words.pop(0)
				words = set(words)
				sets[ id ] = { 'name': name, 'elements': words}
				ids |= words
	return [ sets, len( ids ) ]
(sets, population_size) = load_sets(param.sets)


#Récupération des 4 sets de p-values pour les 4 mesures implémentées

results_binomial = []
query_size = len(query)
g = population_size 
for id in sets:


	elements = sets[ id ][ 'elements' ]

	common_elements = elements.intersection( query )
	c= len(common_elements)
	t= len(elements)
	pval = binom.cdf( query_size - len(common_elements), query_size, 1 - float(len(elements))/g)
	
	r = ComparedSet( id, sets[id]['name'], len(common_elements), len(elements), pval, elements, common_elements)
	results_binomial.append( r )

results_hypergeometric = []
query_size = len(query)
g = population_size 
for id in sets:
	elements = sets[ id ][ 'elements' ]
	common_elements = elements.intersection( query )
	c= len(common_elements)
	t= len(elements)
	pval = hypergeom.sf(len(common_elements)-1, g, len(elements), query_size)
	r = ComparedSet( id, sets[id]['name'], len(common_elements), len(elements), pval, elements, common_elements)
	results_hypergeometric.append( r )

results_coverage = []
query_size = len(query)
g = population_size 
for id in sets:
	elements = sets[ id ][ 'elements' ]
	common_elements = elements.intersection( query )
	c= len(common_elements)
	t= len(elements)
	pval = (1-(c/query_size)*(c/t))
	r = ComparedSet( id, sets[id]['name'], len(common_elements), len(elements), pval, elements, common_elements)
	results_coverage.append( r )

results_chi2 = []
query_size = len(query)
g = population_size 
for id in sets:
	elements = sets[ id ][ 'elements' ]
	common_elements = elements.intersection( query )
	c= len(common_elements)
	t= len(elements)
	obs = np.array([[c, query_size-c,query_size], [t-c, g-query_size-t+c,g-query_size],[t,g-t,g]])
	x, pval, dof, expctd = chi2_contingency(obs)
	r = ComparedSet( id, sets[id]['name'], len(common_elements), len(elements), pval, elements, common_elements)
	results_chi2.append( r )

#Trie des sets en fonction des p-values
results_binomial.sort(key=lambda an_item: an_item.pvalue) #trie par pvaleur croissante
results_hypergeometric.sort(key=lambda an_item: an_item.pvalue) #trie par pvaleur croissante
results_coverage.sort(key=lambda an_item: an_item.pvalue) #trie par pvaleur croissante
results_chi2.sort(key=lambda an_item: an_item.pvalue) #trie par pvaleur croissante


#Récupération dans des listes tous les résultats significatifs pour les 4 mesures (on récupère tout pour coverage dès qu'il y a au moins un élément en commun)
binom = []
hyper = []
cover = []
chi2 = []


for r in results_coverage:
	# print("%s\t%s\t%s/%s\t%s" % ( r.id, r.pvalue, r.common, r.size, r.name))
	if r.common>0: 
		cover.append(r.name)

for r in results_binomial:
	if r.pvalue < param.alpha: 
		binom.append(r.name)

for r in results_hypergeometric:
	if r.pvalue < param.alpha: 
		hyper.append(r.name)

for r in results_chi2:
	if r.pvalue < param.alpha: 
		chi2.append(r.name)


#Vérification de la présence de résultat dans les sets, si il y a des résultats on fait le calcul, sinon c'est 0
if len(set(binom)) != 0:
	binom_s = len(set(binom))/len(set(binom))*100
else:
	binom_s = 0

if len(set(hyper)) != 0:
	hyper_s = len(set(hyper))/len(set(hyper))*100
else:
	hyper_s = 0

if len(set(chi2)) != 0:
	chi2_s = len(set(chi2))/len(set(chi2))*100
else:
	chi2_s = 0

if len(set(cover)) != 0:
	cover_s = len(set(cover))/len(set(cover))*100
else:
	cover_s = 0


#On fait ici le calcul de recouvrement, si l'un des deux sets est null et bien le recouvrement est null, sinon on fait le calcul
if binom_s == 0 or hyper_s == 0:
	binom_hyper = 0
else:
	binom_hyper = (len(set(binom).intersection(hyper))/(len(set(binom).union(hyper))))*100

if binom_s == 0 or cover_s == 0:
	binom_cover = 0
else:
	binom_cover = (len(set(binom).intersection(cover))/(len(set(binom).union(cover))))*100

if binom_s == 0 or chi2_s == 0:
	binom_chi2 = 0
else:
	binom_chi2 = (len(set(binom).intersection(chi2))/(len(set(binom).union(chi2))))*100

if hyper_s == 0 or cover_s == 0:
	hyper_cover = 0
else:
	hyper_cover = (len(set(hyper).intersection(cover))/(len(set(hyper).union(cover))))*100

if chi2_s == 0 or hyper_s == 0:
	hyper_chi2 = 0
else:
	hyper_chi2 = (len(set(hyper).intersection(chi2))/(len(set(hyper).union(chi2))))*100

if cover_s == 0 or chi2_s == 0:
	cover_chi2 = 0
else:
	cover_chi2 = (len(set(cover).intersection(chi2))/(len(set(cover).union(chi2))))*100



#print du dataFrame dans la console

myDataframe = pd.DataFrame({'methode' : ['binomiale', 'hypergéométrique','coverage','chi2'],
							'binomiale' : [binom_s,binom_hyper,binom_cover,binom_chi2], 
							'hypergéométrique' : [binom_hyper,hyper_s,hyper_cover,hyper_chi2],
							'coverage' : [binom_cover,hyper_cover,cover_s,cover_chi2],
							'chi2' : [binom_chi2,hyper_chi2,cover_chi2,chi2_s]})

print('Table de recouvrement (en pourcentage) des méthodes avec tous les résultats significatifs: ')
print(myDataframe)
