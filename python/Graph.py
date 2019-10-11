#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 18:11:10 2018

@author: ckj2181a
"""

from pprint import pprint
import numpy as np # for Floyd-Warshall matrices
import time 
import copy #for transpose

# TP1 functions
###############

 
def create_graph(directed = True, weighted = False, weight_attribute = None): # TP1
    g = { 'nodes': {}, 'edges': {}, 'nb_edges': 0, 'directed': directed, 'weighted': weighted, 'weight_attribute': weight_attribute }
    return g
 
def add_node(g, n, attributes = None): # TP1
    if n not in g['nodes']: # ensure node does not already exist
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['nodes'][n] = attributes
        g['edges'][n] = {} # init outgoing edges
    return g['nodes'][n] # return node attributes 


def add_edge(g, n1, n2, attributes = None, n1_attributes = None, n2_attributes = None): # TP1
    # create nodes if they do not exist
    if n1 not in g['nodes']: add_node(g, n1, n1_attributes) # ensure n1 exists
    if n2 not in g['nodes']: add_node(g, n2, n2_attributes) # ensure n2 exists
    # add edge(s) only if they do not exist
    if n2 not in g['edges'][n1]:
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['edges'][n1][n2] = attributes
        if not g['directed']:
            g['edges'][n2][n1] = g['edges'][n1][n2] # share the same attributes as n1->n2
        g['nb_edges'] += 1
    return g['edges'][n1][n2] # return edge attributes


def add_edgeT(g, n1, n2, attributes = None, n1_attributes = None, n2_attributes = None): # TP1
    # create nodes if they do not exist
    temp = n1
    n1 = n2
    n2 = temp
    if n1 not in g['nodes']: add_node(g, n1, n1_attributes) # ensure n1 exists
    if n2 not in g['nodes']: add_node(g, n2, n2_attributes) # ensure n2 exists
    # add edge(s) only if they do not exist
    if n2 not in g['edges'][n1]:
        if attributes is None: # create empty attributes if not provided
            attributes = {}
        g['edges'][n1][n2] = attributes
        if not g['directed']:
            g['edges'][n2][n1] = g['edges'][n1][n2] # share the same attributes as n1->n2
        g['nb_edges'] += 1
    return g['edges'][n1][n2] # return edge attributes

def load_SIF(filename, directed=True): # TP1
    # line syntax: nodeD <relationship type> nodeE nodeF nodeB
    g = create_graph(directed) # new empty graph
    with open(filename) as f: # OPEN FILE
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip() # read next line and remove ending whitespaces
        while row:
            vals = row.split('\t') # split line on tab
            for i in range(2, len(vals)):
                att = { 'type': vals[1] } # set edge type
                add_edge(g, vals[0], vals[i], att)
            row = f.readline().rstrip() # read next line
    return g # return graph



def BFS(G, s):
    BFS = {'colour':{},'distance':{}, 'predecesseur':{}} #3 dictionnaires: couleur, distance, predecesseur
    for i in G['nodes']:
        BFS['colour'][i] = 'WHITE'
        BFS['distance'][i] = float('inf')
        BFS['predecesseur'][i] = None
    BFS['colour'][s] = 'GRAY'
    BFS['distance'][s] = 0
    Q = [s]
    while (len(Q) > 0):
        u = Q.pop(0)
        for v in G['edges'][u]:
            if BFS['colour'][v] == 'WHITE':
                BFS['colour'][v] = 'GRAY'
                BFS['distance'][v] = BFS['distance'][u] + 1
                BFS['predecesseur'][v] = u
                Q.append(v)
        BFS['colour'][u] = 'BLACK'
    return(BFS)
    
    
    
def DFS(G):
    DFS = {'colour':{},'debut':{}, 'fin':{}, 'predecesseur':{}, 'class':{}, 'time':0} #4 dictionnaires: couleur, distance, predecesseur
    G['dfs'] = DFS
    for u in G['nodes']:
        DFS['colour'][u] = 'WHITE'
        DFS['predecesseur'][u] = None
    for u in G['nodes']:
        if DFS['colour'][u] == 'WHITE':
            DFS_VISIT(G, u)
    del G['dfs']
    return(DFS)
    
    
    
def DFS_VISIT(G,u):
    DFS = G['dfs']
    DFS['colour'][u] = 'GRAY'
    DFS['time'] = DFS['time'] + 1
    DFS['debut'][u] = DFS['time']
    for v in G['edges'][u]:
        if DFS['colour'][v] == 'WHITE':
            DFS['predecesseur'][v] = u
            DFS_VISIT(G, v)
            DFS['class'][(u,v)] = 'TREE EDGE'
        elif DFS['colour'][v] == 'GRAY':
            DFS['class'][(u,v)] = 'BACK EDGE'
        elif DFS['debut'][u] > DFS['debut'][v]:
            DFS['class'][(u,v)] = 'CROSS EDGE'
        else:
            DFS['class'][(u,v)] = 'FORWARD EDGE'
    DFS['colour'][u] = 'BLACK'
    DFS['time'] = DFS['time'] + 1
    DFS['fin'][u] = DFS['time']
    
def DFS_node(G, node):
    DFS = {'colour':{},'debut':{}, 'fin':{}, 'predecesseur':{}, 'class':{}, 'time':0} #4 dictionnaires: couleur, distance, predecesseur
    G['dfs'] = DFS
    for u in G['nodes']:
        DFS['colour'][u] = 'WHITE' 
        DFS['predecesseur'][u] = None
    DFS_VISIT(G, node) # on sort le DFS de la boucle pour le faire une seule fois
    del G['dfs']
    return(DFS)
    
def is_acyclic(G):
    dfs = DFS(G)
    for i in dfs['class']:
        if dfs['class'][i] == 'BACK EDGE': # si detection back edge alors ce n'est pas acyclic
            return(False)
    return(True)
    
    
    
    
def topological_sort(G):
    dfs = DFS(G)
    sorted_dfs = sorted(dfs['fin'], key=(dfs['fin']).get, reverse=True) #on sort les noeuds en fonction de le fin (décroissant)
    return(sorted_dfs)
    
    
    
def load_TAB(filename, directed=True, weighted=False, weight_attribute=None): 
    """
    parse a TAB file (as cytoscape format) and returns a graph.
 
    line syntax: id1    id2    att1    att2    att3    ...
    """
    g = create_graph(directed, weighted, weight_attribute)
    with open(filename) as f: 
        # GET COLUMNS NAMES
        tmp = f.readline().rstrip()
        attNames= tmp.split('\t')
        # REMOVES FIRST TWO COLUMNS WHICH CORRESPONDS TO THE LABELS OF THE CONNECTED VERTICES
        attNames.pop(0)
        attNames.pop(0)
        # PROCESS THE REMAINING LINES
        row = f.readline().rstrip()
        while row:
            vals = row.split('\t')
            v1 = vals.pop(0)
            v2 = vals.pop(0)
            att = {}
            for i in range(len(attNames)):
                att[ attNames[i] ] = vals[i]
            add_edge(g, v1, v2, att)
            row = f.readline().rstrip() # NEXT LINE
    return g


def INTIALIZE_SINGLE_SOURCE(G, s):
    BF = {'distance':{}, 'predecessor': {}}
    for v in G['nodes']:
        BF['distance'][v] = float('inf')
        BF['predecessor'][v] = None
    BF['distance'][s] = 0
    return(BF)
    

def RELAX(u, v, BF, G):
    w = G['weight_attribute']
    if (BF['distance'][v] > (BF['distance'][u] + float(G['edges'][u][v][w]))):
        BF['distance'][v] = BF['distance'][u] + float(G['edges'][u][v][w])
        BF['predecessor'][v] = u

def BELLMAN_FORD(G, s):
    BF = INTIALIZE_SINGLE_SOURCE(G, s)
    for i in range(1,len(G['nodes'])-1):
        for u in G['edges']:
            for v in G['edges'][u]:
                RELAX(u,v, BF, G)
    return(BF)
    
def adjacency_matrix(G):
    matrix = np.full( (len(G['nodes']),len(G['nodes'])), np.inf )
    iden = sorted(G['nodes'])
    
    for u in G['edges']:
        matrix[iden.index(u),iden.index(u)] = 0
        for v in G['edges'][u]:
            matrix[iden.index(u),iden.index(v)] = G['edges'][u][v]['weight']
    return(matrix,iden)
    
    
    
def FloydWarshall(G):
    W, iden = adjacency_matrix(G)
    D = W
    N = np.full((len(D),len(D)), -1,dtype = np.int)
    for x in range(len(W)):
        for y in range(len(W)):
            D[x][y] = W[x][y]
            N[x][y] = y
    for k in range(len(W)):
        for i in range(len(W)):
            for j in range(len(W)):
                if D[i][k] + D[k][j] < D[i][j]:
                    D[i][j] = D[i][k] + D[k][j]
                    N[i][j] = N[i][k]
    return(D,N)

def FloydWarshallPath(G, i, j):
    W, iden = adjacency_matrix(G)
    D, N = FloydWarshall(G)
    if D[iden.index(i)][iden.index(j)] == float('inf') :
        raise Exception('il ny a pas de chemin entre i et j')
    chemin = [i]
    k = N[iden.index(i)][iden.index(j)]
    while k != iden.index(j) :
        chemin.append(iden[k])
        k= N[k][iden.index(j)]
    weight= D[iden.index(i)][iden.index(j)]
    chemin.append(j)
    return(chemin,weight)
    
    
def diametreGraph(G):
    adj, iden = adjacency_matrix(G)
    D, N = FloydWarshall(G)
    chemin = []
    for i in range(len(D)):
        for j in range(len(D)):
            if D[i,j] == float('inf'):
                D[i,j] = float('-inf')
    diametre = np.amax(D)
    coordonnées = list(np.where(D == diametre))
    coordonnées_depart = iden[coordonnées[0][0]]
    coordonnées_fin = iden[coordonnées[1][0]]
    FW_path, length_path = FloydWarshallPath(G, coordonnées_depart, coordonnées_fin)
    chemin.extend(FW_path)
    return(chemin,length_path)

# projet:
    
def source(G):
    source = []
    pas_source=[]
    for i in G['edges']:
        for u in G['edges'][i]:
            pas_source.append(u) #si le sommet est dans les successeurs ce n'est pas une source
    for i in G['nodes'] :
        if i not in pas_source:
            source.append(i) # les autres sommets sont une source
    return(source)
    
def sinks(G):
    puit=[]
    for u in G['nodes']:
        if G['edges'][u] == {}: #si pas de successeur alors c'est un puit
            puit.append(u)
    return(puit)
    
def transpose_g(G):
    temp = create_graph()
    gT = copy.deepcopy(G) # créer une vraie copie du graph pour ne pas modifier le graph iunitial
    gT['edges'] = {}
    for u in G['edges']:
        for v in G['edges'][u]:
            add_edgeT(temp, u,v) #utilisation de la fonction edgeT pour mettre des edges inversés
    gT['edges'] = temp['edges']
    return(gT)
    
    
## Dijkstra
   
def Dijkstra(G,w,s):
    INTIALIZE_SINGLE_SOURCE(G, s)
    Q = list(G['nodes'])
    while Q:
        u = min()


    

def RELAX(u, v, BF, G):
    w = G['weight_attribute']
    if (BF['distance'][v] > (BF['distance'][u] + float(G['edges'][u][v][w]))):
        BF['distance'][v] = BF['distance'][u] + float(G['edges'][u][v][w])
        BF['predecessor'][v] = u

def BELLMAN_FORD(G, s):
    BF = INTIALIZE_SINGLE_SOURCE(G, s)
    for i in range(1,len(G['nodes'])-1):
        for u in G['edges']:
            for v in G['edges'][u]:
                RELAX(u,v, BF, G)

#######################     fonction utilisées pour le max_depth     ##################
def DFS_VISIT_path(G,u):
    G2 =copy.deepcopy(G)
    DFS = G['dfs']
    DFS['colour'][u] = 'GRAY'
    DFS['time'] = DFS['time'] + 1
    DFS['timeD'] = DFS['timeD'] + 1
    DFS['debut'][u] = DFS['time']
    DFS['distance'][u] = DFS['timeD']
    for v in G['edges'][u]:
        if DFS['colour'][v] == 'WHITE':
            DFS['predecesseur'][v] = u
            DFS_VISIT_path(G, v)
            DFS['class'][(u,v)] = 'TREE EDGE'
        elif DFS['colour'][v] == 'GRAY':
            DFS['class'][(u,v)] = 'BACK EDGE'
        elif DFS['debut'][u] > DFS['debut'][v]:
            i=0
            G2['dfs'] = DFS_node(G2,v)
            for key in DFS['colour'].keys():
                if G2['dfs']['colour'][key] == 'BLACK': 
                    i = i+1
            if i > 1:
                for key in G2['dfs']['colour'].keys():
                    if G2['dfs']['colour'][key] == 'BLACK':
                        DFS['colour'][key] = 'WHITE'
        else:
            DFS['class'][(u,v)] = 'FORWARD EDGE'
    DFS['colour'][u] = 'BLACK'
    DFS['time'] = DFS['time'] + 1 
    DFS['timeD'] = DFS['timeD'] - 1
    DFS['fin'][u] = DFS['time']
    DFS['intervalle'][u] = DFS['fin'][u] - DFS['debut'][u]
    
def DFS_node_path(G, node):
    DFS = {'colour':{},'debut':{}, 'fin':{}, 'predecesseur':{}, 'class':{}, 'distance':{}, 'intervalle':{}, 'time':0,'timeD':0} #4 dictionnaires: couleur, distance, predecesseur
    G['dfs'] = DFS
    for u in G['nodes']:
        DFS['colour'][u] = 'WHITE'
        DFS['predecesseur'][u] = None
    DFS_VISIT_path(G, node)
    del G['dfs']
    return(DFS)
    

def DFS_VISIT_p(G,u):
    DFS = G['dfs']
    DFS['colour'][u] = 'GRAY'
    DFS['time'] = DFS['time'] + 1
    DFS['timeD'] = DFS['timeD'] + 1
    DFS['debut'][u] = DFS['time']
    DFS['distance'][u] = DFS['timeD']
    for v in G['edges'][u]:
        if DFS['colour'][v] == 'WHITE':
            DFS['predecesseur'][v] = u
            DFS_VISIT_p(G, v)
            DFS['class'][(u,v)] = 'TREE EDGE'
        elif DFS['colour'][v] == 'GRAY':
            DFS['class'][(u,v)] = 'BACK EDGE'
        elif DFS['debut'][u] > DFS['debut'][v]:
            
            DFS['class'][(u,v)] = 'CROSS EDGE'
        else:
            DFS['class'][(u,v)] = 'FORWARD EDGE'
    DFS['colour'][u] = 'BLACK'
    DFS['time'] = DFS['time'] + 1 
    DFS['timeD'] = DFS['timeD'] - 1
    DFS['fin'][u] = DFS['time']
    DFS['intervalle'][u] = DFS['fin'][u] - DFS['debut'][u]
    
def DFS_node_p(G, node):
    DFS = {'colour':{},'debut':{}, 'fin':{}, 'predecesseur':{}, 'class':{}, 'distance':{}, 'intervalle':{}, 'time':0,'timeD':0} #4 dictionnaires: couleur, distance, predecesseur
    G['dfs'] = DFS
    for u in G['nodes']:
        DFS['colour'][u] = 'WHITE'
        DFS['predecesseur'][u] = None
    DFS_VISIT_p(G, node)
    del G['dfs']
    return(DFS)

def BFS2(G, s):
    BFS = {'colour':{},'distance':{}, 'predecesseur':{}} #3 dictionnaires: couleur, distance, predecesseur
    for i in G['nodes']:
        BFS['colour'][i] = 'WHITE'
        BFS['distance'][i] = 0
        BFS['predecesseur'][i] = None
    BFS['colour'][s] = 'GRAY'
    BFS['distance'][s] = 1
    Q = [s]
    while (len(Q) > 0):
        u = Q.pop(0)
        for v in G['edges'][u]:
            if BFS['colour'][v] == 'WHITE':
                BFS['colour'][v] = 'GRAY'
                BFS['predecesseur'][v] = u
            Q.append(v)
            if BFS['distance'][u] > BFS['distance'][v]:
                BFS['distance'][v] = BFS['distance'][u] +1
                BFS['predecesseur'][v] = u
        BFS['colour'][u] = 'BLACK'
    return(BFS)
    
    #☺ max(list(BFS2['distance'].values())) != max(list(BFS3['distance'].values()))

#for i in h['edges'].values() :
#    for j in i.values():
#        print(j)
#############
# TP1 Tests

def TP1():
    
    print('----Test fonction create_graph ----')
    
    g = create_graph()
    pprint(g)
    
    
    print('----Test fonction Add node----')
       
    add_node(g,'Z')
    pprint(g)
    
    
    print('----Test fonction add_edge Z -> B----')
   
    add_edge(g, 'Z', 'B')
    pprint(g)
    
    
    print('----Test fonction Load_SIF à partir de M1BBS_Graphe_dressing.sif et Directed.graph.with.cycle.slide29.sif----')
       
    graph_dress = load_SIF('M1BBS_Graphe_dressing.sif')
    pprint(graph_dress)
    
    di_graph_with_cycle = load_SIF('Directed.graph.with.cycle.slide29.sif')
    pprint(di_graph_with_cycle)
    
    
    print('---- Test fonction Load_TAB à partir de M1BBS_Graphe_Bellman-Ford.tab et  M1BBS_Graphe_Floyd-Warshall.tab----')
       
    B_FORD_graph = load_TAB('M1BBS_Graphe_Bellman-Ford.tab', weighted = True, weight_attribute = 'weight')
    pprint(B_FORD_graph)
    
    FWarshall_graph = load_TAB('M1BBS_Graphe_Floyd-Warshall.tab', weighted = True, weight_attribute = 'weight')
    pprint(FWarshall_graph)
    
    
    print('----Test fonction BFS à partir du graph orienté dressing_graph à partir de chemise----')
    
    start_time = time.time()
    bfs_dress = BFS(graph_dress, 'chemise')
    
    print('Temps d exécution de BFS:',time.time()-start_time,'s')
    print('Le BFS à partir de chemise est:')
    pprint(bfs_dress)
    
    print('----Test fonction DFS à partir du graph orienté dressing_graph----')
    
    start_time = time.time()
    dfs_dress = DFS(graph_dress)
    print('Temps d exécution de DFS:',time.time()-start_time,'s')
    print('Le DFS est:')
    pprint(dfs_dress)
    
    print('----Test fonction is_acyclic à partir de directed_graph_with_cycle_29----')
    
    start_time = time.time()
    acyclic = is_acyclic(di_graph_with_cycle)
    print('Temps d exécution de is_acyclic:',time.time()-start_time,'s')
    print('Le graph est acyclic:', acyclic)

    print('----Test fonction topological_sort à partir de directed_graph_with_cycle_29----')
    
    start_time = time.time()
    topo_sort = topological_sort(di_graph_with_cycle)
    print('Temps d exécution de topological_sort:',time.time()-start_time,'s')
    print('Sommets du graph après un tri topologique:', topo_sort)
    
#############
# TP2 Tests
def TP2() :
    
    print('----importer les graphs utilisés----')

    graph_dress = load_SIF('M1BBS_Graphe_dressing.sif')
    pprint(graph_dress)
    B_FORD_graph = load_TAB('M1BBS_Graphe_Bellman-Ford.tab', weighted = True, weight_attribute = 'weight')
    pprint(B_FORD_graph)
    
    FWarshall_graph = load_TAB('M1BBS_Graphe_Floyd-Warshall.tab', weighted = True, weight_attribute = 'weight')
    pprint(FWarshall_graph)
    
    print('----Test fonction BELLMAN_FORD à partir de B_FORD_graph----')
    
    start_time = time.time()
    sommet_A_BF = BELLMAN_FORD(B_FORD_graph, 'A')
    print('Temps d exécution de BELLMAN_FORD:',time.time()-start_time,'s')
    print('Le Bellman Ford à partir du sommet A est:', sommet_A_BF)
    
    print('----test fonction FloydWarshall à partir de FWarshall_graph----')
    
    adj_FW, iden_FW = adjacency_matrix(FWarshall_graph) #1ere étape
    start_time = time.time()
    D_FW, N_FW = FloydWarshall(FWarshall_graph)
    print('Temps d exécution de FloydWarshall:',time.time()-start_time,'s')
    print('Les distances et les successeurs de FloydWarshall sont:')
    print('D:')
    pprint(D_FW)
    print('N:')
    pprint(N_FW)
    
    print('----Test fonction FloydWarshallPath grace aux deux matrices ainsi créées----')
    
    start_time = time.time()
    FW_path_A_to_C, length_path_A_to_C = FloydWarshallPath(FWarshall_graph, 'A', 'C')
    print('Temps d exécution de FloydWarshallPath:',time.time()-start_time,'s')
    print('Le chemin le plus court entre A et C est', FW_path_A_to_C, ' et sa longueur est ', length_path_A_to_C)
    
    print('----Test fonction diametreGraph----')
    
    start_time = time.time()
    chemin,diametre = diametreGraph(FWarshall_graph)
    print('Temps d exécution de diametreGraph:',time.time()-start_time,'s')
    print('Le diametre du graphe est', diametre,'et le chemin correspondant est', chemin)
    
    print('----Test fonction Loading sources----')
    
    start_time = time.time()
    liste_sources = source(graph_dress)
    print('Temps d exécution de sources:',time.time()-start_time,'s')
    print('Les sommets sources sont', liste_sources)
    
    print('----Test fonction sinks----')
    
    start_time = time.time()
    liste_sinks = sinks(graph_dress)
    print('Temps d exécution de sinks:',time.time()-start_time,'s')
    print('Les sommets puits sont', liste_sinks)



#############
# projet Tests
def projet() :
    #trois fonctions on été créées pour le max_depth
    
   print('----Test fonction projet----')
   graph_dress = load_SIF('M1BBS_Graphe_dressing.sif')
   add_edge(graph_dress,'pantalon','chemise')
   add_edge(graph_dress,'chemise','chaussures')
   add_edge(graph_dress,'chaussures','cravate')
   add_edge(graph_dress,'ceinture','chaussures')
   add_edge(graph_dress,'sous-vetements','chemise')
   pprint(graph_dress)
   
   
   print('----Test fonction max_depth 2 donc sans prise en compte des cross edges ----')
   start_time = time.time()
   dfs2 = DFS_node_p(graph_dress, 'sous-vetements') 
   print('Temps d exécution de sinks:',time.time()-start_time,'s')
   print('test sur le nouveau graph dress (profondeur):',max(dfs2['distance'].values()))
   
   
   print('----Test fonction max_depth 1 avec prise en compte des cross edges----')
   start_time = time.time()
   dfs = DFS_node_path(graph_dress, 'sous-vetements') 
   print('Temps d exécution de sinks:',time.time()-start_time,'s')
   print('test sur le nouveau graph dress (profondeur):',max(dfs['distance'].values()))
   
   
   print('----Test fonction max_depth final en utilisant le BFS ----')
   start_time = time.time()
   BFS = BFS2(graph_dress, 'sous-vetements') 
   print('Temps d exécution de sinks:',time.time()-start_time,'s')
   print('test sur le nouveau graph dress (profondeur):',max(BFS['distance'].values()))

#############
 
# Perform tests if not imported by another script
if __name__ == "__main__":
    TP1()
    TP2()
    projet()











    
#G2 =copy.deepcopy(graph_dress)

























