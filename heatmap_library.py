#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from taxonomy_library import *
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys 
from collections import defaultdict
import seaborn as sns

def filter_sw(G):
	toremove = list()
	for n1,n2,d in G.edges(data=True):
		if len(d["sw"]) == 1:
			toremove.append((n1,n2))
	G.remove_edges_from(toremove)
	G.remove_nodes_from(list(nx.isolates(G)))
	return(G)

def open_sw(sw):
	with open(sw) as f:
		files = f.read().splitlines()

	graphs=list()
	for f in files:
	 	graphs.append(nx.read_gml(f))

	for g in range(0,len(graphs)):
		nx.set_node_attributes(graphs[g],list(), "sw")
		nx.set_edge_attributes(graphs[g],list(), "sw")

	Gsw = nx.Graph()
	for i in range(0,len(graphs)):
		Gsw = nx.compose(Gsw,graphs[i])

	dswnodes = defaultdict(list)
	dswedges = defaultdict(list)
	dswweights = defaultdict(list)
	for n in Gsw.nodes:
		dswnodes[n] = []
	for n1,n2,d in Gsw.edges(data=True):
		dswedges[(n1,n2)] = []
		dswweights[(n1,n2)] = []

	cpt_g = 1
	for g in graphs:
		for n in g.nodes:
			dswnodes[n].append(cpt_g)
		for n1,n2,d in g.edges(data=True):
			if d["weight"] > 0:
				dswedges[(n1,n2)].append(cpt_g)
				dswweights[(n1,n2)].append(d["weight"])
		cpt_g += 1

	for k,v in dswnodes.items():
		Gsw.nodes[k]["sw"] = v
	for k,v in dswedges.items():
		Gsw.edges[k]["sw"] = v

	for k,v in dswweights.items():
		Gsw.edges[k]["weights"] = v

	Gsw = filter_sw(Gsw) #remove edges in only one sw

	return(Gsw)


def crossing_seasons(Gs,seasons):
	# with open(seasons) as f:
	# 	files = f.read().splitlines()

	# graphs=list()
	# for f in files:
	#  	graphs.append(nx.read_gml(f))

	for g in range(0,len(seasons)):
		nx.set_node_attributes(seasons[g],list(), "seasons")
		nx.set_edge_attributes(seasons[g],list(), "seasons")

	Gs = nx.Graph()
	for i in range(0,len(seasons)):
		Gs = nx.compose(Gs,seasons[i])

	dsnodes = defaultdict(list)
	dsedges = defaultdict(list)
	# dsweights = defaultdict(list)
	for n in Gs.nodes:
		dsnodes[n] = []
	for n1,n2,d in Gs.edges(data=True):
		dsedges[(n1,n2)] = []
		# dsweights[(n1,n2)] = []

	cpt_g = 1
	for g in seasons:
		for n in g.nodes:
			dsnodes[n].append(cpt_g)
		for n1,n2,d in g.edges(data=True):
			# if d["weight"] > 0:
				dsedges[(n1,n2)].append(cpt_g)
				# dsweights[(n1,n2)].append(d["weight"])
		cpt_g += 1

	for k,v in dsnodes.items():
		Gs.nodes[k]["seasons"] = v
	for k,v in dsedges.items():
		Gs.edges[k]["seasons"] = v

	# for k,v in dswweights.items():
	# 	Gsw.edges[k]["weights"] = v

	# Gsw = filter_sw(Gsw) #remove edges in only one sw

	return(Gs)



def startend_dates(str):
	fd = re.findall(r"[0-9]{4}-[0-9]{2}-[0-9]{2}",str)[0]
	ld = re.findall(r"[0-9]{4}-[0-9]{2}-[0-9]{2}",str)[1]
	return(fd,ld)


def make_heatmap_node(G, sw, node):

	with open(sw) as f:
		files = f.read().splitlines()

	ldates = list()
	for i in range(0, len(files)):
	 	fd,ld=startend_dates(files[i])
	 	ldates.append(fd+" - "+ld)

	edges = G.edges(node,data=True)
	neighbors = list(nx.neighbors(G,node))
	df = pd.DataFrame(0, index = neighbors, columns = ldates )

	for e in edges:
		neighbor = e[1]
		dates = [ldates[i-1] for i in e[2]["sw"] ]
		weights = e[2]["weights"]
		comb = list(zip(dates,weights))
		for c in comb:
			df.loc[neighbor, c[0]] = round(c[1],2)

	df = df.loc[df.gt(0).sum(axis=1).sort_values(ascending=False).index]

	fig,ax = plt.subplots(figsize = (15,10)) #frameon = False for no background
	ax = sns.heatmap(df,cmap="cubehelix_r",linewidths = 0.1 ,linecolor="black",vmin = 0, vmax=1)
	ax.set_title('Interactions of '+node)

	plt.ion()
	plt.tight_layout()
	plt.draw()
	plt.pause(0.001)
	input("\n  Press [enter] to continue.")

########################

def make_heatmap_edge(G, files, graphs, node1, node2):

	fig,ax = plt.subplots(figsize = (15,5)) #frameon = False for no background

	list_sw = G.edges[node1,node2]["sw"]
	dates = list()
	for sw in files:
		fd,ld=startend_dates(sw)
		dates.append(fd+" - "+ld)

	df = pd.DataFrame(np.nan, index = [node1+"-"+node2], columns = dates )

	for sw in list_sw:
		df.at[node1+"-"+node2,dates[sw]] = graphs[sw].edges[node1,node2]["weight"]


	list_weights = list(df.loc[node1+"-"+node2])
	mean_weight = np.nansum(list_weights)/np.count_nonzero(~np.isnan(list_weights))
	# print("mean of weight is :",mean_weight)

	if mean_weight>=0:#positive weighted edges:rouge
		ax = sns.heatmap(df,cmap="Reds",linewidths = 0.1 ,linecolor="black",vmin=float(df.T.min())-0.01,vmax=float(df.T.max())+0.01)#,annot=True)
	else: #negative weighted edge:bleu
		ax = sns.heatmap(df,cmap="Blues_r",linewidths = 0.1 ,linecolor="black",vmin=float(df.T.min())-0.01,vmax=float(df.T.max())+0.01)#,annot=True)

	ax.set_title(node1+"-"+node2+' interaction')
	plt.tight_layout()
	plt.ion()
	plt.show()
	plt.draw()
	plt.pause(0.001)
	input("\n  Press [enter] to continue.")

#############################

def make_heatmap_patterns(G, files, patterns , pattern):

	sg = G.subgraph(pattern.split("-"))
	simple_draw_network(sg)

	fig,ax = plt.subplots(figsize = (15,5)) #frameon = False for no background
	
	dates = list()
	for sw in files:
		fd,ld=startend_dates(sw)
		dates.append(fd+" - "+ld)

	df = pd.DataFrame(np.nan, index = [pattern], columns = dates )

	sw = dict(patterns)[pattern]

	for i in sw:
		df.at[pattern,dates[i]] = 1

	ax = sns.heatmap(df,cmap="Reds",linewidths = 0.1 ,linecolor="black")
	ax.set_title(pattern+' pattern')
	plt.tight_layout()
	plt.ion()
	plt.show()
	plt.draw()
	plt.pause(0.001)
	input("\n  Press [enter] to continue.")



