#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from taxonomy_library import *
from heatmap_library import *
import sys, getopt

# Get full command-line arguments
full_cmd_arguments = sys.argv

# Keep all but the first
argument_list = full_cmd_arguments[1:]

short_options = ""
long_options = ["graph=","pos","neg","seasons=","min_cor=","graphml","out=","test","gml","annotate","clustering=","merge","sw=","cross_database"] 

try :
	arguments, values = getopt.getopt(argument_list, short_options,long_options)
except getopt.error as err:
	# Output error, and return with an error code
    print (str(err))
    sys.exit(2)


# by default
graph = None
sign = 1 #positive by default
min_cor = 0
fseasons = None
graphml = 0
gml = 0
out = "out"
test = 0
clustering = None
annotate = False
merge = False
internodes = None #keep this
interedges = None 
fsw = None
cross_database = 0

 # Evaluate given options
for current_argument, current_value in arguments:
	if current_argument in ("--graph"):
		graph = current_value
	if current_argument in ("--min_cor"):
		min_cor = float(current_value)
	if current_argument in ("--seasons"):
		fseasons = current_value
		print("* Seasons mode *")
	if current_argument in ("--pos"):
		sign = 1
	if current_argument in ("--neg"):
		sign = 0
	if current_argument in ("--graphml"):
		graphml = 1
	if current_argument in ("--gml"):
		gml = 1
	if current_argument in ("--out"):
		out = current_value
	if current_argument in ("--sw"):
		fsw = current_value
		print("* Sliding windows mode *")
	if current_argument in ("--test"):
		test = 1
	if current_argument in ("--annotate"):
		annotate = True
	if current_argument in ("--merge"):
		merge = True
	if current_argument in ("--cross_database"):
		cross_database = True
	if current_argument in ("--clustering"):
		clustering = current_value
		if clustering not in ["greedy","eigen","betweenness",None]:
			print("Invalid clustering method. Please select 'greedy' or 'eigen' or 'betweenness' or not at all.")
			exit(1)


## abundance_seasons parameters
# it is normalized tss abundances
# df = pd.read_csv("useful_tables/df_norm_ab_seasons.tsv",sep="\t",index_col=0)
# it is normalized tss abundances + normalization of each ASvs abundance
norm_df = pd.read_csv("useful_tables/norm_abundances.tsv",sep="\t",index_col=0)
norm_df.index = pd.DatetimeIndex(norm_df.index)
	
# ab_columns = df.columns[0:24833]
# ab_seasons= pd.DataFrame(df.groupby(["Season"])[ab_columns].sum()).T[["Spring","Summer","Autumn","Winter"]]

if fseasons != None:
	#seasons crossing 
	print("Seasons path :",fseasons)
	with open(fseasons) as f:
		files_seasons = f.read().splitlines()
	#open seasons
	dseasons = dict()
	for i in range(0,len(files_seasons)):
		try:
			dseasons[str(i+1)] = nx.read_gml(files_seasons[i])
		except nx.exception.NetworkXError:
			dseasons[str(i+1)] = nx.read_graphml(files_seasons[i])
		seasons = list(dseasons.values())

	globalSeasons = nx.compose(seasons[0],seasons[1])
	for i in range(2,len(seasons)):
		globalSeasons = nx.compose(globalSeasons,seasons[i])

	dexclunodes = dict()
	dexcluedges = dict()
	saveseasons = seasons.copy()

	#exclusives
	for i in range(0,len(seasons)):
	    g = seasons.pop(i)
	    dexclunodes[i+1] = exclusive_nodes(g,seasons)
	    dexcluedges[i+1] = exclusive_edges(g,seasons)
	    seasons = saveseasons.copy()

	    #########

		##########

	# for k,v in dexcluedges.items():
	#     print(k,":",len(v))

	# if graph == None and annotate:
	# 	globalSeasons = get_nodes_taxonomy(globalSeasons)
	# 	globalSeasons = get_temp_niches(globalSeasons) 
	# 	globalSeasons = get_metavars(globalSeasons)
	# 	globalSeasons = make_cross_databases(globalSeasons)



#############################################################################################################
G = nx.Graph()
if graph != None and fsw == None and fseasons == None:
	try:
		G = nx.read_gml(graph)
	except nx.exception.NetworkXError:
		G = nx.read_graphml(graph)
	if fseasons:
		G = seasons_crossing(G, dseasons, verbose=False)
elif fseasons != None:
	G = globalSeasons
	G = crossing_seasons(G, seasons)
if fsw!= None:
	Gsw = open_sw(fsw)
	if graph != None:
		try:
			G = nx.read_gml(graph)
		except nx.exception.NetworkXError:
			G = nx.read_graphml(graph)
	Frozen = nx.subgraph(Gsw,G.nodes())
	G = Frozen.copy()


G = filtering_edges(G,min_cor,verbose=True,color=True,sign=sign)
if merge:
	G = contract_network(G)
# G = filter_identitity_taxonomy(G, 90)
G = get_closeness_network(G)
G = get_degrees_network(G)
G = get_betweenness_network(G)
G = get_degrees_centrality_network(G)
G = get_subgraph_centrality_network(G)

if annotate or fsw != None:
	G = get_nodes_taxonomy(G)
	G = get_temp_niches(G) 
	G = get_metavars(G)
	
	if cross_database:
		G = get_graph_lineages(G)
		G = make_cross_databases(G)

    
#estrada keystone index
attributes = {n:[d["degreescentrality"],d["betweennesscentrality"],d["closenesscentrality"],d["subgraphcentrality"],d["degrees"]] for (n,d) in G.nodes(data=True)}

attributes = dict(sorted( attributes.items(), key= lambda k: (k[1][0],k[1][1],k[1][2],k[1][3]), reverse=False ))

rank = {key: rank for rank, key in enumerate(sorted(attributes, key=attributes.get, reverse=True), 1)}
nx.set_node_attributes(G, rank, "keystoneindexrank")

#find corresponding subgraph
lsubgraphs = connected_component_subgraphs(G)
cpt = len(lsubgraphs)
for sg in lsubgraphs[::-1]: #in reverse order
	for n in sg.nodes():
		G.nodes[n]["subgraph"] = cpt
	cpt+=1
	

if annotate:
	for n1,n2,d in G.edges(data=True):
		if cross_database:
		    d["pidaInteractionType"] = str(d["pidaInteractionType"])
		    d["pidaInteractionLink"] = str(d["pidaInteractionLink"])
		    d["pidaInteractionRef"] = str(d["pidaInteractionRef"])
		    d["globiInteractionType"] = str(d["globiInteractionType"])
		    d["globiInteractionLink"] = str(d["globiInteractionLink"])
		    d["globiInteractionRef"] = str(d["globiInteractionRef"])
		if fsw != None:
			d["weights"] = str(d["weights"])
			d["sw"] = str(d["sw"])
		if fseasons != None:
			d["seasons"] = str(d["seasons"])
	for n,d in G.nodes(data=True):
		G.nodes[n]["key"]=n
		if cross_database:
			d["taxid"] = str(d["taxid"])
		if fsw != None:
			d["sw"] = str(d["sw"])
		if fseasons != None:
			d["seasons"] = str(d["seasons"])
		if "contraction" in d.keys(): #if merging
			del d['contraction']

	#if graphml:	
	nx.write_graphml(G,out+".graphml")
	print(out,".graphml was saved !",sep="")
	#	exit()
	#if gml:
	#nx.write_gml(G,out+".gml")
	#print("\n",out,".gml was saved !\n",sep="")
	exit()



merged = []
# #find patterns
# print("\nDetection of patterns ...")
# print("root is :",list(rank.keys())[0])
# G, cli, tri, cir = find_circles_and_complete_patterns(G,list(rank.keys())[0],max_circle=4)
# cli.extend(tri)
# cli.extend(cir)
# merged, notmerged = merge_patterns(cli)
# merged.extend(notmerged)
# merged.sort(key=len) # sort by len of patterns
# merged = [ "-".join(p) for p in merged ]

make_graph_topology(G)
print("\nProportions of the interactions : ")
cptP, cptE = print_proportions_interactions(G)

# print("\nZoom on Prokaryotes : ")
# myProk = print_proportions_prokaryotes(G,cptP)
# # study_metavar_by_taxo_group(G,myProk)
# print("\nZoom on Eukaryotes : ")
# myEuk = print_proportions_eukaryotes(G,cptE)

if test:
	simple_draw_network(lsubgraphs[5], myax= None, clustering =clustering)
	plt.show()
else:
	while True:
		
		while True:
			print("\n * Write 0 to leave *")
			print(" * Write -1 to see all the nodes *")
			print(" * Write -2 to zoom in an interaction *")
			print(" * Write -3 to see all the connected components *")
			print(" * Write -4 to see the detected patterns *")
			print(" * Write -5 to make a NODE taxonomy research *")
			print(" * Write -6 to make an EDGE taxonomy research *")
			print(" * Write -7 to have a deeper view on the databases cross edges *")
			if fseasons != None:
				print(" * Write -8 to have a view on the seasons *")
			try:
				n = input("\n\t***\n Input : ")
				#whole network
				if n in G.nodes :

					print("Found in subgraph",G.nodes[n]["subgraph"])

					print_a_node(G, n, rank, taxonomy = False)
					print_node_interactions(G, n, taxonomy=False)
					
					# #the whole subgraph
					# if lsubgraphs[cpt-1].order() < 500:
					# 	simple_draw_network(lsubgraphs[cpt-1], nodes = n, myax= None, clustering =clustering)

					if fsw != None:
						make_heatmap_node(G, fsw, n)
					if fseasons != None:
						make_heatmap_node(G, fseasons, n)


				elif n in merged :

					nodes = n.split("-")

					#find corresponding subgraph
					cpt = len(lsubgraphs)
					for sg in lsubgraphs[::-1]: #in reverse order
						if nodes[0] in sg.nodes():
							print(" ",n,"found in subgraph",cpt,":",sg.order(),"nodes,",sg.size(),"edges")
							break
						else:
							cpt -= 1

					# for node in nodes:
					# 	print_taxonomy(G,node)

					sg = G.subgraph(nodes)

					#show weights
					print("\n Weights :")
					for n1,n2,d in sg.edges(data=True):
						print("",n1,"-",n2,":",round(d["weight"],2))
					
					fig = plt.figure(figsize = (20,15)) #frameon = False for no background
					ax1 = fig.add_subplot(221)
					print_metaniches(sorted(list(sg.nodes)),ax=ax1)


					ax3 = fig.add_subplot(223)
					norm_df[sorted(list(sg.nodes))].plot(kind="line",ax=ax3,xticks = list(norm_df.index)[::2],rot=90,fontsize=6)

					ax = fig.add_subplot(122)
					# ax inside the function
					simple_draw_network(sg,myax=ax, clustering =clustering)
					plt.ion()
					fig.subplots_adjust(bottom=0.08)
					plt.show()
					plt.draw()
					plt.pause(0.001)
					input("\n  Press [enter] to continue.")

					#the whole subgraph
					if lsubgraphs[cpt-1].order() < 500:
						simple_draw_network(lsubgraphs[cpt-1], nodes, clustering =clustering)


				else:
					
					n = int(n)
					
					if n == 0:
						print("\tYou exit")
						exit(0)

					elif n == -1:
						print_list_nodes(G,taxonomy=True)

					elif n == -2:
						n1 = input("\n Give the first node label : ")
						print(n1)
						if n1 in G.nodes():
							print_node_interactions(G,n1, taxonomy=False, clustering = None, figure=True, allInteractions = False)
						else:
							(n1,"not found in the graph")


					elif n == -3:
						cpt = len(lsubgraphs)
						for sg in lsubgraphs[::-1]: #in reverse order
							print(cpt,":",sg.order(),"nodes,",sg.size(),"edges")
							cpt -= 1

					elif n == -4:
						print("\n  List of patterns :")
						for p in merged:
							print(" -",p,":",p.count("-")+1,"nodes")

					elif n == -5:
						cpt = 0

						search = input("\n Search : ")
						list_nodes = list()

						for n,d in G.nodes(data=True):
							tax = d["taxonomy"]
							# print(n,tax)
							if re.findall(search,tax):
								# print(" -",n,":",tax,"- degree(s) :",d["degrees"],"- closeness value :",d["closeness"])
								list_nodes.append(n)
								cpt += 1
						sG = G.subgraph(list_nodes)
						print_list_nodes(G,subG = sG,taxonomy=True)
						print("\n  => Number of occurrences :",cpt)

					elif n == -6:
						cpt = 0

						search1 = input("\n Search 1 : ")
						search2 = input(" Search 2 : ")

						nodes = G.nodes()

						matches = list()

						for n1 in nodes:
							d1 = G.nodes[n1]
							tax1 = d1["taxonomy"]
							# print(n,tax)
							if re.findall(search1,tax1):
							# 	print("MATCH")
								connected = [ e[1] for e in G.edges(n1)]
								for n2 in connected:
									tax2 = G.nodes[n2]["taxonomy"]
									if re.findall(search2,tax2):
										d2 = G.nodes[n2]
										print("\n***")
										print("weigth :",G.edges[n1,n2]["weight"])
										print_a_node(G,n1,rank,taxonomy=True)
										print_a_node(G,n2,rank,taxonomy=True)
										cpt += 1
										matches.append([n1,n2])
						print("\n  => Number of occurrences :",cpt)

						answer = int(input("Want to see all interactions (1:yes/0:no)"))
						if answer:
							for m in matches:
								print_node_interactions(G,m[0], n2=m[1], taxonomy=False, clustering = None, figure=True, allInteractions = False)

					elif n == -7 :
						deeper_view_cross_databases(G)

					elif n == -8 :
						
						for i in range(1,len(seasons)+1):
							print("-",i)
						s = int(input("Which season (0 to erase the crossing season): "))

						if s != 0:
							internodes = intersection_nodes(season, G)
							interedges = intersection_edges(season,G)
							crossSeason = nx.edge_subgraph(season, interedges)
							excluedges = [ e for e in dexcluedges[s] if e in crossSeason.edges() ]

							for e in excluedges:
								crossSeason.edges[e]["color"] = "#4888FF"
							simple_draw_network(crossSeason,myax=None,clustering=clustering)
							plt.show()
						
						else:
							internodes = None
							interedges = None

					elif n > 0:

						sg = lsubgraphs[n-1]

						if sign=="both":
							notfrozenpos = sg.copy()
							notfrozenneg = sg.copy()

						if sg.order() <= 15:

							fig = plt.figure(figsize = (20,15)) #frameon = False for no background

							ax1 = fig.add_subplot(221)
							print_metaniches(sorted(list(sg.nodes())),ax = ax1)

							ax3 = fig.add_subplot(223)
							norm_df[sorted(list(sg.nodes))].plot(kind="line",ax=ax3,xticks = list(norm_df.index)[::2],rot=90,fontsize=6)

							ax = fig.add_subplot(122)
							# ax inside the function
							simple_draw_network(sg,myax=ax, clustering =clustering,nodes=internodes, edges=interedges)

							if sign=="both":
								sgpos = filtering_edges(notfrozenpos,0,sign=1)
								sgneg = filtering_edges(notfrozenneg,0,sign=0)
								fig2 = plt.figure(figsize = (15,15))
								ax1 = fig2.add_subplot(211)
								simple_draw_network(sgpos,myax=ax1, clustering =clustering)
								ax2 = fig2.add_subplot(212)
								simple_draw_network(sgneg,myax=ax2, clustering =clustering)
								
							plt.ion()
							fig.subplots_adjust(bottom=0.08)
							plt.show()
							plt.draw()
							plt.pause(0.001)
							input("\n  Press [enter] to continue.")

						else:
							simple_draw_network(sg, clustering =clustering, nodes = internodes, edges=interedges)
							plt.ion()
							plt.show()
							plt.draw()
							plt.pause(0.001)
							input("\n  Press [enter] to continue.")
							
			except ValueError as e:
				print("/!\\",e)
				print("Value error has occurred, please retry")

			except IndexError as e:
				print("/!\\",e)
				print("Index error has occurred, please retry")



