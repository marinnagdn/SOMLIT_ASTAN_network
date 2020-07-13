#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from useful import *
warnings.filterwarnings("ignore", category=Warning)


####################

#data

norm_df = pd.read_csv("useful_tables/norm_abundances.tsv",sep="\t",index_col=0)
norm_df.index = pd.DatetimeIndex(norm_df.index)

clr_df = pd.read_csv("useful_tables/ASV_CLR_abundances.tsv",sep="\t",index_col=0)
clr_df.index = pd.DatetimeIndex(clr_df.index)

dfDatesSeasons = pd.read_csv("/home/marinna/Bureau/scan_networks/useful_tables/dates_and_seasons.tsv",sep="\t",index_col=0)



#####################

def get_font_size(G):
	# for network display #
	
	nbr_nodes = len(G.nodes)
	if nbr_nodes > 1000:
		return(2)
	elif nbr_nodes > 100:
		return(3)
	elif nbr_nodes > 50:
		return(4)
	elif nbr_nodes > 10:
		return(5)
	else :
		return(6)


######################

def get_sizeMaxNodes(G):
	# for network display #
	
	nbr_nodes = len(G.nodes)
	if nbr_nodes > 1000:
		return(50)
	elif nbr_nodes > 100:
		return(100)
	elif nbr_nodes > 50:
		return(200)
	elif nbr_nodes > 10:
		return(300)
	else :
		return(500)




##########

#def histogram(mydict):
#	values = list(mydict.values())
#	# mini = floor(min(values))
#	# maxi = ceil(max(values))
#	mini = 8
#	maxi = 17
#	ar = np.array(values)

#	fig, ax = plt.subplots()
#	N, bins, patches = plt.hist(ar, edgecolor = "white", bins = maxi-mini)

#	colorsRed = list(Color("#FFC1C1").range_to(Color("#B80000"),int(len(bins)/2)))
#	colors = list(Color("#3333A7").range_to(Color("#AEFFF7"),int(len(bins)/2)))
#	colors.extend(colorsRed)
#	colors = [str(c) for c in colors]
#	
#	for i in range(0,maxi-mini):
#		patches[i].set_facecolor(colors[i])

#	ax.set_facecolor('white')

#	plt.show()

#########

def make_graph_topology(G):

		print("\n*** Graph topologic informations ***\n")
		print("Number of nodes : ",G.order())
		print("Number of edges : ",G.size())
		print("Connectivity :",round(G.size()/G.order(),2))
		print("Number of connected components :",nx.number_connected_components(G))
		print("Degree distribution :",nx.degree_histogram(G))
		mean_degree = sum([degree * node for degree, node in enumerate(nx.degree_histogram(G))]) / G.order()
		print("Mean degrees :", round(mean_degree,2))
		weights = nx.get_edge_attributes(G,"weight")
		print("Mean weight :",round(np.mean(list(weights.values())),2))
		is_connected = nx.is_connected(G)
		if is_connected:

			excentricity = nx.eccentricity(G)
			# print("Excentricity :", excentricity)
			print("Radius :",nx.radius(G, excentricity))
			print("Diameter :", nx.diameter(G, excentricity))
			print("Center : ",nx.center(G, excentricity))
			print("Periphery :",nx.periphery(G, excentricity))

		else:
		    print("LCC : ",len(max(nx.connected_components(G), key=len)))

		# More complex infos
		clustering = nx.clustering(G)
		print("Clustering coefficient :",round(np.mean(list(clustering.values())),3))
		print("Degree assortativity :",round(nx.degree_assortativity_coefficient(G),3))


#########

def find_circle_patterns(G,root,max_circle=5):
    patterns = list(nx.cycle_basis(G,root=root))
    # print("len patterns :",len(patterns))
    circle_patterns=list()
    r = list(range(4,max_circle+1))

    for p in patterns:
        sg = G.subgraph(p)
        if len(p) >= sg.size() and sg.size() in r:
            circle_patterns.append(p)

    num_patterns=list()
    for i in range(0,len(circle_patterns)):
    	num_patterns.append( [ int(p.split("_")[0]) for p in circle_patterns[i] ])
    circle_patterns = [x for _,x in sorted(zip(num_patterns,circle_patterns))]
#    print("circle patterns : ",circle_patterns)
    return(G,circle_patterns)


#########

def find_circles_and_complete_patterns(G,root,max_circle=5):
	G, circles = find_circle_patterns(G,root,max_circle)
	cli = list(nx.enumerate_all_cliques(G))

#remove 1 and 2 and 3 nodes cliques
	while cli != [] and len(cli[0]) != 3 :
		cli.pop(0)
 
	tri = list()
	while cli != [] and len(cli[0]) != 4 :
		tri.append(cli[0])
		cli.pop(0)
        
        
    #remove triangles in 4-plets for example
	if cli != 0:
		for c in cli[::-1]:
			for i in range(len(c)-1,2,-1):
				for comb in itertools.combinations(c, i):
					comb = list(comb)
					for comb in cli:  
						cli.pop(cli.index(comb))
#	print("complete patterns (>= 4 nodes) :",cli)
#	print("triangles :",tri)
	return(G, cli, tri, circles)
	
#####################

#def find_complete_patterns(G):
#	cli = list(nx.enumerate_all_cliques(G))
#	cpt_c = 0

#	#remove 1 and 2 nodes cliques
#	while cli != [] and len(cli[0]) != 3 :
#		cli.pop(0)

#	for c in cli[::-1]:
#	    for i in range(len(c)-1,2,-1):
#	        for comb in itertools.combinations(c, i):
#	            comb = list(comb)
#	            if comb in cli:  
#	                cli.pop(cli.index(comb))

#	for p in cli:
#		for n in p:
#			G.nodes[n]["color"] = "#FF4E00"
#			for e in G.edges(n):
#				G.edges[e]['color'] = "#FF4E00"

#	return(G, cli)
    

###########
    
def merge_patterns(pat):
    #if two triangles are side by side 
    patcopy = [p.copy() for p in pat]
    rindex = list()
    for i in range(0,len(pat)-1):
        	for j in range(i+1,len(pat)):
        		if len(set(pat[i]).intersection(set(pat[j]))) > 0:
        			pat[i].extend(pat[j])
        			rindex.append(j)
    rindex=list(set(rindex))
    rpatterns = [pat[r] for r in rindex ]
    for r in rpatterns:
        	pat.remove(r)
    common = [ p for p in pat if p in patcopy ]
    merged = [ list(set(m)) for m in pat if m not in common ]
#    print("merged :",merged)
    
    return(merged,common)

#####################

# def seasons_crossing(G,seasons,verbose=True):

# 	interactions_seasons = list() #juste pour le print du nombre d'interactions trouvées en commun
# 	#find seasons
# 	# dseasonsnodes = defaultdict(list)
# 	dseasonsnodes = dict() #dictionnaire simple avec des chaînes de caractères et pas des lites pour pouvoir enregistrer le fichier au format graphml ou gml
# 	for n in G.nodes():
# 		dseasonsnodes[n]="" #initialisation string vide pour chaque noeud

# 	for n1,n2,d in G.edges(data=True):
# 		d["seasons"] = "" #on initialise à chaîne vide pour concaténer
# 		sn1 = [k for k in seasons.keys() if n1 in seasons[k].nodes() ] #le noeud est il dans une saison ?
# 		# dseasonsnodes[n1].extend(s)
# 		for s in sn1:
# 			if s not in dseasonsnodes[n1]: #pour ne pas avoir plus fois la même saison
# 				dseasonsnodes[n1]+=str(s) #on paste

# 		sn2 = [k for k in seasons.keys() if n2 in seasons[k].nodes() ]
# 		# dseasonsnodes[n2].extend(s)
# 		for s in sn2:
# 			if s not in dseasonsnodes[n2]: #pour ne pas avoir plus fois la même saison
# 				dseasonsnodes[n2]+=str(s) #on paste

# 		sedge = [k for k in seasons.keys() if (n1,n2) in seasons[k].edges() ]
# 		n = [(n1,n2) for k in seasons.keys() if (n1,n2) in seasons[k].edges() ]
# 		interactions_seasons.extend(n)
# 		# d["seasons"] = s
# 		if sedge != []:
# 			for s in sedge:
# 				if s not in d["seasons"]:
# 					d["seasons"]+=str(s)

# 	# for k in dseasonsnodes.keys():
# 	# 	dseasonsnodes[k] = list(set(dseasonsnodes[k])) #plus besoin avec les chaînes de caractères

# 	nx.set_node_attributes(G,dseasonsnodes,"seasons")

# 	interactions_seasons = list(set(interactions_seasons))
# 	# print(interactions_seasons)
# 	print(len(interactions_seasons),"interactions seasons found over",G.size())

# 	if verbose:
	
# 	#print seasons crossing #########
# 		print("\ncross seasons and nodes :")
# 		dseasons = {"":0,"1":0,"2":0,"3":0,"4":0,"12":0,"13":0,"14":0,"23":0,"24":0,"34":0,"123":0,"234":0,"124":0,"134":0,"1234":0}
# 		for n in G.nodes:
# 			dseasons["".join(sorted(G.nodes[n]["seasons"]))]+=1	
# 		for k,v in dseasons.items():
# 			if v != 0:
# 				print(k,":",v)
# 		print("sum :",G.order())

# 		print("\ncross seasons and edges :")
# 		dseasons = {"":0,"1":0,"2":0,"3":0,"4":0,"12":0,"13":0,"14":0,"23":0,"24":0,"34":0,"123":0,"234":0,"124":0,"134":0,"1234":0}
# 		for n1,n2 in G.edges:
# 			dseasons["".join(sorted(G.edges[n1,n2]["seasons"]))]+=1	
# 		for k,v in dseasons.items():
# 			if v != 0:
# 				print(k,":",v)
# 		print("sum :",G.size())

# 	return(G)


#print 
def print_metaniches(list_nodes,fig=None,ax=None):
	# metaniches representation #
	
	df_niches_asv = pd.read_csv("useful_tables/ASVs_metavars.tsv",sep="\t",index_col = 0)
	df_niches_seasons = pd.read_csv("useful_tables/seasons_metavars.tsv",sep="\t",index_col = 0)
	vars = ["1_Q2","2_Q2","3_Q2","4_Q2"]
	df_niches_asv = df_niches_asv[vars]
	df_niches_asv.columns = ["1","2","3","4"]

	df_niches_asv = df_niches_seasons.append(df_niches_asv.loc[list_nodes])
	
	if fig!=None:
		ax = fig.add_subplot(1, 1, 1)

	# Move left y-axis and bottim x-axis to centre, passing through (0,0)
	ax.spines['left'].set_position('center')
	ax.spines['bottom'].set_position('center')
	# Eliminate upper and right axes
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	# Show ticks in the left and lower axes only
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')

	#minimum and maximum
	ax.set_xlim(-1, 1)
	ax.set_ylim(-1, 1)

	#graduation
	ax.xaxis.set_ticks(np.arange(-1,1.25,0.5))
	ax.yaxis.set_ticks(np.arange(-1,1.25,0.5))

	#data
	xall = list()
	yall = list()
	for i in range(0,df_niches_asv.shape[0]):
	    xall.append([0,-df_niches_asv.at[df_niches_asv.index[i],df_niches_asv.columns[1]],0,df_niches_asv.at[df_niches_asv.index[i],df_niches_asv.columns[3]],0])

	    yall.append( [df_niches_asv.at[df_niches_asv.index[i],df_niches_asv.columns[0]],0,-df_niches_asv.at[df_niches_asv.index[i],df_niches_asv.columns[2]],0,df_niches_asv.at[df_niches_asv.index[i],df_niches_asv.columns[0]]])

	#plot
	colors=["green","magenta","orange","turquoise"]
	for i in range(0,4):
	    plt.fill(xall[i],yall[i],color=colors[i],alpha=0.15)

	for i in range(0+4,len(xall)):
	    plt.plot(xall[i],yall[i],'-o',color=mycolormap(i-4))

	#legend
	legend_dict = { 'Spring' : 'green', 'Summer' : 'magenta', 'Autumn' : 'orange', "Winter": "turquoise" }
	patchList = []
	for key in legend_dict:
	        data_key = mpatches.Patch(color=legend_dict[key], label=key,alpha=0.15)
	        patchList.append(data_key)
	legend1 = plt.legend(handles=patchList)
	# ax.add_artist(legend1)
	# plt.legend(df_niches_asv.index[4:],loc="lower right")


#######

def get_temp_niches(G):
	# get optimal temperature for each node #
	
	df_niches_all = pd.read_csv("useful_tables/niches_RO.tsv",sep="\t",index_col=0)
	niches_temp = df_niches_all[["T (Température)_Q1","T (Température)_Q2","T (Température)_Q3"]] #interval temp of each ASV
	for n in G.nodes:
		temprange = list(niches_temp.loc[n])
		G.nodes[n]["optimatemp"] = temprange[1]
	return(G)

####################

def make_venn_diagram_edges(graphs,names,colors,title):
	cmap = ["Set 1"]
	dgraphs = dict()
	cpt_g = 0
	for g in graphs:
		edges = list(graphs[cpt_g].edges)
		dgraphs[names[cpt_g]] = set(edges)
		cpt_g += 1
	print(len(intersection_edges(graphs)))
	venn(dgraphs, cmap = colors)
	plt.title(title,fontsize = 15,fontweight="bold")
	plt.show()


####################

def thresholding_algo(y, lag, filtre, influence):

    signals = [-1]*len(y)
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    # stdFilter[lag - 1] = np.std(y[0:lag])

    # if lag == 1:
    #     stdFilter[lag-1] = np.std(y[0:2])
    # else:
    stdFilter[lag-1] = np.std(y[0:lag])

    for i in range(lag, len(y)):
        #sommes nous dans un pic
        if abs(y[i]) > avgFilter[i-1]: #> filtre * stdFilter [i-1]:
            # if y[i] > threshold and y[i] > avgFilter[i-1]-diffavg:
            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
        else:
            # signals[i] = -1 #0 normalement
            filteredY[i] = y[i]
        avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
#         if lag == 1:
# 	        stdFilter[lag-1] = np.std(y[(i-lag+1):i+1])
#         else:
        #stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        stdFilter = np.std([avgFilter[i-1],filteredY[i]])
        if abs(y[i] - avgFilter[i-1]) > filtre * stdFilter:
            # if y[i] != 0:
            if abs(y[i]) > avgFilter[i-1] and y[i] != 0:
                signals[i] = 1
            else:
                signals[i] = -1


    list(zip(list(range(0,len(y))),list(zip(y,signals))))
    # # correction : un -1 tout seul entre deux 1 est considéré comme pb comptage
    # s = signals.copy()
    # for i in range(lag, len(y)-1):
    #     if signals[i] == -1 and s[i-1] == 1 and s[i+1] == 1:
    #         # print(i)
    #         signals[i] = 1

    y = y[::-1]
    rev_signals = [-1]*len(y)
    rev_avgFilter = [0]*len(y)
    rev_stdFilter = [0]*len(y)
    rev_avgFilter[lag-1] = np.mean(y[0:lag])
    # rev_stdFilter[lag - 1] = np.std(y[0:lag])
    # if lag == 1:
    #     rev_stdFilter[lag-1] = np.std(y[0:2])
    # else:
    rev_stdFilter[lag-1] = y[i]

    for i in range(len(y)-lag-1, len(y)):
        #sommes nous dans un pic
        if abs(y[i]) > rev_avgFilter[i-1]: #> filtre * stdFilter [i-1]:
            # if y[i] > threshold and y[i] > avgFilter[i-1]-diffavg:
            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
        else:
            # signals[i] = -1 #0 normalement
            filteredY[i] = y[i]
        rev_avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
#         if lag == 1:
# 	        stdFilter[lag-1] = np.std(y[(i-lag+1):i+1])
#         else:
        #stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        stdFilter = np.std([rev_avgFilter[i-1],filteredY[i]])
        if abs(y[i] - rev_avgFilter[i-1]) > filtre * stdFilter:
            if abs(y[i]) > rev_avgFilter[i-1] and y[i] != 0:
                rev_signals[i] = 1
            else:
                rev_signals[i] = -1

    # for i in range(lag, len(y)):
    #     if abs(y[i] - rev_avgFilter[i-1]) > filtre * rev_stdFilter [i-1]:
    #         # if y[i] > threshold and y[i] > rev_avgFilter[i-1]+diffavg:
    #         if abs(y[i]) > avgFilter[i] and y[i] != 0:
    #             rev_signals[i] = 1
    #         else:
    #             rev_signals[i] = -1

    #         filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
    #     else:
    #         rev_signals[i] = -1 #0 normalement
    #         filteredY[i] = y[i]
    #     rev_avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
    #     # if lag == 1:
    #     #     rev_stdFilter[lag-1] = np.std(y[(i-lag+1):i+1])
    #     # else:
    #     rev_stdFilter[i] = np.std(y[(i-lag+1):i+1])


    signals[0:lag] = rev_signals[-lag:] #il semblerait qu'il faille inverser les signaux ???

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))

################

def detect_regular_associations(G):
	lag, filtre, influence = 2, 0.2, 0.5
	lseasons = list(set(dfDatesSeasons["Ab4Season"]))
	lyears = list(set(dfDatesSeasons["Year"]))
	for asv1,asv2,d in G.edges(data=True):
		y1 = clr_df[asv1].to_list()
		y2 = clr_df[asv2].to_list()
		res1 = thresholding_algo(y1,lag,filtre,influence)
		res2 = thresholding_algo(y2,lag,filtre,influence)
		comp = [ 1 if (res1["signals"][i]==1 and res2["signals"][i]==1) else -1 if (res1["signals"][i]==-1 and res2["signals"][i]==-1) else 0 for i in range(0,len(res1["signals"]))]
		match_index = [i for i,c in enumerate(comp) if c == 1]
		match_dates = [ list(df_clr.Date)[i] for i in match_index ]
		match_seasons = dfDatesSeasons.loc[match_dates,["Ab4Season","Year"]]
		match_seasons["Year"] = match_seasons["Year"].astype(str)
		match = list(match_seasons.agg("-".join,axis=1))
		for s in lseasons:
			ls = [s+"-"+str(y) for y in lyears]
			check =  all(item in match for item in ls)
			if check:
				d["regular"]=True
				break

##################

def show_regular_associations(asv1,asv2):
	lag, filtre, influence = 2, 0.2, 0.5
	y1 = clr_df[asv1].to_list()
	y2 = clr_df[asv2].to_list()
	res1 = thresholding_algo(y1,lag,filtre,influence)
	res2 = thresholding_algo(y2,lag,filtre,influence)
	comp = [ 1 if (res1["signals"][i]==1 and res2["signals"][i]==1) else -1 if (res1["signals"][i]==-1 and res2["signals"][i]==-1) else 0 for i in range(0,len(res1["signals"]))  ]

	fig = plt.figure(figsize=(20,15))

	ax1 = fig.add_subplot(411)
	plt.xlim((0, 69))
	plt.ylim((-1, 1.1))
	plt.step(range(0,69),comp,color="black",where="mid")

	ax2 = fig.add_subplot(412)
	pd.DataFrame([res1["signals"],res2["signals"]]).T.plot(drawstyle="steps-mid",legend=False,alpha = 0.7, ax=ax2,xticks=[],xlim = (0,69))
													   
	ax3 = fig.add_subplot(413)
	clr_df[[asv1,asv2]].plot(kind="line",ax=ax3,
						 xticks = (0, 69),
						 xlim = (0,69),
# 						 xticks=list(df_clr.index),
						 rot=90,
						 legend=False
 						 # ,
 						 # xlim = (df_clr.index[0],df_clr.index[-2])
						 )

	ax4 = fig.add_subplot(414)
	norm_df[[asv1,asv2]].plot(kind="line",ax=ax4,
 						 xticks = (0, 69),
 						 xlim = (0,69),
					     legend=False,
# 						 xticks=list(df_clr.index),
 						 rot=90
 						 # ,
 						 # xlim = (df_clr.index[0],df_clr.index[-2])
 						 )
	
	plt.ion()
	plt.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	plt.show()
	plt.draw()
	plt.pause(0.001)
	input("\n  Press [enter] to continue.")






# https://github.com/rkistner/chinese-postman/issues/21
def connected_component_subgraphs(G):
	# list of connected components #
	subgraphs = list()
	for c in sorted(nx.connected_components(G), key=len, reverse=True):
		subgraphs.append(G.subgraph(c))
	lsubgraphs = list()
	for i, sg in enumerate(subgraphs):
		lsubgraphs.append(sg)
	return(lsubgraphs)


#####################

def get_width_edges(G):
# for network display #
	nbr_nodes = len(G.nodes)
	#size of edges
	weightEdges = list(nx.get_edge_attributes(G,'weight').values())
	weightEdges = [abs(w) for w in weightEdges]
	if nbr_nodes > 100:
		max_weightEdges = 0.5
	elif nbr_nodes > 30 :
		max_weightEdges = 2.5
	else :
		max_weightEdges = 8
	try:
		return([ (i - 0)/(1) * (max_weightEdges - 0.05) + 0.5  for i in weightEdges] )
	except : #if max = min to avoid division by 0 error
		return(1)

#####################

def get_width_nodes(G):
# for network display #
	if G.order() > 70:
		max_widthNodes = 2.5
		min_widthNodes = 1
	else :
		max_widthNodes = 5
		min_widthNodes = 1
	return(max_widthNodes,min_widthNodes)


######################

def filtering_edges(G,min_cor,verbose=True,color=True,sign=None):
# filter edges according weight of correlation #

	#filtering edges
	removeEdges = []
	for N1,N2,E in G.edges(data=True):
		weight = E["weight"]
		E["style"] = "solid"
		E["color"] = "black"
		# dashed for negative edges
		if weight < 0 :
			if color :
				E["color"] = "#CF2900" 
			if sign == 1 : #if we want positive correlations only
				removeEdges.append([N1,N2])
				# print(N1,N2)
		else: # negative correlations only
			if sign == 0 :
			 	removeEdges.append([N1,N2])
			# print("Negative edge :",weight)
		if abs(weight) < min_cor and [N1,N2] not in removeEdges:
			removeEdges.append([N1,N2])

	if verbose:
		print("Removing ",len(removeEdges)," edges under minimum correlation threshold or with the wrong desired sign of correlation ... (",min_cor,"). ",G.size()," edges remaining.",sep="")

	# remove under threshold edges
	# must be done after the prec bloc, otherwise modification of G within it
	for N1,N2 in removeEdges:
		G.remove_edge(N1,N2)

	# remove isolate nodes
	isolateNodes = list(nx.isolates(G))
	G.remove_nodes_from(isolateNodes)
	if verbose:
		print("Removing ",len(isolateNodes)," isolate nodes ... (",G.order()," nodes remaining)\n",sep="")
	return(G)

	###########

def get_closeness_network(G):
	# get closeness centrality of each node #
	closeness = nx.closeness_centrality(G)
	closeness = {k:round(v,4) for k,v in closeness.items()}
	nx.set_node_attributes(G, closeness, "closenesscentrality")
	return(G)

	###########


def get_degrees_network(G):
	# get degrees of each node #
	degrees = dict(nx.degree(G))
	nx.set_node_attributes(G, degrees, "degrees")
	return(G)

	###########

def get_betweenness_network(G):
	# get betweenness centrality of each node #
	betweenness = nx.betweenness_centrality(G)
	betweenness = {k:round(v,4) for k,v in betweenness.items()}
	nx.set_node_attributes(G, betweenness, "betweennesscentrality")
	return(G)

	###########

def get_degrees_centrality_network(G):
	# get degrees centrality of each node #
	degrees = nx.degree_centrality(G)
	degrees = {k:round(v,4) for k,v in degrees.items()}
	nx.set_node_attributes(G, degrees, "degreescentrality")
	return(G)

	###########

def get_subgraph_centrality_network(G):
	# get subgraph centrality of each node #
	subgraph = nx.subgraph_centrality(G)
	subgraph = {k:round(v,4) for k,v in subgraph.items()}
	nx.set_node_attributes(G, subgraph, "subgraphcentrality")
	return(G)

	###########

def leading_eigen_clustering(G):
	# posG = G

	# attributes_nodes = ['mv', 'amplicon', 'occurrence', 'sequence', 'silvaidentity', 'silvareferences', 'silvataxonomy', 'IDTAXAtaxonomy', 'IDTAXAconfidence', 'pr2identity', 'pr2references', 'pr2taxonomy', 'COIdbidentity', 'COIdbreferences', 'COIdbtaxonomy', 'fraction']

	# attributes_edges = list(list(G.edges(data=True))[0][2].keys())[1:]

	# for N in posG.nodes:
	# 	for an in attributes_nodes:
	# 		del(posG.nodes[N][an])

	# for N1,N2,E in posG.edges(data=True):
		# for ae in attributes_edges:
		# 	del[E[ae]]
	# for N1,N2,E in posG.edges(data=True):
	# 	if E["weight"] < 0:
	# 		E["weight"] = abs(E["weight"])
	nx.write_graphml(G,"jeter.gml")

	iG = ig.read("jeter.gml",format="graphml")

	clusters = iG.community_leading_eigenvector(weights="weight")#,clusters=30)

	nbr_clusters = len(clusters)
	print(nbr_clusters,"cluster(s) was/were detected with Leading Eigenvector Algorithm !")

	cpt_c = 0
	for c in clusters.subgraphs():
		vs = ig.VertexSeq(c)
		for v in c.vs:

			# print("cluster",cpt_c+1)
			id = v["id"]

			G.nodes(data=True)[id]["clustereigen"] = cpt_c
		cpt_c += 1
	return(G, nbr_clusters)

	###########

def edge_betweenness_clustering(G):

	nx.write_graphml(G,"jeter.graphml")
	iG = ig.read("jeter.graphml",format="graphml")

	clusters = iG.community_edge_betweenness()
	clusters = clusters.as_clustering()

	nbr_clusters = len(clusters)
	print(nbr_clusters,"cluster(s) was/were detected with Edge Betweenness Clustering Algorithm !")

	cpt_c = 0
	for c in clusters.subgraphs():
		vs = ig.VertexSeq(c)
		for v in c.vs:
			id = v["id"]
			G.nodes(data=True)[id]["clusterbetweenness"] = cpt_c
		cpt_c += 1

	return(G, nbr_clusters)

#########


def intersection_nodes(graph,reference_graph):
	return list(set(graph.nodes).intersection(set(reference_graph.nodes)))

def intersection_nodes_multigraphs(graphs): 
	inter = set(graphs[0].nodes)
	for i in range(1, len(graphs)):
		inter = inter.intersection(set(graphs[i].nodes))
	return inter

def exclusive_nodes(graph,graphs):
    return ([n for n in graph.nodes if n not in set(sum([list(g.nodes) for g in graphs],[]))])

def intersection_edges(graph,reference_graph):
    inter = [(n1,n2) for (n1,n2) in reference_graph.edges if (n1,n2) in graph.edges() or (n2,n1) in graph.edges()]
    return inter

def intersection_edges_multigraphs(graphs): 
	inter = set(graphs[0].edges)
	for i in range(1, len(graphs)):
		edges =  list(graphs[i].edges)
		reverse_edges = [ e[::-1] for e in edges] #to get identical edges whatever the sense they are written
		edges.extend(reverse_edges)
		inter = inter.intersection(edges)
	return inter

def exclusive_edges(graph,graphs):
    return ([e for e in graph.edges if e not in set(sum([list(g.edges) for g in graphs],[]))])


#########

#no division in subgraphs (already subgraph is preferred)
#no clustering
def simple_draw_network(G, nodes = None, edges = None, myax=None, clustering = None):
	# basic display of the graph#

	#make_graph_topology(G)

	if myax == None:
		f, myax = plt.subplots(1, 1,figsize=(15,15))

	plt.axis('off')
	plt.tight_layout(pad=2, w_pad=5, h_pad=3.0)
	
	datasets = {"18SV1V2" : ["o","pr2taxonomy"],
		"18SV4" : ["d","pr2taxonomy"],
		"COI" : ["^","COIdbtaxonomy"],
		"16SV4V5" : ["p","silvataxonomy"]}

	degree = dict(nx.degree(G))
	nx.set_node_attributes(G, degree, "degree")

	closeness = nx.closeness_centrality(G).values()
	closeness = [ round(c,2) for c in closeness ]
	closenessRank = [int(i)-1 for i in list(pd.DataFrame(closeness).rank(method="dense",ascending=False)[0])]
	Dcloseness = dict(zip(list(G.nodes), closeness)) 
	DclosenessRank = dict(zip(list(G.nodes), closenessRank)) 
	nx.set_node_attributes(G, Dcloseness, "rankcloseness")

	font_size = get_font_size(G)*2

	# colorsLittle = list(Color("#D8F0FF").range_to(Color("#5959FF"),int(max(degree.values())+1)))
	# colorsBig = list(Color("#EEFFD5").range_to(Color("green"),int(max(degree.values())+1)))
	labelNodes = dict()
	for n,d in G.nodes(data=True):
		labelNodes[n] = n
	# 	if d["fraction"] == "3":
	# 		d["color"] = str(colorsBig[d["degree"]])
	# 	else:
	# 		d["color"] = str(colorsLittle[d["degree"]])
	# G = optima_projection(G)


	pos = nx.nx_agraph.graphviz_layout(G, prog='sfdp', args="-start=1")

	sizeMaxNodes = get_sizeMaxNodes(G)*4
	widthEdges = get_width_edges(G)
	colorEdges = list(nx.get_edge_attributes(G,'color').values())
	styleEdges = list(nx.get_edge_attributes(G,'style').values())
	labelEdges = nx.get_edge_attributes(G,'label')

	nx.draw_networkx_nodes(G,
				pos=pos,
				node_color = "#636363",
				node_size = 1,
				linewidths = 0,
				alpha = 0.3
				)


	#interactive visualisation
	dataset = list(nx.get_node_attributes(G,"dataset").values())
	tax = [l.split("|") for l in list(nx.get_node_attributes(G,"taxonomy").values())]
	tax = [ "\n".join(l) for l in  tax]
	idtax = list(nx.get_node_attributes(G,"idtaxonomy").values())
	idtax = [ str(l)+"%" for l in  idtax]
	rank = list(nx.get_node_attributes(G,"keystoneindexrank").values())
	degrees = list(nx.get_node_attributes(G,"degrees").values())
	closeness = list(nx.get_node_attributes(G,"closenesscentrality").values())
	betweenness = list(nx.get_node_attributes(G,"betweennesscentrality").values())
#	try:
#		seasons = list(nx.get_node_attributes(G,"seasons").values())
#		ziplist = list(zip(dataset,tax,idtax,rank,degrees,closeness,betweenness,seasons))
#		attributes = ["dataset : ","","id taxonomy :","\nrank : ","degrees : ","closeness : ","betweenness : ","\nseasons : "]
#	except :
	ziplist = list(zip(dataset,tax,idtax,rank,degrees,closeness,betweenness))
	attributes = ["dataset : ","","id taxonomy :","\nrank : ","degrees : ","closeness : ","betweenness : "]
	
	labels = list()
	for l in ziplist:
		cpt_a = 0
		object = ""
		for ll in l:
			attribute = attributes[cpt_a]+str(ll)+"\n"
			object += attribute
			cpt_a += 1
		labels.append(object[:-1])
	
	cursor = mplcursors.cursor(myax, hover=True)
	cursor.connect(
	"add", lambda sel: sel.annotation.set_text(labels[sel.target.index])) #when clicking on the nodes, display its full name

	widthNodesG = list(nx.get_node_attributes(G,"closenesscentrality").values()) #get min and max of the whole subgraph and not only within the dataset !

	if nodes != None or edges != None:
		cmap = None
		vmin = None
		vmax = None
		nx.set_node_attributes(G,"#C3C3C3","color")
		if nodes != None:
			if isinstance(nodes, list):
				for node in nodes:
					if node in G.nodes:
						G.nodes[node]["color"] = "#FFBC00"
			else:
				G.nodes[nodes]["color"] = "#FFBC00"
		if edges!=None:# and myax == None:
			if isinstance(edges, list):
				for edge in nodes:
					if edge in G.edges:
						G.edges[edge]["color"] = "#FFBC00"
			else:
				G.edges[edge]["color"] = "#FFBC00"
	elif clustering == None:
		# colorsG = list(nx.get_node_attributes(G,"optimatemp").values())
		nx.set_node_attributes(G,nx.get_node_attributes(G,"optimatemp"),"color")
		cmap = cmocean.cm.thermal
		vmin = 8
		vmax = 17
	elif clustering == "eigen":
		G,nbr_clusters = leading_eigen_clustering(G)
		vmin = min(nx.get_node_attributes(G,"clustereigen"))
		vmax = max(nx.get_node_attributes(G,"clustereigen"))
		cmap = None
		for n,d in G.nodes(data=True):
			d["color"] = mycolormap(d["clustereigen"])
	elif clustering == "betweenness":
		G,nbr_clusters = edge_betweenness_clustering(G)
		vmin = min(nx.get_node_attributes(G,"clusterbetweenness"))
		vmax = max(nx.get_node_attributes(G,"clusterbetweenness"))
		cmap = None
		for n,d in G.nodes(data=True):
			d["color"] = mycolormap(d["clusterbetweenness"])


	

	# dcolors={"":"#DADADA","1":"green","2":"#CE4E9F","3":"#F68230","4":"#5AEFFC","12":"grey","13":"grey","14":"grey","23":"grey","24":"grey","34":"grey","123":"grey","234":"grey","124":"grey","134":"grey","1234":"grey"}
	
	# cmap =  plt.cm.Purples


	for dataset in datasets.keys():

		nodelist = [N for N,D in G.nodes(data=True) if D['dataset']==dataset] 
		sG = G.subgraph(nodelist)

		colors = list(nx.get_node_attributes(sG,"color").values())

		# colors = list(nx.get_node_attributes(sG,"seasons").values())
		# for i in range(0,len(colors)):
		# 	colors[i] = dcolors[colors[i]]

		widthNodes = list(nx.get_node_attributes(sG,"closenesscentrality").values())

		max_widthNodes,min_widthNodes = get_width_nodes(sG)
		# to see the difference of the contour of the nodes
		widthNodes = [ (w - min(widthNodesG))/(max(widthNodesG) * (max_widthNodes - min_widthNodes)) + min_widthNodes for w in widthNodes ]

		# sizeNodes = nx.get_node_attributes(sG,"abundance").values()
		# sizeNodes = [int(float(i)) for i in sizeNodes]
		# sizeNodes = [(A / max(sizeNodes)) * (sizeMaxNodes - 50 ) + 500 for A in sizeNodes]


		nx.draw_networkx_nodes(sG,
				pos=pos,
				node_color = "#636363",
				node_shape = datasets[dataset][0],
				node_size = 1,
				linewidths = 0,
				alpha = 0.3,
				label = dataset,
				)

			# draw nodes
		nx.draw_networkx_nodes(sG,
			pos=pos,
			node_color = colors,
			node_shape = datasets[dataset][0],
			# node_size = sizeNodes,
			alpha=0.9,
			linewidths= widthNodes,
			edgecolors = "#636363",
			cmap = cmap,
			vmin = vmin,
			vmax = vmax
			)


	# draw edges
	nx.draw_networkx_edges(G,
		pos = pos,
		edge_color = colorEdges,
		style = styleEdges,
		width = widthEdges,
		)

	# draw labels of the nodes (taxonomic data)
	nx.draw_networkx_labels(G,
		pos = pos,
		labels = labelNodes,
		font_size = font_size
	)

	# if G.size() <  100 :
	# 	nx.draw_networkx_edge_labels(G,
	# 		pos = pos,
	# 		font_size = 8,
	# 		edge_labels = labelEdges
	# 	)

	plt.legend(numpoints = 1,fontsize = 10,markerscale = 10)

	if myax == None:
		plt.savefig("out.png", bbox_inches='tight',figsize=(20,20), dpi = 300)
		# plt.ion()
		# plt.show()
		# plt.draw()
		# plt.pause(0.001)
		# input("\n  Press [enter] to continue.")


	# plt.savefig(output+"_"+str(cpt)+".png", bbox_inches='tight',figsize=(20,30), dpi = 500)
	# plt.clf() #clean plot
	# print("-> Saved ",output,"_",cpt,".png !",sep="")
