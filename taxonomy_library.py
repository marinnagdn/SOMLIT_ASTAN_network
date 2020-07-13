from network_library import *



def get_metavars(G):
# compute metavariables for the metaniches #

	df_niches = pd.read_csv("useful_tables/ASVs_metavars.tsv",sep="\t",index_col = 0)
	for i in range(1,5):
		metavar = str(i)+"_Q2"
		name = "metavar"+str(i)
		for n in G.nodes:
			G.nodes[n][name] = round(df_niches.at[n,metavar],3)
	return(G)

#####

# def study_metavar_by_taxo_group(G,mydict):
# 	vars = ["metavar1","metavar2","metavar3","metavar4"]
# 	for k,v in mydict.items():
# 		print(k,":")
# 		metavars = defaultdict(list)
# 		for node in v:
# 			for var in vars:
# 				metavars[var].append(G.nodes[node][var])
# 		for km,vm in metavars.items():
# 			print(km,":",vm)
# 				# print(var,":",G.nodes[node][var])

######

def get_nodes_taxonomy(G):
	# get taxonomy for each node #	
	datasets = {"18SV1V2" : "silva",
		"18SV4" : "silva",
		"COI" : "COIdb",
		"16SV4V5" : "silva"}
	for n,d in G.nodes(data=True):
		d["taxonomy"] = d[datasets[d["dataset"]]+"taxonomy"]
		d["idtaxonomy"] = d[datasets[d["dataset"]]+"identity"].replace(",",".")
		d["reftaxonomy"] = datasets[d["dataset"]]
	return(G)


#####
#def get_network_taxonomy(G, maxranktax = 3):

#	# for n,d in G.nodes(data=True):
#	# 	d["taxonomy"] = d[datasets[d["dataset"]][1]]

#	for n1,n2,d in G.edges(data=True):
#		n1tax = G.nodes[n1]["taxonomy"].split("|")
#		n2tax = G.nodes[n2]["taxonomy"].split("|")
#		n1seq = G.nodes[n1]["sequence"]
#		n2seq = G.nodes[n2]["sequence"]
#		n1dataset = G.nodes[n1]["dataset"]
#		n2dataset = G.nodes[n2]["dataset"]
#		d["tax"] = ""

#		if "Unknown" in n1tax or "Unknown" in n2tax or "None" in n1tax or "None" in n2tax:
#				# print("\nUnknown : ")
#				# print(n1tax)
#				# print(n2tax)
#				d["tax"] = "?"
#		else:
#			minlen = min(len(n1tax),len(n2tax)) 
#			if n1tax[0:minlen] == n2tax[0:minlen]:
#				if minlen > maxranktax:
#					# print("\nProbably the same : ")
#					# print(n1tax)
#					# print(n2tax)
#					if n1dataset == n2dataset:
#						if n1seq == n2seq:
#							d["tax"] = "==" #même sequence
#						else :
#							d["tax"] = "=/=" #pas même sequence
#					else:
#						d["tax"] = "==?" #almost sure they're the same : merge them? #different dataset : we cannot compare
#						# print("Almost sure :")
#						# print(n1tax,n1dataset)
#						# print(n2tax,n2dataset)
#				else:
#					# print("\n Not sure the same : ")
#					# print(n1tax)
#					# print(n2tax)
#					if n1dataset == n2dataset:
#						if n1seq == n2seq:
#							d["tax"] = "==" #même sequence
#						else :
#							d["tax"] = "=/=" #pas même sequence
#					else:
#						d["tax"] = "=?" #different dataset : we cannot compare
#			else:
#				same = False
#				i = minlen-1
#				# print(minlen,n1tax,n2tax)
#				while same == False and i >= 0: #on s'arrete avant le regne parce que dire que ce sont des eucaryotes = genial
#					
#					if n1tax[i] == n2tax[i]:
#						same=True
#						d["tax"] = "="+str(i)
#						# print("\nSame tax up to",i,":")
#						# print(n1tax)
#						# print(n2tax)
#			
#					else:
#						i -= 1
#				
#				if same == False:
#					# print("\nDifferent : ")
#					# print(n1tax)
#					# print(n2tax)
#					d["tax"] = "=/="
#					# d["style"] = "dashed"

#	# cpt_same = 0
#	print("\n*** Description of the taxonomies ***")
#	cpt_same = 0
#	cpt_probably_same = 0
#	cpt_little_probably_same = 0
#	cpt_different = 0
#	cpt_inconnu = 0
#	for n1,n2,d in G.edges(data=True):
#		tax = d["tax"]
#		if tax == "==":
#			cpt_same+=1
#		elif tax == "==?":
#			cpt_probably_same+=1 #different dataset we can't compare
#		elif tax == "=?":
#			cpt_little_probably_same+=1 #different dataset we can't compare and < maxtaxrank
#		elif tax == "?":
#			cpt_inconnu +=1
#		# elif tax == "==":
#		# 	cpt_same +=1
#		else :
#			cpt_different +=1

#	# print("Very probably the same :",cpt_same,"interactions")
#	print("Same sequences :",cpt_same,"interactions ->",round(cpt_same/G.size()*100,2),"%")
#	print("Probably the same :",cpt_probably_same,"->",round(cpt_probably_same/G.size()*100,2),"%")
#	print("Little probably the same :",cpt_little_probably_same,"->",round(cpt_little_probably_same/G.size()*100,2),"%")
#	print("Different :",cpt_different,"->",round(cpt_different/G.size()*100,2),"%")
#	print("With an unknown node :",cpt_inconnu,"->",round(cpt_inconnu/G.size()*100,2),"%")

#	return(G)

	###########

#def clean_ambiguous_interactions(G,maxranktax = 2):

#	ambiguous_edges = list()
#	cpt_id_nodes = 0
#	for n1,n2,d in G.edges(data=True):
#		if d["tax"] == "=?":
#			# print("\n***\n",G.nodes[n1]["taxonomy"],"\n",G.nodes[n1]["taxonomy"])
#			cpt_id_nodes += 1
#		if d["tax"] == "=?" or d["tax"] == "?" or d["tax"] == "==?" :
#			ambiguous_edges.append([n1,n2])
#		else:
#			try:
#				ranktax =  int(d["tax"][1])
#				if ranktax > maxranktax:
#					ambiguous_edges.append([n1,n2])
#			except:
#				pass

#	print("Removing ",len(ambiguous_edges)," ambiguous edges ... (",round((len(ambiguous_edges)/G.size()*100),2),"%) with ",cpt_id_nodes," possible identical species nodes interacting (",round((cpt_id_nodes/G.size()*100),2),"%)",sep="")
#	G.remove_edges_from(ambiguous_edges)
#	G.remove_nodes_from(list(nx.isolates(G)))
#	return (G)
def contract_neighrbors_network(G):
	cpt_contraction = 1
	while cpt_contraction != 0:
		dcontracted = dict()
		cpt_contraction = 0
		for n1,n2,d in G.edges(data=True):
			while n1 not in G.nodes():
				n1 = dcontracted[n1]
			while n2 not in G.nodes():
				n2 = dcontracted[n2]
			if n1==n2:
				continue
			d1 = G.nodes(data=True)[n1]
			d2 = G.nodes(data=True)[n2]
			if d1["dataset"] == d2["dataset"]:
				if d1["sequence"] == d2["sequence"]:
					G = nx.contracted_nodes(G, n1, n2, self_loops=False)
					dcontracted[n2] = n1
					cpt_contraction += 1
				elif d1["taxonomy"] == d2["taxonomy"]:
					G = nx.contracted_nodes(G, n1, n2, self_loops=False)
					dcontracted[n2] = n1
					cpt_contraction += 1
			elif d1["taxonomy"] == d2["taxonomy"]:
				G = nx.contracted_nodes(G, n1, n2, self_loops=False)
				cpt_contraction += 1
				dcontracted[n2] = n1
		print("nombre de contractions :",cpt_contraction)
	return(G)

##############

def contract_network(G):
	cpt_contraction = 1
	while cpt_contraction != 0:
		dcontracted = dict()
		cpt_contraction = 0
		for n,d in G.nodes(data=True):
			if n in dcontracted.keys():
				n = dcontracted[n]
			d = G.nodes(data=True)[n]
# 			print(n)
			tax = d["taxonomy"]
			seq = d["sequence"]
			revG = G.copy()
			revG.remove_node(n)
			for revn,revd in revG.nodes(data=True):
				revseq = revd["sequence"]
				revtax = revd["taxonomy"]
				if revseq == seq or revtax == tax:
					G = nx.contracted_nodes(G, n, revn, self_loops=False)
					dcontracted[revn] = n
					cpt_contraction += 1
		print("nombre contractions :",cpt_contraction)
	return(G)

#############

def filter_identitity_taxonomy(G, thr):
	#remove nodes under a threshold of identity of taxonomy#
	
	removeNodes = list()
	idtax = nx.get_node_attributes(G,"idtaxonomy")
	# print(idtax)
	for k,v in idtax.items():
		try:
			if float(v) < thr:
				# print(k,v)
				removeNodes.append(k)
		except ValueError:
			pass

	print("Removing",len(removeNodes),"under threshold taxonomy identity of",thr,"% ...")

	G.remove_nodes_from(removeNodes)
	G.remove_nodes_from(list(nx.isolates(G)))

	return(G)

#############

ranks = ["domain","phylum","class","order","family","genus"]
def print_a_node(G,k,rank,taxonomy=False,cross_database=False):
	if taxonomy:
			print("\n * rank ",rank[k],":",k)
			print_taxonomy(G,k)
			print(" ",G.nodes[k]["degrees"]," degree(s) - closeness : ",G.nodes[k]["closenesscentrality"]," - betweenness : ", G.nodes[k]["betweennesscentrality"]," - subgraph : ",G.nodes[k]["subgraphcentrality"],sep="")
			try:

				print(" season(s) :",G.nodes[k]["seasons"])
			except:
				pass
			# print(" temperature niche : [",G.nodes[k][""])
			# print(" interaction season(s) :",G.edges[node_pivot,k]["seasons"])
	else:
			print("\n * rank ",rank[k]," : ",k," (",G.nodes[k]["dataset"],") - ",G.nodes[k]["degrees"]," degree(s) - closeness : ",G.nodes[k]["closenesscentrality"]," - betweenness : ", G.nodes[k]["betweennesscentrality"]," - subgraph : ",G.nodes[k]["subgraphcentrality"],sep="")
			try:
				print(" season(s) :",G.nodes[k]["seasons"])
			except:
				pass
	
	if cross_database:
		for e in G.edges(k,data=True):
			if e[2]["pidalevel"] != -1 and e[2]["globilevel"] != -1:
				print("\t*",e[1])
				print("\t pida : ",ranks[e[2]["pidalevel"]])
				print("\t globi : ",ranks[e[2]["globilevel"]])
			elif e[2]["pidalevel"] != -1:
				print("\t*",e[1])
				print("\t pida : ",ranks[e[2]["pidalevel"]])
			elif e[2]["globilevel"] != -1:
				print("\t*",e[1])
				print("\t globi : ",ranks[e[2]["globilevel"]])



	#############

def print_list_nodes(G, subG = None,taxonomy = False):

	rank = nx.get_node_attributes(G,"keystoneindexrank")
	rank = {k: v for k, v in sorted(rank.items(), key=lambda item: item[1],reverse=True)}

	if subG != None:
		rank = {k:v for k,v in rank.items() if k in subG.nodes()}

	if taxonomy:
		for k,v in rank.items():
			print_a_node(G,k,rank,taxonomy=True)
	else:
		for k,v in rank.items():
			print_a_node(G,k,rank, taxonomy=False)

	

#####

# def print_node_interactions(G,n, taxonomy=False):
# 	edges = G.edges(n)
# 	print("\n * Interacting with : ")
# 	for e in edges:
# 		tax = G.nodes[e[1]]["taxonomy"].split("|")
# 		try:
# 			print("\n  - ",e[1]," (weight : ",round(G.edges[n,e[1]]["weight"],2),", Interaction seasons : ",G.edges[n,e[1]]["seasons"],")",sep="")
# 			print("  Temperature niche overlap :",G.edges[n,e[1]]["tempoverlap"],"%")
# 			print("  Node seasons : ",G.nodes[e[1]]["seasons"])
# 			common_seasons = [ s for s in G.nodes[n]["seasons"] if s in G.nodes[e[1]]["seasons"] ]
# 			print("  Common seasons :",common_seasons)
# 		except:
# 			print("\n  - ",e[1]," (weight : ",round(G.edges[n,e[1]]["weight"],2),")",sep="")
		
# 		print("  Dataset : ",G.nodes[e[1]]["dataset"])
# 		# print("  Niche : ",G.nodes[e[1]]["temprange"])
		
# 		if taxonomy :
# 			for i in range(0,len(tax)):
# 				print(" "*(i+3),i,":",tax[i])

def print_node_interactions(G,n, n2=None, taxonomy=False, clustering = None, figure=True,allInteractions=True):
	edges = G.edges(n,data=True)
	ranks = {(e[0],e[1]):e[2]["weight"] for e in edges}
	ranks = dict(sorted( ranks.items(), key= lambda k: k[1], reverse=False ))

	if n2 == None:

		print("\n * Interacting with : ")
		cpt = 1
		for e in ranks:
			tax = G.nodes[e[1]]["taxonomy"].split("|")
			try:
				print("\n",cpt,") ",e[1]," (weight : ",round(G.edges[n,e[1]]["weight"],2),", Interaction seasons : ",G.edges[n,e[1]]["seasons"],")",sep="")
	# 			print("  Temperature niche overlap :",G.edges[n,e[1]]["tempoverlap"],"%")
				print("  Node seasons : ",G.nodes[e[1]]["seasons"])
				common_seasons = [ s for s in G.nodes[n]["seasons"] if s in G.nodes[e[1]]["seasons"] ]
				print("  Common seasons :",common_seasons)
			except:
				print("\n",cpt,") ",e[1]," (weight : ",round(G.edges[n,e[1]]["weight"],2),")",sep="")
			
			print("  Dataset : ",G.nodes[e[1]]["dataset"])
			# print("  Niche : ",G.nodes[e[1]]["temprange"])
			if taxonomy :
				for i in range(0,len(tax)):
					print(" "*(i+3),i,":",tax[i])
			cpt += 1


	if figure :
		if allInteractions:
			list_nodes = [r[1] for r in ranks]
			list_nodes.append(n)
			fig = plt.figure(figsize = (20,15)) #frameon = False for no background
			ax1 = fig.add_subplot(221)	
			print_metaniches(sorted(list_nodes),ax=ax1)
			ax3 = fig.add_subplot(223)
			norm_df[sorted(list_nodes)].plot(kind="line",ax=ax3,xticks = list(norm_df.index)[::2],rot=90,fontsize=6,color=mycolormap())
			ax = fig.add_subplot(122)
			sg = G.subgraph(list_nodes)
			simple_draw_network(sg, myax=ax, clustering =clustering)
			plt.ion()
			fig.subplots_adjust(bottom=0.1)
			plt.show()
			plt.draw()
			plt.pause(0.001)
			input("\n  Press [enter] to continue.")
		else:

			if n2 == None:
				nbr_int= 1
				nbr_int = int(input("\n Give the number of the interaction (0 to exit) : "))
				while nbr_int != 0:
					list_nodes = list(list(ranks.keys())[nbr_int-1])
					fig = plt.figure(figsize = (20,15)) #frameon = False for no background
					ax1 = fig.add_subplot(221)	
					print_metaniches(sorted(list_nodes),ax=ax1)
					ax3 = fig.add_subplot(223)
					norm_df[sorted(list_nodes)].plot(kind="line",ax=ax3,xticks = list(norm_df.index)[::2],rot=90,fontsize=6,color=mycolormap())
					ax = fig.add_subplot(122)
					sg = G.subgraph(list_nodes)
					simple_draw_network(sg, myax=ax, clustering =clustering)
					plt.ion()
					fig.subplots_adjust(bottom=0.1)
					plt.show()
					plt.draw()
					plt.pause(0.001)
					input("\n  Press [enter] to continue.")
					# show_regular_associations(list_nodes[0],list_nodes[1])
					nbr_int = int(input("\n Give the number of the interaction (0 to exit) : "))
			else:
				list_nodes = [n,n2]
				fig = plt.figure(figsize = (20,15)) #frameon = False for no background
				ax1 = fig.add_subplot(221)	
				print_metaniches(sorted(list_nodes),ax=ax1)
				ax3 = fig.add_subplot(223)
				norm_df[sorted(list_nodes)].plot(kind="line",ax=ax3,xticks = list(norm_df.index)[::2],rot=90,fontsize=6,color=mycolormap())
				ax = fig.add_subplot(122)
				sg = G.subgraph(list_nodes)
				simple_draw_network(sg, myax=ax, clustering =clustering)
				# plt.ion()
				fig.subplots_adjust(bottom=0.1)
				# plt.draw()
				# plt.pause(0.001)
				# input("\n  Press [enter] to continue.")
			plt.show()



#####
def print_taxonomy(G,n):

	datasets = {"18SV1V2" : "pr2",
		"18SV4" : "pr2",
		"COI" : "COIdb",
		"16SV4V5" : "silva"}
	dataset = G.nodes[n]["dataset"]
	print(" -",n,":")
	try:
		print(" * seasons :",G.nodes[n]["seasons"])
	except:
		pass
	print(" * Dataset :",dataset)
	print(" * Taxonomy (",datasets[dataset],", ",G.nodes[n][datasets[dataset]+"identity"],"%) :",sep="")
	tax = G.nodes[n]["taxonomy"].split("|")
	for i in range(0,len(tax)):
		print(" "*(i+1),i,":",tax[i])

#######
def print_proportions_interactions(G,verbose=True):

	cptPP = 0
	cptPE = 0 
	cptEE = 0
	cptP = 0
	cptE = 0

	for n1,n2,e in G.edges(data=True):
		d1 = G.nodes[n1]["dataset"]
		d2 = G.nodes[n2]["dataset"]
		if d1[0:2] == "16" and d2[0:2] == "16":
			cptPP+= 1
		elif (d1[0:2] == "18" or d1[0:2] == "CO") and (d2[0:2] == "18" or d1[0:2] == "CO"):
			cptEE += 1

	cptPE = G.size() - cptPP - cptEE

	for n,d in G.nodes(data=True):
		if d["dataset"][0:2] == "16":
			cptP += 1
		elif d["dataset"][0:2] == "18" or d["dataset"][0:2] == "CO":
			cptE += 1

	if verbose:
		print("\tNumber of Prokaryotes : ",cptP, " (",round(cptP/G.order()*100,1),"%)",sep="")
		print("\tNumber of Eukaryotes : ",cptE, " (",round(cptE/G.order()*100,1),"%)",sep="")
		
		print("\tNumber of Prokaryotes - Prokaryotes interactions : ", cptPP, " (",round(cptPP/G.size()*100,1),"%)",sep="")
		print("\tNumber of Eukaryotes - Eukaryotes interactions : ", cptEE, " (",round(cptEE/G.size()*100,1),"%)",sep="")
		print("\tNumber of Prokaryotes - Eukaryotes interactions : ", cptPE, " (",round(cptPE/G.size()*100,1),"%)",sep="")

	return(cptP,cptE)

#######

def print_proportions_prokaryotes(G,cptP):

	# groups = ["Alphaproteobacteria","Gammaprotebacteria","Bacteroidia","Verrucomicrobiae","Flavobacteriia","Actinobacteria","Cyanobacteriia","Deltaproteobacteria","Opitutae","Planctomycetes","Betaproteobacteria","Acidimicrobiia","Thermoplasmata"]
	groups = ['Verrucomicrobiota', 'Alphaproteobacteria', 'Gammaproteobacteria','Cyanobacteria', "Betaproteobacteria",'Thermoplasmatota', 'Actinobacteriota', 'Planctomycetota', 'Bacteroidota', 'Margulisbacteria', 'Crenarchaeota', 'Marinimicrobia', 'Desulfobacterota', 'Nitrospinota', 'Myxococcota', 'Bdellovibrionota', 'Fusobacteriota', 'Chloroflexi', 'Dadabacteria', 'SAR324_clade', 'Gemmatimonadota']
	# groups = list(set([ d["taxonomy"].split("|")[1] for n,d in G.nodes(data=True) if len(d["taxonomy"].split("|"))>=2 and d["dataset"][0:2] == "16"]))
	# ids = [] * len(groups)
	dgroups = defaultdict(list)

	for n,d in G.nodes(data=True):
		try:
			tax = d["taxonomy"]

			for g in groups:
				if re.findall(g,tax):
					neighbors = nx.neighbors(G,n)
					dgroups[g].extend(neighbors)
		except IndexError:
			pass
		# tax = d["taxonomy"]
		# for g in groups:
		# 	if re.findall(g,tax):
		# 		dgroups[g].append(n)

	for k,v in dgroups.items():
		print("\t",k," : ",len(v)," (",round(len(v)/cptP*100,1),"%)",sep="")

	return dgroups

####

def print_proportions_eukaryotes(G,cptE,level):

	groups = list(set([ d["lineage"][level] for n,d in G.nodes(data=True) if (d["dataset"][0:2] == "18" or d["dataset"][0:2] == "CO")]))

	dgroups = defaultdict(list)

	for n,d in G.nodes(data=True):
		try:
			tax = d["taxonomy"]
			for g in groups:
				if re.findall(g,tax):
					neighbors = nx.neighbors(G,n)
					dgroups[g].extend(neighbors)
		except IndexError:
			pass

	for k,v in dgroups.items():
		print("\t",k," : ",len(v)," (",round(len(v)/cptE*100,1),"%)",sep="")

	return dgroups

#######

def print_connected_component(G, n):
	#print a node, its neighbors and the whole connected component associated#
	
	print("\n * Subnetwork :")
	nodes = nx.shortest_path(G,n).keys()
	sg = G.subgraph(nodes)
	simple_draw_network(sg, n)

	
 #######

def get_index_match_last(ranktax,tax):
	# lastax = [ t.split(";")[-2].split(" ")[0] for t in tax["taxonomy"].values] #deuxième split pour exemple Richelia
	lastax = [ t.split(";")[-2] for t in tax["taxonomy"].values] #deuxième split pour exemple Richelia
	
	found = 0
	indexmatchlastax = 0
	while not found and ranktax != "":
		try:
			# indexmatchlastax = lastax.index(ranktax.split(";")[-1])
			indexmatchlastax = int([i for i, s in enumerate(lastax) if ranktax.split(";")[-1].capitalize() == s][-1])
			found = 1
		except IndexError:
			ranktax=(";").join(ranktax.split(";")[:-1])

	return indexmatchlastax
 
 
def get_silva_lineage_from_taxid(taxid,tax):	
	taxo = tax.loc[taxid]["taxonomy"]
	rank = tax.loc[taxid]["rank"]
	if taxo.split(";")[0] =="Eukaryota":
		ranks = ["domain","major_clade","kingdom","phylum","class","order","family","genus"]
	elif taxo.split(";")[0] == "Bacteria" or taxo.split(";")[0] == "Archaea":
		ranks = ["domain","phylum","class","order","family","genus"]
	lineage = [0] * len(ranks)
	indexrank=0
	while True:
		try:
			indexrank = ranks.index(rank)
			lineage[indexrank] = taxid
		except ValueError:
			taxo = (";").join(taxo.split(";")[0:len(taxo.split(";"))-2])
			indexmatchlastax = get_index_match_last(taxo,tax)
			taxid = tax.iloc[indexmatchlastax].name
		break
	cutaxo = taxo 
	for i in range(len(cutaxo.split(";"))-1,1,-1):
		try:
			cutaxo = cutaxo.split(";")[0:len(cutaxo.split(";"))-2]
			cutaxo = (";").join(cutaxo)
			if cutaxo[-1] != ";":
				cutaxo+=";"
			matchrow = tax[tax["taxonomy"] == cutaxo]		   
			cutindexrank = ranks.index(matchrow["rank"].values[0])#on prend l'index du rank
			lineage[cutindexrank] = matchrow.index.values[0]   
		except ValueError:
			pass #rank non désiré 
	return(lineage)

###########


def get_silva_lineage_from_ranktax(ranktax,tax):
	indexmatchlastax = get_index_match_last(ranktax,tax)
	if indexmatchlastax == 0:
		return 0
	else:
		try:
			taxid = tax.iloc[indexmatchlastax].name
			lineage = get_silva_lineage_from_taxid(taxid,tax)
			return(lineage)
		except KeyError:
			return 0
			

def get_graph_lineages(G):
	for n,d in G.nodes(data=True):
		taxonomy = d["taxonomy"]
		taxonomy = (";").join([t.split(":")[1] if ":" in t else t for t in taxonomy.split("|")])
		taxonomy = taxonomy.replace("_"," ")
		lin = get_silva_lineage_from_ranktax(taxonomy,tax)
		if lin != 0 and lin[0] == 4:
			lineage = [lin[0]]
			lineage.extend(lin[3:])
			lin = lineage
		elif lin == 0 or lin ==  [0]*8:
			lin = [0,0,0,0,0,0]
		G.nodes[n]["lineage"] = str(lin)
		for i in range(1,len(lin)+1):
			try:
				d["lineage"+str(i)]=tax.loc[lin[i-1],"taxonomy"].split(";")[-2]
			except KeyError:
				d["lineage"+str(i)]="0"
		try:
			taxid = lin[[lin.index(i) for i in lin if i != 0][-1]]
			G.nodes[n]["taxid"] = taxid
			G.nodes[n]["silvataxidtaxonomy"] = tax.loc[taxid]["taxonomy"]
		except IndexError:
			G.nodes[n]["taxid"] = 0
			G.nodes[n]["silvataxidtaxonomy"] = "None"
		
	return(G)



####

def build_reference_globi_graph(level,globi):
	globi_array = globi[["lineage_org1","lineage_org2","interactionTypeName","referenceDoi","referenceCitation"]].to_numpy()
	
	for row in globi_array:
		first = ast.literal_eval(row[0])
		second = ast.literal_eval(row[1])
		if first != 0 and first[0]==4:
			lineage = [first[0]]
			lineage.extend(first[3:])
			row[0] = str(lineage)
		if second != 0 and second[0]==4:
			lineage = [second[0]]
			lineage.extend(second[3:])
			row[1] = str(lineage)
		row[0] = ast.literal_eval(row[0])[level]
		row[1] = ast.literal_eval(row[1])[level]

	globi_graph = nx.Graph()
	globi_graph.add_edges_from(globi_array[:,[0,1]])
	
	for e in globi_graph.edges():
		globi_graph.edges[e]["InteractionType"] = list()
		globi_graph.edges[e]["InteractionLink"] = list()
		globi_graph.edges[e]["InteractionRef"] = list()
	
	for row in globi_array:
		globi_graph.edges[row[0],row[1]]["InteractionType"].append(row[2])
		globi_graph.edges[row[0],row[1]]["InteractionLink"].append(row[3])
		globi_graph.edges[row[0],row[1]]["InteractionRef"].append(row[4])
		
	for e in globi_graph.edges():
		globi_graph.edges[e]["InteractionType"] = list(set(globi_graph.edges[e]["InteractionType"]))
		globi_graph.edges[e]["InteractionLink"] = list(set(globi_graph.edges[e]["InteractionLink"]))
		globi_graph.edges[e]["InteractionRef"] = list(set(globi_graph.edges[e]["InteractionRef"]))
		
	if 0 in globi_graph.nodes():
		globi_graph.remove_node(0)
		
	return(globi_graph)
			
####

def build_reference_pida_graph(level,pida):
	pida_array = pida[["lineage_org1","lineage_org2","Ecological interaction","Link","Reference"]].to_numpy()
	
	for row in pida_array:
		first = ast.literal_eval(row[0])
		second = ast.literal_eval(row[1])
		if first != 0 and first[0]==4:
			lineage = [first[0]]
			lineage.extend(first[3:])
			row[0] = str(lineage)
		if second != 0 and second[0]==4:
			lineage = [second[0]]
			lineage.extend(second[3:])
			row[1] = str(lineage)
		row[0] = ast.literal_eval(row[0])[level]
		row[1] = ast.literal_eval(row[1])[level]

	pida_graph = nx.Graph()
	pida_graph.add_edges_from(pida_array[:,[0,1]])
	
	for e in pida_graph.edges():
		pida_graph.edges[e]["InteractionType"] = list()
		pida_graph.edges[e]["InteractionLink"] = list()
		pida_graph.edges[e]["InteractionRef"] = list()
	
	for row in pida_array:
		pida_graph.edges[row[0],row[1]]["InteractionType"].append(row[2])
		pida_graph.edges[row[0],row[1]]["InteractionLink"].append(row[3])
		pida_graph.edges[row[0],row[1]]["InteractionRef"].append(row[4])
		
	for e in pida_graph.edges():
		pida_graph.edges[e]["InteractionType"] = list(set(pida_graph.edges[e]["InteractionType"]))
		pida_graph.edges[e]["InteractionLink"] = list(set(pida_graph.edges[e]["InteractionLink"]))
		pida_graph.edges[e]["InteractionRef"] = list(set(pida_graph.edges[e]["InteractionRef"]))
		
	if 0 in pida_graph.nodes():
		pida_graph.remove_node(0)
		
	return(pida_graph)

######


def cross_graphs_with_pida(G,ref,intere,level):
	mydict = defaultdict(list)
	for n in G.nodes():
		mydict[int(ast.literal_eval(str(G.nodes[n]["lineage"]))[level])].append(n)
		
	interG = nx.Graph()
	interG.add_edges_from(intere) 
	
	for n,d in interG.nodes(data=True):
		d["keys"] = (";").join(mydict[n])
		d["lineage"] = n

	for e in interG.edges():
		interG.edges[e]["InteractionType"] = ref.edges[e]["InteractionType"] #pour parasitism mettre dans bon ordre
		interG.edges[e]["InteractionLink"] = ref.edges[e]["InteractionLink"]
		interG.edges[e]["InteractionRef"] = ref.edges[e]["InteractionRef"]
		
	interG = nx.relabel_nodes(interG, lambda x: tax.loc[x]["taxonomy"]) #b pour bis
	
	#now annotate main graph
	for n1,n2,d in interG.edges(data=True):
		keys1 = interG.nodes[n1]["keys"].split(";")
		keys2 = interG.nodes[n2]["keys"].split(";")
		combi = list(itertools.product(keys1,keys2))
		edges = list(G.edges)
		common_edges = [comb for comb in combi if comb in edges]
		for e in common_edges:
			if G.edges[e]["pidalevel"] == -1:
				G.edges[e]["pidalevel"] = level
				G.edges[e]["pidaInteractionType"] =  d["InteractionType"]
				G.edges[e]["pidaInteractionLink"] =  d["InteractionLink"]
				G.edges[e]["pidaInteractionRef"] =  d["InteractionRef"]
	return G,interG #common_edges are edges to remove from Gcopy !!! to not annotate them twice
	

def cross_graphs_with_globi(G,ref,intere,level):
	mydict = defaultdict(list)
	for n in G.nodes():
		mydict[int(ast.literal_eval(str(G.nodes[n]["lineage"]))[level])].append(n)
		
	interG = nx.Graph()
	interG.add_edges_from(intere) 
	
	for n,d in interG.nodes(data=True):
		d["keys"] = (";").join(mydict[n])
		d["lineage"] = n

	for e in interG.edges():
		interG.edges[e]["InteractionType"] = ref.edges[e]["InteractionType"] #pour parasitism mettre dans bon ordre
		interG.edges[e]["InteractionLink"] = ref.edges[e]["InteractionLink"]
		interG.edges[e]["InteractionRef"] = ref.edges[e]["InteractionRef"]
		
	interG = nx.relabel_nodes(interG, lambda x: tax.loc[x]["taxonomy"]) #b pour bis
	
	#now annotate main graph
	for n1,n2,d in interG.edges(data=True):
		keys1 = interG.nodes[n1]["keys"].split(";")
		keys2 = interG.nodes[n2]["keys"].split(";")
		combi = list(itertools.product(keys1,keys2))
		edges = list(G.edges)
		common_edges = [comb for comb in combi if comb in edges]
		for e in common_edges:
			if G.edges[e]["globilevel"] == -1:
				G.edges[e]["globilevel"] = level
				G.edges[e]["globiInteractionType"] =  d["InteractionType"]
				G.edges[e]["globiInteractionLink"] =  d["InteractionLink"]
				G.edges[e]["globiInteractionRef"] =  d["InteractionRef"]
	return G,interG #common_edges are edges to remove from Gcopy !!! to not annotate them twice
			

def cross_pida_database(G):
	print("\n*** PIDA DATABASE***")
	pida = pd.read_csv("useful_tables/pidaAssociationsSILVA_V3.tsv",sep="\t",index_col=0)
	pida = pida[pida.lineage_org1 != "0"]
	pida = pida[pida.lineage_org2 != "0"]

	ranks = ["domain","phylum","class","order","family","genus"]
	interGraphs = list()
	for level in range(len(ranks)-1,0,-1):
		print("rank :",ranks[level])
		ref = build_reference_pida_graph(level,pida)
		Gb = nx.relabel_nodes(G, lambda x: int(ast.literal_eval(str(G.nodes[x]["lineage"]))[level])) #b pour bis
		intere = intersection_edges(ref,Gb)
		G,interGraph=cross_graphs_with_pida(G,ref,intere,level) 
		interGraphs.append(interGraph)
		print("   -> number interactions in common :",len(intere))
	return(G,interGraphs)

	
def cross_globi_database(G):

	print("\n*** GLOBI DATABASE***")
	globi = pd.read_csv("useful_tables/globiAssociationsSILVA_V2.tsv",sep="\t",index_col=0)
	globi = globi[globi.lineage_org1 != "0"]
	globi = globi[globi.lineage_org2 != "0"]
	globi = globi.dropna(subset=["lineage_org1"])
	globi = globi.dropna(subset=["lineage_org2"])
	ranks = ["domain","phylum","class","order","family","genus"]
	interGraphs = list()
	for level in range(len(ranks)-1,0,-1):
		print("rank :",ranks[level])
		ref = build_reference_globi_graph(level,globi)
		Gb = nx.relabel_nodes(G, lambda x: int(ast.literal_eval(str(G.nodes[x]["lineage"]))[level])) #b pour bis
		intere = intersection_edges(ref,Gb)
		G,interGraph=cross_graphs_with_globi(G,ref,intere,level) 
		interGraphs.append(interGraph)
		print("   -> number interactions in common :",len(intere))
	return(G,interGraphs)

tax = pd.read_csv("useful_tables/tax_slv_ssu_138.txt",
					 sep="\t",
					 header=None,
					 index_col=1,
					 names = ["taxonomy","rank","bug","jesaispas"])
tax = tax[["taxonomy","rank"]]

def make_cross_databases(G):
	print("\nCrossing databases ...")
	ranks = ["domain","phylum","class","order","family","genus"]
	# G = get_graph_lineages(G)
	nx.set_edge_attributes(G,-1,"pidalevel") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"pidaInteractionType") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"pidaInteractionLink") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"pidaInteractionRef") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,-1,"globilevel") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"globiInteractionType") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"globiInteractionLink") #level -1 par défaut tant que c'est pas trouvé
	nx.set_edge_attributes(G,None,"globiInteractionRef") #level -1 par défaut tant que c'est pas trouvé
	G,interGraphsPida = cross_pida_database(G)
	G,interGraphsGlobi = cross_globi_database(G)
	#annotate edges
	for n1,n2,d in G.edges(data=True):
		if d["pidalevel"] != -1 and d["globilevel"] != -1:
			d["color"] = "#526FF3"
			# d["label"] = "pida : "+ranks[d["pidalevel"]]+" - globi : "+ranks[d["globilevel"]]
			d["label"] = "pida : "+str(d["pidaInteractionType"])+"\nglobi : "+str(d["globiInteractionType"])
		elif d["pidalevel"] != -1:
			# d["color"] = "#FF44AC"
			d["label"] = "pida : "+str(d["pidaInteractionType"])
		elif d["globilevel"] != -1:
			# d["color"] = "#0D876E"
			d["label"] = "globi : "+str(d["globiInteractionType"])
	return(G)
	
	
def deeper_view_cross_databases(G):
	while True:
		try:
			print("1 : phylum")
			print("2 : class")
			print("3 : order")
			print("4 : family")
			print("5 : genre")
			level = int(input("Which level : "))
			indexrank = level -1
			rank = ranks[indexrank]
			if level == 0:
				break
			else:
				sg = nx.Graph()
				for n1,n2,d in G.edges(data=True):
					if d["pidalevel"] == level or d["globilevel"] == level :
						sg.add_edge(n1,n2)
						sg.edges[n1,n2].update(d)

				for n,d in sg.nodes(data=True):
					sg.nodes[n].update(G.nodes[n])
				print_list_nodes(G, subG = sg,taxonomy=True)
				simple_draw_network(sg)
				plt.show()
				return
				
		except IndexError:
			print("Index error has occurred, please retry")
	
	 	
