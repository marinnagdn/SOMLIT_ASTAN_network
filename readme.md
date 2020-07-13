## Preparation

0) Download jupyter_notebook_abundance_metaniches.zip and uncompress it
1) Download and install Anaconda for Windows (https://docs.anaconda.com/anaconda/install/windows/)
2) Open the Anaconda Navigator
3) Launch Jupyter Notebook app in your web navigator
4) Nagigate in your folders until you are in jupyter_notebook_abundance_metaniches

## Abundances and niches 

1) Open the file ploting_abundance_and_metaniches.ipynb
=> You are on the Jupyter Notebook
I didn't succeed to hide the code so far. 
2) Input data :
	- excel_path : A csv file such as example.csv.
		(One node is to find a node and all its neighbors, using the id of the node or a given taxonomic keyword such as Rhizaria)
		
		OR
		
	- fasta_path : A fasta file such as sequences_mariarita.fasta or sequences_nathalie.fasta 
		Using this option you will have a table with the following informations : label of the matched node, the id fasta, the taxonomy, and the reference. 
	
	(When you use one of the option, you have to let only "" with nothing inside to the unused variable. 
	
	- graph_path : The graph path. So far you have the choice between :
		- the global graph of FlashWeave "flashweaveclr+1.graphml" based on CLR + Pseudocount 1 abundances
		- the global graph of fLSA "flsaclr+1.graphml" based on CLR + Pseudocount 1 abundances
		
	- output_path : path of output at the pdf format
	
#### Example :
		graph_path = "flsaclr+1.graphml" #the path of the graph
		fasta_path = "sequences_mariarita.fasta" #fasta file. if no just let the "" empty
		excel_path = "" #excel file. if no just let the "" empty
		output_path = "outputs/mariarita_lsa.pdf" #path of your output with the extension.pd
	
	
3) When all is ready, just click on the little double triangle button at the top of the page to execute at once all the cells.
4) Once it's computed, you can check the results in the pdf file. 

## Interaction Heatmaps

1) Open the file in interaction_heatmaps.ipynb
2) Input Data :
	- graph_path : enter the path of the graph
	- output_path : the path of the output at the pdf format
	
#### Example :
		graph_path = "flsaclr+1.graphml" #the path of the graph
		output_path = "outputs/test.pdf" #path of your output with the extension.pdf
	
3) Keywords
	- kw1 : first keyword 
	- kw2 : second keyword
	
You can either give two keywords or one keyword. 

#### Example : 
		kw1="Rhizaria"
		kw2=""
	
	=> It will generate a heatmap with Rhizaria versus all other organisms, including Rhizaria
		
	OR

		kw1="Cercozoa"
		kw2="Mollusca"
		
	=> It will generate a heatmap with Cercozoa versus Mollusca
	
4) Taxonomic levels.

	You can submit different taxonomic levels depending your request.
	
	level1 corresponds to kw1, and level2 to kw2.
	
	1 - Phylum
	2 - Class
	3 - Order
	4 - Family
	5 - Genus
	
#### Example :

		level1 = 1 #Phylum
		level2 = 5 #Genus 

If an error occurs during the execution, maybe try with a lower taxonomic level. I didn't check each situation yet. 
		


