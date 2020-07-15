#!/usr/bin/Rscript

library(igraph)
library(data.table)
library(magrittr)

args <- commandArgs(T)
graph = args[1]
output = args[2]


####################################

graph <- read_graph(graph,format="graphml")
graph_table <- as_long_data_frame(graph) %>% data.table

colnames(graph_table)

weights = graph_table[,weight]
pida = graph_table[,.SD,.SDcols=grep("^pida",colnames(graph_table))]
globi = graph_table[,.SD,.SDcols=grep("^globi",colnames(graph_table))]

tmp1 <- graph_table[,.SD,.SDcols=grep("^from",colnames(graph_table))]
setnames(tmp1,"from","graphid")
setnames(tmp1,colnames(tmp1),sub("^from_","",colnames(tmp1)))

tmp2 <- graph_table[,.SD,.SDcols=grep("^to",colnames(graph_table))]
setnames(tmp2,"to","graphid")
setnames(tmp2,colnames(tmp2),sub("^to_","",colnames(tmp2)))


node_data <- rbind(tmp1,tmp2) %>% unique
colnames(node_data)

#manque le rank, reftaxonomy, seasons???, numero du subgraph

filter.node_data = node_data[,c("keystoneindexrank","key","dataset","degrees","degreescentrality","closenesscentrality","betweennesscentrality","subgraphcentrality",
                         "abundance","occurrence","reftaxonomy","idtaxonomy","taxonomy","optimatemp","amplicon","sequence","lineage","taxid","silvataxidtaxonomy")]

sorted.filter.node_data = filter.node_data[order(node_data[,keystoneindexrank]),]

write.table(sorted.filter.node_data,paste(output,"_nodes.tsv",sep=""),sep="\t",row.names=F,quote=F)

    ####edges

setnames(tmp1,colnames(tmp1),paste(sub("^to_","",colnames(tmp1)),"_1",sep=""))
setnames(tmp2,colnames(tmp2),paste(sub("^to_","",colnames(tmp2)),"_2",sep=""))

edges_data = data.table(cbind(tmp1,tmp2,weight=weights))
edges_data = cbind(edges_data,pida)
edges_data = cbind(edges_data,globi)

colnames(edges_data)

key1 = paste(c("keystoneindexrank","key","dataset","idtaxonomy","taxonomy","optimatemp","amplicon","lineage","taxid","silvataxidtaxonomy"),"_1",sep="")
key2 = paste(c("keystoneindexrank","key","dataset","idtaxonomy","taxonomy","optimatemp","amplicon","lineage","taxid","silvataxidtaxonomy"),"_2",sep="")
# filter.edge_data=edges_data[,c("weight",
#                                paste(c("keystoneindexrank","label","dataset","idtaxonomy","taxonomy","optimatemp","amplicon"),"_1",sep=""),
#                                paste(c("keystoneindexrank","label","dataset","idtaxonomy","taxonomy","optimatemp","amplicon"),"_2",sep=""),
#                                colnames(pida),
#                                colnames(globi)
#                                )]
vec = c("weight",sapply(seq_along(key1), function(i) append(key1[i], key2[i], i)),colnames(pida),colnames(globi))
filter.edge_data = data.frame(edges_data)[,vec]
colnames(filter.edge_data)
# filter.edge_data = 

# ranks.keyindex = rank(rowMeans(filter.edge_data[,c(1,2)]))
# ranks.weight = rev(order(filter.edge_data[,c(13)]))
# 
# head(ranks.keyindex)
# head(ranks.weight)

#sor nodes with first rank node on the left
rev.vec = c("weight",sapply(seq_along(key2), function(i) append(key2[i], key1[i], i)),colnames(pida),colnames(globi))

filter.edge_data[apply(filter.edge_data[,c(2,3)],1,is.unsorted),][,vec]=filter.edge_data[apply(filter.edge_data[,c(2,3)],1,is.unsorted),][,rev.vec]

ranks = order(filter.edge_data[,"keystoneindexrank_1"],
              filter.edge_data[,"weight"],decreasing=TRUE)
tmp = filter.edge_data[ranks,]
sorted.filter.edge_data = tmp[order(tmp[,"keystoneindexrank_1"]),]

head(sorted.filter.edge_data,5)


write.table(sorted.filter.edge_data,paste(output,"_edges.tsv",sep=""),sep="\t",row.names=F,quote=F)
