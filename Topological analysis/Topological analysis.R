#### Topological analysis

## Start clean
rm(list = ls())

## Make sure the current working directory is right
getwd()

## Load functions
source("fun.R")
n <- 10

## Download and install the package
#install.packages("igraph")

## Load package
library(igraph)

## ---- TOB1 ----
# lcc1 TOB1 48h
edges_lcc1_TOB1_48h <- read.csv("Network/eLCC1_TOB1_48h.csv", stringsAsFactors=F)
lcc1_TOB1_48h <- graph.data.frame(edges_lcc1_TOB1_48h, directed=F) # note to assign undirected graph
topn_lcc1_TOB1_48h <- topn_list(edges_lcc1_TOB1_48h, n)
topn_lcc1_TOB1_48h

# lcc9 TOB1 48h
edges_lcc9_TOB1_48h <- read.csv("Network/eLCC9_TOB1_48h.csv", stringsAsFactors=F)
lcc9_TOB1_48h <- graph.data.frame(edges_lcc9_TOB1_48h, directed=F)
topn_lcc9_TOB1_48h <- topn_list(edges_lcc9_TOB1_48h, n)
topn_lcc9_TOB1_48h 

# mcf7 TOB1 48h
edges_mcf7_TOB1_48h <- read.csv("Network/eMCF7_TOB1_48h.csv", stringsAsFactors=F)
mcf7_TOB1_48h <- graph.data.frame(edges_mcf7_TOB1_48h, directed=F)
topn_mcf7_TOB1_48h <- topn_list(edges_mcf7_TOB1_48h, n)
topn_mcf7_TOB1_48h 

# lcc1 TOB1 96h
edges_lcc1_TOB1_96h <- read.csv("Network/eLCC1_TOB1_96h.csv", stringsAsFactors=F)
lcc1_TOB1_96h <- graph.data.frame(edges_lcc1_TOB1_96h, directed=F) 
topn_lcc1_TOB1_96h <- topn_list(edges_lcc1_TOB1_96h, n)
topn_lcc1_TOB1_96h

# lcc9 TOB1 96h
edges_lcc9_TOB1_96h <- read.csv("Network/eLCC9_TOB1_96h.csv", stringsAsFactors=F)
lcc9_TOB1_96h <- graph.data.frame(edges_lcc9_TOB1_96h, directed=F) 
topn_lcc9_TOB1_96h <- topn_list(edges_lcc9_TOB1_96h, n)
topn_lcc9_TOB1_96h

# mcf7 TOB1 96h
edges_mcf7_TOB1_96h <- read.csv("Network/eMCF7_TOB1_96h.csv", stringsAsFactors=F)
mcf7_TOB1_96h <- graph.data.frame(edges_mcf7_TOB1_96h, directed=F)
topn_mcf7_TOB1_96h <- topn_list(edges_mcf7_TOB1_96h, n)
topn_mcf7_TOB1_96h 

# lcc1 TOB1 144h
edges_lcc1_TOB1_144h <- read.csv("Network/eLCC1_TOB1_144h.csv", stringsAsFactors=F)
lcc1_TOB1_144h <- graph.data.frame(edges_lcc1_TOB1_144h, directed=F) 
topn_lcc1_TOB1_144h <- topn_list(edges_lcc1_TOB1_144h, n)
topn_lcc1_TOB1_144h

# lcc9 TOB1 144h
edges_lcc9_TOB1_144h <- read.csv("Network/eLCC9_TOB1_144h.csv", stringsAsFactors=F)
lcc9_TOB1_144h <- graph.data.frame(edges_lcc9_TOB1_144h, directed=F) 
topn_lcc9_TOB1_144h <- topn_list(edges_lcc9_TOB1_144h, n)
topn_lcc9_TOB1_144h

# mcf7 TOB1 144h
edges_mcf7_TOB1_144h <- read.csv("Network/eMCF7_TOB1_144h.csv", stringsAsFactors=F)
mcf7_TOB1_144h <- graph.data.frame(edges_mcf7_TOB1_144h, directed=F)
topn_mcf7_TOB1_144h <- topn_list(edges_mcf7_TOB1_144h, n)
topn_mcf7_TOB1_144h 

## ---- POLR2B ----
# lcc1 POLR2B 48h
edges_lcc1_POLR2B_48h <- read.csv("Network/eLCC1_POLR2B_48h.csv", stringsAsFactors=F)
lcc1_POLR2B_48h <- graph.data.frame(edges_lcc1_POLR2B_48h, directed=F) # note to assign undirected graph
topn_lcc1_POLR2B_48h <- topn_list(edges_lcc1_POLR2B_48h, n)
topn_lcc1_POLR2B_48h

# lcc1 POLR2B 96h
edges_lcc1_POLR2B_96h <- read.csv("Network/eLCC1_POLR2B_96h.csv", stringsAsFactors=F)
lcc1_POLR2B_96h <- graph.data.frame(edges_lcc1_POLR2B_96h, directed=F) # note to assign undirected graph
topn_lcc1_POLR2B_96h <- topn_list(edges_lcc1_POLR2B_96h, n)
topn_lcc1_POLR2B_96h

# lcc1 POLR2B 144h
edges_lcc1_POLR2B_144h <- read.csv("Network/eLCC1_POLR2B_144h.csv", stringsAsFactors=F)
lcc1_POLR2B_144h <- graph.data.frame(edges_lcc1_POLR2B_144h, directed=F) # note to assign undirected graph
topn_lcc1_POLR2B_144h <- topn_list(edges_lcc1_POLR2B_144h, n)
topn_lcc1_POLR2B_144h

# lcc9 POLR2B 48h
edges_lcc9_POLR2B_48h <- read.csv("Network/eLCC9_POLR2B_48h.csv", stringsAsFactors=F)
lcc9_POLR2B_48h <- graph.data.frame(edges_lcc9_POLR2B_48h, directed=F) # note to assign undirected graph
topn_lcc9_POLR2B_48h <- topn_list(edges_lcc9_POLR2B_48h, n)
topn_lcc9_POLR2B_48h

# lcc9 POLR2B 96h
edges_lcc9_POLR2B_96h <- read.csv("Network/eLCC9_POLR2B_96h.csv", stringsAsFactors=F)
lcc9_POLR2B_96h <- graph.data.frame(edges_lcc9_POLR2B_96h, directed=F) # note to assign undirected graph
topn_lcc9_POLR2B_96h <- topn_list(edges_lcc9_POLR2B_96h, n)
topn_lcc9_POLR2B_96h

# lcc9 POLR2B 144h
edges_lcc9_POLR2B_144h <- read.csv("Network/eLCC9_POLR2B_144h.csv", stringsAsFactors=F)
lcc9_POLR2B_144h <- graph.data.frame(edges_lcc9_POLR2B_144h, directed=F) # note to assign undirected graph
topn_lcc9_POLR2B_144h <- topn_list(edges_lcc9_POLR2B_144h, n)
topn_lcc9_POLR2B_144h

# mcf7 POLR2B 48h
edges_mcf7_POLR2B_48h <- read.csv("Network/eMCF7_POLR2B_48h.csv", stringsAsFactors=F)
mcf7_POLR2B_48h <- graph.data.frame(edges_mcf7_POLR2B_48h, directed=F) # note to assign undirected graph
topn_mcf7_POLR2B_48h <- topn_list(edges_mcf7_POLR2B_48h, n)
topn_mcf7_POLR2B_48h

# mcf7 POLR2B 96h
edges_mcf7_POLR2B_96h <- read.csv("Network/eMCF7_POLR2B_96h.csv", stringsAsFactors=F)
mcf7_POLR2B_96h <- graph.data.frame(edges_mcf7_POLR2B_96h, directed=F) # note to assign undirected graph
topn_mcf7_POLR2B_96h <- topn_list(edges_mcf7_POLR2B_96h, n)
topn_mcf7_POLR2B_96h

# mcf7 POLR2B 144h
edges_mcf7_POLR2B_144h <- read.csv("Network/eMCF7_POLR2B_144h.csv", stringsAsFactors=F)
mcf7_POLR2B_144h <- graph.data.frame(edges_mcf7_POLR2B_144h, directed=F) # note to assign undirected graph
topn_mcf7_POLR2B_144h <- topn_list(edges_mcf7_POLR2B_144h, n)
topn_mcf7_POLR2B_144h

## ---- CYR61 ----
# lcc1 CYR61 48h
edges_lcc1_CYR61_48h <- read.csv("Network/eLCC1_CYR61_48h.csv", stringsAsFactors=F)
lcc1_CYR61_48h <- graph.data.frame(edges_lcc1_CYR61_48h, directed=F) 
topn_lcc1_CYR61_48h <- topn_list(edges_lcc1_CYR61_48h, n)
topn_lcc1_CYR61_48h

# lcc1 CYR61 96h
edges_lcc1_CYR61_96h <- read.csv("Network/eLCC1_CYR61_96h.csv", stringsAsFactors=F)
lcc1_CYR61_96h <- graph.data.frame(edges_lcc1_CYR61_96h, directed=F) 
topn_lcc1_CYR61_96h <- topn_list(edges_lcc1_CYR61_96h, n)
topn_lcc1_CYR61_96h

# lcc1 CYR61 144h
edges_lcc1_CYR61_144h <- read.csv("Network/eLCC1_CYR61_144h.csv", stringsAsFactors=F)
lcc1_CYR61_144h <- graph.data.frame(edges_lcc1_CYR61_144h, directed=F) 
topn_lcc1_CYR61_144h <- topn_list(edges_lcc1_CYR61_144h, n)
topn_lcc1_CYR61_144h

# lcc9 CYR61 48h
edges_lcc9_CYR61_48h <- read.csv("Network/elcc9_CYR61_48h.csv", stringsAsFactors=F)
lcc9_CYR61_48h <- graph.data.frame(edges_lcc9_CYR61_48h, directed=F) 
topn_lcc9_CYR61_48h <- topn_list(edges_lcc9_CYR61_48h, n)
topn_lcc9_CYR61_48h

# lcc9 CYR61 96h
edges_lcc9_CYR61_96h <- read.csv("Network/elcc9_CYR61_96h.csv", stringsAsFactors=F)
lcc9_CYR61_96h <- graph.data.frame(edges_lcc9_CYR61_96h, directed=F) 
topn_lcc9_CYR61_96h <- topn_list(edges_lcc9_CYR61_96h, n)
topn_lcc9_CYR61_96h

# lcc9 CYR61 144h
edges_lcc9_CYR61_144h <- read.csv("Network/elcc9_CYR61_144h.csv", stringsAsFactors=F)
lcc9_CYR61_144h <- graph.data.frame(edges_lcc9_CYR61_144h, directed=F) 
topn_lcc9_CYR61_144h <- topn_list(edges_lcc9_CYR61_144h, n)
topn_lcc9_CYR61_144h

# mcf7 CYR61 48h
edges_mcf7_CYR61_48h <- read.csv("Network/eMCF7_CYR61_48h.csv", stringsAsFactors=F)
mcf7_CYR61_48h <- graph.data.frame(edges_mcf7_CYR61_48h, directed=F) 
topn_mcf7_CYR61_48h <- topn_list(edges_mcf7_CYR61_48h, n)
topn_mcf7_CYR61_48h

# mcf7 CYR61 96h
edges_mcf7_CYR61_96h <- read.csv("Network/eMCF7_CYR61_96h.csv", stringsAsFactors=F)
mcf7_CYR61_96h <- graph.data.frame(edges_mcf7_CYR61_96h, directed=F) 
topn_mcf7_CYR61_96h <- topn_list(edges_mcf7_CYR61_96h, n)
topn_mcf7_CYR61_96h

# mcf7 CYR61 144h
edges_mcf7_CYR61_144h <- read.csv("Network/eMCF7_CYR61_144h.csv", stringsAsFactors=F)
mcf7_CYR61_144h <- graph.data.frame(edges_mcf7_CYR61_144h, directed=F) 
topn_mcf7_CYR61_144h <- topn_list(edges_mcf7_CYR61_144h, n)
topn_mcf7_CYR61_144h

## ---- PSMC5 ----
# lcc1 PSMC5 48h
edges_lcc1_PSMC5_48h <- read.csv("Network/eLCC1_PSMC5_48h.csv", stringsAsFactors=F)
lcc1_PSMC5_48h <- graph.data.frame(edges_lcc1_PSMC5_48h, directed=F) 
topn_lcc1_PSMC5_48h <- topn_list(edges_lcc1_PSMC5_48h, n)
topn_lcc1_PSMC5_48h

# lcc1 PSMC5 96h
edges_lcc1_PSMC5_96h <- read.csv("Network/eLCC1_PSMC5_96h.csv", stringsAsFactors=F)
lcc1_PSMC5_96h <- graph.data.frame(edges_lcc1_PSMC5_96h, directed=F) 
topn_lcc1_PSMC5_96h <- topn_list(edges_lcc1_PSMC5_96h, n)
topn_lcc1_PSMC5_96h

# lcc1 PSMC5 144h
edges_lcc1_PSMC5_144h <- read.csv("Network/eLCC1_PSMC5_144h.csv", stringsAsFactors=F)
lcc1_PSMC5_144h <- graph.data.frame(edges_lcc1_PSMC5_144h, directed=F) 
topn_lcc1_PSMC5_144h <- topn_list(edges_lcc1_PSMC5_144h, n)
topn_lcc1_PSMC5_144h

# lcc9 PSMC5 48h
edges_lcc9_PSMC5_48h <- read.csv("Network/elcc9_PSMC5_48h.csv", stringsAsFactors=F)
lcc9_PSMC5_48h <- graph.data.frame(edges_lcc9_PSMC5_48h, directed=F) 
topn_lcc9_PSMC5_48h <- topn_list(edges_lcc9_PSMC5_48h, n)
topn_lcc9_PSMC5_48h

# lcc9 PSMC5 96h
edges_lcc9_PSMC5_96h <- read.csv("Network/elcc9_PSMC5_96h.csv", stringsAsFactors=F)
lcc9_PSMC5_96h <- graph.data.frame(edges_lcc9_PSMC5_96h, directed=F) 
topn_lcc9_PSMC5_96h <- topn_list(edges_lcc9_PSMC5_96h, n)
topn_lcc9_PSMC5_96h

# lcc9 PSMC5 144h
# no network available

# mcf7 PSMC5 48h
edges_mcf7_PSMC5_48h <- read.csv("Network/eMCF7_PSMC5_48h.csv", stringsAsFactors=F)
mcf7_PSMC5_48h <- graph.data.frame(edges_mcf7_PSMC5_48h, directed=F) 
topn_mcf7_PSMC5_48h <- topn_list(edges_mcf7_PSMC5_48h, n)
topn_mcf7_PSMC5_48h

# mcf7 PSMC5 96h
edges_mcf7_PSMC5_96h <- read.csv("Network/eMCF7_PSMC5_96h.csv", stringsAsFactors=F)
mcf7_PSMC5_96h <- graph.data.frame(edges_mcf7_PSMC5_96h, directed=F) 
topn_mcf7_PSMC5_96h <- topn_list(edges_mcf7_PSMC5_96h, n)
topn_mcf7_PSMC5_96h

# mcf7 PSMC5 144h
edges_mcf7_PSMC5_144h <- read.csv("Network/eMCF7_PSMC5_144h.csv", stringsAsFactors=F)
mcf7_PSMC5_144h <- graph.data.frame(edges_mcf7_PSMC5_144h, directed=F) 
topn_mcf7_PSMC5_144h <- topn_list(edges_mcf7_PSMC5_144h, n)
topn_mcf7_PSMC5_144h

## Take a closer look using mcf7_TOB1_48h as an example
# node names
graph <- mcf7_TOB1_48h
V(graph)$name 
# node degrees
degree(graph) 
sort(degree(graph), decreasing=T)
# betweenness
round(betweenness(graph), 2)
sort(betweenness(graph), decreasing=T)
# closeness
round(closeness(graph), 4)
sort(closeness(graph), decreasing=T)
# eigenvector centrality
round(eigen_centrality(graph)$vector, 2)
sort(eigen_centrality(graph)$vector, decreasing=T)
# page rank
round(page_rank(graph)$vector, 3)
sort(page_rank(graph)$vector, decreasing=T)
# correlation map
cent_graph <- list(`Degree`=degree(graph),
                   `Betweenness`=betweenness(graph),
                   `Closeness`=closeness(graph),
                   `Eigenvector`=eigen_centrality(graph)$vector,
                   `PageRank`=page_rank(graph)$vector)
pairs(cent_graph, lower.panel=function(x,y) {
    usr <- par("usr")
    text(mean(usr[1:2]), mean(usr[3:4]), round(cor(x,y), 3), cex=2, col="blue")
} )

