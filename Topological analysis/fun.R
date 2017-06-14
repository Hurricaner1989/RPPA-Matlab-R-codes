topn_list <- function(data, n){
    graph <- graph.data.frame(data, directed=F) # note to assign undirected graph
    topn_list <- list('Node'=V(graph)$name, 
                      'Degree'=sort(degree(graph), decreasing=T)[1:n],
                      'Betweenness'=sort(betweenness(graph), decreasing=T)[1:n],
                      'Closeness'=sort(closeness(graph), decreasing=T)[1:n],
                      'Eigenvector'=sort(eigen_centrality(graph)$vector, decreasing=T)[1:n],
                      'PageRank'=sort(page_rank(graph)$vector, decreasing=T)[1:n])
}



