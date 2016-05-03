createGODAG <-
function(sigNodes, ontology = "BP") {

  nodeLabel <- new.env(hash = T, parent = emptyenv())
 
  isNodeInDAG <- function(node) {
    return(exists(node, envir = nodeLabel, mode = 'logical', inherits = FALSE))
  }
  setNodeInDAG <- function(node) {
    assign(node, TRUE, envir = nodeLabel)
  }

  GOParents <- get(paste('GO', ontology, 'PARENTS', sep = ''))

  ROOT <- as.character(revmap(GOParents)$all)

  adjList <- as.list(GOParents)  
  edgeEnv <- new.env(hash = T, parent = emptyenv())  
  envAddEdge <- function(u, v, type) {
    assign(v, switch(type, is_a = 0, part_of = 1, -1), envir = get(u, envir = edgeEnv))
  }
  
  createNodesRelationship <- function(node) {
    if(isNodeInDAG(node))
      return(1)
    
    ## we put the node in the graph and we get his parents
    setNodeInDAG(node)    # we visit the node
    assign(node, new.env(hash = T, parent = emptyenv()), envir = edgeEnv) # adj list
  
    if(node == ROOT) 
      return(2)

    adjNodes <- adjList[[node]]

    if(length(adjNodes) == 0)
      cat('\n There are no adj nodes for node: ', node, '\n')
          
    for(i in 1:length(adjNodes)) {
      x <- as.character(adjNodes[i])
      envAddEdge(node, x, names(adjNodes[i]))
      createNodesRelationship(x)
    }

    return(0)
  }

  ## we start from the most specific nodes
  lapply(sigNodes, createNodesRelationship)
  
  .graphNodes <- ls(edgeEnv)
  .edgeList <- eapply(edgeEnv,
                      function(adjEnv) {
                        aux <- as.list(adjEnv)
                        return(list(edges = match(names(aux), .graphNodes),
                                    weights = as.numeric(aux)))
                      })
  
  ## now we can build the graphNEL object
  GOgraph.topo <- new('graphNEL',
                      nodes = .graphNodes,
                      edgeL = .edgeList,
                      edgemode = 'directed')
  
  require("SparseM") || stop("package SparseM is required")
  
  GOgraph.topo <- sparseM2Graph(t(graph2SparseM(GOgraph.topo, TRUE)),
                        .graphNodes, edgemode = "directed")
  return(GOgraph.topo)
}

