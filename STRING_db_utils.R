require (data.table)

#' @param identifiers uniprot, gene names, etc. Anything found in the string aliases file. character
#' @stringAliasFile  path to file with aliases. Download from string at: 
#' @returns A data table with two columns `string` and `alias`. Alias is the identifiers

GetStringIDMapping <- function (identifiers, 
                                stringAliasFile = "/Users/ben/Downloads/9606.protein.aliases.v11.5.txt.gz"){
  aliases <- fread (stringAliasFile)
  message ("There are ", sum (identifiers %in% aliases$alias), " proteins found in string, and ",
           sum (!identifiers %in% aliases$alias), " not in string")
  # attempt to get only the best alias mapping per identifier
  # loop through the sources in the alias table, favoring the source that identifies the most aliases first
  # after selecting each source, remove the identified alias from the list of aliases we need to identify.
  # take the next best source based on number of proteins it contains an alias for
  # continue until no sources are left for any alias.
  goodAliases <- data.table(string = character(), alias = character())
  while(TRUE)
  {
    topSource <- aliases[!alias %in% goodAliases$alias & alias %in% identifiers, .(N = length(unique(alias))), by = source][which.max(N), source]
    if(length(topSource) == 0){
      break
    }
    newAliases <- aliases[source == topSource & !alias %in% goodAliases$alias & alias %in% identifiers, .(string = `#string_protein_id`, alias)]
    numAliases <- nrow(newAliases)
    newAliases <- newAliases[, .SD[1], by = alias]
    if (nrow(newAliases) != numAliases)
      message ("1 alias to multiple string IDs found. Taking the first found in the alias table")
    
    goodAliases <- rbind (newAliases,
                          goodAliases)
    
    message("Using aliases from ", topSource, " and now have ", nrow(goodAliases), " aliases mapped")
  }
  
  #stringsOI <- unique(goodAliases$string)
  goodAliases

  }


#' @param ids uniprot, gene names, etc. Anything found in the string aliases file. character
#' @stringAliasFile  path to file with aliases. Download from string at: 
#' @returns a vector same length as `ids`, but replaced with stringIDs where possible, or NA when not
GetStringIDMapping.inOrder <- function(ids, stringAliasFile = "/Users/ben/Downloads/9606.protein.aliases.v11.5.txt.gz" ){
  GetStringIDMapping(ids, stringAliasFile)[ids, string, on = "alias"]
}




require (igraph)

# namedThresholds = c(experimental = 900)

subsetSTRINGasGraph <- function(stringIDs, stringFile = "/Users/ben/Downloads/9606.protein.links.detailed.v11.5.txt.gz", threshold = 600, namedThresholds = c()){
  
  string <- fread (stringFile)
  
  if (length(stringIDs) == 1 & is.numeric(stringIDs)){
    stringIDs <- sample(unique(c(string$protein1, string$protein2)), size = stringIDs)
  }
  
  if ("data.frame" %in% class(stringIDs)){
    vTable <- as.data.table(stringIDs)
    setnames(vTable, old = colnames(vTable)[1], new = "name")
    stringIDs <- vTable$name
  } else{
    vTable <- data.table(name = stringIDs)
  }
  
  
  edges <- string[combined_score > threshold & protein1 %in% stringIDs & protein2 %in% stringIDs]
  
  for (name in names(namedThresholds)){
    edges <- edges[ edges[[name]] >= namedThresholds[name],]
      
  }
  
  # de-duplicate
  edges[protein1 > protein2, c("protein1", "protein2") := .(protein2, protein1)]
  edges <- edges[, .SD[1], by = .(protein1, protein2)]
  
  
  ss <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vTable)
  
  return (ss)
}

plotSTRINGgraph <- function (stringIDs, stringFile = NULL, threshold = 600){
  if ('igraph' %in% class(stringIDs)){
    stringGraph <- stringIDs
  }else{
    stringGraph <- subsetSTRINGasGraph(stringIDs, stringFile, threshold = threshold)
  }
  
  com <- igraph::walktrap.community(stringGraph)
  V(stringGraph)$walkTrap <- com$membership
  
  p <- ggraph(stringGraph, layout = "fr")+
    geom_edge_fan(color = "gray") +
    geom_node_point(aes (color = as.factor(walkTrap)), show.legend = FALSE)  +
    geom_node_text(aes(label = gene, color = as.factor(walkTrap)), repel = TRUE, max.overlaps = 15, show.legend = FALSE) +
    theme(panel.background = NULL)
  
  return(p)
}




# find one hop edges in a network (not string specific)

findConnectorNodeEdges <- function (edgeTable, stringsOI, minScore = 600){
  
  # subset to get all edges from one of the stringsOI nodes
  allEdges <- edgeTable[combined_score > minScore & 
                          (protein1 %in% stringsOI |  # include single edges to outside group
                             protein2 %in% stringsOI )]
  
  # any node in the above edges that is not in stringsOI is a candidate singleHop node
  candidateSteppingNodes <- setdiff(unique(c(allEdges$protein1, allEdges$protein2)), stringsOI)
  
  
  # just the edges out from (or in to) 
  outEdges <- rbind(
    allEdges[candidateSteppingNodes,,on = "protein1"][, cand := protein1][protein1 > protein2, c("protein1", "protein2") := .(protein2, protein1)][, .(protein1, protein2, cand, weight = combined_score/1000)],
    allEdges[candidateSteppingNodes,, on = "protein2"][, cand := protein2][protein1 > protein2, c("protein1", "protein2") := .(protein2, protein1)][, .(protein1, protein2, cand, weight = combined_score/1000)]
  )
  
  # collapse over protein1, protein2, cand, taking best weight if there are multiple different
  outEdges <- outEdges[, .(weight = max(weight)), by = .(protein1, protein2, cand)]
  goodCandidates <- outEdges[, .N, by = cand][N > 1,cand]
  outEdges[cand %in% goodCandidates, .(protein1, protein2, weight)]
}

# subnetwork string to a table, include one hop edges


GetStringSubNetwork <- function (stringsOI,
                                 stringFile = "/Users/ben/Downloads/9606.protein.links.detailed.v11.5.txt.gz",
                                 minString = 600,
                                 oneHopConnections = FALSE){
  
  if ("data.table" %in% class(stringFile)){
    string <- stringFile    
  } else{
    string <- fread (stringFile)
  }
  
  
  if (length(stringsOI) == 1 & is.numeric(stringsOI)){
    stringsOI <- sample(unique(c(string$protein1, string$protein2)), size = stringsOI)
    print (sort(stringsOI))
  }
  
  # edges between stringsOI
  sigNetworkSet <- string[combined_score > minString & 
                            (protein1 %in% stringsOI &  
                               protein2 %in% stringsOI ),
                          .(protein1, protein2, weight  = combined_score/1000)]
  
  if (oneHopConnections == TRUE){
    # edges to 1-hop connectors
    singleHopEdges <- findConnectorNodeEdges (string, stringsOI, minScore = minString)
    return (rbind(sigNetworkSet, singleHopEdges))
  }
  return (sigNetworkSet)
}