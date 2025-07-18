utils::globalVariables(c("Source"))

#' @title get_network_matches
#' @description This function calculates network matches for a given graph and gene list using a transcription network and an aggregation method.
#' @param ENSG_Genelist A character vector of gene identifiers.
#' @param organism A character string specifying the organism. Options are 'human' or 'mouse'. Default is 'human'.
#' @param aggregator A character string specifying the aggregation method. Options are 'mean' or 'degree_importance'. Default is 'mean'.
#' @param network_data A data frame containing the transcription network data. If NULL, the function will download the data from a predefined URL.The columns of the data frame should be 'source', 'target', 'edge_type', and 'evidence'.
#' @return A data frame with columns 'ENSG', 'repression_matches', and 'activation_matches' containing the vertex names and their respective feature values.
#' @export
get_network_matches = function(ENSG_Genelist, organism, aggregator = 'mean', network_data = NULL){

  #Loading the transcription network data from web link https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
  # and converting it to an igraph object

  if (!organism %in% c("human", "mouse")) {
    stop("Organism not supported")
  }

  if(is.null(network_data)){
    human_path = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
    mouse_path = "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv"

    if (organism == 'human') {
      trascription_network <- utils::read.table(human_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else if (organism == 'mouse') {
      trascription_network <- utils::read.table(mouse_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }

    colnames(trascription_network) = c('source','target','edge_type','evidence')
  }else{
    trascription_network = network_data
  }

  g <- igraph::graph_from_data_frame(trascription_network[,1:2], directed = TRUE)
  X = data.frame('ENSG' = names(igraph::V(g)),'value' = 0)
  X$value[which(X$ENSG %in% ENSG_Genelist)] = 1
  igraph::V(g)$feature =  X$value
  g = message_passing(g, transcription_network = trascription_network, agg=aggregator)
  res_md = as_dataframe_mp_matches(g)
  return(res_md)
}

#' @title message_passing
#' @description This function performs message passing on a graph based on a transcription network and an aggregation method.
#' @param graph An igraph object representing the graph.
#' @param transcription_network A data frame containing the transcription network with columns 'target' and 'edge_type'.
#' @param agg A character string specifying the aggregation method. Options are 'mean' or 'degree_importance'. Default is 'mean'.
#' @return The updated graph with new vertex attributes 'feature_rep' and 'features_act'.
#' @export
message_passing <- function(graph, transcription_network, agg = 'mean') {
  features_rep <- igraph::V(graph)$feature
  names(features_rep) = names(igraph::V(graph))
  features_act = features_rep
  data_seeds = data.frame('Node' = names(igraph::V(graph)),'Seed' = igraph::V(graph)$feature)


  for (v in names(igraph::V(graph))) {
    neighbors <- igraph::neighbors(graph, v)

    if(length(neighbors != 0)){

      messages <- igraph::V(graph)$feature[neighbors]

      if (agg == 'mean'){
        # Mean Aggregation method
        aggregated_message_rep = 0
        aggregated_message_act = 0
        for(k in 1:length(neighbors)){
          n = names(neighbors)[k]
          if(names(neighbors)[k] %in% transcription_network$target){

            if(c('Repression') %in% transcription_network$edge_type[transcription_network$target == names(neighbors)[k]]){

              aggregated_message_rep = aggregated_message_rep + data_seeds$Seed[which(data_seeds$Node == names(neighbors)[k])]

            }else if(c('Activation') %in% transcription_network$edge_type[transcription_network$target == names(neighbors)[k]]){

              aggregated_message_act = aggregated_message_act + data_seeds$Seed[which(data_seeds$Node == names(neighbors)[k])]

            }}
        }

        features_rep[v] <- mean(aggregated_message_rep,as.numeric(features_rep[v]))
        features_act[v] <- mean(aggregated_message_act,as.numeric(features_rep[v]))

      }else if (agg == 'degree_importance'){
        # Degree based importance Aggregation method
        aggregated_message_rep = 0
        aggregated_message_act = 0

        for(k in 1:length(neighbors)){
          n = names(neighbors)[k]

          if(names(neighbors)[k] %in% transcription_network$target){

            if(c('Repression') %in% transcription_network$edge_type[transcription_network$target == names(neighbors)[k]]){

              aggregated_message_rep = aggregated_message_rep + data_seeds$Seed[which(data_seeds$Node == names(neighbors)[k])]/igraph::degree(graph)[n]

            }else if(c('Activation') %in% transcription_network$edge_type[transcription_network$target == names(neighbors)[k]]){

              aggregated_message_act = aggregated_message_act + data_seeds$Seed[which(data_seeds$Node == names(neighbors)[k])]/igraph::degree(graph)[n]

            }}

        }

        features_rep[v] <- mean(aggregated_message_rep) + as.numeric(features_rep[v])/igraph::degree(graph)[v]
        features_act[v] <- mean(aggregated_message_act) + as.numeric(features_rep[v])/igraph::degree(graph)[v]
      }

    }


  }

  igraph::V(graph)$feature_rep <- features_rep
  igraph::V(graph)$features_act <- features_act

  return(graph)
}

#' @title as_dataframe_mp_matches
#' @description This function converts the graph's vertex attributes 'feature_rep' and 'features_act' into a data frame.
#' @param g An igraph object representing the graph.
#' @return A data frame with columns 'ENSG', 'repression_matches', and 'activation_matches' containing the vertex names and their respective feature values.
#' @export
as_dataframe_mp_matches = function(g){

  rep_feature = igraph::V(g)$feature_rep
  act_feature = igraph::V(g)$features_act


  res_md = data.frame('ENSG' = names(igraph::V(g)),
                      'repression_matches' = rep_feature,
                      'activation_matches' = act_feature)
  return(res_md)
}
