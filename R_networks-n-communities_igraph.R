# Load necessary libraries
library(igraph)
library(ggplot2)

#set working directory

#create function of extracting data to csv
write_to_csv <-  function (data, filename) {
  write.csv( data, file= filename, row.names=FALSE)
}

#--------------------------
# 1- Create igraph objects
#--------------------------

  #create Folder path of my edges
  files_path <- "Networks/edges_P/edges_P_wthVerticesInfo"
  
  # List CSV files in the directory in chronological order
  file_names <- list.files(path = files_path, pattern = "\\.csv$", full.names = TRUE)


#1.1- convert edges to igraph objects
  # Initialise a list to store graphs
  graphs_list <- list()
  
  # Read CSV files, create igraph objects
  for (file in file_names) {
    
    edge_list <- read.csv(file)
    
    # Create vertices data frame directly from edge list
    vertices_info <- unique(rbind(
      data.frame(name=edge_list$gene1, type=edge_list$genetype1),
      data.frame(name=edge_list$gene2, type=edge_list$genetype2)
    ))
    
    # Create graph from edge list - Select only necessary columns: genes & avg_correlation
    g <- graph_from_data_frame(
      d = edge_list[, c("gene1", "gene2", "avg_correlation")],
      vertices = vertices_info, 
      directed = FALSE)
    
    # store the graph with basename of the file - without the path!
    graphs_list[[basename(file)]] <- g
  }

  #check that both vertices and edges attributes have been added
  print(graphs_list[["edges_E1.csv"]], e=TRUE, v=TRUE)

#------------------------------
#2- Perform community detection
#------------------------------

# 2.1 Create Function to perform community detection and plot results and plots save as pdf
      #Note - I couldnt color the nodes by their gene type

    # Define the function to perform community detection and plot
    perform_community_detection <- function(graph, method, graph_name) {
      
      # Set up the PDF file for output, naming the file after the graph and method
      pdf(paste0(graph_name, "_", method, ".pdf"), width = 7, height = 5)  # Adjust size as needed
      
      # Perform community detection based on the specified method
      communities <- switch(method,
                            "walktrap" = cluster_walktrap(graph, weights = E(graph)$weight),
                            "louvain" = cluster_louvain(graph, weights = E(graph)$weight),
                            "eigenvector" = cluster_leading_eigen(graph),
                            NULL)
      
      if (is.null(communities)) {
        stop("Invalid method specified")
      }
  
  
      # Define colors for each type of gene
      genetype_colors <- setNames(c("blue", "red3"), c("protein_coding", "lncRNA"))
      
      # Plot the graph with community detection results
      plot(communities, graph, 
           main = paste(graph_name, "-", method),  # Include graph_name in the plot title
           vertex.size = 5,                       # Adjust vertex size
           vertex.label.cex = 0.8,                # Adjust vertex label text size
           vertex.label.color = "black", 
           vertex.color = genetype_colors[V(graph)$type],  # Color nodes by their type
           edge.color = "darkgrey",                 # Use light grey for edges
           edge.arrow.size = 0.5)                 # Adjust arrow size for directed graphs if necessary
      
      # Add a legend to the plot
      legend("bottomright", 
             title = "Gene Type",
             legend = names(genetype_colors), 
             fill = genetype_colors,
             cex = 0.8,
             bty = "n")  # 'n' for no box around the legend
      
      # Close the PDF device
      dev.off()
    }

    # Example of usage within a loop over named graphs
    for (name in names(graphs_list)) {
      g <- graphs_list[[name]]
      perform_community_detection(g, "walktrap", name)
      perform_community_detection(g, "louvain", name)
      perform_community_detection(g, "eigenvector", name)
    }


#----------------------------------------
# 3- Now extract communities information
#----------------------------------------

#3.1 create function to extract community information
    extract_community_info <- function(graph, communities, filename) {
      # Extract community membership
      membership <- membership(communities)
      
      # Create a data frame to store community data
      community_data <- data.frame(gene = names(V(graph)), 
                                   community = membership, 
                                   stringsAsFactors = FALSE)
      write_to_csv(community_data, filename)
      
      return(community_data)
    }


#3.2 create function to calculate and save modularity scores
    calculate_modularity <- function(graph, communities, filename) {
      #calculate modularity
      modularity_score <- modularity(communities)
      
      #create dataframe
      modularity_data <- data.frame(Modularity= modularity_score)
      
      write_to_csv(modularity_data, filename)
      
      return(modularity_score)
    }


#3.3 create function to calculate and save centrality measures
    calculate_centrality <- function(graph, filename) {
      # Calculate betweenness
      betweenness_centrality <- betweenness(graph, v = V(graph), directed = FALSE)
      
      # Calculate other centralities as needed
      eigenvector_centrality <- evcent(graph)$vector
      
      centrality_data <- data.frame(gene = names(V(graph)), 
                                    betweenness = betweenness_centrality, 
                                    eigenvector = eigenvector_centrality, 
                                    stringsAsFactors = FALSE)
      
      return(centrality_data)
    }


#3.4 combine all above fucntions in one function to perform analysis

    community_detection <- function(graph, method, graph_name) {
      # Perform community detection based on the specified method
      communities <- switch(method,
                            "walktrap" = cluster_walktrap(graph, weights = E(graph)$weight),
                            "louvain" = cluster_louvain(graph, weights = E(graph)$weight),
                            "eigenvector" = cluster_leading_eigen(graph),
                            stop("Invalid method specified"))
      
      # Extract and return community information
      community_data <- extract_community_info(graph, communities)
      
      # Optionally, calculate modularity and save or use it
      modularity_score <- calculate_modularity(graph, communities)
      cat("Modularity Score for", method, "on", graph_name, ":", modularity_score, "\n")
      
      # Calculate and print centrality measures if needed
      centrality_data <- calculate_centrality(graph)
      
      # Function can be extended to handle these results as needed
      
      return(list(community_data = community_data, 
                  modularity_score = modularity_score, 
                  centrality_data = centrality_data))
    }

# Usage example within a loop over named graphs
community_results <- list()

for (name in names(graphs_list)) {
  g <- graphs_list[[name]]
  community_results[[paste(name, "walktrap")]] <- community_detection(g, "walktrap", name)
  community_results[[paste(name, "louvain")]]   <- community_detection(g, "louvain", name)
  community_results[[paste(name, "eigenvector")]] <- community_detection(g, "eigenvector", name)
}


#*********************************************************************************

library(igraph)


# Function to perform community detection and extract metrics
perform_community_detection <- function(graph, method, graph_name) {
  # Perform community detection based on the specified method
  communities <- switch(method,
                        "walktrap" = cluster_walktrap(graph, weights = E(graph)$weight),
                        "louvain" = cluster_louvain(graph, weights = E(graph)$weight),
                        "eigenvector" = cluster_leading_eigen(graph),
                        stop("Invalid method specified"))
  
  # Extract community information
  membership <- membership(communities)
  community_data <- data.frame(gene = names(V(graph)), community = membership, stringsAsFactors = FALSE)
  
  # Calculate modularity
  modularity_score <- modularity(communities)
  
  # Calculate centrality measures
  betweenness_centrality <- betweenness(graph, v = V(graph), directed = FALSE)
  eigenvector_centrality <- evcent(graph)$vector
  
  # Combine all results into one data frame
  results <- data.frame(
    GraphName = graph_name,
    Method = method,
    Gene = names(V(graph)),
    Community = membership,
    Modularity = modularity_score,
    Betweenness = betweenness_centrality,
    Eigenvector = eigenvector_centrality
  )
  
  return(results)
}

# Example usage within a loop over named graphs, save combined results
all_results <- list()

for (name in names(graphs_list)) {
  g <- graphs_list[[name]]
  results_walktrap <- perform_community_detection(g, "walktrap", name)
  results_louvain <- perform_community_detection(g, "louvain", name)
  results_eigenvector <- perform_community_detection(g, "eigenvector", name)
  
  # Combine results from different methods
  combined_results <- rbind(results_walktrap, results_louvain, results_eigenvector)
  all_results[[name]] <- combined_results



  #hello everyonE
  
  # Save combined results to a CSV file
  write_to_csv(combined_results, paste0(name, "_combined_results.csv"))
}

