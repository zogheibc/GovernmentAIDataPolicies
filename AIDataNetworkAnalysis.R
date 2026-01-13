##load libraries
library(tidyverse)
library(igraph)
library(ggnetwork)
library(intergraph)
library(readxl)
library(data.table)
library(tidygraph)
library(ggraph)
library(viridis)

##set seed
set.seed(613)


##load data
networkData <- read_excel("networkanalysis.xlsx", 
                          sheet = "network", col_types = c("text", 
                                                           "text", "numeric", "numeric"))

nodeDefinitions <- read_excel("networkanalysis.xlsx", 
                              sheet = "definitions")


##make table of nodes
##rename column in our nodeDefinitions table (nodes and info)
nodeDefinitions <- nodeDefinitions |>
  rename(node = document)
##make lowercase
nodeDefinitions <- nodeDefinitions |>
  mutate(across(c(engtitle, frtitle), tolower))
##reorder columns
move_before <- function(df, col_to_move, before_col) {
  cols <- names(df)
  cols <- cols[cols != col_to_move]  # remove the column first
  before_idx <- match(before_col, cols)
  new_order <- append(cols, col_to_move, after = before_idx - 1)
  df[new_order]
}

nodeDefinitions <- move_before(nodeDefinitions, "node", "level")


##make table of edges
networkData <- networkData |>
  rename(from = node1)
networkData <- networkData |>
  rename(to = node2)
##remove rows with 0 interactions
networkData <- networkData |>
  filter(connection != 0)




##### MAKE OUR NETWORK


ig <- igraph::graph_from_data_frame(d = networkData, vertices = nodeDefinitions, directed = FALSE)

tg <- tidygraph::as_tbl_graph(ig) |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(label = shorttitle,
                topic = topic,
                level = level)

##make 'level' a factor so we can order the legend
tg <- tg |>
  activate(nodes) |>
  mutate(level = factor(level, levels = c("Federal", "Provincial", "Municipal")))


##### MAKE OUR GRAPH

##create a graph using the Fruchterman-Reingold layout
tg |>
  ggraph::ggraph(layout = "stress") +
  
  ##modify our edge appearance
  geom_edge_arc(
    colour = "gray60",
    lineend = "round",
    strength = .1, ##how curved the arcs are
    aes(
      edge_width = networkData$numberofconnections
    )
  ) +
  
  ##change the appearance of our nodes
  geom_node_point(
    aes(fill = level,
        shape = topic),
    size = 6,
    color = "black",
    stroke = 0.5
  ) +
  
  #add text labels to our nodes
  geom_node_text(aes(label = shorttitle),
                 repel = TRUE,
                 point.padding = unit(0.2, "lines"),
                 colour = "gray10",
                 size = 3,
                 fontface = "bold"
  ) +
  
  ##adjust edge width a
  scale_edge_width(range = c(0, 3.5)) +

  ##set colours for fills of nodes; set shapes for nodes
  scale_shape_manual(values = 21:25) +
  scale_fill_viridis_d(option = "C") +
  
  ##set graph background color as white
  theme_graph(background = "white") +
  
  ##modify our legend appearance and remove some features from the legend
  theme(
    legend.position = "right",
    # suppress legend title
    legend.title = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(override.aes = list(size = 6, color = "black")),
    edge_width = FALSE,
    edge_alpha = FALSE
  )




#####ANALYSIS OF THE NETWORK


##some properties of the network
cat("Number of nodes:", gorder(ig), "\n")
cat("Number of edges:", gsize(ig), "\n")
cat("Network density:", edge_density(ig), "\n\n")

##get degree and betweenness
node_names <- V(ig)$name
##degree (number of connections per node)
degree_vals <- degree(ig)
##betweenness (how much is a node the bridge between other nodes)
betweenness_vals <- betweenness(ig, weights = 1/E(ig)$weight)
##closeness (how close is a node to other nodes in the network)
closeness_vals <- closeness(ig, weights = 1/E(ig)$weight)
##eigenvector centrality (nodes are more important if they're connected to other nodes that are also important)
eigen_vals <- eigen_centrality(ig, weights = E(ig)$weight)$vector
##collect centrality data (above) into a dataframe and print
centrality_df <- tibble(
  node = node_names,
  betweenness = betweenness_vals,
  closeness = closeness_vals,
  eigen = eigen_vals
)
print(centrality_df) 
##export as csv
#write_csv(centrality_df, "centrality_results.csv")

##ASSORTATIVITY BY LEVEL
##check assortativity by attribute to see if docs at the same level of govt are more likely to be connected to each other than by chance
##convert 'level' to numeric factors for assortativity_nominal
V(ig)$level_factor <- as.factor(V(ig)$level)
##nominal assortativity by level
assortativity_level <- assortativity_nominal(
  graph = ig,
  types = V(ig)$level_factor,
  directed = FALSE
)
cat("Assortativity by level:", assortativity_level, "\n")
##value = 0.781; this network shows HOMOPHILY; docs are more likely to reference other docs at the same level

##ASSORTATIVITY BY DOC TYPE
##check assortativity by attribute to see if docs at the same level of govt are more likely to be connected to each other than by chance
##convert 'topic' to numeric factors for assortativity_nominal
V(ig)$topic_factor <- as.factor(V(ig)$topic)
##nominal assortativity by topic
assortativity_topic <- assortativity_nominal(
  graph = ig,
  types = V(ig)$topic_factor,
  directed = FALSE
)
cat("Assortativity by topic:", assortativity_topic, "\n")
##value = 0.238; this network shows HOMOPHILY; documents are more likely to reference other docs on the same topic (AI/Data/Info)



###calculate odds of same topic connections
##convert to adjacency matrix
adj <- as_adjacency_matrix(ig, sparse = FALSE)
topics_vec <- V(ig)$topic

###make logical matrix where value is TRUE if nodes have the same topic
same_topic <- outer(topics_vec, topics_vec, FUN = "==")

##count actual connections
same_conn <- sum(adj[same_topic & upper.tri(adj)])       # same-topic connections
diff_conn <- sum(adj[(!same_topic) & upper.tri(adj)])    # different-topic connections

##count possible connections
possible_same <- sum(same_topic[upper.tri(same_topic)])
possible_diff <- sum((!same_topic)[upper.tri(same_topic)])

##do ratio of actual to possible to get odds ratio
odds_ratio_topic <- (same_conn / (possible_same - same_conn)) /
  (diff_conn / (possible_diff - diff_conn))

cat("Odds ratio of same-topic connections:", odds_ratio_topic, "\n")
##########nodes are 3x as likely to connect to nodes of the same topic than to nodes of diff topics



###calculate odds of same level connections
##convert to adjacency matrix
adj <- as_adjacency_matrix(ig, sparse = FALSE)
levels_vec <- V(ig)$level

###make logical matrix where value is TRUE if nodes have the same topic
same_level <- outer(levels_vec, levels_vec, FUN = "==")

##count actual connections
same_conn <- sum(adj[same_level & upper.tri(adj)])       # same-topic connections
diff_conn <- sum(adj[(!same_level) & upper.tri(adj)])    # different-topic connections

##count possible connections
possible_same <- sum(same_level[upper.tri(same_level)])
possible_diff <- sum((!same_level)[upper.tri(same_level)])

##do ratio of actual to possible to get odds ratio
odds_ratio_level <- (same_conn / (possible_same - same_conn)) /
  (diff_conn / (possible_diff - diff_conn))

cat("Odds ratio of same-topic connections:", odds_ratio_level, "\n")
##########nodes are 44x as likely to connect to nodes of the same level than to nodes of diff levels





