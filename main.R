# clean research code

library(bio3d)
library(reshape2)
library(ggplot2)
library(igraph)
library(ThemePark)

####################
#    FUNCTIONS     #
####################

euclidean_distance = function(v1, v2) {
  # calculates the euclidean distance between two 3D points
  # inputs:
  # v1: the first point
  # v2: the second point
  if (length(v1) != length(v2)) {
    stop("Vectors v1 and v2 must have the same length")
  } else {
    return(sqrt(sum((v1 - v2)^2)))
  }
}

calc_dist_mat = function(file, anchor = "CA") {
  # calculates the distance matrix of a protein pdb file
  # inputs:
  # file: path to the pdb file that is being analyed
  # anchor: anchor atom 
  # returns: the distance matrix of the inputted protein
  pdb = read.pdb(file)
  selected_idx = which(pdb$atom$elety == anchor)
  filtered_table = pdb$atom[selected_idx, ]
  coords = filtered_table[, c("x", "y", "z")] # [n, 3]
  n = nrow(coords)
  dist_mat = matrix(nrow = n, ncol = n) 
  
  for (i in 1:n) {
    a1 = coords[i,]
    for (j in 1:n) {
      a2 = coords[j,]
      dist = euclidean_distance(a1, a2)
      dist_mat[i, j] = dist
      dist_mat[j, i] = dist
    }
  }
  return(dist_mat)
}

make_contact_map = function(dist_mat, threshold = 8) {
  # produces a binary matrix of amino acid contacts from a distance matrix and
  # specified threshold, with 1 representing contact and 0 representing no contact
  # input/dist_mat: the raw, numeric distance matrix of a protein
  # return: calculated matrix
  dist_vec = as.integer(dist_mat < threshold)
  binary_mat = matrix(dist_vec, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
  return(binary_mat)
}

map_to_graph = function(binary_mat) {
  # turns a binary protein contact matrix into a network-style graph
  # input/binary_mat: the binary matrix of protein contacts
  # return: constructed graph object
  network_graph = graph_from_adjacency_matrix(binary_mat, mode = "undirected")
  return(network_graph)
}

plot_distance_map = function(dist_mat) {
  # plots a heatmap of a calulated distance matrix
  # input/dist_mat: a distance matrix of all amino acids in a protein
  m = melt(dist_mat)
  ggplot(data = m, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(256)) +
    theme_minimal()
}

plot_contact_map = function(binary_mat) {
  # creates a binary heatmap of amino acid contacts
  # input/binary_mat: a binary contact matrix of all amino acids in a protein
  m = melt(binary_mat)
  ggplot(data = m, aes(x = Var1, y = Var2, fill = factor(value))) +
    geom_tile() +
    scale_fill_manual(values = c("white", "black")) +
    theme_minimal()
}

plot_network = function(network_graph) {
  # creates a network-style visualization of all contacts in a protein
  # input/network_graph: the calculated graph from an adjacency matrix
  plot(network_graph, layout = layout.kamada.kawai)
}

##########
#  MAIN  #
##########

d1 = calc_dist_mat("~/Downloads/statstart_R/research_project/data/6obi.pdb")
d2 = calc_dist_mat("~/Downloads/statstart_R/research_project/data/2jst.pdb")
d3 = calc_dist_mat("~/Downloads/statstart_R/research_project/data/8qdn.pdb")
d4 = calc_dist_mat("~/Downloads/statstart_R/research_project/data/2nog.pdb")

# create a function for the histograms of the degrees of nodes