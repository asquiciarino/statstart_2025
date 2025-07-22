# research project day 1

a = c(1,2,3)
b = c(2,3,4)

euclidean_distance = function(v1, v2) {
  if (length(v1) != length(v2)) {
    stop("Vectors v1 and v2 must have the same length")
  } else {
    return(sqrt(sum((v1 - v2)^2)))
  }
}

euclidean_distance(a,b)

library(bio3d)

my_pdb = read.pdb("~/Downloads/statstart_R/research_project/data/6obi.pdb")
my_pdb$seqres
my_pdb$xyz
my_pdb$atom$x

calc_dist_mat = function(file, anchor = "CA") {
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

c = calc_dist_mat("~/Downloads/statstart_R/research_project/data/6obi.pdb")

plot_distance_heatmap = function(dist_mat) {
  heatmap(dist_matrix, Rowv = NA, Colv = NA,
          col = heat.colors(256), scale = "none",
          xlab = "Residue Index", ylab = "Residue Index",
          main = "Distance Matrix Heatmap")
}

library(reshape2)
library(ggplot2)

ggplot_distance_heatmap = function(dist_mat) {
  m = melt(dist_mat)
  ggplot(data = m, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = heat.colors(256)) +
    theme_minimal()
}

ggplot_distance_heatmap(c)

d = calc_dist_mat("~/Downloads/statstart_R/research_project/data/6xb8.pdb")
ggplot_distance_heatmap(d)

e = calc_dist_mat("~/Downloads/statstart_R/research_project/data/2jst.pdb")
ggplot_distance_heatmap(e)

c[upper.tri(c)] == c[lower.tri(c)]

matrix(nrow = n, ncol = n)
m = (matrix(NA, 5, 5))
m[1:25] = 1:25
t(m)
all((m[upper.tri(m)]) == t(m[lower.tri(m)]))

# try to figure out how to figure out distance between three points (pythagorean theorem-esque, not right triangles)
# distance between 2 pts can give you the distance to the third

make_contact_map = function(dist_mat) {
  dist_vec = as.integer(dist_mat < 8)
  return(matrix(dist_vec, nrow = nrow(dist_mat), ncol = ncol(dist_mat)))
}

ggplot_distance_heatmap(make_contact_map(c))

library(igraph)

network_graph_c = graph_from_adjacency_matrix(make_contact_map(c), mode = "undirected")
plot(network_graph, layout = layout.kamada.kawai)
degree(network_graph_c)
hist(degree(network_graph_c))

make_contact_map(e)

network_graph_e = graph_from_adjacency_matrix(make_contact_map(e), mode = "undirected")
plot(network_graph, layout = layout.kamada.kawai)

degree(network_graph_e)
hist(degree(network_graph_e))


deg_c <- degree(network_graph_c)
deg_e <- degree(network_graph_e)

common_breaks <- hist(deg_e, plot = FALSE)$breaks

hist(deg_c, 
     col = rgb(1, 0, 0, 0.5),
     breaks = common_breaks,
     xlim = range(common_breaks),
     main = "Degree Distributions",
     xlab = "Degree",
     ylab = "Frequency")

hist(deg_e, 
     col = rgb(0, 0, 1, 0.5),
     breaks = common_breaks,
     add = TRUE)

legend("topright", 
       legend = c("6obi", "2jst"),
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))
