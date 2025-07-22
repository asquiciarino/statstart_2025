[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Protein Contact Maps
Our code aims to construct various visualizations of pariwise amino acid distances within proteins.
We first wrote a function to find the continuous, pairwise, amino acid distance matrices of each protein, using another function we created that finds the euclidean distance between two points.
From this distance matrix, we produced a binary matrix of amino acid contacts from a distance matrix and specified threshold.
We then turned the binary protein contact matrix into a network-style graph.
We have created heatmaps from continuous distance matrices, contact maps from binary adjacency matrices, and networks from those same binary matrices. 

### Acknowedgements
StatStart 2025, Mentor 