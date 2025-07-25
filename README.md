# 2D Visualizations of Protein Interactions

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This project uses 2D matrix-based techniques to explore the structural biology of proteins. By visualizing pairwise amino acid distances, contact maps, and interaction networks, we aim to uncover structural patterns and functional features from PDB files.

> Avery Squiciarino and Maliha Parmar  
> StatStart 2025 @ Harvard T.H. Chan School of Public Health

---

## Background

Proteins are long chains of amino acids that can fold into complex 3D structures. These structures determine how proteins function and interact with other molecules, playing an important role in drug discovery.

Protein structure can be categorized into:
- **Primary**: Amino acid sequence
- **Secondary**: Local folding (e.g., α-helices, β-sheets)
- **Tertiary**: Overall 3D structure of one chain
- **Quaternary**: Assembly of multiple polypeptide chains

Rather than working directly with 3D models, we converted structural data into **2D matrices** to simplify analysis and highlight biologically meaningful relationships.

---

## What Does the Code Do?

Our R code performs the following steps:

1. **Reads PDB files** using the `bio3d` package.
2. **Calculates pairwise Euclidean distances** between amino acids.
3. **Generates a distance matrix** and converts it to a binary **contact map** using a distance threshold.
4. **Builds a network graph** from contact data.
5. **Visualizes**:
   - Continuous distance heatmaps
   - Binary contact maps
   - Contact networks
   - Degree distributions

---

## Functions

| Function            | Purpose |
|---------------------|---------|
| `euclidean_distance()`| Calculates 3D Euclidean distance between two coordinate vectors |
| `calc_dist_mat()`     | Creates a pairwise distance matrix for a protein (using α-carbon atoms by default) |
| `make_contact_map()`  | Converts the distance matrix to a binary contact matrix with a defined threshold (default: 8Å) |
| `map_to_graph()`      | Turns the contact matrix into a network graph object using `igraph` |
| `plot_distance_map()` | Plots a continuous heatmap of amino acid distances |
| `plot_contact_map()`  | Plots a binary contact map |
| `plot_network()`      | Visualizes the binary contact graph object with a network-like layout |

---

## Graph Significance

- **Distance Matrices** and **Heatmaps** reveal intra- and inter-chain interactions.
- **Contact Maps** highlight which amino acids are in close proximity, uncovering important structural insights.
- **Network Graphs** provide an intuitive visualization of residue contacts.
- **Degree Distributions** compare contacts within and across different proteins, giving an overview of their connectivity.

---

## Acknowledgements

This project was created during the [StatStart 2025](https://hsph.harvard.edu/fellowship-special-program/statstart/) summer program.  
Special thanks to our mentor Mohammad Haddadnia for his guidance and feedback.

> Collaborators: [@malihaparmar](https://github.com/malihaparmar), [@asquiciarino](https://github.com/asquiciarino)

