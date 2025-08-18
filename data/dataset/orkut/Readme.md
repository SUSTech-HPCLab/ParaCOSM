# Orkut dataset

The Orkut dataset, provided by the SNAP (Stanford Network Analysis Project) laboratory and originally curated by Alan Mislove et al., represents a large-scale online social network. Orkut was a free social networking platform where users could establish friendships and form or join user-defined groups, which are considered ground-truth communities in this dataset. 

The dataset captures the Orkut friendship network, consisting of 3,072,441 nodes (users) and 117,185,083 edges (friendship connections), with a fully connected largest weakly and strongly connected component. It exhibits an average clustering coefficient of 0.1666, contains 627,584,181 triangles, and has a diameter of 9, with a 90-percentile effective diameter of 4.8. This dataset is particularly valuable for studying social network structures, community detection, and graph-based algorithms due to its scale and the presence of well-defined community structures.

## Dataset information

| Statistic                          | Value            |
|------------------------------------|------------------|
| Nodes                              | 3,072,441        |
| Edges                              | 117,185,083      |
| Nodes in largest WCC               | 3,072,441 (1.000)|
| Edges in largest WCC               | 117,185,083 (1.000)|
| Nodes in largest SCC               | 3,072,441 (1.000)|
| Edges in largest SCC               | 117,185,083 (1.000)|
| Average clustering coefficient      | 0.1666           |
| Number of triangles                | 627,584,181      |
| Fraction of closed triangles       | 0.01414          |
| Diameter (longest shortest path)   | 9                |
| 90-percentile effective diameter   | 4.8              |


## Download link

https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz

## Reference

https://snap.stanford.edu/data/com-Orkut.html

J. Yang and J. Leskovec. Defining and Evaluating Network Communities based on Ground-truth. ICDM, 2012.