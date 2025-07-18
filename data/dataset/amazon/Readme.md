# Amazon dataset


The Amazon dataset was collected by crawling the Amazon website, based on the "Customers Who Bought This Item Also Bought" feature, where edges represent frequent co-purchases between products. Product categories define ground-truth communities, with connected components as communities (excluding those with fewer than 3 nodes). It is useful for studying co-purchase patterns and community detection in e-commerce networks.

## Dataset information

| Statistic                          | Value            |
|------------------------------------|------------------|
| Nodes                              | 334,863          |
| Edges                              | 925,872          |
| Nodes in largest WCC               | 334,863 (1.000)  |
| Edges in largest WCC               | 925,872 (1.000)  |
| Nodes in largest SCC               | 334,863 (1.000)  |
| Edges in largest SCC               | 925,872 (1.000)  |
| Average clustering coefficient      | 0.3967           |
| Number of triangles                | 667,129          |
| Fraction of closed triangles       | 0.07925          |
| Diameter (longest shortest path)   | 44               |
| 90-percentile effective diameter   | 15               |


## Download link

https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz


## Reference

https://snap.stanford.edu/data/com-Amazon.html

J. Yang and J. Leskovec. Defining and Evaluating Network Communities based on Ground-truth. ICDM, 2012.