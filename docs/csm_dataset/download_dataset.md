# CSM Dataset 汇总


## CSM related passage

目前CSM算法有这些文章：

#### Algorithm

Tc-match: Fast time-constrained continuous subgraph matching **VLDB 24**

CSMN: An Efficient Continuous Subgraph Matching Network for Dynamic Graphs **MLBDBI 24**

NewSP: A New Search Process for Continuous Subgraph Matching over Dynamic Graphs **ICDE 24**

CSM-TopK: Continuous Subgraph Matching with TopK Density Constraints **ICDE 24**

Efficient Multi-Query Oriented Continuous Subgraph Matching **ICDE 24**

Rapidflow: An efficient approach to continuous subgraph matching  **VLDB 22**

Symmetric continuous subgraph matching with bidirectional dynamic programming **VLDB 21**

Turboflux: A fast continuous subgraph matching system for streaming graph data   **SIGMOD 18**

Fast continuous subgraph matching over streaming graphs via backtracking reduction  **Proceedings of the ACM on Management of Data 23**

#### Parallelization

GPU-Accelerated Batch-Dynamic Subgraph Matching **ICDE 24**

Accelerating **Continuous Subgraph Matching** on Heterogeneous Processors  **学位论文**

GCSM: GPU-Accelerated Continuous Subgraph Matching for Large Graphs **IPDPS 24**

#### Survey
In-depth Analysis of **Continuous Subgraph Matching** in a Common Delta Query Compilation Framework    **期刊 Proceedings of the ACM on Management of Data**

A survey of continuous subgraph matching for dynamic graphs  **期刊 Knowledge and Information Systems** 

An in-depth study of continuous subgraph matching **VLDB22**



# Dataset from source

## Stanford SNAP Lab

link: https://snap.stanford.edu/data/index.html

| [com-Orkut](https://snap.stanford.edu/data/com-Orkut.html)               | Undirected, Communities | 3,072,441 | 117,185,083 | 6,288,363                            | Orkut online social network   |
| ------------------------------------------------------------------------ | ----------------------- | --------- | ----------- | ------------------------------------ | ----------------------------- |
| [soc-LiveJournal1](https://snap.stanford.edu/data/soc-LiveJournal1.html) | Directed                | 4,847,571 | 68,993,773  | LiveJournal online social network    |                               |
| [soc-Pokec](https://snap.stanford.edu/data/soc-Pokec.html)               | Directed                | 1,632,803 | 30,622,564  | Pokec online social network          |                               |
| [gemsec-Facebook](https://snap.stanford.edu/data/gemsec-Facebook.html)   | Undirected              | 134,833   | 1,380,293   | Gemsec Facebook dataset              |                               |
| [musae-github](https://snap.stanford.edu/data/github-social.html)        | Undirected              | 37,700    | 289,003     | Social network of Github developers. |                               |
| [com-Youtube](https://snap.stanford.edu/data/com-Youtube.html)           | Undirected, Communities | 1,134,890 | 2,987,624   | 8,385                                | Youtube online social network |
| [com-DBLP](https://snap.stanford.edu/data/com-DBLP.html)                 | Undirected, Communities | 317,080   | 1,049,866   | 13,477                               | DBLP collaboration network    |
| [com-Amazon](https://snap.stanford.edu/data/com-Amazon.html)             | Undirected, Communities | 334,863   | 925,872     | 75,149                               | Amazon product network        |
| [email-Eu-core](https://snap.stanford.edu/data/email-Eu-core.html)       | Directed, Communities   | 1,005     | 25,571      | 42                                   | E-mail network                |

## DC-2012 And UK-2007

### Introduction to the WebDataCommons Hyperlink Graph Dataset

The WebDataCommons Hyperlink Graph dataset, referenced in the paper *Tesseract: Distributed, General Graph Pattern Mining on Evolving Graphs* [EuroSys '21], was extracted from the 2012 Common Crawl web corpus. It captures hyperlinks between web pages, where nodes represent web pages and directed edges represent hyperlinks. The dataset is used for research in search algorithms, spam detection, and graph analysis scalability.[](http://webdatacommons.org/)

| Statistic                          | Value                     |
|------------------------------------|---------------------------|
| Nodes                              | 3,500,000,000             |
| Edges                              | 128,000,000,000           |


**Source**: Web Data Commons - Hyperlink Graphs. Available at: [https://webdatacommons.org/hyperlinkgraph/](https://webdatacommons.org/hyperlinkgraph/)[](http://webdatacommons.org/)

### WebGraph Format
The dataset is provided in the WebGraph format, an early compressed graph format developed by the [Laboratory for Web Algorithmics](https://law.di.unimi.it/index.php). Another dataset in this format is [uk-2007-05](https://law.di.unimi.it/webdata/uk-2007-05/). Due to the outdated Maven package for WebGraph, pre-converted data is available, as detailed in [Laboratory for Web Algorithmics Dataset Format Conversion | Chuckie's Blog](https://chuckiewill.github.io/2022/06/29/Graph/LWA-dataset-format-conversion/).

### Decompression Instructions
To process the WebGraph format dataset, use the following commands:

 Convert ClueWeb12 Graph
```bash
java -cp "lib/*" it.unimi.dsi.webgraph.ASCIIGraph clueweb12 clueweb12
```

 Convert Host Graph
```bash
java -cp "lib/*" it.unimi.dsi.webgraph.ASCIIGraph hostgraph hostgraph
```

 Convert to Edge List
```bash
java -cp "lib/*" it.unimi.dsi.webgraph.ArcListASCIIGraph hostgraph hostgraph-edgelist.txt
```

These commands leverage the WebGraph library to convert the compressed graph data into usable formats, such as ASCII or edge lists, for further analysis.

## GunRock data

Gunrock: GPU Graph Analytics

[gunrock/datasets/soc-twitter-2010/Makefile at main · gunrock/gunrock](https://github.com/gunrock/gunrock/blob/main/datasets/soc-twitter-2010/Makefile)


---

# Dataset from paper & format

##  Symbi & TurboISO & TurboFlux

[[2104.00886] Symmetric Continuous Subgraph Matching with Bidirectional Dynamic Programming](https://arxiv.org/abs/2104.00886)

[TurboFlux | Proceedings of the 2018 International Conference on Management of Data](https://dl.acm.org/doi/10.1145/3183713.3196917)

[Turboiso | Proceedings of the 2013 ACM SIGMOD International Conference on Management of Data](https://dl.acm.org/doi/10.1145/2463676.2465300)

In-depth Analysis of Continuous Subgraph Matching in a CommonDeltaQueryCompilation Framework


link: [SNUCSE-CTA/SymBi: Symmetric Continuous Subgraph Matching with Bidirectional Dynamic Programming](https://github.com/SNUCSE-CTA/SymBi)

Dataset
- Netflow ([download](https://drive.google.com/file/d/1g1NVK29k27V76JrSAVhKOZ2qvyP_kj--/view?usp=sharing))
- LSBench_x1 ([download](https://drive.google.com/file/d/1MmfMTE2SKPOJaaSibPC-5OPXh7brDulf/view?usp=sharing))
- LSBench_x5 ([download](https://drive.google.com/file/d/1k0eQdZYElAsBqk5NfhW414FsZuhDNsqG/view?usp=sharing))
- LSBench_x25 ([download](https://drive.google.com/file/d/1c_pGxT18H-BOL0cirdTLJ-0qbqcaDADx/view?usp=sharing))
- stacksboverflow 1.32GB
- yahoo 70.2 MB

Format
```txt
t # s 1
v 0 0 -1
v 1 0 -1
v 2 0 -1
v 3 0 -1
v 4 0 -1
v 5 0 -1
e 1 7 41
e 7 4 0
e 6 1 41
e 0 6 15
e 0 1 33
```



---
## TC-Match 数据格式

vldb 2024
该图问题的主要特点为，每个边有特定timestamp，跟传统的CSM问题不太一致

https://github.com/Sh-Fang/TCMatch.git

Both the stream graph and query graph are edge-based.

Each stream edge is represented by: `(v_source_id,v_target_id,edge_label,v_source_label,v_target_label,timestamp)`.

Each query edge is represented by: `(qid,u_source_id,u_target_id,edge_label,u_source_label,u_target_label)`.

### Query Graph
 Each line in the query graph file represent an edge or a time-constraint order.

1. An edge is represented by `e <qid> <u_source_id> <u_target_id> <edge_label> <u_source_label> <u_target_label>`.
```txt
e 0 1 2 0 80 80
e 1 1 3 1 443 21 
e 2 2 3 0 8080 22
b 0 1 2
```

### Stream Graph

Each line in the stream graph file represent a time-edge.

An edge is represented by `<v_source_id> <v_target_id> <edge_label> <v_source_label> <v_target_label> <timestamp>`.

For example:
```
10 20 0 8080 80 1
10 30 1 443 21 2
20 30 0 8080 22 3
```


---
## SIGMOD 2023 Journal CaLiG 


Fast Continuous Subgraph Matching over Streaming Graphs via Backtracking Reduction

https://github.com/JackChuengQAQ/CaLiG.git

### data graph / query graph
```
v <id> <label>
...
e <id1> <id2>
...
```
- Lines starting with "v" represent vertices;
- Lines starting with "e" represent edges.

### stream
```
e <id1> <id2>
...
e <-id1-1> <-id2-1>
...
```
- `e <id1> <id2>` represents the addition of the edge `(id1, id2)`;
- `e <-id1-1> <-id2-1>` represents the deletion of the edge `(id1, id2)`.

---
## ICDE 24 NewSP   

https://github.com/hnuGraph/NewSP.git
### Query Graph

Each line in the query graph file represent a vertex or an edge.
1. A vertex is represented by `v <vertex-id> <vertex-label>`;
2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`.

The two endpoints of an edge must appear before the edge. For example,
```
v 0 0
v 1 0
e 0 1 0
v 2 1
e 0 2 1
e 2 1 2
```

### Initial Data Graph

The initial data graph file has the same format as the query graph file.

### Graph Update Stream

Graph update stream is a collection of insertions and deletions of a vertex or an edge.

1. A vertex insertion is represented by `v <vertex-id> <vertex-label>`;
2. A vertex deletion is represented by `-v <vertex-id> <vertex-label>`;
3. An edge insertion is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`;
4. An edge deletion is represented by `-e <vertex-id-1> <vertex-id-2> <edge-label>`;

The vertex or edge to be deleted must exist in the graph, and the label must be the same as that in the graph. If an edge is inserted to the data graph, both its endpoints must exist. For example,

```
v 3 1
e 2 3 2
-v 2 1
-e 0 1 0
```

## Datasets

We provide 4 datasets in our experiment

1. Amazon dataset [[link](https://snap.stanford.edu/data/com-Amazon.html)]
2. Livejournal dataset [[link](https://snap.stanford.edu/data/soc-LiveJournal1.html)]
3. LSBench dataset [[link](https://code.google.com/archive/p/lsbench/)]
4. Netflow dataset [[link](https://catalog.caida.org/dataset/passive_2013_pcap)]

| **Datasets** |         **Type**         | **Vertexes** | **Edges**  | **Average Degree** |
| :----------: | :----------------------: | :----------: | :--------: | :----------------: |
|    Amazon    |     Product network      |   403,394    | 2,433,408  |       12.06        |
| Livejournal  |    Community network     |  4,847,571   | 42,841,237 |       17.68        |
|   LSBench    | Benchmark data generator |  5,210,099   | 20,270,676 |        7.78        |
|   Netflow    |     Network traffic      |  3,114,895   | 2,849,732  |        1.83        |


## VLDB22 Summary on CSM

https://github.com/RapidsAtHKUST/ContinuousSubgraphMatching.git
1. iedyn
2. TurboFlux (2018)
3. Graphflow (2017)
4. SymBi
5. SJ-Tree

## Input File Format
Both the input query graph and data graph are vertex- and edge-labeled. Each vertex is represented by a distinct unsigned integer (from 0 to 4294967295). There is at most one edge between two arbitrary vertices.

### Query Graph
Each line in the query graph file represent a vertex or an edge.
1. A vertex is represented by `v <vertex-id> <vertex-label>`;
2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`.

The two endpoints of an edge must appear before the edge. For example,
```
v 0 0
v 1 0
e 0 1 0
v 2 1
e 0 2 1
e 2 1 2
```

### Initial Data Graph
The initial data graph file has the same format as the query graph file.

### Graph Update Stream

Graph update stream is a collection of insertions and deletions of a vertex or an edge.
1. A vertex insertion is represented by `v <vertex-id> <vertex-label>`;
2. A vertex deletion is represented by `-v <vertex-id> <vertex-label>`;
3. An edge insertion is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`;
4. An edge deletion is represented by `-e <vertex-id-1> <vertex-id-2> <edge-label>`;
The vertex or edge to be deleted must exist in the graph, and the label must be the same as that in the graph. If an edge is inserted to the data graph, both its endpoints must exist. For example,

```
v 3 1
e 2 3 2
-v 2 1
-e 0 1 0
```

Bro的数据集：

| **Datasets** |         **Type**         | **Vertexes** | **Edges**  | **Average Degree** |
| :----------: | :----------------------: | :----------: | :--------: | :----------------: |
|    Amazon    |     Product network      |   403,394    | 2,433,408  |       12.06        |
| Livejournal  |    Community network     |  4,847,571   | 42,841,237 |       17.68        |
|   LSBench    | Benchmark data generator |  5,210,099   | 20,270,676 |        7.78        |
|   Netflow    |     Network traffic      |  3,114,895   | 2,849,732  |        1.83        |

incresement数据：自己在原始图里选取的
query 数据：自己在原始图里sample的

## VLDB22 RapidFlow

vldb 2022

https://github.com/shixuansun/RapidFlow.git
The datasets and baseline methods in our experiments can be found in this [repository](https://github.com/RapidsAtHKUST/ContinuousSubgraphMatching), which is an in-depth study of existing continuous subgraph matching methods by our research team.




