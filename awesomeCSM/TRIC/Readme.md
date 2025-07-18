# Efficient Continuous Multi-Query Processing over Graph Streams

The paper titled **"Efficient Continuous Multi-Query Processing over Graph Streams"** introduces a novel method for efficiently handling *multiple continuous subgraph queries* on *evolving graph streams*, a scenario common in applications like social networks, biological networks, and network security.



### Link

https://arxiv.org/abs/1902.05134
https://github.com/Sh-Fang/TRIC-reproduction.git

### 🌐 Motivation

Graphs in many real-world applications evolve rapidly. Applications often want to detect patterns (e.g., spam, attack behavior, protein interactions) by issuing *continuous subgraph queries*. But:

* Most existing methods only support **single queries**.
* Evaluating **many concurrent queries** is inefficient, especially when their patterns overlap.

---

### 🧠 Key Idea

Introduce **TRIC (TRIe-based Clustering)**:

* Observes that **many queries share common subgraph patterns**.
* **Indexes and clusters** the queries using **tries** built from shared *query paths*.
* **Shares computation and materialized views** across queries to reduce redundant work.

---

### 🔧 Method Overview

1. **Query Decomposition**:

   * Each query graph is broken into a minimal set of **covering paths** (sequences of edges that represent the whole query).
   * Shared paths across queries are identified.

2. **Trie-based Indexing (TRIC)**:

   * The paths are stored in **tries**, where **shared prefixes** are naturally merged.
   * This clusters queries that share structure.

3. **Query Answering on Updates**:

   * When a new graph update (e.g., an edge is added) arrives:

     * TRIC identifies **which queries are potentially affected**.
     * It **incrementally updates materialized views** and joins only what’s necessary.
     * Efficiently determines which queries match after each update.

4. **Optimized Variant (TRIC+)**:

   * Adds **caching** of intermediate join results and hash tables for further speedups.

---

### 📊 Baselines and Evaluation

Compared against:

* **INV / INV+**: Inverted index-based methods.
* **INC / INC+**: Inverted index + incremental join strategies.
* **Neo4j**: A commercial graph database.

On datasets from **social networks (LDBC SNB)**, **NYC taxi rides**, and **BioGRID** (protein interactions), the results show:

* **TRIC outperforms** all baselines by **up to two orders of magnitude** in query processing time.
* TRIC+ offers even better performance with caching.
* TRIC scales well with increasing query count, query overlap, and graph size.

---

### 📌 Contributions

* First formalization of **continuous multi-query answering over graph streams**.
* Proposal of **TRIC**: a novel, efficient method for clustering and answering large sets of continuous graph queries.
* Extensive evaluation demonstrating significant speedups and scalability over alternatives.


