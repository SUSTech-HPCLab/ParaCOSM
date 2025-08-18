# TC-Match: Fast Time-Constrained Continuous Subgraph Matching
The paper titled "TC-Match: Fast Time-Constrained Continuous Subgraph Matching" proposes an innovative approach for efficiently detecting structural patterns in streaming graphs under timing constraints, addressing critical needs in applications like cyberattack detection and credit card fraud detection.

## Link

https://www.vldb.org/pvldb/vol17/p2791-yang.pdf

https://github.com/Sh-Fang/TCMatch.git


## 🌐 Motivation
Streaming graphs, where edges and vertices are continuously updated, are prevalent in real-time applications such as network security and financial transaction monitoring. Key challenges include:

Time-constrained subgraph matching: Existing methods struggle to incorporate timing order constraints, crucial for applications like detecting sequential cyberattack patterns or fraudulent transactions.
High computational overhead: State-of-the-art solutions, such as Timing, incur significant index space and intermediate result maintenance costs, limiting scalability.


## 🧠 Key Idea
The paper introduces TC-Match, a novel method for Time-Constrained Continuous Subgraph Matching (TCSM) that:

Utilizes a space- and time-efficient index structure called Candidate Storage Structure (CSS), designed as a k-partite graph to store partial matches and timing order information.
Supports incremental updates for edge insertions and deletions, minimizing redundant computations.
Ensures isomorphism and timing constraints are met for accurate and efficient matching.


## 🔧 Method Overview

### Candidate Storage Structure (CSS):

A k-partite graph where each node represents an edge in the data graph ( G ), grouped by matching edges in the query graph ( Q ).
Edges in CSS connect nodes if their corresponding edges in ( G ) are adjacent and satisfy the timing constraints of ( Q ).
Incorporates node states to track valid matches, reducing redundant traversals.


### Incremental Matching:

For edge insertions, TC-Match updates the CSS index and identifies new matches by traversing only the affected neighborhood.
Edge deletions are handled symmetrically by removing invalid matches and updating the index.
The approach leverages depth-first traversal to generate matching orders efficiently.


### Efficient Query Processing:

Filters out irrelevant edges in ( G ) that do not match any edge in ( Q ).
Uses the CSS to prune unpromising edges or substructures, minimizing computational overhead.




## 📊 Baselines and Evaluation
TC-Match is compared against:

Timing: The state-of-the-art algorithm for TCSM, which relies on materializing intermediate results.
CaliG: A baseline method for subgraph matching, limited by its inability to handle multi-labeled graphs efficiently.

Evaluations on datasets such as WikiTalk, CAIDA, LiveJournal, and LSBench demonstrate:

Superior performance: TC-Match is significantly faster than Timing and CaliG, particularly on complex datasets like LiveJournal, where Timing fails due to high materialization costs.
Low memory usage: The CSS index is space-efficient, consuming less memory than Timing and other baselines.
Scalability: TC-Match handles large query sets with minimal unsolved queries, excelling on datasets like WikiTalk and CAIDA.
Query time breakdown: The CSS update phase dominates the processing time, highlighting its effectiveness in filtering unpromising structures, though enumeration time varies with match volume on datasets like CAIDA.


## 📌 Contributions

Novel TCSM framework: First comprehensive solution for time-constrained continuous subgraph matching in streaming graphs.
Efficient CSS index: A lightweight, edge-centric structure that encapsulates partial matches and timing constraints, reducing space and time costs.
Symmetric update handling: Efficient processing of both edge insertions and deletions, ensuring robust performance in dynamic graph environments.
Extensive evaluation: Demonstrates significant improvements in speed, memory efficiency, and scalability over existing methods like Timing and CaliG.
