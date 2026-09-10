# ParaCOSM Dataset Inventory

Inventory date: 2026-09-10
Archive directory: `/home/haibin/tpds/datasets`

The archives were inspected without extracting them. ZIP integrity tests and
the bzip2 integrity test completed successfully.

## Summary

| Dataset | Archive | Compressed size | Vertex records | `Graph::NumVertices()` | Initial edges | Insertions | Queries | Integrity |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Amazon | `amazon.zip` | 17,513,877 B | 403,394 | 403,394 | 2,199,067 | 244,341 | 511 | pass |
| LiveJournal | `livejournal.zip` | 340,509,622 B | 4,846,609 | 4,847,571 | 38,566,113 | 4,285,124 | 511 | pass |
| LSBench | `lsbench.zip` | 100,739,802 B | 5,210,099 | 5,210,099 | 17,982,522 | 2,288,154 | 500 | pass |
| Orkut | `orkut.tar.bz2` | 616,123,674 B | 3,072,441 | 3,072,627 | 117,185,083 | not included | 1,500 | pass |

`Graph::NumVertices()` is the largest vertex ID plus one, because the graph
storage uses vertex IDs as vector indices. LiveJournal and Orkut contain gaps
in their vertex ID spaces, so their vertex-record counts are slightly lower.
This distinction should be preserved in the experimental setup table.

All three supplied update streams contain edge insertions only. No edge
deletions or vertex updates were found. They therefore match the current
insertion-only, fixed-vertex scope of the CPU/GPU versioned executors.

## Archive layout

### Amazon

- Initial graph: `AZ/data_graph/data.graph`
- Update stream: `AZ/data_graph/insertion.graph`
- Queries: `AZ/{6v,7v,8v,9v,10v}/Q_1...Q_100`
- Additional queries: 11 files under `AZ/11v`
- Total uncompressed archive content: 46,110,375 B

### LiveJournal

- Initial graph: `LJ/data_graph/data.graph`
- Update stream: `LJ/data_graph/insertion.graph`
- Queries: `LJ/{6v,7v,8v,9v,10v}/Q_1...Q_100`
- Additional queries: 11 files under `LJ/11v`
- Total uncompressed archive content: 859,450,095 B

### LSBench

- Initial graph: `lsbench_x1/data_graph/data.graph`
- Update stream: `lsbench_x1/data_graph/insertion.graph`
- Queries: `lsbench_x1/random_walk/<size>_self/<shape>/Q_*`
- Total uncompressed archive content: 465,360,160 B

Query distribution:

| Size | Tree | Sparse | Dense | Total |
|---:|---:|---:|---:|---:|
| 6 | 15 | 85 | 0 | 100 |
| 7 | 15 | 85 | 0 | 100 |
| 8 | 15 | 83 | 2 | 100 |
| 9 | 15 | 83 | 2 | 100 |
| 10 | 15 | 80 | 5 | 100 |

### Orkut

- Initial graph: `./data.graph`
- Queries: `./query_graphs/{tree,sparse,dense}_{6,7,8,9,10}/Q_0...Q_99`
- Each size/shape group contains 100 queries (15 groups, 1,500 total)
- Plot images are also present under `query_graphs/plot`
- No update stream was found in the archive

Orkut cannot enter the same E0/E1 update-stream experiment until an insertion
stream is supplied or generated with a documented method and random seed.

## Checksums

```text
028e1e30a54bcc06a349e3f29a03db2532b73eee5b9e0196e485fb3f1edbd859  amazon.zip
13b6f770ee6cec4938685ccc2e8692d33e9c35281102a264874b2efbd3d654d3  livejournal.zip
3827cb7b250a0cc1f4e0dec9945e20da365ac3f74d93d94d2917fcf8ce2e0da0  lsbench.zip
2b0ac91e11926751b7d99d61ee04ecec52ffab81a3923e517c8cce0a8c61f3e1  orkut.tar.bz2
```

## Recommended next action

Extract each archive into a separate directory rather than extracting Orkut at
the dataset root, because its archive members begin with `./data.graph` and
`./query_graphs`:

```bash
mkdir -p amazon livejournal lsbench orkut
unzip amazon.zip -d amazon
unzip livejournal.zip -d livejournal
unzip lsbench.zip -d lsbench
tar -xjf orkut.tar.bz2 -C orkut
```

After extraction, run a small E0 pilot on Amazon 6v before launching the full
correctness matrix. LiveJournal's 4.28 million-update stream should first be
tested using a reproducibly selected prefix so that failures are inexpensive
to diagnose.
