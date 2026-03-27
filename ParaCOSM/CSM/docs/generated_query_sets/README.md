# Generated Query Sets

This directory records stable generated query-set locations and the commands used to build and evaluate them.

The current active formal graphflow/newsp experiment dataset is intended to be created at:

- `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20`

Expected contents:

- `accepted_queries/`
- `accepted_summary.csv`
- `generated_queries/`
- `logs/`

The dataset is generated from the stable seed query:

- `/home/v-haibinlai/haibin/paracosm_data/8_self/tree/Q_21`

with the following constraints:

- 8 vertices
- runtime at least 30 seconds
- no timeout under a 120-second limit
- validated with 8 threads during generation

## Formal Graphflow Dataset

The stable 20-query formal dataset was generated successfully at:

- `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20`

It contains:

- `accepted_queries/`
- `accepted_summary.csv`
- `generated_queries/`
- `logs/`

Generation command:

```bash
cd /home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM
/usr/bin/python3 scripts/generate_slow_queries.py \
	--seed-dir /home/v-haibinlai/haibin/paracosm_data/8_self/tree \
	--seed-query-ids 21 \
	--target-count 20 \
	--max-attempts 40 \
	--min-seconds 30 \
	--metric max \
	--threads 8 \
	--time-limit 120 \
	--initial-time-limit 120 \
	--run-timeout 300 \
	--random-seed 20260401 \
	--output-dir /home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20 \
	--keep-all-logs \
	--require-no-timeout
```

## Full Graphflow Experiment

Full batch experiment output directory:

- `/home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20_graphflow_compare_v2`

Batch command:

```bash
cd /home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM
/usr/bin/python3 scripts/compare_csm_algorithms_batch.py \
	--query-dir /home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20/accepted_queries \
	--serial-algorithm graphflow \
	--parallel-algorithm parallel_graphflow \
	--serial-threads 1 \
	--parallel-threads 8 \
	--serial-auto-tuning 0 \
	--parallel-auto-tuning 0 \
	--time-limit 180 \
	--initial-time-limit 180 \
	--run-timeout 900 \
	--output-dir /home/v-haibinlai/haibin/ParaCOSM/ParaCOSM/CSM/formal_graphflow_queryset_20_graphflow_compare_v2
```

Final status summary:

- total = 20
- initial-mismatch = 13
- pass = 2
- timeout-mismatch = 2
- timeout-mismatch-progress = 3

Interpretation:

- `parallel_graphflow` is often faster on these queries
- but it frequently diverges from `graphflow` in initial match count and/or processed edge-update count
- only 2 of the 20 formal queries produced fully aligned statistics in this run
