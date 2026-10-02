# Stronger CPU algorithm A/B/A

All methods include exact distances and canonical-equivalent positive-cost DAG edge sets. CPU batches use 10 workers. Eligibility is cached outside timing for both devices. Seven measured samples and two checked warmups. These remain exploratory, non-isolated desktop measurements. GPU choice is an oracle across 12 tested variants.

| Case / N / E / Q | Cached heap + DAG A1 / A2 ms | Dial + DAG A1 / A2 ms | Best GPU B ms | GPU method |
|---|---:|---:|---:|---|
| random / 1024 / 8192 / 1 | 0.053 / 0.059 | 0.028 / 0.030 | 0.634 | pull32_interleaved |
| random / 4096 / 32768 / 256 | 13.030 / 13.024 | 7.056 / 6.925 | 4.001 | pull32_gpu_dag |
| random / 16384 / 131072 / 64 | 17.500 / 17.268 | 7.696 / 7.677 | 3.851 | pull32_gpu_dag |
| random / 65536 / 524288 / 1 | 11.474 / 11.530 | 4.384 / 4.372 | 1.244 | pull32_gpu_dag |
| random / 4096 / 131072 / 64 | 6.672 / 6.662 | 3.854 / 3.865 | 2.215 | pull32_gpu_dag |
| clos / 800 / 320000 / 64 | 5.125 / 5.391 | 5.421 / 5.537 | 8.372 | pull32_gpu_dag |
| grid / 10000 / 39600 / 1 | 0.201 / 0.194 | 0.126 / 0.129 | 3.251 | pull32_gpu_dag |
