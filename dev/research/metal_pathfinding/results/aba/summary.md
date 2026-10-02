Exploratory observations, not a production performance acceptance result.
Inspect load-*.txt and gpu-load-*.txt; speedup claims require a quiet machine.
The prototype matches supported DAG edge sets, not CPU predecessor order.

All times are medians in milliseconds per batch, with nine timed samples.
CPU A columns use serial SPF for Q=1 and the 10-worker pool otherwise.
GPU includes host distance copies and serial CPU ECMP DAG reconstruction.
Upload and pipeline setup are excluded and recorded separately in raw data.
The table picks the fastest measured GPU method for each case.

| Case | V / E | Q | CPU A1 | GPU B | CPU A2 | GPU method | Distance-only B | CPU drift |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| random | 128 / 1,024 | 1 | 0.031 | 0.242 | 0.010 | metal_shared | 0.239 | -67.4% |
| random | 128 / 1,024 | 64 | 0.234 | 0.527 | 0.162 | metal_shared | 0.291 | -30.7% |
| random | 1,024 / 8,192 | 1 | 0.145 | 0.707 | 0.115 | metal_pull | 0.673 | -21.0% |
| random | 1,024 / 8,192 | 64 | 1.273 | 3.384 | 1.330 | metal_group | 1.054 | +4.5% |
| fabric | 288 / 16,384 | 1 | 0.107 | 0.585 | 0.110 | metal_pull | 0.524 | +3.0% |
| fabric | 288 / 16,384 | 64 | 1.011 | 4.538 | 0.992 | metal_group | 1.189 | -1.9% |
| grid | 1,024 / 3,968 | 1 | 0.037 | 0.374 | 0.045 | metal_shared | 0.359 | +22.8% |
| grid | 1,024 / 3,968 | 64 | 0.498 | 1.724 | 0.556 | metal_group | 0.854 | +11.6% |
| chain | 256 / 255 | 1 | 0.002 | 0.225 | 0.003 | metal_group | 0.224 | +13.0% |
| random_masked_int64_reverse | 1,024 / 8,192 | 8 | 0.150 | 0.523 | 0.149 | metal_pull | 0.135 | -0.2% |
| random | 4,096 / 32,768 | 1 | 0.641 | 0.458 | 0.665 | metal_pull | 0.326 | +3.8% |
| random | 4,096 / 32,768 | 64 | 5.984 | 12.983 | 5.952 | metal_group | 3.123 | -0.5% |
| random | 4,096 / 32,768 | 256 | 23.461 | 47.235 | 23.387 | metal_group | 8.932 | -0.3% |
| random | 16,384 / 131,072 | 1 | 4.095 | 1.529 | 3.915 | metal_pull | 0.892 | -4.4% |
| random | 16,384 / 131,072 | 64 | 28.290 | 54.004 | 28.429 | metal_group | 9.176 | +0.5% |
| random | 65,536 / 524,288 | 1 | 20.762 | 4.922 | 22.345 | metal_pull | 1.921 | +7.6% |
| random | 65,536 / 524,288 | 64 | 174.082 | 277.240 | 189.990 | metal_group | 78.233 | +9.1% |
| chain | 2,048 / 2,047 | 1 | 0.020 | 4.804 | 0.020 | metal_group | 4.793 | -1.1% |
| random | 4,096 / 131,072 | 1 | 1.762 | 0.843 | 1.763 | metal_pull | 0.452 | +0.0% |
| random | 4,096 / 131,072 | 64 | 12.774 | 35.626 | 12.923 | metal_group | 7.957 | +1.2% |
| clos | 200 / 20,000 | 1 | 0.146 | 0.254 | 0.144 | metal_shared | 0.190 | -1.4% |
| clos | 200 / 20,000 | 64 | 1.181 | 5.015 | 1.699 | metal_shared | 0.885 | +43.9% |
| clos | 800 / 320,000 | 1 | 2.330 | 2.590 | 2.368 | metal_shared | 1.607 | +1.6% |
| clos | 800 / 320,000 | 64 | 17.462 | 72.452 | 17.634 | metal_shared | 5.681 | +1.0% |
| grid | 10,000 / 39,600 | 1 | 0.614 | 5.349 | 0.616 | metal_pull | 5.217 | +0.4% |
| grid | 10,000 / 39,600 | 8 | 0.829 | 6.016 | 0.818 | metal_pull | 4.933 | -1.3% |
| grid | 40,000 / 159,200 | 1 | 2.913 | 8.408 | 2.937 | metal_pull | 7.867 | +0.9% |

Distance-only timings stop before DAG reconstruction; they are not SPF API timings.
CPU distance-only and serial full-batch controls are also available in a1/a2.csv.

Distances-only comparison (both sides omit DAG construction):

| Case | V / E | Q | CPU A1 | GPU B | CPU A2 | GPU method |
|---|---:|---:|---:|---:|---:|---|
| random | 128 / 1,024 | 1 | 0.011 | 0.239 | 0.005 | metal_shared |
| random | 128 / 1,024 | 64 | 0.126 | 0.291 | 0.102 | metal_shared |
| random | 1,024 / 8,192 | 1 | 0.080 | 0.673 | 0.056 | metal_pull |
| random | 1,024 / 8,192 | 64 | 0.634 | 1.054 | 0.657 | metal_group |
| fabric | 288 / 16,384 | 1 | 0.022 | 0.524 | 0.022 | metal_pull |
| fabric | 288 / 16,384 | 64 | 0.194 | 1.168 | 0.219 | metal_shared |
| grid | 1,024 / 3,968 | 1 | 0.010 | 0.359 | 0.012 | metal_shared |
| grid | 1,024 / 3,968 | 64 | 0.146 | 0.854 | 0.157 | metal_group |
| chain | 256 / 255 | 1 | 0.001 | 0.224 | 0.001 | metal_group |
| random_masked_int64_reverse | 1,024 / 8,192 | 8 | 0.093 | 0.135 | 0.095 | metal_pull |
| random | 4,096 / 32,768 | 1 | 0.305 | 0.326 | 0.322 | metal_pull |
| random | 4,096 / 32,768 | 64 | 2.968 | 3.123 | 2.946 | metal_group |
| random | 4,096 / 32,768 | 256 | 12.019 | 8.932 | 11.368 | metal_group |
| random | 16,384 / 131,072 | 1 | 2.076 | 0.892 | 1.960 | metal_pull |
| random | 16,384 / 131,072 | 64 | 14.898 | 9.176 | 15.059 | metal_group |
| random | 65,536 / 524,288 | 1 | 10.450 | 1.921 | 11.777 | metal_pull |
| random | 65,536 / 524,288 | 64 | 85.846 | 78.233 | 93.938 | metal_group |
| chain | 2,048 / 2,047 | 1 | 0.004 | 4.793 | 0.004 | metal_group |
| random | 4,096 / 131,072 | 1 | 0.665 | 0.452 | 0.681 | metal_pull |
| random | 4,096 / 131,072 | 64 | 5.610 | 7.957 | 5.638 | metal_group |
| clos | 200 / 20,000 | 1 | 0.024 | 0.190 | 0.024 | metal_shared |
| clos | 200 / 20,000 | 64 | 0.211 | 0.832 | 0.208 | metal_group |
| clos | 800 / 320,000 | 1 | 0.356 | 1.607 | 0.355 | metal_shared |
| clos | 800 / 320,000 | 64 | 2.672 | 5.681 | 2.661 | metal_shared |
| grid | 10,000 / 39,600 | 1 | 0.175 | 5.217 | 0.183 | metal_pull |
| grid | 10,000 / 39,600 | 8 | 0.247 | 4.933 | 0.246 | metal_pull |
| grid | 40,000 / 159,200 | 1 | 0.863 | 7.867 | 0.878 | metal_pull |
