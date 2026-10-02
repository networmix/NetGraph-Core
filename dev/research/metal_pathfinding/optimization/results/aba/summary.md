# Exploratory A/B/A measurements

Warm full-result latency in ms, including host distances and positive-cost predecessor DAG. Seven measured samples plus two checked warmups per row. CPU batches use 10 reusable workers. Best variant is selected after measurement; this is an oracle comparison, not an implemented dispatcher. Old GPU columns cover the original pull and group variants; the shared-memory variant was not rerun. The machine was not isolated; these are not accepted performance claims.

| Case | N / E / Q | CPU full A1 / A2 | Old GPU A1 / A2 | Best tested optimized GPU | ms | CPU BFS full A1 / A2 |
|---|---:|---:|---:|---|---:|---:|
| random | 1024 / 8192 / 1 | 0.086 / 0.097 | 0.915 / 0.972 | pull32_interleaved | 0.599 | — |
| random | 1024 / 8192 / 64 | 1.248 / 1.340 | 3.818 / 3.769 | old_group_parallel_dag | 0.853 | — |
| clos | 200 / 20000 / 64 | 1.072 / 1.076 | 4.768 / 4.870 | old_group_parallel_dag | 1.461 | 0.669 / 0.633 |
| grid | 1024 / 3968 / 1 | 0.039 / 0.039 | 1.177 / 1.166 | old_group_parallel_dag | 0.366 | 0.009 / 0.009 |
| hub | 1024 / 3070 / 1 | 0.052 / 0.052 | 4.908 / 4.938 | frontier32_simd | 0.784 | — |
| random_masked_int64_reverse | 1024 / 8192 / 8 | 0.160 / 0.153 | 0.518 / 0.551 | bucket32_d16 (64-bit pull fallback) | 0.224 | — |
| random | 4096 / 32768 / 256 | 21.391 / 21.281 | 47.812 / 49.310 | pull64_interleaved | 7.003 | — |
| random | 16384 / 131072 / 1 | 3.662 / 3.886 | 2.098 / 1.870 | pull32_interleaved | 0.705 | — |
| random | 16384 / 131072 / 64 | 27.693 / 27.721 | 57.748 / 53.106 | pull64_interleaved | 7.473 | — |
| random | 65536 / 524288 / 1 | 19.316 / 19.240 | 5.184 / 5.854 | pull32_gpu_dag | 1.373 | — |
| clos | 800 / 320000 / 1 | 2.178 / 2.227 | 1.915 / 2.377 | frontier32_simd | 1.604 | 0.908 / 0.938 |
| clos | 800 / 320000 / 64 | 17.042 / 16.525 | 69.896 / 70.228 | pull32_gpu_dag | 8.185 | 8.296 / 7.885 |
| grid | 10000 / 39600 / 1 | 0.570 / 0.582 | 3.836 / 5.874 | pull32_interleaved | 3.378 | 0.089 / 0.088 |
| grid | 10000 / 39600 / 8 | 0.834 / 0.823 | 5.953 / 5.979 | frontier32_gpu_dag | 4.150 | 0.286 / 0.312 |
| chain | 2048 / 2047 / 1 | 0.017 / 0.017 | 4.834 / 4.908 | old_group_parallel_dag | 4.847 | — |
| hub | 4096 / 12286 / 64 | 2.438 / 2.338 | 28.559 / 28.801 | frontier32_simd | 3.680 | — |
| random | 4096 / 131072 / 1 | 1.590 / 1.518 | 0.906 / 1.040 | frontier32_simd | 0.547 | — |
| random | 4096 / 131072 / 64 | 12.302 / 12.817 | 31.930 / 37.791 | pull32_gpu_dag | 2.051 | — |

## All optimized variants

No variants omitted. Columns: full wall / distance wall / device work, ms. Device work includes DAG kernels where selected. Edge visits count SSSP adjacency entries, excluding bucket vertex scans, queue management and DAG extraction.

### ('random', '1024', '8192', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 1.423 / 1.376 / 0.424 | 131072 | 16 | 2 | 0 |
| old_group_parallel_dag | 0.834 / 0.819 / 0.680 | 106496 | 13 | 1 | 0 |
| pull64_barrier | 0.659 / 0.646 / 0.391 | 131072 | 16 | 2 | 0 |
| pull64_interleaved | 0.656 / 0.643 / 0.375 | 131072 | 16 | 2 | 0 |
| pull32_interleaved | 0.599 / 0.586 / 0.338 | 131072 | 16 | 2 | 0 |
| pull32_gpu_dag | 0.960 / 0.604 / 0.406 | 131072 | 16 | 4 | 0 |
| frontier32_vertex | 0.924 / 0.910 / 0.634 | 21088 | 16 | 2 | 0 |
| frontier32_simd | 0.661 / 0.647 / 0.380 | 20272 | 16 | 2 | 0 |
| frontier32_gpu_dag | 1.003 / 0.652 / 0.452 | 20384 | 16 | 4 | 0 |
| bucket32_d4 | 1.895 / 1.881 / 1.276 | 8232 | 32 | 4 | 0 |
| bucket32_d16 | 1.422 / 1.407 / 0.920 | 9336 | 24 | 3 | 0 |
| bucket32_d64 | 0.954 / 0.939 / 0.654 | 19264 | 16 | 2 | 0 |
### ('random', '1024', '8192', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 1.973 / 1.703 / 1.236 | 12582912 | 24 | 3 | 0 |
| old_group_parallel_dag | 0.853 / 0.566 / 0.380 | 7217152 | 17 | 1 | 0 |
| pull64_barrier | 1.223 / 0.949 / 0.470 | 12582912 | 24 | 3 | 0 |
| pull64_interleaved | 0.940 / 0.671 / 0.268 | 12582912 | 24 | 3 | 0 |
| pull32_interleaved | 0.863 / 0.595 / 0.197 | 12582912 | 24 | 3 | 0 |
| pull32_gpu_dag | 0.998 / 0.579 / 0.199 | 12582912 | 24 | 5 | 0 |
| frontier32_vertex | 0.870 / 0.629 / 0.239 | 1289824 | 16 | 2 | 0 |
| frontier32_simd | 1.028 / 0.796 / 0.456 | 1202328 | 16 | 2 | 0 |
| frontier32_gpu_dag | 1.197 / 0.890 / 0.509 | 1201208 | 16 | 4 | 0 |
| bucket32_d4 | 5.929 / 5.638 / 4.211 | 528528 | 56 | 7 | 0 |
| bucket32_d16 | 1.624 / 1.347 / 0.712 | 601456 | 32 | 4 | 0 |
| bucket32_d64 | 1.403 / 1.126 / 0.671 | 1154720 | 24 | 3 | 0 |
### ('clos', '200', '20000', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 2.137 / 1.514 / 1.016 | 10240000 | 8 | 1 | 0 |
| old_group_parallel_dag | 1.461 / 0.833 / 0.303 | 3840000 | 3 | 1 | 0 |
| pull64_barrier | 2.072 / 1.524 / 1.000 | 10240000 | 8 | 1 | 0 |
| pull64_interleaved | 1.903 / 1.332 / 0.757 | 10240000 | 8 | 1 | 0 |
| pull32_interleaved | 1.893 / 1.311 / 0.688 | 10240000 | 8 | 1 | 0 |
| pull32_gpu_dag | 2.092 / 1.226 / 1.013 | 10240000 | 8 | 3 | 0 |
| frontier32_vertex | 2.382 / 1.557 / 0.961 | 1280000 | 8 | 1 | 0 |
| frontier32_simd | 1.552 / 0.924 / 0.318 | 1280000 | 8 | 1 | 0 |
| frontier32_gpu_dag | 1.979 / 0.982 / 0.837 | 1280000 | 8 | 3 | 0 |
| bucket32_d4 | 1.572 / 1.025 / 0.481 | 1280000 | 8 | 1 | 0 |
| bucket32_d16 | 1.739 / 1.166 / 0.492 | 1280000 | 8 | 1 | 0 |
| bucket32_d64 | 1.693 / 1.072 / 0.490 | 1280000 | 8 | 1 | 0 |
### ('grid', '1024', '3968', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 1.747 / 1.739 / 0.651 | 253952 | 64 | 8 | 0 |
| old_group_parallel_dag | 0.366 / 0.359 / 0.227 | 249984 | 63 | 1 | 0 |
| pull64_barrier | 1.119 / 1.111 / 0.177 | 253952 | 64 | 8 | 0 |
| pull64_interleaved | 1.144 / 1.136 / 0.177 | 253952 | 64 | 8 | 0 |
| pull32_interleaved | 1.089 / 1.081 / 0.164 | 253952 | 64 | 8 | 0 |
| pull32_gpu_dag | 1.271 / 0.918 / 0.173 | 253952 | 64 | 10 | 0 |
| frontier32_vertex | 1.313 / 1.305 / 0.528 | 3968 | 64 | 8 | 0 |
| frontier32_simd | 1.180 / 1.172 / 0.396 | 3968 | 64 | 8 | 0 |
| frontier32_gpu_dag | 1.394 / 1.178 / 0.407 | 3968 | 64 | 10 | 0 |
| bucket32_d4 | 1.450 / 1.441 / 0.637 | 3968 | 64 | 8 | 0 |
| bucket32_d16 | 1.451 / 1.443 / 0.636 | 3968 | 64 | 8 | 0 |
| bucket32_d64 | 1.466 / 1.459 / 0.637 | 3968 | 64 | 8 | 0 |
### ('hub', '1024', '3070', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 6.046 / 6.030 / 5.544 | 98240 | 32 | 4 | 0 |
| old_group_parallel_dag | 4.888 / 4.866 / 4.721 | 95170 | 31 | 1 | 0 |
| pull64_barrier | 5.016 / 4.988 / 4.313 | 98240 | 32 | 4 | 0 |
| pull64_interleaved | 4.986 / 4.962 / 4.321 | 98240 | 32 | 4 | 0 |
| pull32_interleaved | 4.397 / 4.374 / 3.684 | 98240 | 32 | 4 | 0 |
| pull32_gpu_dag | 5.116 / 4.358 / 4.101 | 98240 | 32 | 6 | 0 |
| frontier32_vertex | 1.436 / 1.428 / 0.852 | 3128 | 32 | 4 | 0 |
| frontier32_simd | 0.784 / 0.776 / 0.233 | 3128 | 32 | 4 | 0 |
| frontier32_gpu_dag | 1.467 / 0.770 / 0.649 | 3128 | 32 | 6 | 0 |
| bucket32_d4 | 0.869 / 0.860 / 0.345 | 3100 | 32 | 4 | 0 |
| bucket32_d16 | 0.851 / 0.843 / 0.337 | 3124 | 32 | 4 | 0 |
| bucket32_d64 | 0.886 / 0.878 / 0.327 | 3128 | 32 | 4 | 0 |
### ('random_masked_int64_reverse', '1024', '8192', '8')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 0.232 / 0.167 / 0.032 | 524288 | 8 | 1 | 0 |
| old_group_parallel_dag | 0.253 / 0.194 / 0.053 | 368640 | 8 | 1 | 0 |
| pull64_barrier | 0.248 / 0.190 / 0.029 | 524288 | 8 | 1 | 0 |
| pull64_interleaved | 0.264 / 0.197 / 0.029 | 524288 | 8 | 1 | 0 |
| pull32_interleaved | 0.252 / 0.198 / 0.029 | 524288 | 8 | 1 | 1 |
| pull32_gpu_dag | 0.499 / 0.174 / 0.049 | 524288 | 8 | 3 | 1 |
| frontier32_vertex | 0.238 / 0.188 / 0.029 | 524288 | 8 | 1 | 1 |
| frontier32_simd | 0.256 / 0.192 / 0.029 | 524288 | 8 | 1 | 1 |
| frontier32_gpu_dag | 0.504 / 0.184 / 0.047 | 524288 | 8 | 3 | 1 |
| bucket32_d4 | 0.243 / 0.198 / 0.029 | 524288 | 8 | 1 | 1 |
| bucket32_d16 | 0.224 / 0.178 / 0.029 | 524288 | 8 | 1 | 1 |
| bucket32_d64 | 0.261 / 0.199 / 0.029 | 524288 | 8 | 1 | 1 |
### ('random', '4096', '32768', '256')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 14.955 / 11.144 / 9.933 | 201326592 | 24 | 3 | 0 |
| old_group_parallel_dag | 11.270 / 7.608 / 6.915 | 136544256 | 20 | 1 | 0 |
| pull64_barrier | 13.565 / 9.762 / 8.501 | 201326592 | 24 | 3 | 0 |
| pull64_interleaved | 7.003 / 3.340 / 2.081 | 201326592 | 24 | 3 | 0 |
| pull32_interleaved | 9.000 / 5.200 / 3.525 | 201326592 | 24 | 3 | 0 |
| pull32_gpu_dag | 8.790 / 6.950 / 7.134 | 201326592 | 24 | 5 | 0 |
| frontier32_vertex | 7.936 / 4.094 / 2.380 | 19766184 | 24 | 3 | 0 |
| frontier32_simd | 11.324 / 7.659 / 6.414 | 20599968 | 24 | 3 | 0 |
| frontier32_gpu_dag | 9.247 / 7.874 / 7.025 | 20610776 | 24 | 5 | 0 |
| bucket32_d4 | 16.328 / 12.614 / 10.131 | 8454072 | 72 | 9 | 0 |
| bucket32_d16 | 13.121 / 9.301 / 7.673 | 9728592 | 48 | 6 | 0 |
| bucket32_d64 | 14.196 / 9.623 / 8.264 | 18345256 | 24 | 3 | 0 |
### ('random', '16384', '131072', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 0.964 / 0.659 / 0.244 | 3145728 | 24 | 3 | 0 |
| old_group_parallel_dag | 5.290 / 4.968 / 4.819 | 2359296 | 18 | 1 | 0 |
| pull64_barrier | 0.761 / 0.571 / 0.164 | 3145728 | 24 | 3 | 0 |
| pull64_interleaved | 0.828 / 0.564 / 0.164 | 3145728 | 24 | 3 | 0 |
| pull32_interleaved | 0.705 / 0.521 / 0.137 | 3145728 | 24 | 3 | 0 |
| pull32_gpu_dag | 0.880 / 0.513 / 0.158 | 3145728 | 24 | 5 | 0 |
| frontier32_vertex | 0.832 / 0.619 / 0.258 | 377056 | 24 | 3 | 0 |
| frontier32_simd | 0.806 / 0.593 / 0.237 | 335680 | 24 | 3 | 0 |
| frontier32_gpu_dag | 0.983 / 0.612 / 0.269 | 334824 | 24 | 5 | 0 |
| bucket32_d4 | 1.833 / 1.546 / 0.654 | 132016 | 56 | 7 | 0 |
| bucket32_d16 | 1.148 / 0.917 / 0.398 | 148256 | 32 | 4 | 0 |
| bucket32_d64 | 0.937 / 0.712 / 0.340 | 268728 | 24 | 3 | 0 |
### ('random', '16384', '131072', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 16.913 / 12.651 / 11.412 | 201326592 | 24 | 3 | 0 |
| old_group_parallel_dag | 13.239 / 9.052 / 8.276 | 158203904 | 22 | 1 | 0 |
| pull64_barrier | 15.959 / 11.386 / 10.305 | 201326592 | 24 | 3 | 0 |
| pull64_interleaved | 7.473 / 3.322 / 1.950 | 201326592 | 24 | 3 | 0 |
| pull32_interleaved | 8.710 / 3.948 / 2.593 | 201326592 | 24 | 3 | 0 |
| pull32_gpu_dag | 8.810 / 7.334 / 7.227 | 201326592 | 24 | 5 | 0 |
| frontier32_vertex | 10.144 / 5.760 / 3.640 | 21296680 | 24 | 3 | 0 |
| frontier32_simd | 12.189 / 7.941 / 6.799 | 22006672 | 24 | 3 | 0 |
| frontier32_gpu_dag | 9.203 / 7.809 / 7.471 | 21985376 | 24 | 5 | 0 |
| bucket32_d4 | 16.825 / 12.607 / 10.280 | 8451400 | 72 | 9 | 0 |
| bucket32_d16 | 13.809 / 9.687 / 7.926 | 9572312 | 48 | 6 | 0 |
| bucket32_d64 | 12.602 / 8.439 / 7.400 | 17413440 | 24 | 3 | 0 |
### ('random', '65536', '524288', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 3.545 / 2.038 / 1.061 | 12582912 | 24 | 3 | 0 |
| old_group_parallel_dag | 30.358 / 28.978 / 28.314 | 11534336 | 22 | 1 | 0 |
| pull64_barrier | 3.163 / 1.729 / 0.758 | 12582912 | 24 | 3 | 0 |
| pull64_interleaved | 3.048 / 1.519 / 0.703 | 12582912 | 24 | 3 | 0 |
| pull32_interleaved | 2.913 / 1.338 / 0.502 | 12582912 | 24 | 3 | 0 |
| pull32_gpu_dag | 1.373 / 0.882 / 0.571 | 12582912 | 24 | 5 | 0 |
| frontier32_vertex | 2.737 / 1.259 / 0.435 | 1734632 | 24 | 3 | 0 |
| frontier32_simd | 3.045 / 1.576 / 0.640 | 1561768 | 24 | 3 | 0 |
| frontier32_gpu_dag | 1.461 / 1.045 / 0.686 | 1558160 | 24 | 5 | 0 |
| bucket32_d4 | 4.357 / 2.761 / 1.174 | 528256 | 64 | 8 | 0 |
| bucket32_d16 | 3.481 / 2.047 / 0.836 | 599256 | 40 | 5 | 0 |
| bucket32_d64 | 3.187 / 1.681 / 0.677 | 1050376 | 24 | 3 | 0 |
### ('clos', '800', '320000', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 2.142 / 1.366 / 0.744 | 2560000 | 8 | 1 | 0 |
| old_group_parallel_dag | 2.326 / 1.577 / 0.905 | 960000 | 3 | 1 | 0 |
| pull64_barrier | 2.017 / 1.253 / 0.699 | 2560000 | 8 | 1 | 0 |
| pull64_interleaved | 2.265 / 1.493 / 0.878 | 2560000 | 8 | 1 | 0 |
| pull32_interleaved | 2.801 / 1.980 / 1.363 | 2560000 | 8 | 1 | 0 |
| pull32_gpu_dag | 2.133 / 1.256 / 0.994 | 2560000 | 8 | 3 | 0 |
| frontier32_vertex | 2.376 / 1.575 / 0.933 | 320000 | 8 | 1 | 0 |
| frontier32_simd | 1.604 / 0.852 / 0.361 | 320000 | 8 | 1 | 0 |
| frontier32_gpu_dag | 1.714 / 0.824 / 0.606 | 320000 | 8 | 3 | 0 |
| bucket32_d4 | 1.660 / 0.873 / 0.446 | 320000 | 8 | 1 | 0 |
| bucket32_d16 | 1.716 / 0.904 / 0.444 | 320000 | 8 | 1 | 0 |
| bucket32_d64 | 1.698 / 0.939 / 0.453 | 320000 | 8 | 1 | 0 |
### ('clos', '800', '320000', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 26.735 / 19.418 / 18.854 | 163840000 | 8 | 1 | 0 |
| old_group_parallel_dag | 12.637 / 5.915 / 5.412 | 61440000 | 3 | 1 | 0 |
| pull64_barrier | 24.737 / 18.165 / 17.725 | 163840000 | 8 | 1 | 0 |
| pull64_interleaved | 13.270 / 6.122 / 5.405 | 163840000 | 8 | 1 | 0 |
| pull32_interleaved | 10.983 / 4.274 / 3.647 | 163840000 | 8 | 1 | 0 |
| pull32_gpu_dag | 8.185 / 4.308 / 6.164 | 163840000 | 8 | 3 | 0 |
| frontier32_vertex | 14.231 / 7.148 / 6.877 | 20480000 | 8 | 1 | 0 |
| frontier32_simd | 9.324 / 2.065 / 1.458 | 20480000 | 8 | 1 | 0 |
| frontier32_gpu_dag | 15.944 / 2.089 / 14.063 | 20480000 | 8 | 3 | 0 |
| bucket32_d4 | 10.000 / 2.529 / 2.040 | 20480000 | 8 | 1 | 0 |
| bucket32_d16 | 9.610 / 2.183 / 1.863 | 20480000 | 8 | 1 | 0 |
| bucket32_d64 | 9.988 / 2.684 / 2.106 | 20480000 | 8 | 1 | 0 |
### ('grid', '10000', '39600', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 4.013 / 3.947 / 0.706 | 7920000 | 200 | 25 | 0 |
| old_group_parallel_dag | 8.014 / 7.948 / 7.804 | 7880400 | 199 | 1 | 0 |
| pull64_barrier | 4.120 / 3.996 / 0.682 | 7920000 | 200 | 25 | 0 |
| pull64_interleaved | 3.963 / 3.894 / 0.656 | 7920000 | 200 | 25 | 0 |
| pull32_interleaved | 3.378 / 3.314 / 0.604 | 7920000 | 200 | 25 | 0 |
| pull32_gpu_dag | 3.492 / 3.086 / 0.616 | 7920000 | 200 | 27 | 0 |
| frontier32_vertex | 4.642 / 4.568 / 1.859 | 39600 | 200 | 25 | 0 |
| frontier32_simd | 4.040 / 3.931 / 1.337 | 39600 | 200 | 25 | 0 |
| frontier32_gpu_dag | 4.283 / 3.993 / 1.302 | 39600 | 200 | 27 | 0 |
| bucket32_d4 | 4.910 / 4.826 / 2.156 | 39600 | 200 | 25 | 0 |
| bucket32_d16 | 5.336 / 5.270 / 2.375 | 39600 | 200 | 25 | 0 |
| bucket32_d64 | 5.193 / 5.129 / 2.313 | 39600 | 200 | 25 | 0 |
### ('grid', '10000', '39600', '8')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 4.347 / 4.145 / 1.436 | 63360000 | 200 | 25 | 0 |
| old_group_parallel_dag | 8.292 / 8.082 / 7.572 | 51480000 | 199 | 1 | 0 |
| pull64_barrier | 4.480 / 4.280 / 1.391 | 63360000 | 200 | 25 | 0 |
| pull64_interleaved | 4.328 / 4.167 / 1.330 | 63360000 | 200 | 25 | 0 |
| pull32_interleaved | 4.287 / 4.032 / 1.368 | 63360000 | 200 | 25 | 0 |
| pull32_gpu_dag | 4.436 / 4.109 / 1.562 | 63360000 | 200 | 27 | 0 |
| frontier32_vertex | 5.906 / 5.657 / 2.502 | 316800 | 200 | 25 | 0 |
| frontier32_simd | 4.526 / 4.257 / 1.461 | 316800 | 200 | 25 | 0 |
| frontier32_gpu_dag | 4.150 / 3.853 / 1.370 | 316800 | 200 | 27 | 0 |
| bucket32_d4 | 7.676 / 7.433 / 3.707 | 316800 | 200 | 25 | 0 |
| bucket32_d16 | 7.596 / 7.402 / 3.678 | 316800 | 200 | 25 | 0 |
| bucket32_d64 | 7.603 / 7.391 / 3.721 | 316800 | 200 | 25 | 0 |
### ('chain', '2048', '2047', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 30.065 / 30.058 / 4.364 | 4192256 | 2048 | 256 | 0 |
| old_group_parallel_dag | 4.847 / 4.825 / 4.696 | 4192256 | 2048 | 1 | 0 |
| pull64_barrier | 28.539 / 28.532 / 4.712 | 4192256 | 2048 | 256 | 0 |
| pull64_interleaved | 32.552 / 32.545 / 6.166 | 4192256 | 2048 | 256 | 0 |
| pull32_interleaved | 33.313 / 33.304 / 7.668 | 4192256 | 2048 | 256 | 0 |
| pull32_gpu_dag | 30.619 / 30.342 / 4.891 | 4192256 | 2048 | 258 | 0 |
| frontier32_vertex | 37.679 / 37.671 / 12.262 | 2047 | 2048 | 256 | 0 |
| frontier32_simd | 41.050 / 41.041 / 15.764 | 2047 | 2048 | 256 | 0 |
| frontier32_gpu_dag | 40.148 / 39.758 / 13.087 | 2047 | 2048 | 258 | 0 |
| bucket32_d4 | 49.001 / 48.994 / 21.208 | 2047 | 2056 | 257 | 0 |
| bucket32_d16 | 48.755 / 48.747 / 21.149 | 2047 | 2056 | 257 | 0 |
| bucket32_d64 | 48.401 / 48.395 / 21.022 | 2047 | 2056 | 257 | 0 |
### ('hub', '4096', '12286', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 32.652 / 32.266 / 30.925 | 31452160 | 40 | 5 | 0 |
| old_group_parallel_dag | 26.147 / 25.761 / 25.032 | 27176632 | 38 | 1 | 0 |
| pull64_barrier | 25.410 / 25.029 / 23.756 | 31452160 | 40 | 5 | 0 |
| pull64_interleaved | 33.830 / 33.457 / 31.908 | 31452160 | 40 | 5 | 0 |
| pull32_interleaved | 29.085 / 28.621 / 27.390 | 31452160 | 40 | 5 | 0 |
| pull32_gpu_dag | 31.331 / 28.778 / 29.382 | 31452160 | 40 | 7 | 0 |
| frontier32_vertex | 20.164 / 19.779 / 18.439 | 1847854 | 40 | 5 | 0 |
| frontier32_simd | 3.680 / 3.294 / 1.977 | 1683322 | 40 | 5 | 0 |
| frontier32_gpu_dag | 5.438 / 3.358 / 3.575 | 1683332 | 40 | 7 | 0 |
| bucket32_d4 | 5.120 / 4.740 / 3.072 | 788342 | 48 | 6 | 0 |
| bucket32_d16 | 4.911 / 4.544 / 3.014 | 1124798 | 48 | 6 | 0 |
| bucket32_d64 | 4.974 / 4.580 / 3.351 | 1929942 | 40 | 5 | 0 |
### ('random', '4096', '131072', '1')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 0.700 / 0.521 / 0.269 | 2097152 | 16 | 2 | 0 |
| old_group_parallel_dag | 2.531 / 2.326 / 2.178 | 1572864 | 12 | 1 | 0 |
| pull64_barrier | 0.670 / 0.484 / 0.244 | 2097152 | 16 | 2 | 0 |
| pull64_interleaved | 0.708 / 0.538 / 0.266 | 2097152 | 16 | 2 | 0 |
| pull32_interleaved | 0.603 / 0.418 / 0.180 | 2097152 | 16 | 2 | 0 |
| pull32_gpu_dag | 0.828 / 0.458 / 0.237 | 2097152 | 16 | 4 | 0 |
| frontier32_vertex | 0.827 / 0.638 / 0.358 | 421728 | 16 | 2 | 0 |
| frontier32_simd | 0.547 / 0.381 / 0.130 | 385248 | 16 | 2 | 0 |
| frontier32_gpu_dag | 0.790 / 0.404 / 0.188 | 385888 | 16 | 4 | 0 |
| bucket32_d4 | 0.878 / 0.693 / 0.279 | 138720 | 24 | 3 | 0 |
| bucket32_d16 | 0.602 / 0.418 / 0.176 | 177024 | 16 | 2 | 0 |
| bucket32_d64 | 0.676 / 0.489 / 0.213 | 379328 | 16 | 2 | 0 |
### ('random', '4096', '131072', '64')

| Variant | Full / distance / device ms | Edge visits | Rounds | Host waits | 64-bit fallback |
|---|---:|---:|---:|---:|---|
| old_pull_parallel_dag | 9.211 / 7.184 / 6.254 | 134217728 | 16 | 2 | 0 |
| old_group_parallel_dag | 6.357 / 3.995 / 3.300 | 105644032 | 15 | 1 | 0 |
| pull64_barrier | 11.066 / 9.020 / 8.132 | 134217728 | 16 | 2 | 0 |
| pull64_interleaved | 4.039 / 1.981 / 1.106 | 134217728 | 16 | 2 | 0 |
| pull32_interleaved | 3.814 / 1.788 / 1.031 | 134217728 | 16 | 2 | 0 |
| pull32_gpu_dag | 2.051 / 1.482 / 0.861 | 134217728 | 16 | 4 | 0 |
| frontier32_vertex | 5.418 / 3.315 / 2.452 | 22442464 | 16 | 2 | 0 |
| frontier32_simd | 5.078 / 2.991 / 1.970 | 22729248 | 16 | 2 | 0 |
| frontier32_gpu_dag | 4.086 / 2.828 / 2.683 | 22826080 | 16 | 4 | 0 |
| bucket32_d4 | 4.731 / 2.704 / 1.615 | 8805184 | 24 | 3 | 0 |
| bucket32_d16 | 4.926 / 2.956 / 1.859 | 12113536 | 24 | 3 | 0 |
| bucket32_d64 | 5.330 / 3.310 / 2.396 | 23278880 | 16 | 2 | 0 |
