All times are medians in milliseconds per batch, with nine timed samples.
CPU A columns use serial SPF for Q=1 and the 10-worker pool otherwise.
GPU includes host distance copies and serial CPU ECMP DAG reconstruction.
Upload and pipeline setup are excluded and recorded separately in raw data.
The table picks the fastest measured GPU method for each case.

| Case | V / E | Q | CPU A1 | GPU B | CPU A2 | GPU method | Distance-only B | CPU drift |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| random | 128 / 1,024 | 1 | 0.010 | 0.265 | 0.011 | metal_shared | 0.260 | +12.6% |
| random | 128 / 1,024 | 64 | 0.182 | 0.598 | 0.157 | metal_group | 0.323 | -14.1% |
| random | 1,024 / 8,192 | 1 | 0.104 | 0.752 | 0.109 | metal_pull | 0.709 | +5.1% |
| random | 1,024 / 8,192 | 64 | 1.342 | 3.774 | 1.277 | metal_shared | 1.285 | -4.9% |
| fabric | 288 / 16,384 | 1 | 0.112 | 0.733 | 0.110 | metal_shared | 0.674 | -2.5% |
| fabric | 288 / 16,384 | 64 | 0.972 | 4.629 | 0.995 | metal_group | 1.190 | +2.3% |
| grid | 1,024 / 3,968 | 1 | 0.044 | 0.488 | 0.043 | metal_shared | 0.471 | -2.2% |
| grid | 1,024 / 3,968 | 64 | 0.543 | 1.864 | 0.546 | metal_group | 0.949 | +0.5% |
| chain | 256 / 255 | 1 | 0.003 | 0.659 | 0.003 | metal_shared | 0.657 | -0.0% |
| random_masked_int64_reverse | 1,024 / 8,192 | 8 | 0.156 | 0.682 | 0.203 | metal_pull | 0.246 | +30.6% |
| random | 4,096 / 32,768 | 1 | 0.657 | 0.813 | 0.713 | metal_pull | 0.661 | +8.6% |
| random | 4,096 / 32,768 | 64 | 5.992 | 14.066 | 5.860 | metal_group | 4.463 | -2.2% |
| random | 4,096 / 32,768 | 256 | 23.001 | 49.850 | 22.651 | metal_group | 10.486 | -1.5% |
| random | 16,384 / 131,072 | 1 | 3.887 | 1.386 | 3.929 | metal_pull | 0.696 | +1.1% |
| random | 16,384 / 131,072 | 64 | 28.586 | 54.493 | 28.359 | metal_group | 10.271 | -0.8% |
| random | 65,536 / 524,288 | 1 | 22.864 | 5.859 | 22.647 | metal_pull | 2.820 | -0.9% |
| random | 65,536 / 524,288 | 64 | 181.874 | 269.300 | 190.623 | metal_pull | 71.201 | +4.8% |
| chain | 2,048 / 2,047 | 1 | 0.020 | 4.816 | 0.019 | metal_group | 4.806 | -1.3% |
| random | 4,096 / 131,072 | 1 | 1.540 | 0.897 | 1.680 | metal_pull | 0.484 | +9.1% |
| random | 4,096 / 131,072 | 64 | 12.857 | 36.557 | 12.771 | metal_group | 7.881 | -0.7% |
| clos | 200 / 20,000 | 1 | 0.139 | 0.353 | 0.144 | metal_shared | 0.281 | +4.1% |
| clos | 200 / 20,000 | 64 | 1.167 | 5.050 | 1.280 | metal_group | 0.890 | +9.7% |
| clos | 800 / 320,000 | 1 | 2.314 | 2.640 | 2.457 | metal_shared | 1.565 | +6.2% |
| clos | 800 / 320,000 | 64 | 17.328 | 73.488 | 17.443 | metal_shared | 5.716 | +0.7% |
| grid | 10,000 / 39,600 | 1 | 0.631 | 4.067 | 0.617 | metal_pull | 3.935 | -2.2% |
| grid | 10,000 / 39,600 | 8 | 0.852 | 6.016 | 0.874 | metal_pull | 5.013 | +2.6% |
| grid | 40,000 / 159,200 | 1 | 2.919 | 8.366 | 2.918 | metal_pull | 7.871 | -0.0% |

Distance-only timings stop before DAG reconstruction; they are not SPF API timings.
CPU distance-only and serial full-batch controls are also available in a1/a2.csv.
