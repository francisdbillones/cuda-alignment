YouTube link for video pitch: https://youtu.be/tV8orLYLfi4

YouTube link for midterm update: https://youtu.be/7j1b1YhDKds

Midterm update PDF: https://github.com/francisdbillones/cuda-alignment/blob/main/proposal.pdf


## Building


```sh
cmake -S . -B build
cmake --build build --config Release
```

## Baseline benchmark (original TARA CUDA code)

### Setup

- **Hardware:** NVIDIA Tesla T4 (driver 580.82.07), 2 vCPU Intel Xeon @ 2.00 GHz
- **Toolchain:** CUDA 13.0, CMake 3.31, Ninja; `-DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75`
  - **FF:** `Fruit_Fly_Protein_Final` vs. `Fruit_Fly_DNA_Final` read #10 (6,531 bp)
  - **MC:** `Mouse_Crest_Protein_Final` vs. `Mouse_Crest_DNA_Final` read #100 (3,249 bp)

### Results

| # | Dataset | Frames | Mode       | Kernel time | Wall time |
|---|---------|--------|------------|------------:|----------:|
| 1 | FF      | 3      | CUDA       |      59.7 s |     535 s |
| 2 | FF      | 3      | Sequential |     439.1 s |     495 s |
| 3 | FF      | 6      | CUDA       |      83.5 s |     748 s |
| 4 | FF      | 6      | Sequential |     914.7 s |   1,070 s |
| 5 | MC      | 3      | CUDA       |      27.7 s |     203 s |
| 6 | MC      | 3      | Sequential |     165.1 s |     222 s |
| 7 | MC      | 6      | CUDA       |      40.8 s |     285 s |
| 8 | MC      | 6      | Sequential |     330.2 s |     450 s |

### Minimum to beat

| Dataset | Frames | Kernel time | Wall time              |
|---------|--------|------------:|-----------------------:|
| FF      | 3      |      59.7 s | 495 s (Sequential beats CUDA) |
| FF      | 6      |      83.5 s | 748 s                  |
| MC      | 3      |      27.7 s | 203 s                  |
| MC      | 6      |      40.8 s | 285 s                  |