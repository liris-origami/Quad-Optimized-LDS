# Quad-Optimized-LDS
[![CMake on multiple platforms](https://github.com/liris-origami/Quad-Optimized-LDS/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/liris-origami/Quad-Optimized-LDS/actions/workflows/c-cpp.yml)

Source code of the paper: "Quad-Optimized Low-Discrepancy Sequences",  Victor Ostromoukhov, Nicolas
Bonneel, David Coeurjolly, Jean-Claude Iehl , ACM SIGGRAPH Conference Proceeedings, San Diego, 2024.

<img width="991" alt="Capture d’écran 2024-04-22 à 16 07 52" src="https://github.com/liris-origami/Quad-Optimized-LDS/assets/700165/6397a31a-61fc-44f0-ae1a-b31ace45a8f5">

# Build Instructions

The project uses a classical [cmake](http://cmake.org) worklow. For instance on linux/macos systems:
```
mkdir build
cd build
cmake ..
make
```

# Generating Samples

The main tool is `samplerGF3`:
```
Usage: samplerGF3 [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -d INT                      number of dimensions, default: 8
  -m INT                      size of the matrices, to produce up to 3^m points, default: 5
  -i,--init_file TEXT REQUIRED
                              filename which contains initialization data, default:
```
 
For instance, to generate matrices, and the  $3^5$ Owen scrambled samples (on 5 digits) in dimension 4:

```
samplerGF3 -i ../initIrreducibleGF3.dat -m 5 -d 4
```

which outputs the four matrices in GF(3) and the scrambled samples (cropped output):

```
d= 0
1 0 0 0 0
0 1 0 0 0
0 0 1 0 0
0 0 0 1 0
0 0 0 0 1

d= 1
1 1 2 2 1
0 1 0 2 0
0 0 1 1 1
0 0 0 1 0
0 0 0 0 1

d= 2
1 2 1 2 1
0 1 1 0 2
0 0 1 0 0
0 0 0 1 2
0 0 0 0 1

d= 3
2 2 2 2 2
0 2 1 0 2
0 0 2 0 0
0 0 0 2 2
0 0 0 0 2

4 dimensions.
243 points.
0.296296 0.115226 0.748971 0.835391
0.399177 0.63786 0.641975 0.1893
0.777778 0.876543 0.263374 0.514403
0.0864198 0.444444 0.213992 0.0164609
0.584362 0.711934 0.921811 0.563786
...
```