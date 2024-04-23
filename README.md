# Quad-Optimized-LDS

Source code of the paper: "Quad-Optimized Low-Discrepancy Sequences", ACM SIGGRAPH Conference Proceeedings, San Diego, 2024.

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
 
For instance, to generate matrices for up to $3^4$ samples in dimension 4:

```
samplerGF3 -i ../initIrreducibleGF3.dat -m 4 -d 4
```

which outputs the four matrices in GF(3) and the samples:

```
d= 0
1 0 0 0
0 1 0 0
0 0 1 0
0 0 0 1

d= 1
1 1 2 2
0 1 0 2
0 0 1 1
0 0 0 1

d= 2
1 2 1 2
0 1 1 0
0 0 1 0
0 0 0 1

d= 3
2 2 2 2
0 2 1 0
0 0 2 0
0 0 0 2

4 dimensions.
243 points.
0.691358 0.209877 0.0123457 0.0864198
0.283951 0.481481 0.37037 0.580247
0.54321 0.938272 0.851852 0.938272
0.901235 0.382716 0.753086 0.444444
...
```