# Phantom FHE with CKKS Bootstrapping

This repository contains codes that enables CKKS bootstrapping for [Phantom-FHE](https://github.com/encryptorion-lab/phantom-fhe). 

The codes are borrowed from [NEXUS-CUDA](https://github.com/zju-abclab/NEXUS), a secure inference framework. I reorganized the files to create an installable package so that I can use it as a backend for my other projects. 

### Dependencies

* NVIDIA GPU (Volta or newer)
* CUDA >= 11.0
* CMake >= 3.20
* GCC >= 9.0
* NTL library

> [!WARNING]  
> The NTL library has to be installed from [source](https://libntl.org/doc/tour.html) and reinstalled for each GPU type. 

### Installation

To configure the library: 
```bash
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=true\
 -DCMAKE_INSTALL_PREFIX=/path/to/install
```
if `CMAKE_INSTALL_PREFIX` is not specified, then cmake will install to the default path. 

To build and install the library: 
```bash
cmake --build build --parallel
cmake --install build
```

To test the build: 
```bash
./build/bin/example_context
```