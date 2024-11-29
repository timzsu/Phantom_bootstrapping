# Phantom FHE with CKKS Bootstrapping

This repository contains codes that enables CKKS bootstrapping for [Phantom-FHE](https://github.com/encryptorion-lab/phantom-fhe). 

The codes are borrowed from [NEXUS-CUDA](https://github.com/zju-abclab/NEXUS), a secure inference framework. I reorganized the files to create an installable package so that I can use that as a backend for my other projects. 

```bash
cmake -S . -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=true
cmake --build build --parallel
cmake --install build
```