# ProteinSimulatedAnnealing

This project implements an algorithm for predicting the tertiary structure of proteins using simulated annealing optimization. The implementation focuses on performance optimization through hardware-aware programming techniques.

## Overview

The project aims to predict how a protein's amino acid chain folds in three-dimensional space by finding the structure with minimum energy. This is achieved by calculating various energy components and using simulated annealing to find optimal phi (φ) and psi (ψ) angles that define the protein's 3D structure.

### Energy Components
The algorithm considers four main energy components:
- Ramachandran energy (potential energy associated with different angles)
- Hydrophobic energy (based on 3D distances)
- Electrostatic energy (based on 3D distances)
- Packing energy (based on volume and 3D distances)

## Implementation Requirements

### Architectures
The project requires implementation in multiple versions:
- x86-32+SSE (C and Assembly)
- x86-64+AVX (C and Assembly)
- x86-64+AVX with OpenMP

### Program Usage
Check commands file

Where:
- `<arch>`: Architecture identifier (32c, 64c, 32ompc, 64ompc)
- `<SEQ>`: Input file containing the amino acid sequence
- `<t0>`: Initial temperature for simulated annealing
- `<alpha>`: Cooling rate
- `<k>`: Constant for the algorithm
- `<sd>`: Random seed

### Development Guidelines
1. First implement the algorithm entirely in C as a sequence of function calls
2. Replace performance-critical functions with optimized assembly implementations
3. Compare performance between different optimization stages
4. Implement OpenMP parallel processing for the 64-bit AVX version

### Required Files
- C source files: `pst32c.c`, `pst64c.c`
- Assembly source files: `pst32.nasm`, `pst64.nasm`
- OpenMP version: Files with "_omp" suffix
- Execution scripts: `runpst32`, `runpst64`

## Technical Details

### Core Algorithm Components
1. Backbone reconstruction from angles to 3D coordinates
2. Energy calculation functions:
   - Ramachandran energy calculation
   - Hydrophobic energy calculation
   - Electrostatic energy calculation
   - Packing energy calculation
3. Simulated annealing optimization

### Development Environment
- Programming Languages: C (gcc), x86 Assembly (nasm)
- Assembly Extensions: SSE (32-bit), AVX (64-bit)
- Operating System: Linux (Ubuntu)
- Parallel Processing: OpenMP

## Project Deliverables
1. Source code for all required implementations
2. Technical documentation/report
3. Performance comparison analysis
4. Presentation slides (optional)

## Note
The project should be developed independently. Similar solutions between different groups will receive negative evaluations. Implementation details and file naming conventions will be provided separately before the submission deadline.

[Progetto_Architetture_2024.pdf](https://github.com/user-attachments/files/18458602/Progetto_Architetture_2024.pdf)
