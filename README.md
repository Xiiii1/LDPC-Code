# LDPC Code Simulation

This project implements a simulation of a Low-Density Parity-Check (LDPC) code decoder using different decoding algorithms, including Sum-Product Algorithm (SPA), Min-Sum Algorithm (MSA), and Bit-Flipping (BF) algorithm. The simulation includes the generation of codewords, modulation, AWGN channel noise, and error correction through iterative decoding.

## Overview

The project simulates a communication system where:

1. **Codeword Generation**: LDPC codewords are generated using a Linear Block Code (LBC) generator matrix (`G`).
2. **Transmission Over AWGN Channel**: The generated codewords are passed through an Additive White Gaussian Noise (AWGN) channel, introducing noise to the transmitted signal.
3. **Error Correction**: The received signal is decoded using one of the specified decoding algorithms:
    - **SPA (Sum-Product Algorithm)**
    - **MSA (Min-Sum Algorithm)**
    - **BF (Bit-Flipping Algorithm)**

The goal of the simulation is to evaluate the performance of the LDPC decoder in terms of Bit Error Rate (BER) and Block Error Rate (BLER) under different Signal-to-Noise Ratios (SNRs).

## Features

- **Multiple Decoding Algorithms**: SPA, MSA, and BF algorithms are implemented to demonstrate different error-correction approaches.
- **AWGN Channel Simulation**: The received signal is subject to noise generated based on a Gaussian distribution, with a specified noise level determined by the SNR.
- **LDPC Code Parameters**: The LDPC code's generator and parity-check matrices (`G` and `H`) are loaded from external files.
- **Parallelized Processing**: The simulation uses OpenMP to parallelize parts of the algorithm for faster execution.
- **Error Rate Calculation**: The simulation calculates and prints the Bit Error Rate (BER) and Block Error Rate (BLER) based on the number of errors and the number of decoded blocks.

## Files

- **Sim.txt**: Contains the simulation parameters, such as the number of blocks, iterations, SNR, SEED, and the algorithm to be used.
- **ldpc_G_1023.txt**: LDPC generator matrix (`G`).
- **ldpc_H_1023.txt**: LDPC parity-check matrix (`H`).

## Compilation

To compile the code, use the following command (assuming you're using `g++`):

```bash
g++ -o ldpc_simulation main.cpp -fopenmp
```

The `-fopenmp` flag enables OpenMP for parallel processing.

## Usage
1. Prepare Input Files: Ensure that `Sim.txt`, `ldpc_G_1023.txt`, and `ldpc_H_1023.txt` are correctly placed in the same directory as the program.

2. Run the Simulation: After compiling, run the program using the following command:
```bash
./ldpc_simulation
```
3. View Results: The program will output the simulation results, including the number of decoded blocks, errors, and the final BER/BLER values.

## Parameters in Sim.txt
The `Sim.txt` file should contain the following parameters, each on a separate line:

```php-template
<Number of Blocks>
<Number of Iterations>
<SNR (in dB)>
<Random Seed (SEED)>
<Decoding Algorithm (0 = SPA, 1 = MSA, 2 = BF)>
```
For example:
```
100
50
10.0
12345
0
```
This configuration will simulate 100 blocks, with 50 iterations per block, an SNR of 10 dB, a random seed of 12345, and use the SPA algorithm for decoding.

## Decoding Algorithms
1. **SPA** (Sum-Product Algorithm)
SPA is an iterative algorithm used in LDPC decoding. It uses the sum-product operation to propagate information through the graph of variable and check nodes, eventually converging on the most likely codeword.

2. **MSA** (Min-Sum Algorithm)
MSA is a simplified version of the SPA, where the sum-product operation is replaced by a min-sum operation. It offers lower computational complexity but may result in a higher error rate compared to SPA.

3. **BF** (Bit-Flipping Algorithm)
The Bit-Flipping algorithm is a simpler error correction technique. It flips bits in the codeword based on the syndromes of the received signal, with the goal of driving all syndromes to zero.

## Performance Metrics
The simulation will output the following metrics:

- Number of Decoded Blocks: The total number of blocks successfully decoded.
- Bit Error Rate (BER): The ratio of erroneous bits to the total number of bits transmitted.
- Block Error Rate (BLER): The ratio of erroneous blocks to the total number of blocks transmitted.
Example output:
```makefile
----------------------------------------------------
Number of decoded blocks = 100
Number of iterations = 50
SNR = 10.000000
Algorithm is 0
# total error bits = 23
# total error blocks = 2
BER = 0.00029412
BLER = 0.020000
running time: 12.345678 sec
```

## Dependencies
- C++ Compiler: Ensure you have a C++ compiler that supports the C++11 standard or higher (e.g., `g++`).
- OpenMP: The program uses OpenMP for parallelization, so your compiler should support it.
