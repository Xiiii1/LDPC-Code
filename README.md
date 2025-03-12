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
g++ -o ldpc_simulation ldpc_simulation.cpp -fopenmp
