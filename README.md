# Attractor-based models for sequences and pattern generation in neural circuits

This repository contains MATLAB code to reproduce the main simulations in the manuscript:
**Attractor-based models for sequences and pattern generation in neural circuits**. Preprint available at https://www.biorxiv.org/content/10.1101/2025.03.07.642121v2

## Description

The simulations in this repository will print the rate curves for the following figure panels:
- Figure 1A-C: examples of types of attractors realizable with TLNs
- Figure 2A,D: counter and signed counter
- Figure 3D,E: bound and walk gaits (more gaits are included in the script)
- Figure 4B–D: rate curves for different initial conditions of 5-gait network, and transitions via targeted pulses
- Figure 5D: rate curves for quadruped gaits wired with counter network

All scripts are written in MATLAB and are self-contained. Helper functions are included in the `functions/` folder. The helper functions are originally from: https://github.com/nebneuron/CTLN-Basic-2.0/

## Instructions

No installation or setup is required beyond a working MATLAB installation.

Each figure script can be run independently. Comments at the top of each script indicate which figure it reproduces.

## Citation

If you use this code, please cite the preprint:  
Londono Alvarez J, et al. *Attractor-based models for sequences and pattern generation in neural circuits*. bioRxiv, 2025. [https://doi.org/10.1101/2025.03.07.642121](https://doi.org/10.1101/2025.03.07.642121)
