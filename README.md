# Attractor-based models for sequences and pattern generation in neural circuits

This repository contains MATLAB code to reproduce the main simulations in the manuscript:
**Attractor-based models for sequences and pattern generation in neural circuits**. Preprint available at https://www.biorxiv.org/content/10.1101/2025.03.07.642121v2

## Description

Each script reproduces one or more panels from the figures in the paper. 

- `types_of_attractors.m`: examples of types of attractors realizable with TLNs, reproduces **Fig 1A–C**
- `TLN_counters.m`: counter and signed counter dynamics, reproduces **Fig 2A,D**
- `all_isolated_gaits.m`: simulates isolated quadruped gaits (bound, pace, trot, walk, pronk), reproduces **Fig 3D,E** (script includes more than just these panels)
- `all_quadruped_gaits_and_transitions.m`: dynamics for different initial conditions of 5-gait network, and transitions via targeted pulses, reproduces **Fig 4B,D**
- `quadruped_gaits_wired_counter.m`: dynamics of quadruped gaits wired with counter network, reproduces **Fig 5D,E**

Helper functions are in the `functions/` folder. These are adapted from:  
https://github.com/nebneuron/CTLN-Basic-2.0/

## Citation

If you use this code, please cite the preprint:  
Londono Alvarez J, et al. *Attractor-based models for sequences and pattern generation in neural circuits*. bioRxiv, 2025. [https://doi.org/10.1101/2025.03.07.642121](https://doi.org/10.1101/2025.03.07.642121)
