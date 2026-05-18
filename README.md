
# TF-QKD / FSO Quantum Communication Research Framework

<p align="center">
  <img src="https://img.shields.io/badge/Quantum-Communication-purple?style=for-the-badge">
  <img src="https://img.shields.io/badge/TF--QKD-Research-black?style=for-the-badge">
  <img src="https://img.shields.io/badge/Finite-Decoy-LP-darkgreen?style=for-the-badge">
  <img src="https://img.shields.io/badge/Status-Active%20Development-success?style=for-the-badge">
</p>

---

## Overview

A research-oriented implementation and numerical exploration framework
for Twin-Field-Type Quantum Key Distribution (TF-QKD), finite-decoy
security analysis, and future Free-Space Optical (FSO) quantum
communication systems.

This repository is being developed as part of ongoing research and
exploratory work in:

- Quantum Key Distribution (QKD)
- Twin-Field Quantum Cryptography
- Free-Space Optical (FSO) Communication
- Atmospheric Optical Channels
- Satellite Quantum Communication
- Numerical Security Analysis

The framework currently includes:

- asymptotic TF-QKD simulation,
- finite-decoy LP estimation,
- phase-error analysis,
- key-rate optimization,
- robustness diagnostics,
- convergence validation,
- realistic optical-loss modeling.

The long-term objective of this work is to build a modular
research-grade simulation environment capable of studying both
fiber-based and free-space TF-QKD systems under realistic physical
conditions.

---

# Simulation Pipeline

```text
                ┌──────────────────────────────┐
                │     Physical Parameters      │
                │ α, η, pd, θ, ξ, truncation   │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │     Channel Modeling         │
                │ loss, mismatch, dark counts  │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │     Yield Computation        │
                │ pXX , pZZ , photon yields    │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │  Cat-State Decomposition     │
                │ even / odd parity sectors    │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │   Phase Error Estimation     │
                │ Eq. (20–21), truncation      │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │  Finite-Decoy LP Estimation  │
                │ Eq. (23–26), CVXPY solver    │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │     Secret Key Rate          │
                │ Eq. (17–19), PLOB compare    │
                └──────────────┬───────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │ Numerical Diagnostics & Plots│
                │ convergence, robustness      │
                └──────────────────────────────┘
````

---

# Current Status

The simulator has progressed through multiple implementation phases.

---

# Completed Phases

## Phase 1 — Foundation and Channel Modeling

Implemented:

* Binary entropy functions
* Optical transmittance utilities
* dB ↔ transmittance conversion
* X-basis gain calculations
* Z-basis photon-yield calculations
* Dark-count modeling
* Polarization mismatch modeling
* Phase mismatch handling
* Numerical stability protections
* Physical sanity diagnostics

Includes:

* channel-loss sweeps,
* mismatch sensitivity studies,
* detector-noise analysis.

---

## Phase 2 — Cat States and Phase Error Estimation

Implemented:

* Cat-state decomposition
* Even/odd parity separation
* Infinite-decoy phase-error estimation
* Eq. (20) / Eq. (21) implementation
* Photon-number truncation framework
* Inner-sum convergence analysis
* Numerical truncation diagnostics

Includes:

* parity-resolved contributions,
* convergence studies vs truncation parameter,
* numerical validation of phase-error bounds.

---

## Phase 3 — Secret Key Rate and Numerical Analysis

Implemented:

* Full asymptotic secret-key-rate pipeline
* Eq. (17–19) implementation
* PLOB repeaterless bound comparison
* Alpha-parameter optimization
* Key-rate vs channel-loss simulations
* Dark-count sensitivity analysis
* Polarization mismatch robustness
* Phase mismatch robustness
* Numerical convergence diagnostics
* Figure-generation utilities

Current studies include:

* key rate vs total channel loss,
* dark-count scaling,
* phase mismatch sensitivity,
* polarization mismatch sensitivity,
* alpha optimization behavior,
* truncation convergence,
* comparison against repeaterless bounds.

---

## Phase 4 — Finite-Decoy Linear Programming Framework

Implemented:

* Finite-decoy LP formulation
* Eq. (23–26) implementation
* Poisson coefficient matrix construction
* Residual-tail handling
* LP feasibility diagnostics
* Photon-pair yield upper-bound estimation
* Finite-decoy phase-error estimation
* Finite vs infinite decoy comparison
* CVXPY-based optimization pipeline
* LP consistency validation framework

Implemented diagnostics include:

* LP upper-bound validation,
* yield consistency checks,
* finite ≤ infinite key-rate verification,
* truncation-tail validation,
* convergence diagnostics,
* cutoff-range consistency analysis.

Current LP framework supports:

* arbitrary decoy-intensity sets,
* configurable truncation parameter (`Mcut`),
* batch yield estimation,
* finite-decoy key-rate optimization.

---

# Current Numerical Characteristics

The framework currently uses:

* asymptotic security assumptions,
* explicit photon-number truncation,
* parity-resolved cat-state expansions,
* direct numerical evaluation of Eq. (35)-type expressions,
* LP-based finite-decoy estimation,
* scalar/grid alpha optimization,
* dense parameter sweeps.

The simulator becomes computationally expensive for:

* large truncation parameters,
* high-dimensional LP sweeps,
* dense optimization grids,
* high-loss asymptotic studies.

Fast-analysis modes are therefore included for exploratory studies.

---

# Validation Status

The framework currently includes automated validation checks for:

* Poisson normalization,
* LP feasibility,
* physical yield bounds,
* finite-decoy consistency,
* truncation-tail adequacy,
* phase-error consistency,
* entropy-function correctness,
* numerical positivity constraints.

The implementation is now internally physics-consistent for the
implemented asymptotic and finite-decoy analyses.

---

# Ongoing Work

Current development directions include:

* reproduction of paper-level finite-decoy plots,
* tighter LP formulations,
* decoy-state optimization,
* improved numerical stability,
* performance optimization.

---

# Planned Extensions

## Free-Space Optical (FSO) Channels

Current exploratory work includes early-stage implementation of:

* atmospheric attenuation models,
* turbulence-aware channel studies,
* stochastic optical fading,
* pointing-error analysis,
* atmospheric-loss modeling.

Planned future extensions include:

* Gamma-Gamma turbulence,
* log-normal turbulence,
* scintillation studies,
* adaptive optics,
* aperture averaging,
* beam-wander analysis.

---

## Satellite Quantum Communication

Planned future studies include:

* uplink/downlink asymmetry,
* orbital-loss modeling,
* satellite pointing jitter,
* dynamic atmospheric channels,
* tracking and alignment analysis.

---

## Advanced Numerical Studies

Potential future work:

* asymmetric TF-QKD,
* finite-size effects,
* detector-efficiency mismatch,
* parameter-estimation analysis,
* composable security extensions.

---

# Repository Structure

```text
notebooks/
│
├── Phase1_Foundation_and_Channel_Model.ipynb
├── Phase2_CatStates_and_PhaseError.ipynb
├── Phase3_KeyRate_Optimization_and_Robustness.ipynb
├── Phase4_Finite_Decoy_LP.ipynb
├── Phase5_FSO_Channel_Modeling.ipynb
```

---

# Important Note

This repository is still an actively developing research-oriented
implementation intended for:

* conceptual understanding,
* numerical experimentation,
* protocol exploration,
* reproducibility studies.

Although major physical inconsistencies and LP-structure issues have
been resolved, the framework should not yet be treated as a fully
validated production-grade TF-QKD simulator.

Some finite-decoy regimes remain conservative due to LP looseness and
limited decoy configurations.

---

# Tools and Libraries

Primary libraries:

* Python
* NumPy
* SciPy
* Matplotlib
* CVXPY
* Jupyter Notebook

---

# References

Primary references:

* Marcos Curty, Koji Azuma, and Hoi-Kwong Lo
  *Simple security proof of twin-field type quantum key distribution protocol* (2018)

* Lucamarini et al.
  *Overcoming the rate–distance limit of quantum key distribution without quantum repeaters*
  Nature (2018)

---

# Research Watermark

```text
────────────────────────────────────────────────────────────
TF-QKD / FSO Quantum Communication Research Framework
Developed by Anand Kumar

B.Tech Electronics and Communication Engineering
Guru Gobind Singh Indraprastha University (GGSIPU)

Research Areas:
Quantum Communication • TF-QKD • FSO Systems
Satellite Quantum Communication • Atmospheric Channels

© Anand Kumar — Numerical Research Framework
────────────────────────────────────────────────────────────
```

---

# Author

## Anand Kumar

B.Tech — Electronics and Communication Engineering (2024–2028)
Guru Gobind Singh Indraprastha University (GGSIPU), New Delhi

Research Interests:

* Quantum Communication
* Quantum Key Distribution (QKD)
* Free-Space Optical (FSO) Communication
* Atmospheric Channel Modeling
* Satellite Quantum Communication
* Quantum Networks

Email: [anandambastha72@gmail.com](mailto:anandambastha72@gmail.com)

GitHub: https://github.com/Anand-Ambastha

LinkedIn: https://www.linkedin.com/in/anand-kumar05
