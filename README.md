# Twin-Field-Type QKD Simulation Framework

This repository contains a phased exploratory implementation of a 2018 twin-field-type quantum key distribution (TF-type QKD) framework in Python, based primarily on the simplified security-proof formulation introduced by Curty, Azuma, and Lo.

The primary objective of this work is to develop a practical and
conceptual understanding of TF-QKD system behavior, key-rate scaling,
and channel sensitivity under realistic optical losses.

---

## Current Status

The implementation is currently being developed incrementally in
multiple phases.

### Completed Phases

#### Phase 1 — Foundation and Channel Modeling

Implemented:

- Binary entropy functions
- Channel transmittance utilities
- Optical loss modeling
- X-basis gain calculations
- Z-basis yield calculations
- Dark-count modeling
- Polarisation mismatch modeling
- Phase mismatch handling
- Numerical stability guards
- Physical sanity checks

#### Phase 2 — Cat States and Phase Error Estimation

Implemented:

- Cat-state decomposition framework
- Even/odd photon-parity separation
- Phase-error upper-bound estimation
- Eq. (20) / Eq. (21) implementation
- Photon-number truncation framework
- Infinite-decoy yield estimation
- Numerical truncation diagnostics

#### Phase 3 — Secret Key Rate and Numerical Analysis

Implemented:

- Full asymptotic secret-key-rate pipeline
- Eq. (17–19) implementation
- PLOB repeaterless bound
- Alpha-parameter optimization
- Key-rate vs loss analysis
- Dark-count sensitivity studies
- Photon truncation convergence analysis
- Numerical robustness diagnostics
- Figure-generation framework

Current numerical studies include:

- secret key rate vs channel loss
- dark-count sensitivity
- phase mismatch sensitivity
- polarization mismatch sensitivity
- photon-number truncation convergence
- comparison against the PLOB bound
- alpha optimization behavior

---

## Current Numerical Characteristics

The implementation currently uses:

- asymptotic infinite-decoy assumptions
- explicit photon-number truncation
- nested combinatorial summations
- direct numerical evaluation of Eq. (35)-type expressions
- scalar/grid-based alpha optimization

The simulator is computationally expensive for:

- large truncation parameters
- dense parameter sweeps
- high-loss asymptotic studies

Current optimization and convergence studies are therefore performed
with moderate truncation limits and configurable fast-analysis modes.

---

## Ongoing Work

The next stages involve:

* Reproduction of paper-level plots
* Finite-decoy linear programming methods

Future extensions may also include:

* Free-space optical (FSO) channel effects
* Atmospheric turbulence modeling
* Pointing-error analysis
* Asymmetric channel studies

---

## Repository Structure

```text
notebooks/
│
├── Phase1_Foundation_and_Channel_Model.ipynb
├── Phase2_CatStates_and_PhaseError.ipynb
├── Phase3_KeyRate_Optimization_Robustness_Studies.ipynb
├── Phase4_Finite_Decoy_LP.ipynb
```

---

## Important Note

This repository is currently an exploratory and learning-oriented
implementation intended for research understanding and reproducibility.

The implementation is still under active development and should not yet
be considered a fully validated research-grade simulator.

---

## Tools and Libraries

* Python
* NumPy
* SciPy
* Matplotlib
* CVXPY (planned)

---

## References

Primary references:

* Marcos Curty, Koji Azuma, and Hoi-Kwong Lo,
  "Simple security proof of twin-field type quantum key distribution protocol," 2018.

* Lucamarini et al.,
  "Overcoming the rate–distance limit of quantum key distribution without quantum repeaters,"
  Nature, 2018.

---

## Author

**Anand Kumar** ,

B.Tech Electronics and Communication Engineering (2024–2028),
Guru Gobind Singh Indraprastha University (GGSIPU), New Delhi

Research Interests:

* Quantum Communication
* Free-Space Optical (FSO) Communication
* Quantum Key Distribution (QKD)
* Atmospheric Channel Modeling
* Satellite Quantum Communication

Email: [anandambastha72@gmail.com](mailto:anandambastha72@gmail.com)
GitHub: [https://github.com/Anand-Ambastha](https://github.com/Anand-Ambastha)
LinkedIn: [https://www.linkedin.com/in/anand-kumar05](https://www.linkedin.com/in/anand-kumar05)
