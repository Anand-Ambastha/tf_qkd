# Research Notebooks

This section documents the implementation roadmap of the Twin-Field Quantum Key Distribution (TF-QKD) research framework. The notebooks progress from protocol foundations and security analysis to practical finite-decoy implementations. An additional Free-Space Optical (FSO) extension investigates the deployment of TF-QKD over atmospheric channels.

---

## Phase 1 – Foundations and Channel Modelling

### Foundation and Channel Model

Introduces the theoretical foundations of Twin-Field Quantum Key Distribution (TF-QKD), including protocol assumptions, optical losses, detector models, channel transmittance calculations, and baseline simulation infrastructure. This phase establishes the mathematical and computational framework used throughout the project.

---

## Phase 2 – Quantum States and Security Analysis

### Cat States and Phase Error Analysis

Implements the quantum-state representation employed in TF-QKD and studies phase-error estimation techniques required for security analysis. The notebook investigates the relationship between observable quantities and the security parameters governing secret-key generation.

---

## Phase 3 – Key Rate Optimization and Robustness

### Key Rate Optimization and Robustness

Explores parameter optimization strategies for maximizing secret-key rates under realistic operating conditions. Sensitivity analyses are performed to evaluate protocol robustness against channel loss, detector imperfections, and environmental disturbances.

---

## Phase 4 – Practical Finite-Decoy Implementation

### Finite Decoy Linear Programming

Implements finite-decoy-state estimation using linear-programming methods. This phase transitions from idealized theoretical assumptions to experimentally realistic scenarios, enabling practical secret-key-rate evaluation with finite data samples.

---

## Extension Module – Free-Space Optical (FSO) Quantum Channels

## Phase 1 Extension – Deterministic FSO Channel Modelling

### FSO Channel

This notebook extends the Phase 1 channel modelling framework to free-space optical (FSO) communication links. The implementation focuses on deterministic propagation effects, including geometric beam spreading, atmospheric attenuation, link-budget calculations, and receiver coupling losses.

The objective is to establish a baseline FSO channel model before introducing stochastic effects such as turbulence-induced fading, scintillation, beam wander, or pointing errors in future developments.

Key topics include:

- Free-space path loss
- Atmospheric attenuation
- Beam divergence effects
- Receiver aperture collection efficiency
- End-to-end channel transmittance
- Deterministic link-budget analysis

This extension provides the foundational channel model required for future TF-QKD over FSO investigations and satellite-based quantum communication studies.

---

## Repository Resources

The repository also contains:

- Research notes and derivations
- Literature reviews and survey papers
- TF-QKD security-proof references
- Technical documentation and project reports
- Supporting implementation resources

Together, these materials form a comprehensive research framework for studying Twin-Field Quantum Key Distribution and its potential deployment in advanced quantum communication networks.