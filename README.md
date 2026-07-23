# Sample-and-Hold Control of the Semilinear Heat Equation via $L^2$ Residue Separation

This repository contains the official MATLAB implementation accompanying the paper:

> **Sample-and-hold control of the semilinear heat equation via the $L^2$ residue separation**  
> Anton Selivanov and Emilia Fridman  
> *IEEE Conference on Decision and Control (CDC), 2026.*

---

## 📌 Overview

This repository provides LMI-based tools to analyze and synthesize sampled-data boundary controllers for semilinear reaction–diffusion PDEs. 

By leveraging the **$L^2$ residue separation method**, the approach avoids boundary lifting transformations (dynamic extension), eliminating the need for complex generalized hold mechanisms and enabling physical implementation via a standard **Zero-Order Hold (ZOH)**.

### Key Features
* **Continuous-time Design (`ARE_th1.m`):** Solves coupled Algebraic Riccati Equations (AREs) for output-feedback synthesis while preventing spillover modes.
* **Sampled-Data LMI Feasibility (`LMI_th1.m`):** Constructs and solves linear matrix inequality (LMI) conditions derived via a tailored dominant-mode Lyapunov–Krasovskii functional (LKF).
* **Maximum Sampling Period Calculation (`example.m`):** Evaluates stability and automatically searches for the maximum allowable sampling period $h_{\max}$ via bisection.

---

## 🛠️ Prerequisites & Dependencies

To run the scripts, you need **MATLAB** along with the following free toolboxes/solvers:

1. **[YALMIP](https://yalmip.github.io/)** (Optimization modeling parser)
2. **SDP Solver** (At least one of the following):
   * [SeDuMi](http://sedumi.ie.lehigh.edu/) *(Recommended)*
   * [MOSEK](https://www.mosek.com/)
   * [SDPT3](https://github.com/sqlp/sdpt3)

*Note: Ensure YALMIP and your chosen SDP solver are added to your MATLAB path (`addpath(genpath(...))`).*

---

## 📁 Repository Structure

```text
.
├── example.m       % Main script: Parameter setup, ARE solver, LMI check, and h_max search
├── ARE_th1.m       % Function solving the continuous-time AREs (Theorem 1)
├── LMI_th1.m       % Function constructing and verifying the sampled-data LMIs (Theorem 1)
└── README.md       % Repository documentation****
