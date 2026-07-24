# Sample-and-Hold Control of the Semilinear Heat Equation via $L^2$ Residue Separation

MATLAB code accompanying the paper:

> **Sample-and-hold control of the semilinear heat equation via the $L^2$ residue separation**  
> Anton Selivanov and Emilia Fridman  
> *IEEE Conference on Decision and Control (CDC), 2026.*

---

### 🚀 File Description
* **`ARE_th1.m`:** Solves coupled Algebraic Riccati Equations (AREs) for continuous-time output-feedback synthesis.
* **`LMI_th1.m`:** Solves linear matrix inequalities (LMIs) to verify stability under sample-and-hold actuation for a given sampling period.
* **`example.m`:** Runs the full synthesis and analysis for the paper's numerical example.

---

### 🛠️ Dependencies
Requires **MATLAB** with:
* **[YALMIP](https://yalmip.github.io/)** (Optimization parser)
* An SDP solver (e.g., [SeDuMi](http://sedumi.ie.lehigh.edu/), [MOSEK](https://www.mosek.com/), or [SDPT3](https://github.com/sqlp/sdpt3))

---

### 📁 Repository Structure
```text
.
├── example.m   % Main script: Parameter setup, ARE solver, LMI check, and h_max search
├── ARE_th1.m   % Function solving continuous-time AREs (Theorem 1)
└── LMI_th1.m   % Function verifying sampled-data LMIs (Theorem 1)
