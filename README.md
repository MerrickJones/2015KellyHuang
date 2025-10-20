<!-- Enable MathJax for LaTeX equation rendering -->
<script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
<script id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>

# Bayesian Updating for One-Dimensional Consolidation (Kelly & Huang, 2015)

[![Python](https://img.shields.io/badge/Python-3.12%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Last Commit](https://img.shields.io/github/last-commit/MerrickJones/2015KellyHuang)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.compgeo.2015.06.010-blue.svg)](https://doi.org/10.1016/j.compgeo.2015.06.010)

This repository implements a **Bayesian back-analysis** of **one-dimensional consolidation** following the methodology of **Kelly & Huang (2015)**.  
It employs a Differential Evolution Adaptive Metropolis (**DREAM**) Markov Chain Monte Carlo (**MCMC**) algorithm to update geotechnical parameters governing settlement over time.

---

## 🧩 Overview

The program performs parameter inference for a 1D consolidation problem, progressively updating estimates of:
- **m_v** — coefficient of volume compressibility (1/kPa)  
- **H** — drainage path length (m)  
- **γ_f** — fill unit weight (kN/m³)  
- **c_v** — coefficient of consolidation (m²/year)

Given observed settlement data, the DREAM MCMC sampler explores the posterior distributions of these parameters, updating their estimates as additional field observations become available.

Two variants are included:
1. **example_consolidation.py** — baseline Bayesian updating for the full dataset.  
2. **progressive_DREAM.py** — progressive assimilation experiment, updating posteriors sequentially with each new observation.

---

## 📁 File Structure

| File | Description |
|------|--------------|
| `kellyFORTRAN.py` | Forward deterministic model implementing Eq. (4) from Kelly & Huang (2015). Computes theoretical settlement vs. time for given parameters. |
| `DREAM_Suite.py` | Simplified DREAM MCMC engine implementing Differential Evolution–Metropolis updates for posterior sampling. |
| `example_consolidation.py` | Main single-run Bayesian updating experiment using the full observation dataset. |
| `progressive_DREAM.py` | Progressive (data-assimilation-style) version performing sequential posterior updates as new settlement data become available. |
| `settlement_obs.txt` | Observed settlement measurements (m) used for calibration. |
| `install.py` | Simple environment setup script that installs required dependencies (`numpy`, `scipy`, `matplotlib`). |

---

## ⚙️ Dependencies

To install dependencies manually:

```bash
pip install numpy scipy matplotlib
```

Or simply run:

```bash
python install.py
```

---

## 🚀 Running the Models

### 1. Baseline DREAM MCMC
Runs full Bayesian updating using all available settlement data:

```bash
python example_consolidation.py
```

**Outputs:**
- `back_analysis_results.txt` — summary of prior/posterior/true parameters  
- `back_analysis_summary.png` — graphical summary of posterior distributions  
- Main figure showing settlement prediction with 95 % credible interval

---

### 2. Progressive Bayesian Updating
Runs progressive DREAM MCMC, adding one new observation at a time:

```bash
python progressive_DREAM.py
```

**Outputs:**
- `Figure_1.png` — progressive settlement predictions for each observation  
- `progress_SETT_results.txt` — text summary of convergence accuracy  
- `Figure_2.png` — four-panel summary: Bayesian update, embankment schematic, and posteriors for `m_v` and `c_v`  
- `Figure_3_progress_SETT_results.png` — visual rendering of the results text

---

## 📊 Workflow Summary

1. **Forward Model:**  
   The program computes theoretical settlement \( s(t) \) using the analytical solution from Kelly & Huang (2015):

   \[
   s(t) = m_v \sigma H \left[ 1 - \sum_{m=0}^{\infty} \frac{2}{M^2} e^{-M^2 T_v} \right],
   \quad M = \frac{\pi}{2}(2m+1), \quad T_v = \frac{c_v t}{H^2}
   \]

2. **Bayesian Updating:**  
   DREAM MCMC explores the posterior parameter space, updating the prior estimates by minimising the difference between predicted and observed settlements.

3. **Progressive Assimilation:**  
   The progressive version performs incremental updates, identifying how many field observations are required before settlement predictions are within 5 % of the true value.

---

## 📈 Key Features

- Full Bayesian framework for soil consolidation analysis  
- DREAM(ZS)-style multi-chain adaptive MCMC  
- Progressive assimilation demonstrating model convergence  
- Built-in plotting and result export  
- Modular, extensible codebase for coupling with external solvers (e.g., CAOS, PLAXIS, Settle3)

---

## 📚 Reference

Kelly, R. B., & Huang, J. (2015).  
*Bayesian updating of consolidation parameters from field measurements.*  
**Computers and Geotechnics**, 69, 496 – 507.  
[https://doi.org/10.1016/j.compgeo.2015.06.010](https://doi.org/10.1016/j.compgeo.2015.06.010)

---

## 🧠 Suggested Extensions

- Implement **lognormal priors** for strictly positive parameters (`m_v`, `c_v`).  
- Parallelise the DREAM sampler using `multiprocessing.Pool` for multi-chain acceleration.  
- Integrate with **CAOS** or **PLAXIS** to perform real-data back-analyses of embankment case studies.  
- Add **convergence diagnostics** (Gelman–Rubin, autocorrelation, trace plots).  
- Expand to 2D and 3D consolidation models for comparison with FEM solvers.

---

## 🧾 Citation

If you use this repository or adapt the DREAM-based consolidation workflow in your research, please cite:

> Jones, M. (2025). *Bayesian Updating for One-Dimensional Consolidation (Kelly & Huang, 2015) [Code repository].* GitHub.  
> [https://github.com/MerrickJones/2015KellyHuang](https://github.com/MerrickJones/2015KellyHuang)

---

## 👤 Author

**Merrick Jones (2025)**  
PhD Candidate, University of Newcastle  
*(Bayesian Back-Analysis for Embankments on Soft Soils)*  
📧 merrick.jones@uon.edu.au  

---

## 🪪 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
