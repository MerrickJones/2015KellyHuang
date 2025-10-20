# Bayesian Updating for One‑Dimensional Consolidation (Kelly & Huang, 2015)

[![Python](https://img.shields.io/badge/Python-3.11%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Last Commit](https://img.shields.io/github/last-commit/MerrickJones/2015KellyHuang)
[![Paper DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.compgeo.2015.06.010-blue.svg)](https://doi.org/10.1016/j.compgeo.2015.06.010)

This repository implements a **Bayesian back‑analysis** of **one‑dimensional consolidation** following the methodology of **Kelly & Huang (2015)**.  
It employs a Differential Evolution Adaptive Metropolis (**DREAM**) Markov Chain Monte Carlo (**MCMC**) algorithm to update geotechnical parameters governing settlement through time.

> ℹ️ GitHub’s Markdown does **not** render LaTeX by default.  
> All equations below are provided in a GitHub‑friendly monospaced format, so they display correctly on the repository page.

---

## 🧩 Overview

The program estimates parameters for a 1D consolidation model by calibrating theoretical settlement to observations:

- **m_v** — coefficient of volume compressibility (1/kPa)  
- **H** — drainage path length (m)  
- **γ_f** — fill unit weight (kN/m³)  
- **c_v** — coefficient of consolidation (m²/time; units consistent with *t*)

Two variants are included:
1. `example_consolidation.py` — baseline Bayesian updating using the full dataset.  
2. `progressive_DREAM.py` — progressive (sequential) updating as new observations arrive.

---

## 📁 File Structure

| File | Description |
|------|-------------|
| `kellyFORTRAN.py` | Forward deterministic model implementing the Kelly & Huang (2015) series solution for settlement vs. time. |
| `DREAM_Suite.py` | Minimal DREAM‑style MCMC sampler (differential‑evolution proposals, multi‑chain). |
| `example_consolidation.py` | Single‑run Bayesian calibration using all observations. |
| `progressive_DREAM.py` | Sequential assimilation; re‑estimates posteriors as each new data point is added. |
| `settlement_obs.txt` | Observed settlement measurements used for calibration. |
| `install.py` | Convenience installer for core dependencies (`numpy`, `scipy`, `matplotlib`). |

---

## ⚙️ Dependencies

Install manually:
```bash
pip install numpy scipy matplotlib
```
or run:
```bash
python install.py
```

---

## 🚀 How to Run

### 1) Baseline DREAM MCMC
```bash
python example_consolidation.py
```
**Outputs (typical):**
- `back_analysis_results.txt` — prior/posterior summary
- `back_analysis_summary.png` — posterior plots and fit
- settlement prediction plot with 95% credible interval

### 2) Progressive (Sequential) Updating
```bash
python progressive_DREAM.py
```
**Outputs (typical):**
- `Figure_1.png` — evolving predictions over observation times
- `progress_SETT_results.txt` — text summary of convergence
- `Figure_2.png` — panel view: update curve + marginal posteriors for *m_v* and *c_v*
- `Figure_3_progress_SETT_results.png` — annotated results page

---

## 📊 Forward Model (GitHub‑Friendly Formula)

The forward model computes theoretical settlement *s(t)* using the standard odd‑term series solution (double drainage), as used by Kelly & Huang (2015).

\[
s(t) = m_v \sigma H \left[1 - \sum_{m=0}^{\infty} \frac{2}{M^2} e^{-M^2 T_v}\right],
\quad M = \frac{\pi}{2}(2m + 1), \quad T_v = \frac{c_v t}{H^2}
\]

where:  
- *m_v* = coefficient of volume compressibility,  
- *σ* = applied stress increment,  
- *H* = drainage path length,  
- *c_v* = coefficient of consolidation,  
- *t* = time.  

In practice, the infinite series is truncated to a large integer *N* (e.g., 200–500 terms) such that the tail is negligible.

**Implementation notes**
- In code, the infinite series is truncated to a large integer **N** (e.g., 200–500 terms) such that the tail is negligible for the *t* of interest.  
- Units of **c_v** and **t** must be consistent (e.g., m²/s with seconds, or m²/yr with years).  
- For single‑drainage conditions, adjust **H** accordingly (drainage path length).

---

## 🔄 Bayesian Updating Workflow

1. **Define priors** for (*m_v, c_v, H, γ_f*). Positivity‑constrained parameters are commonly assigned **lognormal** priors.  
2. **Simulate** settlement predictions *s(t | θ)* using the forward model for parameter vector **θ**.  
3. **Evaluate likelihood** from data misfit (e.g., Gaussian errors with σ_obs or a robust alternative).  
4. **Sample posteriors** with DREAM (multiple chains, DE‑style proposals, periodic adaptation).  
5. **Diagnose convergence** (trace, R‑hat, ESS, autocorrelation).  
6. **Summarise & predict** (posterior means/medians, credible intervals, posterior‑predictive checks).

---

## 📈 Features

- DREAM(ZS)‑style **multi‑chain, self‑tuning** proposals (Differential‑Evolution moves)
- **Progressive assimilation** example for field monitoring workflows
- Clear separation of **forward model** and **inference** logic
- Ready to couple with external solvers (CAOS, PLAXIS, Settle3)
- Export of summaries and publication‑quality figures

---

## 🧠 Suggested Extensions

- Switch to **lognormal priors** for strictly positive parameters (*m_v*, *c_v*).  
- Parallelise chains via `multiprocessing` for speed‑up on multi‑core CPUs.  
- Add convergence metrics (Gelman–Rubin R̂, effective sample size).  
- Compare **single vs double drainage** and **radial consolidation** extensions.  
- Integrate **measurement‑error models** (e.g., heteroscedastic noise).

---

## 📚 Reference

Kelly, R. B., & Huang, J. (2015). *Bayesian updating of consolidation parameters from field measurements.*  
**Computers and Geotechnics**, 69, 496–507.  
https://doi.org/10.1016/j.compgeo.2015.06.010

---

## 🧾 Citation

If you use this repository or adapt the workflow, please cite:
> Jones, M. (2025). *Bayesian Updating for One‑Dimensional Consolidation (Kelly & Huang, 2015) [Code repository].* GitHub.  
> https://github.com/MerrickJones/2015KellyHuang

---

## 👤 Author

**Merrick Jones (2025)**  
PhD Candidate, University of Newcastle  
*(Bayesian Back‑Analysis for Embankments on Soft Soils)*  
📧 merrick.jones@uon.edu.au

---

## 🪪 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
