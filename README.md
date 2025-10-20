# Bayesian Updating for One-Dimensional Consolidation (Kelly & Huang, 2015)

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

