## ################################################################################## ##
## PROGRESSIVE DREAM MCMC FOR KELLY & HUANG 2015 - SETTLEMENT PREDICTION EXPERIMENT  ##
## Adapted from example_consolidation.py for progressive data assimilation  
## Implimented by Merrick Jones 18 October 2025                                       ##
## ################################################################################## ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches

# Ensure paths resolve from this script's folder
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Add DREAM_Suite to import path
sys.path.insert(0, script_dir)
from DREAM_Suite import DREAM_Suite

# ===========================================
# CORE FORWARD MODEL (multi-arg)
# ===========================================
def kelly_core(x, n_obs, t_obs_full):
    """
    Kelly & Huang 2015 consolidation model (Eq. 4) with mild numerical stabilization.
    x = [m_v, H, gamma_f, c_v]
    Returns settlement vector corresponding to first n_obs times in t_obs_full.
    """
    H_f = 3.0  # fixed fill thickness from the paper
    m_v, H, gamma_f, c_v = x
    t_obs = np.array(t_obs_full[:n_obs], dtype=float)
    sigma = gamma_f * H_f
    s = np.zeros_like(t_obs, dtype=float)

    for i, t in enumerate(t_obs):
        T_v = c_v * t / (H**2)
        if T_v < 1e-10:
            T_v = 1e-10
        sum_term = 0.0
        for m in range(10):  # truncated series
            M = (np.pi / 2.0) * (2 * m + 1)
            exp_term = np.exp(-M**2 * T_v)
            if exp_term < 1e-15:  # negligible
                break
            sum_term += (2.0 / M**2) * exp_term
        s[i] = m_v * sigma * H * (1.0 - sum_term)
    return s

# ===========================================
# ADAPTER FOR DREAM_Suite (expects 1-arg function)
# ===========================================
# We expose *kellyFORTRAN(theta)* as the callable name DREAM_Suite uses,
# and internally forward to kelly_core using two small adapter settings.
_ADAPTER_N_OBS = None
_T_OBS_FULL = None

def kellyFORTRAN(theta):
    """
    One-argument adapter so DREAM_Suite can call kellyFORTRAN(theta).
    Forward to kelly_core(theta, _ADAPTER_N_OBS, _T_OBS_FULL).
    """
    if _ADAPTER_N_OBS is None or _T_OBS_FULL is None:
        raise RuntimeError("kellyFORTRAN adapter not initialized with n_obs/t_obs_full")
    return kelly_core(theta, _ADAPTER_N_OBS, _T_OBS_FULL)

# ===========================================
# MEASUREMENT DATA
# ===========================================
D_obs = np.array([0.095, 0.125, 0.19, 0.26, 0.37, 0.45, 0.5, 0.51, 0.51, 0.51])  # observed settlements (m)
np.savetxt('settlement_obs.txt', D_obs, fmt='%.4f')
t_obs_full = np.array([0.01, 0.02, 0.05, 0.09, 0.2, 0.4, 0.8, 0.95, 1.1, 1.3])  # times (years)
N_FULL = len(t_obs_full)

# True reference curve (call the core directly)
true_theta = [0.0014, 5.5, 22.0, 80.0]  # m_v, H, γ_f, c_v
S_true = kelly_core(true_theta, N_FULL, t_obs_full)

# ===========================================
# PROBLEM SETTINGS
# ===========================================
DREAMPar = {}
DREAMPar['d'] = 4  # parameters: m_v, H, γ_f, c_v
DREAMPar['lik'] = 11  # Gaussian likelihood

# Priors: Normal (mean from Table 2, std = mean * COV)
prior_means = [0.001, 5.0, 20.0, 40.0]
prior_covs  = [0.3,   0.1,  0.1,  0.5]
prior_stds  = np.array(prior_means) * np.array(prior_covs)

Par_info = {}
Par_info['initial'] = 'latin'
Par_info['boundhandling'] = 'reflect'
Par_info['min'] = np.array(prior_means) - 3 * prior_stds
Par_info['max'] = np.array(prior_means) + 3 * prior_stds

# Convenient named priors
prior_m_v, prior_H, prior_gamma_f, prior_c_v = prior_means

# DREAM configuration: point to the callable *kellyFORTRAN* (adapter name)
Func_name = '__main__.kellyFORTRAN'
options = {'parallel': 'no', 'IO': 'no', 'modout': 'yes'}
method = 'dream'
DREAMPar['N'] = 10
DREAMPar['T'] = 20000  # increase for better convergence if needed

# ===========================================
# EXPERIMENT: PROGRESSIVE MCMC SIMULATIONS
# ===========================================
def run_mcmc(n_obs):
    print(f"\nRunning DREAM MCMC with {n_obs} observation(s)...")

    # Initialize adapter for this run
    global _ADAPTER_N_OBS, _T_OBS_FULL
    _ADAPTER_N_OBS = n_obs
    _T_OBS_FULL = t_obs_full

    # DREAM_Suite expects Meas_info with observed data under 'Y'
    Meas_info = {'Y': D_obs[:n_obs]}

    # Run DREAM
    chain, output, FX, Z, logL = DREAM_Suite(
        method, Func_name, DREAMPar, Par_info, Meas_info, options
    )

    # Post-processing: posterior mean after burn-in
    burn_in = DREAMPar['T'] // 2
    samples = chain[:, burn_in:, :].reshape(-1, DREAMPar['d'])
    post_mean = np.mean(samples, axis=0)

    # Predict full curve for plotting/review (use core directly)
    S_pred_full = kelly_core(post_mean, N_FULL, t_obs_full)

    # Return a thinned sample set for plotting CIs etc. (keeps memory light)
    samples_thin = samples[::10].copy()
    return post_mean, S_pred_full, samples_thin

predictions = []
posteriors = []
samples_by_nobs = []  # store thinned samples for each n_obs
for n_obs in range(1, N_FULL + 1):
    post_mean, S_pred_full, samples_thin = run_mcmc(n_obs)
    posteriors.append(post_mean)
    predictions.append(S_pred_full)
    samples_by_nobs.append(samples_thin)
    print(f"Posterior parameters (m_v, c_v): {post_mean[0]:.4f}, {post_mean[3]:.1f}")

# ===========================================
# FIGURE 1 (unchanged) — and SAVE to disk
# ===========================================
fig, ax = plt.subplots(figsize=(12, 8))

# Plot true settlement
ax.plot(t_obs_full, S_true, 'k--', label='True Settlement', linewidth=2)

# Plot prior prediction (full curve)
prior_params = [prior_m_v, prior_H, prior_gamma_f, prior_c_v]
S_prior = kelly_core(prior_params, N_FULL, t_obs_full)
ax.plot(t_obs_full, S_prior, 'r--', label='Prior Prediction', linewidth=2)

# Plot progressive predictions (each labelled individually)
colors = plt.cm.viridis(np.linspace(0, 1, N_FULL))
for i, (S_pred, n) in enumerate(zip(predictions, range(1, N_FULL + 1))):
    label = f"{n} Observation(s)"  # label each curve
    ax.plot(t_obs_full, S_pred, label=label, linewidth=2, alpha=0.7, color=colors[i])

# Customize plot
ax.set_xlabel('Time (years)', fontsize=12)
ax.set_ylabel('Settlement (m)', fontsize=12)
ax.set_title('PROGRESSIVE SETTLEMENT PREDICTIONS\nKelly & Huang 2015', fontsize=14, weight='bold', pad=20)
ax.invert_yaxis()
ax.set_xscale('log')
ax.legend(fontsize=9, loc='best', ncol=2)
ax.grid(True, alpha=0.3)

plt.tight_layout()
fig1_path = os.path.join(script_dir, 'Figure_1.png')
plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved Figure 1 to: {fig1_path}")

# ===========================================
# ANALYSIS: WHEN IS THE PREDICTION "REASONABLY CLOSE"?
#   -> Find first n_obs with error < 5% and record its index
# ===========================================
conv_index = None  # index in 0-based python terms; n_obs = conv_index + 1
for i, post_mean in enumerate(posteriors):
    n_obs_i = i + 1
    S_pred_full = kelly_core(post_mean, N_FULL, t_obs_full)
    pred_sett_final = S_pred_full[-1]
    true_sett_final = S_true[-1]
    error_percent = abs(pred_sett_final - true_sett_final) / true_sett_final * 100.0
    print(f"With {n_obs_i} observations:")
    print(f"  Posterior m_v: {post_mean[0]:.4f}, c_v: {post_mean[3]:.1f}")
    print(f"  Predicted settlement at t=1.3 years: {pred_sett_final:.3f} m")
    print(f"  True settlement at t=1.3 years: {true_sett_final:.3f} m")
    print(f"  Error: {error_percent:.1f}%")
    if conv_index is None and error_percent < 5.0:
        conv_index = i
        print(f"  ✅ Prediction is reasonably close (<5%) with {n_obs_i} observations!")
        # Do not break; keep printing the rest, but we remember the *first* convergence
# If none met the threshold, fall back to final run
if conv_index is None:
    conv_index = len(posteriors) - 1
conv_n_obs = conv_index + 1

# ===========================================
# SAVE RESULTS
# ===========================================
results_txt_path = os.path.join(script_dir, 'progress_SETT_results.txt')
with open(results_txt_path, 'w', encoding='utf-8') as f:
    f.write("PROGRESSIVE SETTLEMENT PREDICTION RESULTS - KELLY & HUANG 2015\n")
    f.write("=====================================\n")
    for i, post_mean in enumerate(posteriors):
        n_obs_i = i + 1
        S_pred_full = kelly_core(post_mean, N_FULL, t_obs_full)
        pred_sett_final = S_pred_full[-1]
        true_sett_final = S_true[-1]
        error_percent = abs(pred_sett_final - true_sett_final) / true_sett_final * 100.0
        f.write(f"With {n_obs_i} observations:\n")
        f.write(f"  Posterior m_v: {post_mean[0]:.4f}, c_v: {post_mean[3]:.1f}\n")
        f.write(f"  Predicted settlement at t=1.3 years: {pred_sett_final:.3f} m\n")
        f.write(f"  True settlement at t=1.3 years: {true_sett_final:.3f} m\n")
        f.write(f"  Error: {error_percent:.1f}%\n")
        if i == conv_index:
            f.write(f"  ✅ Prediction is reasonably close (<5%) with {n_obs_i} observations!\n")
    f.write("=====================================\n")
print(f"Wrote results to: {results_txt_path}")

# ===========================================
# FIGURE 2 (four-panel summary) — using the *convergence run* (first error < 5%)
#   - Top Left: Bayesian update (Prior → Posterior) with 95% CI
#   - Top Right: Embankment schematic
#   - Bottom Left: Posterior of c_v
#   - Bottom Right: Posterior of m_v
# ===========================================
conv_post = posteriors[conv_index]
post_m_v, post_H, post_gamma_f, post_c_v = conv_post
conv_samples = samples_by_nobs[conv_index]  # thinned samples for that run

# Build credible interval curves for full time vector using posterior samples at convergence
S_samples = np.array([kelly_core(theta, N_FULL, t_obs_full) for theta in conv_samples])
ci_low, ci_high = np.percentile(S_samples, [2.5, 97.5], axis=0)
S_post_full = kelly_core(conv_post, N_FULL, t_obs_full)
S_prior_full = kelly_core([prior_m_v, prior_H, prior_gamma_f, prior_c_v], N_FULL, t_obs_full)

fig2, axs = plt.subplots(2, 2, figsize=(15, 10))

# --- TOP LEFT: Bayesian Update Plot ---
axs[0, 0].plot(t_obs_full, D_obs, 'o', label='Observed', markersize=6, color='black')
axs[0, 0].plot(t_obs_full, S_prior_full, '--', label=f'Prior (c_v={prior_c_v:.1f})', linewidth=2.5, color='red')
axs[0, 0].plot(t_obs_full, S_post_full, '-', label=f'Posterior (c_v={post_c_v:.1f})', linewidth=2.5, color='green')
axs[0, 0].fill_between(t_obs_full, ci_low, ci_high, alpha=0.3, color='blue', label='95% CI')
axs[0, 0].set_xlabel('Time (years)', fontsize=12)
axs[0, 0].set_ylabel('Settlement (m)', fontsize=12)
axs[0, 0].set_title(f'BAYESIAN UPDATE (n_obs={conv_n_obs})\nPrior → Posterior', fontsize=14, weight='bold', pad=20)
axs[0, 0].invert_yaxis()
axs[0, 0].set_xscale('log')
axs[0, 0].legend(fontsize=10)
axs[0, 0].grid(True, alpha=0.3)

# --- TOP RIGHT: Modern Embankment Schematic (clean 2025 style) ---
axschem = axs[0, 1]
soil_cmap = LinearSegmentedColormap.from_list("soil_grad", ["#D68910", "#E59866", "#FAD7A0"])
axschem.imshow(np.linspace(1, 0, 100).reshape(-1, 1),
               cmap=soil_cmap, extent=[0, 10, -5.5, 0],
               aspect='auto', alpha=0.8)  # Gradient soil
axschem.add_patch(
    patches.Polygon([[1.5, 1.5], [8.5, 1.5], [9.5, 0], [0.5, 0]],
                    fc='#85C1E9', ec='#3498DB', lw=2, alpha=0.9)
)
axschem.plot([0, 10], [-5.5, -5.5], color='#17202A', lw=4, linestyle='--')  # base line

# Subtle drainage indicator
for y in np.linspace(-5.5, 0, 10):
    axschem.plot([5, 5], [y, y + 0.3], color='#AED6F1', alpha=0.3, lw=1)
axschem.text(5, -5.75, 'Impermeable', ha='center', fontsize=10, color='#17202A')
axschem.text(5, 0.2, 'Drained', ha='center', fontsize=10, color='#3498DB')

axschem.arrow(5, 2, 0, -0.8, fc='#E74C3C', ec='#E74C3C', lw=3, head_width=0.4, head_length=0.3)
axschem.text(5.5, 1.5, 'σ', color='#E74C3C', ha='center', weight='bold', fontsize=12)

text_params = (
    f"PRIOR → POSTERIOR\n"
    f"m_v: {prior_m_v:.4f} → {post_m_v:.4f}\n"
    f"H: {prior_H:.1f} → {post_H:.1f}\n"
    f"γ_f: {prior_gamma_f:.1f} → {post_gamma_f:.1f}\n"
    f"c_v: {prior_c_v:.1f} → {post_c_v:.1f}"
)
axschem.text(11, 0, text_params, va='top',
             bbox=dict(fc='white', ec='black', pad=1), fontsize=11, weight='bold')
axschem.set_xlim(0, 16)
axschem.set_ylim(-6.5, 2.5)
axschem.set_aspect('equal')
axschem.axis('off')
axschem.set_title('EMBANKMENT SCHEMATIC', fontsize=14, weight='bold', pad=20)

# --- BOTTOM LEFT: c_v Posterior ---
axs[1, 0].hist(conv_samples[:, 3], bins=30, alpha=0.7, color='#3498DB',
               density=True, edgecolor='white')
axs[1, 0].axvline(prior_c_v, color='red', ls='--', lw=2, label=f'Prior: {prior_c_v:.1f}')
axs[1, 0].axvline(post_c_v, color='green', lw=3, label=f'Posterior: {post_c_v:.1f}')
axs[1, 0].axvline(true_theta[3], color='black', ls=':', lw=2, label=f'True: {true_theta[3]:.1f}')
axs[1, 0].set_xlabel('c_v (m²/year)', fontsize=12)
axs[1, 0].set_ylabel('Density', fontsize=12)
axs[1, 0].set_title('POSTERIOR: CONSOLIDATION\nCOEFFICIENT (c_v)', fontsize=14, weight='bold', pad=20)
axs[1, 0].legend(fontsize=10)
axs[1, 0].grid(True, alpha=0.3)

# --- BOTTOM RIGHT: m_v Posterior ---
axs[1, 1].hist(conv_samples[:, 0], bins=30, alpha=0.7, color='#E74C3C',
               density=True, edgecolor='white')
axs[1, 1].axvline(prior_m_v, color='red', ls='--', lw=2, label=f'Prior: {prior_m_v:.4f}')
axs[1, 1].axvline(post_m_v, color='green', lw=3, label=f'Posterior: {post_m_v:.4f}')
axs[1, 1].axvline(true_theta[0], color='black', ls=':', lw=2, label=f'True: {true_theta[0]:.4f}')
axs[1, 1].set_xlabel('m_v (1/kPa)', fontsize=12)
axs[1, 1].set_ylabel('Density', fontsize=12)
axs[1, 1].set_title('POSTERIOR: VOLUME\nCOMPRESSIBILITY (m_v)', fontsize=14, weight='bold', pad=20)
axs[1, 1].legend(fontsize=10)
axs[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.suptitle(f'KELLY & HUANG 2015: BAYESIAN BACK-ANALYSIS (n_obs={conv_n_obs})',
             fontsize=18, weight='bold', y=1.02)
fig2_path = os.path.join(script_dir, 'Figure_2.png')
plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved Figure 2 to: {fig2_path}")

# ===========================================
# FIGURE 3 (save a PNG image rendering of progress_SETT_results.txt)
# ===========================================
try:
    with open(results_txt_path, 'r', encoding='utf-8') as fh:
        txt = fh.read()
        txt = txt.replace('✅', '✔')  # prevents UserWarning about missing glyph
except Exception as e:
    txt = f"Could not open {results_txt_path}\nError: {e}"

fig3 = plt.figure(figsize=(11, 8))  # landscape-ish
plt.text(0.01, 0.99, txt, va='top', ha='left', fontsize=10, family='monospace')
plt.axis('off')
plt.title('Figure 3 — progress_SETT_results.txt', fontsize=14, weight='bold', pad=12)
fig3_path = os.path.join(script_dir, 'Figure_3_progress_SETT_results.png')
plt.savefig(fig3_path, dpi=300, bbox_inches='tight')
plt.show()
print(f"Saved Figure 3 to: {fig3_path}")
