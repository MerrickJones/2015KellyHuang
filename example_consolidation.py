## ################################################################################## ##
## DREAM MCMC FOR KELLY & HUANG 2015 - BAYESIAN UPDATING FOR CONSOLIDATION          ##
## Adapted from Jasper Vrugt's example_5.py for Terzaghi 1D Consolidation            ##
## ################################################################################## ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import norm
from matplotlib.colors import LinearSegmentedColormap

# Ensure paths resolve from this script's folder
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Add DREAM_Suite to import path
sys.path.insert(0, script_dir)
from DREAM_Suite import DREAM_Suite

# ===========================================
# KELLY & HUANG FORWARD MODEL - EQ (4) FOR SETTLEMENT
# ===========================================
def kellyFORTRAN(x):
    """Kelly & Huang 2015 Consolidation Model - Exact Eq (4) with numerical stabilization"""
    # Fixed H_f = 3.0 m (layer thickness from paper)
    H_f = 3.0
    
    # Parameters: x[0]=m_v (1/kPa), x[1]=H (m), x[2]=γ_f (kN/m3), x[3]=c_v (m2/year)
    m_v, H, gamma_f, c_v = x
    
    # Time from Excel (years)
    t_obs = np.array([0.01, 0.02, 0.05, 0.09, 0.2, 0.4, 0.8, 0.95, 1.1, 1.3])
    
    # Load σ = γ_f * H_f
    sigma = gamma_f * H_f
    
    # Settlement s(t) from Eq (4)
    s = np.zeros_like(t_obs)
    for i, t in enumerate(t_obs):
        T_v = c_v * t / H**2
        if T_v < 1e-10:
            T_v = 1e-10
        sum_term = 0
        for m in range(10):  # m from 0 to 9
            M = (np.pi / 2) * (2 * m + 1)
            exp_term = np.exp(-M**2 * T_v)
            if exp_term < 1e-15:
                break
            sum_term += (2 / M**2) * exp_term
        s[i] = m_v * sigma * H * (1 - sum_term)
    
    return s

# ===========================================
# MEASUREMENT DATA FROM EXCEL
# ===========================================
D_obs = np.array([0.095, 0.125, 0.19, 0.26, 0.37, 0.45, 0.5, 0.51, 0.51, 0.51])  # Settlement (m)
np.savetxt('settlement_obs.txt', D_obs, fmt='%.4f')

# True values from Table 1 (for reference)
true_theta = [0.0014, 5.5, 22.0, 80.0]  # m_v, H, γ_f, c_v

# ===========================================
# PROBLEM SETTINGS
# ===========================================
DREAMPar = {}
DREAMPar['d'] = 4          # 4 parameters: m_v, H, γ_f, c_v
DREAMPar['lik'] = 11       # Gaussian likelihood

# Normal distributions for all parameters (mean from Table 2, std = mean * COV)
prior_means = [0.001, 5.0, 20.0, 40.0]
prior_covs = [0.3, 0.1, 0.1, 0.5]  # COVs from Table 2
prior_stds = np.array(prior_means) * np.array(prior_covs)
Par_info = {}
Par_info['initial'] = 'latin'
Par_info['boundhandling'] = 'reflect'
Par_info['min'] = np.array(prior_means) - 3 * prior_stds
Par_info['max'] = np.array(prior_means) + 3 * prior_stds

# Function name
Func_name = '__main__.kellyFORTRAN'

# Measured settlement data
Meas_info = {}
Meas_info['Y'] = np.loadtxt('settlement_obs.txt')

# Options
options = {}
options['parallel'] = 'no'
options['IO'] = 'no'
options['modout'] = 'yes'

# Method
method = 'dream'

# Chain settings
DREAMPar['N'] = 10
DREAMPar['T'] = 20000  # Increased for better convergence

print("🚀 DREAM MCMC: Kelly & Huang 2015 Bayesian Updating")
print(f"True parameters: m_v={true_theta[0]:.4f}, H={true_theta[1]:.1f}, γ_f={true_theta[2]:.1f}, c_v={true_theta[3]:.1f}")
print(f"Observations: {len(D_obs)} settlement points")

# ===========================================
# RUN DREAM
# ===========================================
chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)

# ===========================================
# POST-PROCESSING
# ===========================================
print("\n✅ DREAM COMPLETED! Extracting posterior...")
burn_in = DREAMPar['T'] // 2
samples = chain[:, burn_in:, :].reshape(-1, DREAMPar['d'])
post_mean = np.mean(samples, axis=0)
post_m_v, post_H, post_gamma_f, post_c_v = post_mean

# Prior means from Table 2
prior_m_v = 0.001
prior_H = 5.0
prior_gamma_f = 20.0
prior_c_v = 40.0

# Predictions
t_obs = np.array([0.01, 0.02, 0.05, 0.09, 0.2, 0.4, 0.8, 0.95, 1.1, 1.3])
S_post = kellyFORTRAN(post_mean)
S_prior = kellyFORTRAN([prior_m_v, prior_H, prior_gamma_f, prior_c_v])
S_samples = np.array([kellyFORTRAN(s) for s in samples[::10]])
ci_low, ci_high = np.percentile(S_samples, [2.5, 97.5], axis=0)

# ===========================================
# UPDATED FIGURE LAYOUT
# ===========================================
fig, axs = plt.subplots(2, 2, figsize=(15, 10))

# TOP LEFT: Bayesian Update Plot
axs[0,0].plot(t_obs, D_obs, 'o', label='Observed', markersize=8, color='black')
axs[0,0].plot(t_obs, S_prior, '--', label=f'Prior (c_v={prior_c_v:.1f})', linewidth=3, color='red')
axs[0,0].plot(t_obs, S_post, '-', label=f'Posterior (c_v={post_c_v:.1f})', linewidth=3, color='green')
axs[0,0].fill_between(t_obs, ci_low, ci_high, alpha=0.3, color='blue', label='95% CI')
axs[0,0].set_xlabel('Time (years)', fontsize=12)
axs[0,0].set_ylabel('Settlement (m)', fontsize=12)
axs[0,0].set_title('BAYESIAN UPDATE\nPrior → Posterior', fontsize=14, weight='bold', pad=20)
axs[0,0].invert_yaxis()  # Reverses y-axis for settlement
axs[0,0].set_xscale('log')  # Logarithmic x-axis
axs[0,0].legend(fontsize=10)
axs[0,0].grid(True, alpha=0.3)

# TOP RIGHT: Modern Embankment Schematic (2025 Style)
ax = axs[0,1]
soil_cmap = LinearSegmentedColormap.from_list("soil_grad", ["#D68910", "#E59866", "#FAD7A0"])
ax.imshow(np.linspace(1, 0, 100).reshape(-1, 1), cmap=soil_cmap, extent=[0, 10, -5.5, 0], aspect='auto', alpha=0.8)  # Gradient soil
ax.add_patch(patches.Polygon([[1.5,1.5],[8.5,1.5],[9.5,0],[0.5,0]], fc='#85C1E9', ec='#3498DB', lw=2, alpha=0.9))  # Modern embankment
ax.plot([0,10], [-5.5,-5.5], color='#17202A', lw=4, linestyle='--')  # Modern base line

# Subtle drainage indicator: faint upward gradient lines
for y in np.linspace(-5.5, 0, 10):
    ax.plot([5, 5], [y, y+0.3], color='#AED6F1', alpha=0.3, lw=1)
ax.text(5, -5.75, 'Impermeable', ha='center', fontsize=10, color='#17202A')
ax.text(5, 0.2, 'Drained', ha='center', fontsize=10, color='#3498DB')

ax.arrow(5,2,0,-0.8, fc='#E74C3C', ec='#E74C3C', lw=3, head_width=0.4, head_length=0.3)
ax.text(5.5,1.5,'σ', color='#E74C3C', ha='center', weight='bold', fontsize=12)
text = f"PRIOR → POSTERIOR\nm_v: {prior_m_v:.4f} → {post_m_v:.4f}\nH: {prior_H:.1f} → {post_H:.1f}\nγ_f: {prior_gamma_f:.1f} → {post_gamma_f:.1f}\nc_v: {prior_c_v:.1f} → {post_c_v:.1f}"
ax.text(11,0,text, va='top', bbox=dict(fc='white', ec='black', pad=1), fontsize=11, weight='bold')
ax.set_xlim(0,16)
ax.set_ylim(-6.5,2.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('EMBANKMENT SCHEMATIC', fontsize=14, weight='bold', pad=20)

# BOTTOM LEFT: c_v Posterior
axs[1,0].hist(samples[:,3], bins=30, alpha=0.7, color='#3498DB', density=True, edgecolor='white')
axs[1,0].axvline(prior_c_v, color='red', ls='--', lw=2, label=f'Prior: {prior_c_v:.1f}')
axs[1,0].axvline(post_c_v, color='green', lw=3, label=f'Posterior: {post_c_v:.1f}')
axs[1,0].axvline(true_theta[3], color='black', ls=':', lw=2, label=f'True: {true_theta[3]:.1f}')
axs[1,0].set_xlabel('c_v (m²/year)', fontsize=12)
axs[1,0].set_ylabel('Density', fontsize=12)
axs[1,0].set_title('POSTERIOR: CONSOLIDATION\nCOEFFICIENT (c_v)', fontsize=14, weight='bold', pad=20)
axs[1,0].legend(fontsize=10)
axs[1,0].grid(True, alpha=0.3)

# BOTTOM RIGHT: m_v Posterior
axs[1,1].hist(samples[:,0], bins=30, alpha=0.7, color='#E74C3C', density=True, edgecolor='white')
axs[1,1].axvline(prior_m_v, color='red', ls='--', lw=2, label=f'Prior: {prior_m_v:.4f}')
axs[1,1].axvline(post_m_v, color='green', lw=3, label=f'Posterior: {post_m_v:.4f}')
axs[1,1].axvline(true_theta[0], color='black', ls=':', lw=2, label=f'True: {true_theta[0]:.4f}')
axs[1,1].set_xlabel('m_v (1/kPa)', fontsize=12)
axs[1,1].set_ylabel('Density', fontsize=12)
axs[1,1].set_title('POSTERIOR: VOLUME\nCOMPRESSIBILITY (m_v)', fontsize=14, weight='bold', pad=20)
axs[1,1].legend(fontsize=10)
axs[1,1].grid(True, alpha=0.3)

plt.tight_layout()
plt.suptitle('KELLY & HUANG 2015: BAYESIAN BACK-ANALYSIS', fontsize=18, weight='bold', y=1.02)
plt.show()

# ===========================================
# RESULTS SUMMARY - PRINT TO CONSOLE AND FILE
# ===========================================
summary = f"""=================================================================
🎉 BAYESIAN BACK-ANALYSIS RESULTS - KELLY & HUANG 2015
=================================================================
Parameter        Prior        Posterior      True      Accuracy
mv (1/kPa)     {prior_m_v:10.4f}   {post_m_v:10.4f}   {true_theta[0]:8.4f}   {(post_m_v/true_theta[0]*100-100):+8.1f}%
H (m)          {prior_H:10.1f}   {post_H:10.1f}   {true_theta[1]:8.1f}   {(post_H/true_theta[1]*100-100):+8.1f}%
γf (kN/m3)     {prior_gamma_f:10.1f}   {post_gamma_f:10.1f}   {true_theta[2]:8.1f}   {(post_gamma_f/true_theta[2]*100-100):+8.1f}%
cv (m2/year)   {prior_c_v:10.1f}   {post_c_v:10.1f}   {true_theta[3]:8.1f}   {(post_c_v/true_theta[3]*100-100):+8.1f}%
================================================================="""
print(summary)

# Write summary to text file with UTF-8 encoding
with open(os.path.join(script_dir, 'back_analysis_results.txt'), 'w', encoding='utf-8') as f:
    f.write(summary)

# ===========================================
# ADDITIONAL SUMMARY FIGURE (POPS UP AFTER MAIN FIGURE)
# ===========================================
fig_summary = plt.figure(figsize=(10, 6))
summary_text = f"""=================================================================
BAYESIAN BACK-ANALYSIS RESULTS - KELLY & HUANG 2015
=================================================================
Parameter        Prior        Posterior      True      Accuracy
mv (1/kPa)     {prior_m_v:10.4f}   {post_m_v:10.4f}   {true_theta[0]:8.4f}   {(post_m_v/true_theta[0]*100-100):+8.1f}%
H (m)          {prior_H:10.1f}   {post_H:10.1f}   {true_theta[1]:8.1f}   {(post_H/true_theta[1]*100-100):+8.1f}%
γf (kN/m3)     {prior_gamma_f:10.1f}   {post_gamma_f:10.1f}   {true_theta[2]:8.1f}   {(post_gamma_f/true_theta[2]*100-100):+8.1f}%
cv (m2/year)   {prior_c_v:10.1f}   {post_c_v:10.1f}   {true_theta[3]:8.1f}   {(post_c_v/true_theta[3]*100-100):+8.1f}%
================================================================="""
plt.text(0.5, 0.5, summary_text, ha='center', va='center', fontsize=12, family='monospace')
plt.axis('off')
plt.title('SUMMARY OF BAYESIAN BACK-ANALYSIS RESULTS', fontsize=14, weight='bold', pad=20)
plt.show(block=True)  # Ensures it pops up after main figure
plt.savefig(os.path.join(script_dir, 'back_analysis_summary.png'), dpi=300, bbox_inches='tight')
plt.close()