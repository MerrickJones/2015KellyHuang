# kellyFORTRAN.py
import numpy as np

def kellyFORTRAN(x):
    """Kelly & Huang 2015 Consolidation Model"""
    # Fixed Hf = 3.0m (adjusted to match max settlement ~0.51m)
    Hf = 3.0
    
    # Parameters: x[0]=mv (1/kPa), x[1]=H (m), x[2]=γf (kN/m3), x[3]=cv (m2/year)
    mv, H, gamma_f, cv = x
    
    # Time from Excel (years)
    t_obs = np.array([0.01, 0.02, 0.05, 0.09, 0.2, 0.4, 0.8, 0.95, 1.1, 1.3])
    
    # Load σ = γf * Hf
    sigma = gamma_f * Hf
    
    # Settlement s(t) from Eq (4)
    s = np.zeros_like(t_obs)
    for i, t in enumerate(t_obs):
        Tv = cv * t / H**2
        sum_term = 0
        for m in range(10):  # Sum limited to 10 as in paper
            M = np.pi * (2*m + 1) / 2
            sum_term += (2 / M**2) * np.exp(-M**2 * Tv)
        s[i] = mv * sigma * H * (1 - sum_term)
    
    return s