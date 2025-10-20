# DREAM_Suite.py - FIXED VERSION
import numpy as np
from scipy.stats import norm
import importlib

def DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options):
    """Main DREAM function - FIXED!"""
    # Import the model function
    module_name, func_name = Func_name.split('.')
    module = importlib.import_module(module_name)
    model_func = getattr(module, func_name)
    
    d, N, T = DREAMPar['d'], DREAMPar['N'], DREAMPar['T']
    min_b, max_b = Par_info['min'], Par_info['max']
    Y = Meas_info['Y']
    sigma = 0.05
    
    # Latin hypercube initialization
    chain = np.zeros((N, T+1, d))
    for i in range(d):
        chain[:, 0, i] = np.random.uniform(min_b[i], max_b[i], N)
    
    print("Running DREAM MCMC...")
    for t in range(1, T+1):
        for i in range(N):
            # Differential Evolution proposal
            others = [k for k in range(N) if k != i]
            j1, j2 = np.random.choice(others, 2, replace=False)
            gamma = 2.38 / np.sqrt(2 * d)
            Xp = chain[i, t-1] + gamma * (chain[j1, t-1] - chain[j2, t-1])
            Xp = np.clip(Xp, min_b, max_b)
            
            # Metropolis
            logL_curr = np.sum(norm.logpdf(Y, model_func(chain[i, t-1]), sigma))
            logL_prop = np.sum(norm.logpdf(Y, model_func(Xp), sigma))
            alpha = min(1, np.exp(logL_prop - logL_curr))
            
            if np.random.rand() < alpha:
                chain[i, t] = Xp
            else:
                chain[i, t] = chain[i, t-1]
        
        if t % 1000 == 0:
            print(f"Progress: {t}/{T}")
    
    # Dummy outputs
    output, FX, Z, logL = None, None, None, None
    return chain, output, FX, Z, logL