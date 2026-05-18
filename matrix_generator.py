import numpy as np
import pandas as pd
from numpy.linalg import cond, norm
from scipy.linalg import eigh
import seaborn as sns
import matplotlib.pyplot as plt 

def generate_matrix_pair(n, condition_number =1e2, scale=1.0,seed=None):
    """
    Generate a pair of symmetric positive definite matrices (A,B) of size n x n.

    Parameters:
        n (int): Size of matrices.
        condition_number (float): Desired condition number.
        scale (float): Scaling factor to multiply the matrices.
        seed (int, optional): Random seed for reproducibility.

    Returns:
        A, B (ndarray, ndarray): Symmetric positive definite matrices.
    """
    if seed is not None:
        np.random.seed(seed)
    
    #random orthogonal matrix
    Q,_ = np.linalg.qr(np.random.randn(n,n))
    eigs = np.geomspace(1, condition_number, n)
    A = Q @ np.diag(eigs) @ Q.T * scale

    R,_= np.linalg.qr(np.random.randn(n,n))
    eigs_B = np.geomspace(1, condition_number, n)
    B = R @ np.diag(eigs_B) @ R.T * scale

    return A, B 

def generate_dataset_entry(n, cond_val, scale_val, seed=None):
    """
    Generate a single row of dataset features + eigenvalue error.
    """

    A, B = generate_matrix_pair(n, cond_val, scale=scale_val, seed=seed)
    eigvals_true = eigh(A, B, eigvals_only=True)

    A_pert = np.round(A, decimals=6)
    B_pert = np.round(B, decimals=6)

    eigvals_pert = eigh(A_pert, B_pert, eigvals_only=True)
    eig_error = norm(eigvals_true - eigvals_pert)
    eigvals_true_norm = norm(eigvals_true)
    eig_error_relative = eig_error / eigvals_true_norm if eigvals_true_norm != 0 else np.nan
    log_cond_A = np.log10(cond(A)) if cond(A) > 0 else np.nan
    
    return{
        "size": n,
        "scale": scale_val,
        "cond_input": cond_val,
        "cond_A_actual": cond(A),
        "cond_B_actual": cond(B),
        "fro_norm_A": norm(A, 'fro'),
        "fro_norm_B": norm(B, 'fro'),
        "max_A": np.max(np.abs(A)),
        "max_B": np.max(np.abs(B)),
        "eig_error": eig_error,
        "eig_error_relative": eig_error_relative,
        "log_cond_A": log_cond_A
    }

# --------- Generating full dataset -------------

sizes = [100, 500, 1000, 1500] 
conds = [1e2, 1e4, 1e6 ,1e8, 1e12]
scales = [1e-2, 1.0, 1e2] # Vary s to see the effect of fixed-decimal rounding
num_trials = 10

dataset = []

for n in sizes: 
    for c in conds:
        for s in scales:
            for trial in range(num_trials):
                row = generate_dataset_entry(n, c, scale_val=s, seed=trial)
                dataset.append(row)

df = pd.DataFrame(dataset)

#saving dataset to csv 
df.to_csv("matrix_dataset.csv", index=False)
print("Dataset saved as 'matrix_dataset.csv'")

