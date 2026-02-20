import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import math
import os

def get_file_path(base_path, strategy, p_str, k, seq_type, alpha=None, phi=None, file_type="results"):
    """
    Costruisce i path seguendo esattamente la struttura delle cartelle e i nomi file del codice C++.
    """
    if strategy == "standard":
        if file_type == "results":
            filename = f"results_credits_{p_str}_{k}.csv"
        elif file_type == "density":
            filename = f"density_credits_{p_str}_{k}.csv"
        else: # times
            filename = f"times_credits_{p_str}_{k}.csv"
        return os.path.join(base_path, filename)
    else:
        # Per mixed, il C++ usa std::to_string per alpha e phi (solitamente 6 decimali)
        # Esempio: results_credits_power_0.2_alpha:0.500000_phi:0.600000_k:128.csv
        a_str = f"{alpha:.6f}" if alpha is not None else "0.500000"
        p_val = f"{phi:.6f}" if phi is not None else "0.600000"
        
        if file_type == "results":
            filename = f"results_credits_{seq_type}_{p_str}_alpha:{a_str}_phi:{p_val}_k:{k}.csv"
        elif file_type == "density":
            # Nota: il C++ per la densità mixed non sembra includere alpha/phi nel nome file
            filename = f"density_credits_{seq_type}_{p_str}_{k}.csv"
        else: # times
            # Nota: il C++ ha un piccolo refuso nel codice originale (_phi invece di _phi:)
            filename = f"times_credits_{seq_type}_{p_str}_alpha:{a_str}_phi{p_val}_k:{k}.csv"
            
        return os.path.join(base_path, filename)

def plot_quasi_clique_results_grid(dataset_name, p_str, k_values, strategy="standard", seq_type="erdos", alpha=0.5, phi=0.6):
    base_dir = f"../Esperimenti/{dataset_name}/{strategy}/"
    
    # 1. Caricamento Baseline Dynamic
    if strategy == "standard":
        file_dyn = f"{base_dir}results_dynamic_{p_str}.csv"
        file_dyn_dens = f"{base_dir}density_dynamic_{p_str}.csv"
    else:
        file_dyn = f"{base_dir}results_dynamic_{seq_type}_{p_str}.csv"
        file_dyn_dens = f"{base_dir}density_dynamic_{seq_type}_{p_str}.csv"
    
    df_dyn = None
    df_dyn_dens = None
    try:
        df_dyn = pd.read_csv(file_dyn, header=None).mean(axis=1)
        if os.path.exists(file_dyn_dens):
            df_dyn_dens = pd.read_csv(file_dyn_dens, header=None).mean(axis=1)
    except Exception as e:
        print(f"Warning Baseline: {e}")

    # 2. Setup Grafico
    n_cols = 3
    n_rows = math.ceil(len(k_values) / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows), sharey=True, sharex=True, constrained_layout=True)
    axes_flat = axes.flatten() if len(k_values) > 1 else [axes]

    cmap_base = plt.get_cmap('RdYlGn')
    norm = mcolors.BoundaryNorm([0.0, 0.5, 0.7, 0.8, 0.9, 1.0], cmap_base.N)
    sc = None

    for idx, k in enumerate(k_values):
        ax = axes_flat[idx]
        
        # Plot Baseline (Nera con triangoli rivolti in basso)
        if df_dyn is not None:
            ax.plot(df_dyn, label='Fully-dynamic', color='black', linestyle='--', linewidth=1, alpha=0.4)
            if df_dyn_dens is not None:
                idx_dyn = np.linspace(0, len(df_dyn) - 1, len(df_dyn_dens), dtype=int)
                ax.scatter(idx_dyn, df_dyn.iloc[idx_dyn], c=df_dyn_dens, cmap=cmap_base, norm=norm, 
                           marker='v', s=35, edgecolors='black', linewidths=0.3, zorder=4)

        # Plot Credits (Blu con triangoli rivolti in alto)
        file_cred = get_file_path(base_dir, strategy, p_str, k, seq_type, alpha, phi, "results")
        file_dens = get_file_path(base_dir, strategy, p_str, k, seq_type, alpha, phi, "density")

        if os.path.exists(file_cred):
            df_cred = pd.read_csv(file_cred, header=None).mean(axis=1)
            ax.plot(df_cred, label=f'Incremental', color='steelblue', linewidth=1.5)
            
            if os.path.exists(file_dens):
                df_dens = pd.read_csv(file_dens, header=None).mean(axis=1)
                idx_cred = np.linspace(0, len(df_cred) - 1, len(df_dens), dtype=int)
                sc = ax.scatter(idx_cred, df_cred.iloc[idx_cred], c=df_dens, cmap=cmap_base, norm=norm, 
                                marker='^', s=45, edgecolors='black', linewidths=0.5, zorder=5)
        else:
            ax.set_title(f"k={k} (Mancante)", color='red')

        ax.set_title(f"k = {k}", fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.5)
        if idx == 0: ax.legend(loc='upper left', fontsize='small')

    # Pulizia subplots vuoti
    for i in range(len(k_values), len(axes_flat)):
        fig.delaxes(axes_flat[i])

    if sc is not None:
        cbar = fig.colorbar(sc, ax=axes, shrink=0.7, aspect=25, pad=0.02)
        cbar.set_label('Density', rotation=270, labelpad=15)

    title_suffix = f"({strategy.capitalize()}" + (f" - {seq_type})" if strategy=="mixed" else ")")
    fig.suptitle(f"{dataset_name} {title_suffix} | p={p_str} | α={alpha}, φ={phi}", fontsize=16, fontweight='bold')
    plt.show()

def plot_speedup(dataset_name, p_str, k_values, strategy="standard", seq_type="erdos", alpha=0.5, phi=0.6):
    base_dir = f"../Esperimenti/{dataset_name}/{strategy}/"
    file_dyn_time = f"{base_dir}times_dynamic_{p_str}.csv" if strategy == "standard" else f"{base_dir}times_dynamic_{seq_type}_{p_str}.csv"
    
    try:
        baseline_time = pd.read_csv(file_dyn_time, header=None).mean().mean()
    except:
        print(f"Errore: Baseline tempi non trovata in {file_dyn_time}")
        return

    speedups = []
    valid_ks = []

    for k in k_values:
        file_time = get_file_path(base_dir, strategy, p_str, k, seq_type, alpha, phi, "times")
        if os.path.exists(file_time):
            cred_time = pd.read_csv(file_time, header=None).mean().mean()
            speedups.append(baseline_time / cred_time)
            valid_ks.append(str(k))

    if not speedups: 
        print("Nessun dato di speedup trovato.")
        return

    plt.figure(figsize=(10, 5))
    bars = plt.bar(valid_ks, speedups, color='steelblue', edgecolor='black', zorder=3)
    plt.axhline(y=1.0, color='red', linestyle='--', label='Baseline Threshold', zorder=2)
    
    plt.ylabel("Speedup Multiplier (x)", fontweight='bold')
    plt.xlabel("k (MinHash Size)", fontweight='bold')
    plt.title(f"Speedup vs Baseline | {dataset_name} ({strategy})", fontweight='bold')
    plt.grid(axis='y', linestyle=':', alpha=0.7)
    
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, height + 0.02, f'{height:.2f}x', ha='center', va='bottom', fontweight='bold')
    
    plt.legend()
    plt.tight_layout()
    plt.show()


KS = [32, 64, 80, 128, 150, 256]
DATASET = "ca-hepph"
PROB = "0.2"
STRATEGIA = "mixed" # or "standard"
SEQ = "rich"       # erdos, power, rich
ALPHA = 0.12
PHI = 0.6


plot_speedup(DATASET, PROB, KS, strategy=STRATEGIA, seq_type=SEQ, alpha=ALPHA, phi=PHI)
plot_quasi_clique_results_grid(DATASET, PROB, KS, strategy=STRATEGIA, seq_type=SEQ, alpha=ALPHA, phi=PHI)