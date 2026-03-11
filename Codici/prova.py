import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import sys
from itertools import product

def format_cpp_float(val):
    """Simula std::to_string(float) del C++ (6 decimali)."""
    return f"{float(val):.6f}"

def get_file_path(base_path, strategy, p_str, k, seq_type, alpha, phi, file_type="results"):
    a_str = format_cpp_float(alpha)
    p_val = format_cpp_float(phi)
    
    if strategy == "standard":
        filename = f"{file_type}_credits_{p_str}_alpha{a_str}_phi{p_val}_k{k}.csv"
    else: # mixed
        filename = f"{file_type}_credits_{seq_type}_{p_str}_alpha{a_str}_phi{p_val}_k{k}.csv"
            
    return os.path.join(base_path, filename)

def create_folder_structure(dataset_name):
    """Crea la gerarchia di cartelle richiesta."""
    base_plot_path = os.path.join("plots", dataset_name)
    subfolders = ["standard", "erdos", "power", "rich"]
    
    for sub in subfolders:
        path = os.path.join(base_plot_path, sub)
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Cartella creata: {path}")
    
    return base_plot_path

def generate_all_plots(dataset_name, p_str, k_list, alpha_list, phi_list, strategy="standard", seq_type="erdos"):
    # Path dei dati C++
    base_data_dir = f"../Esperimenti/{dataset_name}/{strategy}/"
    
    # Path di destinazione dei Plot
    dataset_plot_root = create_folder_structure(dataset_name)
    target_folder = "standard" if strategy == "standard" else seq_type
    save_dir = os.path.join(dataset_plot_root, target_folder)
    
    # Caricamento Baseline
    suffix = f"{p_str}" if strategy == "standard" else f"{seq_type}_{p_str}"
    file_dyn = f"{base_data_dir}results_dynamic_{suffix}.csv"
    file_dyn_dens = f"{base_data_dir}density_dynamic_{suffix}.csv"

    df_dyn = pd.read_csv(file_dyn, header=None).mean(axis=1) if os.path.exists(file_dyn) else None
    df_dyn_dens = pd.read_csv(file_dyn_dens, header=None).mean(axis=1) if os.path.exists(file_dyn_dens) else None

    # Configurazione Colormap per le densità (usata per entrambi i grafici)
    cmap = plt.get_cmap('RdYlGn')
    norm = mcolors.BoundaryNorm([0.0, 0.5, 0.7, 0.8, 0.9, 1.0], cmap.N)

    for alpha, phi, k in product(alpha_list, phi_list, k_list):
        f_cred = get_file_path(base_data_dir, strategy, p_str, k, seq_type, alpha, phi, "results")
        f_dens = get_file_path(base_data_dir, strategy, p_str, k, seq_type, alpha, phi, "density")

        if not os.path.exists(f_cred):
            continue

        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 1. Plot Baseline (Grigio)
        if df_dyn is not None:
            ax.plot(df_dyn, color='gray', linestyle='--', alpha=0.7, label='Baseline Dynamic')
            
            # Overlay Densità per Baseline
            if df_dyn_dens is not None:
                idx_dyn = np.linspace(0, len(df_dyn)-1, len(df_dyn_dens), dtype=int)
                ax.scatter(idx_dyn, df_dyn.iloc[idx_dyn], c=df_dyn_dens, cmap=cmap, norm=norm, 
                           marker='^', s=60, edgecolors='black', zorder=5)
        
        # 2. Plot Credits (Verde)
        df_c = pd.read_csv(f_cred, header=None).mean(axis=1)
        ax.plot(df_c, color='steelblue', linewidth=2, label='Credits Algorithm')
        
        # Overlay Densità per Credits
        if os.path.exists(f_dens):
            df_d = pd.read_csv(f_dens, header=None).mean(axis=1)
            idx_c = np.linspace(0, len(df_c)-1, len(df_d), dtype=int)
            sc = ax.scatter(idx_c, df_c.iloc[idx_c], c=df_d, cmap=cmap, norm=norm, 
                            marker='^', s=60, edgecolors='black', zorder=5)
            
            # Aggiungiamo la colorbar una sola volta se abbiamo i dati della densità
            plt.colorbar(sc, ax=ax, label='Density')

        # --- Modifiche applicate per Testo e Assi ---
        ax.set_title(dataset_name, fontweight='bold') # Solo il nome del dataset
        ax.set_xlabel("Number of operations")         # Nuovo asse X
        ax.set_ylabel("Quasi-clique size")            # Nuovo asse Y
        
        ax.legend(loc='lower right')
        ax.grid(True, alpha=0.2)

        # Salvataggio nella sottocartella specifica
        filename = f"a{alpha}_p{phi}_k{k}_prob{p_str}.pdf"
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

# --- Main ---
if __name__ == "__main__":
    ALPHAS = [0.12, 0.3, 0.5]
    PHIS = [0.4, 0.6, 0.8]
    KS = [32, 64, 128, 256, 512]
    
    DATASET = "com-lj.ungraph_clean"
    P_STR = "10"
    
    generate_all_plots(DATASET, P_STR, KS, ALPHAS, PHIS, strategy="standard")
    for mixed_type in ["erdos","rich"]:
        generate_all_plots(DATASET, P_STR, KS, ALPHAS, PHIS, strategy="mixed", seq_type=mixed_type)

    print("\nProcesso completato! Trovi i grafici ordinati in 'plots/{}'".format(DATASET))