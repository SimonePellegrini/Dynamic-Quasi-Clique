import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
from itertools import product

def format_cpp_float(val):
    """Simula std::to_string(float) del C++ (6 decimali)."""
    return f"{float(val):.6f}"

def get_file_path(base_path, strategy, p_str, k, seq_type, alpha, phi, file_type="results"):
    a_str = format_cpp_float(alpha)
    p_val = format_cpp_float(phi)
    
    if strategy == "standard":
        filename = f"{file_type}_credits_{p_str}_alpha:{a_str}_phi:{p_val}_k:{k}.csv"
    else: # mixed
        filename = f"{file_type}_credits_{seq_type}_{p_str}_alpha:{a_str}_phi:{p_val}_k:{k}.csv"
            
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
    # Se la strategia è 'mixed', usiamo seq_type come cartella, altrimenti 'standard'
    target_folder = "standard" if strategy == "standard" else seq_type
    save_dir = os.path.join(dataset_plot_root, target_folder)
    
    # Caricamento Baseline
    suffix = f"{p_str}" if strategy == "standard" else f"{seq_type}_{p_str}"
    file_dyn = f"{base_data_dir}results_dynamic_{suffix}.csv"
    file_dyn_dens = f"{base_data_dir}density_dynamic_{suffix}.csv"
    file_dyn_time = f"{base_data_dir}times_dynamic_{suffix}.csv"

    df_dyn = pd.read_csv(file_dyn, header=None).mean(axis=1) if os.path.exists(file_dyn) else None
    df_dyn_dens = pd.read_csv(file_dyn_dens, header=None).mean(axis=1) if os.path.exists(file_dyn_dens) else None
    baseline_time = pd.read_csv(file_dyn_time, header=None).mean().mean() if os.path.exists(file_dyn_time) else None

    for alpha, phi, k in product(alpha_list, phi_list, k_list):
        f_cred = get_file_path(base_data_dir, strategy, p_str, k, seq_type, alpha, phi, "results")
        f_dens = get_file_path(base_data_dir, strategy, p_str, k, seq_type, alpha, phi, "density")
        f_time = get_file_path(base_data_dir, strategy, p_str, k, seq_type, alpha, phi, "times")

        if not os.path.exists(f_cred):
            continue

        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot Baseline
        if df_dyn is not None:
            ax.plot(df_dyn, color='black', linestyle='--', alpha=0.3, label='Baseline Dynamic')
        
        # Plot Credits
        df_c = pd.read_csv(f_cred, header=None).mean(axis=1)
        ax.plot(df_c, color='steelblue', linewidth=2, label=f'Credits Algorithm')
        
        # Overlay Densità
        if os.path.exists(f_dens):
            df_d = pd.read_csv(f_dens, header=None).mean(axis=1)
            cmap = plt.get_cmap('RdYlGn')
            norm = mcolors.BoundaryNorm([0.0, 0.5, 0.7, 0.8, 0.9, 1.0], cmap.N)
            idx_c = np.linspace(0, len(df_c)-1, len(df_d), dtype=int)
            sc = ax.scatter(idx_c, df_c.iloc[idx_c], c=df_d, cmap=cmap, norm=norm, 
                            marker='^', s=60, edgecolors='black', label='Density Marker', zorder=5)
            plt.colorbar(sc, ax=ax, label='Density')

        # Speedup nel titolo
        title_str = f"{dataset_name} | {target_folder.upper()}\nα={alpha}, φ={phi}, k={k}"
        if baseline_time and os.path.exists(f_time):
            c_time = pd.read_csv(f_time, header=None).mean().mean()
            speedup = baseline_time / c_time
            title_str += f" | Speedup: {speedup:.2f}x"

        ax.set_title(title_str, fontweight='bold')
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Quasi-Clique Size")
        ax.legend(loc='lower right')
        ax.grid(True, alpha=0.2)

        # Salvataggio nella sottocartella specifica
        filename = f"a{alpha}_p{phi}_k{k}.png"
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()

# --- Main ---
if __name__ == "__main__":
    ALPHAS = [0.12, 0.3, 0.5]
    PHIS = [0.4, 0.6, 0.8]
    KS = [32, 64, 128,256,512]
    
    DATASET = "link-dynamic-plwiki"
    P_STR = "0.0"
    
    # 1. Esegui Standard
    generate_all_plots(DATASET, P_STR, KS, ALPHAS, PHIS, strategy="standard")
    
    # 2. Esegui i Mixed
    for mixed_type in []:
        generate_all_plots(DATASET, P_STR, KS, ALPHAS, PHIS, strategy="mixed", seq_type=mixed_type)

    print("\nProcesso completato! Trovi i grafici ordinati in 'plots/{}'".format(DATASET))