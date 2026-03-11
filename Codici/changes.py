import pandas as pd
import os
import sys
from itertools import product

def format_cpp_float(val):
    """Simula std::to_string(float) del C++ (6 decimali)."""
    return f"{float(val):.6f}"

def get_times_dynamic_path(base_dir, strategy, p_str, seq_type):
    """Restituisce il path del file dei tempi per la baseline dinamica."""
    folder = "standard" if strategy == "standard" else "mixed"
    if strategy == "standard":
        filename = f"times_dynamic_{p_str}.csv"
    else:
        filename = f"times_dynamic_{seq_type}_{p_str}.csv"
    return os.path.join(base_dir, folder, filename)

def get_times_credits_path(base_dir, strategy, p_str, seq_type, alpha, phi, k):
    """Restituisce il path del file dei tempi per l'algoritmo credits."""
    folder = "standard" if strategy == "standard" else "mixed"
    a_str = format_cpp_float(alpha)
    p_val = format_cpp_float(phi)
    
    if strategy == "standard":
        filename = f"times_credits_{p_str}_alpha:{a_str}_phi:{p_val}_k:{k}.csv"
    else:
        filename = f"times_credits_{seq_type}_{p_str}_alpha:{a_str}_phi:{p_val}_k:{k}.csv"
        
    return os.path.join(base_dir, folder, filename)

def create_speedup_folder_structure(dataset_name):
    """Crea la gerarchia di cartelle per salvare i risultati dello speedup."""
    base_out_dir = os.path.join("speedup", dataset_name)
    subfolders = ["standard", "Erdos", "rich", "power"]
    
    for sub in subfolders:
        path = os.path.join(base_out_dir, sub)
        os.makedirs(path, exist_ok=True) 
        
    return base_out_dir

def compute_and_save_speedup(dataset_name, p_str, k_list, alpha_list, phi_list):
    """Calcola lo speedup e salva solo i parametri e i tempi."""
    base_data_dir = f"../Esperimenti/{dataset_name}/"
    base_out_dir = create_speedup_folder_structure(dataset_name)
    
    total_results = [] 

    strategies = [
        ("standard", "erdos"),
        ("mixed", "erdos"),
        ("mixed", "rich"),
        ("mixed", "power")
    ]

    for strategy, seq_type in strategies:
        category_results = []
        
        file_dyn = get_times_dynamic_path(base_data_dir, strategy, p_str, seq_type)
        
        if not os.path.exists(file_dyn):
            continue
            
        df_dyn = pd.read_csv(file_dyn, header=None)
        mean_time_dyn = df_dyn.mean(axis=1).values[0]

        for alpha, phi, k in product(alpha_list, phi_list, k_list):
            file_cred = get_times_credits_path(base_data_dir, strategy, p_str, seq_type, alpha, phi, k)
            
            if not os.path.exists(file_cred):
                continue
                
            df_cred = pd.read_csv(file_cred, header=None)
            mean_time_cred = df_cred.mean(axis=1).values[0]
            
            speedup = (mean_time_dyn / mean_time_cred) if mean_time_cred > 0 else float('inf')
            
            # --- Modifica: Salviamo SOLO le info richieste ---
            row = {
                "Alpha": alpha,
                "Phi": phi,
                "K": k,
                "Time_Dynamic_ms": round(mean_time_dyn, 2),
                "Time_Credits_ms": round(mean_time_cred, 2),
                "Speedup": round(speedup, 4)
            }
            category_results.append(row)
            total_results.append(row)

        if category_results:
            if strategy == "standard":
                target_folder = "standard"
            elif seq_type == "erdos":
                target_folder = "Erdos"
            else:
                target_folder = seq_type

            output_file = os.path.join(base_out_dir, target_folder, f"speedup_summary_prob{p_str}.csv")
            df_cat = pd.DataFrame(category_results)
            df_cat.to_csv(output_file, index=False)
            print(f"Salvato: {output_file}")

    if total_results:
        df_total = pd.DataFrame(total_results)
        print("\nAnteprima risultati totali:")
        print(df_total.head(10).to_string())
    else:
        print("\nNessun dato trovato per calcolare lo speedup. Verifica i percorsi dei file generati dal C++.")

# --- Main ---
if __name__ == "__main__":
    ALPHAS = [0.12, 0.3, 0.5]
    PHIS = [0.4, 0.6, 0.8]
    KS = [32, 64, 128, 256, 512]
    
    if len(sys.argv) > 2:
        DATASET = sys.argv[1]
        probability = float(sys.argv[2])
    else:
        DATASET = "link-dynamic-simplewiki"
        prob = ["0.0","0.05","0.1"]
    for probability in prob: 
        P_STR = prob

        print(f"Avvio calcolo speedup per {DATASET} (prob: {probability})...")
        compute_and_save_speedup(DATASET, P_STR, KS, ALPHAS, PHIS)