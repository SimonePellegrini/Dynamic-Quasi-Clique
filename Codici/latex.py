import pandas as pd
import numpy as np
import os

def load_and_mean(filepath):
    """Carica il CSV e restituisce la media su tutte le colonne (esperimenti)."""
    if not os.path.exists(filepath):
        return None
    # Il file CSV ha una riga per ogni step temporale e una colonna per ogni esperimento
    df = pd.read_csv(filepath, header=None)
    return df.mean(axis=1).values

def get_table_rows(dataset, prob, alpha, phi, k, strategy, net_type):
    # Converte prob in formato stringa come fa il C++ (es. 0.05 -> "5")
    p_str = prob
    base_path = f"../Esperimenti/{dataset}/{strategy}"
    
    # Costruisce i percorsi dei file in base alla strategia (standard o mixed)
    if strategy == "standard":
        # Credits
        res_cred_path = f"{base_path}/results_credits_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        den_cred_path = f"{base_path}/density_credits_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        time_cred_path = f"{base_path}/times_credits_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        # Dynamic (Baseline)
        res_dyn_path = f"{base_path}/results_dynamic_{p_str}.csv"
        den_dyn_path = f"{base_path}/density_dynamic_{p_str}.csv"
        time_dyn_path = f"{base_path}/times_dynamic_{p_str}.csv"
    else:
        # Credits
        res_cred_path = f"{base_path}/results_credits_{net_type}_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        den_cred_path = f"{base_path}/density_credits_{net_type}_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        time_cred_path = f"{base_path}/times_credits_{net_type}_{p_str}_alpha{alpha}_phi{phi}_k{k}.csv"
        # Dynamic (Baseline)
        res_dyn_path = f"{base_path}/results_dynamic_{net_type}_{p_str}.csv"
        den_dyn_path = f"{base_path}/density_dynamic_{net_type}_{p_str}.csv"
        time_dyn_path = f"{base_path}/times_dynamic_{net_type}_{p_str}.csv"

    # --- 1. Calcolo Size e Densità per CREDITS ---
    sizes_cred = load_and_mean(res_cred_path)
    if sizes_cred is not None:
        total_ops = len(sizes_cred)
        size_c_60 = sizes_cred[int(0.60 * total_ops) - 1]
        size_c_80 = sizes_cred[int(0.80 * total_ops) - 1]
        size_c_100 = sizes_cred[total_ops - 1]
    else:
        size_c_60 = size_c_80 = size_c_100 = "?"

    densities_cred = load_and_mean(den_cred_path)
    if densities_cred is not None and len(densities_cred) >= 9:
        den_c_60 = densities_cred[-5]
        den_c_80 = densities_cred[-3]
        den_c_100 = densities_cred[-1]
    else:
        den_c_60 = den_c_80 = den_c_100 = "?"

    # --- 2. Calcolo Size e Densità per DYNAMIC (Baseline) ---
    sizes_dyn = load_and_mean(res_dyn_path)
    if sizes_dyn is not None:
        total_ops_dyn = len(sizes_dyn)
        size_d_60 = sizes_dyn[int(0.60 * total_ops_dyn) - 1]
        size_d_80 = sizes_dyn[int(0.80 * total_ops_dyn) - 1]
        size_d_100 = sizes_dyn[total_ops_dyn - 1]
    else:
        size_d_60 = size_d_80 = size_d_100 = "?"

    densities_dyn = load_and_mean(den_dyn_path)
    if densities_dyn is not None and len(densities_dyn) >= 9:
        den_d_60 = densities_dyn[-5]
        den_d_80 = densities_dyn[-3]
        den_d_100 = densities_dyn[-1]
    else:
        den_d_60 = den_d_80 = den_d_100 = "?"

    # --- 3. Calcolo Speedup ---
    time_cred = load_and_mean(time_cred_path)
    time_dyn = load_and_mean(time_dyn_path)
    
    if time_cred is not None and time_dyn is not None:
        avg_time_cred = time_cred[0]
        avg_time_dyn = time_dyn[0]
        speedup = avg_time_dyn / avg_time_cred
        speedup_str = f"{speedup:.2f}\\times"
    else:
        speedup_str = "?\\times"

    # --- Formattazione dell'output LaTeX ---
    def f_size(val):
        return f"$\\approx {int(round(val))}$" if isinstance(val, (int, float, np.float64)) else "$\\approx ?$"
        
    def f_den(val):
        return f"{val:.2f}" if isinstance(val, (int, float, np.float64)) else "?"

    # Riga per la prima tabella (Credits + Speedup)
    row_credits = (f"{dataset:<15} & {f_size(size_c_60):<12} & {f_den(den_c_60):<10} "
                   f"& {f_size(size_c_80):<12} & {f_den(den_c_80):<10} "
                   f"& {f_size(size_c_100):<12} & {f_den(den_c_100):<10} & ${speedup_str}$ \\\\")
                   
    # Riga per la seconda tabella (Baseline Dinamica)
    row_dynamic = (f"{dataset:<15} & {f_size(size_d_60):<12} & {f_den(den_d_60):<10} "
                   f"& {f_size(size_d_80):<12} & {f_den(den_d_80):<10} "
                   f"& {f_size(size_d_100):<12} & {f_den(den_d_100):<10} \\\\")
    
    return row_credits, row_dynamic


if __name__ == "__main__":
    datasets = [
        "linux",
        "ca-cit-hepph",
        "ego-facebook",
        "email-Enron",
        "web-stanford",
        "loc-gowalla",
        "ca-hepph",
        "ca-condmat"
    ]
    
    # PARAMETRI DEGLI ESPERIMENTI
    P_PROB = "0"        
    ALPHA = "0.300000"       
    PHI = "0.800000"         
    K = 64                  
    STRATEGY = "standard"    
    NET_TYPE = "erdos"        

    credits_rows = []
    dynamic_rows = []

    # Calcola le righe per ogni dataset
    for ds in datasets:
        r_cred, r_dyn = get_table_rows(dataset=ds, 
                                       prob=P_PROB, 
                                       alpha=ALPHA, 
                                       phi=PHI, 
                                       k=K, 
                                       strategy=STRATEGY, 
                                       net_type=NET_TYPE)
        credits_rows.append(r_cred)
        dynamic_rows.append(r_dyn)

    # --- Stampa dell'intera tabella LaTeX ---
    print(r"""\begin{table}[h]
\centering

\begin{tabular}{@{}lccccccc@{}}
\toprule
\textbf{Dataset} & \multicolumn{2}{c}{\textbf{60\%}} & \multicolumn{2}{c}{\textbf{80\%}} & \multicolumn{2}{c}{\textbf{100\%}} & \textbf{Speedup} \\ 
\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}
                 & \textbf{size} & \textbf{density} & \textbf{size} & \textbf{density} & \textbf{size} & \textbf{density} & \\ \midrule""")

    for row in credits_rows:
        print(row)

    print(r"""\bottomrule \\
\end{tabular}
\label{tab:incremental_performance_credits}

\vspace{0.5cm} % Uno spazio opzionale tra le due tabelle per estetica

\begin{tabular}{@{}lcccccc@{}}
\toprule
\textbf{Dataset} & \multicolumn{2}{c}{\textbf{60\%}} & \multicolumn{2}{c}{\textbf{80\%}} & \multicolumn{2}{c}{\textbf{100\%}} \\ 
\cmidrule(lr){2-3} \cmidrule(lr){4-5} \cmidrule(lr){6-7}
                 & \textbf{size} & \textbf{density} & \textbf{size} & \textbf{density} & \textbf{size} & \textbf{density} \\ \midrule""")

    for row in dynamic_rows:
        print(row)

    print(r"""\bottomrule
\end{tabular}
\label{tab:incremental_performance_baseline}
\caption{Performance Comparison of the Credits-Based algorithm (first table) and the Baseline (second table) across different Datasets (using $\alpha=0.3, \phi=0.8, k=64$)}
\end{table}""")