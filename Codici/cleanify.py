import networkx as nx

def trasforma_grafo_con_header(file_input, file_output):
    # Usiamo direttamente nx.Graph() per gestire la non-orientazione
    # NetworkX ignorerà automaticamente gli archi inversi (es. se esiste 0-1, ignora 1-0)
    G = nx.Graph()
    
    print(f"Lettura di {file_input}...")
    
    try:
        with open(file_input, 'r') as f:
            # Salta la prima riga del file originale (l'header vecchio)
            next(f)
            
            for linea in f:
                parti = linea.split()
                if len(parti) == 2:
                    # Aggiunge l'arco (gestisce automaticamente i duplicati)
                    G.add_edge(parti[0], parti[1])
        
        # Calcoliamo i nuovi totali dopo la pulizia
        num_nodi = G.number_of_nodes()
        num_archi = G.number_of_edges()
        
        print(f"Elaborazione completata.")
        print(f"Nuovi totali -> Nodi: {num_nodi}, Archi: {num_archi}")

        # Scrittura del file di output
        with open(file_output, 'w') as f_out:
            # Scrive la prima riga con i nuovi metadati
            f_out.write(f"{num_nodi} {num_archi}\n")
            
            # Scrive tutti gli archi unici
            for u, v in G.edges():
                f_out.write(f"{u} {v}\n")
                
        print(f"File salvato correttamente in: {file_output}")

    except FileNotFoundError:
        print("Errore: Il file di input non è stato trovato.")
    except Exception as e:
        print(f"Si è verificato un errore: {e}")

def pulisci_dataset(file_input, file_output):
    try:
        with open(file_input, 'r') as f_in, open(file_output, 'w') as f_out:
            for riga in f_in:
                # Rimuove spazi bianchi all'inizio/fine e divide la riga in colonne
                colonne = riga.split()
                
                # Se la riga non è vuota, prendiamo solo le prime due colonne
                if colonne:
                    nuova_riga = f"{colonne[0]} {colonne[1]}\n"
                    f_out.write(nuova_riga)
        
        print(f"Successo! Il nuovo file è stato salvato come: {file_output}")
    
    except FileNotFoundError:
        print("Errore: Il file di input non è stato trovato.")
    except Exception as e:
        print(f"Si è verificato un errore: {e}")

import os

def clean_graph_dataset(input_path, output_path):
    if not os.path.exists(input_path):
        print(f"Errore: Il file '{input_path}' non esiste.")
        return

    seen_edges = set()
    unique_nodes = set()
    ordered_edges = []

    try:
        with open(input_path, 'r') as f:
            for line in f:
                # Pulizia riga e skip se vuota
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                
                # Conversione in int per normalizzazione
                u, v = int(parts[0]), int(parts[1])
                
                # Creiamo una chiave unica non ordinata (es. {2, 28})
                # Usiamo frozenset o una tupla ordinata per il set dei duplicati
                edge_key = tuple(sorted((u, v)))
                
                if edge_key not in seen_edges:
                    seen_edges.add(edge_key)
                    unique_nodes.add(u)
                    unique_nodes.add(v)
                    # Manteniamo il formato originale della riga
                    ordered_edges.append(f"{u} {v}")

        # Scrittura sul file di output
        num_nodes = len(unique_nodes)
        num_edges = len(ordered_edges)

        with open(output_path, 'w') as f_out:
            # Scriviamo la riga di intestazione con i conteggi
            f_out.write(f"{num_nodes} {num_edges}\n")
            # Scriviamo gli archi preservando l'ordine
            f_out.write("\n".join(ordered_edges))

        print(f"Operazione completata con successo!")
        print(f"Nodi unici: {num_nodes} | Archi unici: {num_edges}")
        print(f"File salvato in: {output_path}")

    except Exception as e:
        print(f"Si è verificato un errore durante l'elaborazione: {e}")


def clean_graph_dataset(input_path, output_path):
    edges = set()
    nodes = set()
    
    try:
        with open(input_path, 'r') as f:
            for line in f:
                # Salta righe vuote o commenti
                if not line.strip() or line.startswith('#'):
                    continue
                
                try:
                    u, v = map(int, line.split())
                    
                    # 1. Rimuove i self-loop
                    if u == v:
                        continue
                    
                    # 2. Aggiunge l'arco (usiamo una tupla ordinata per grafi non orientati)
                    # Se il tuo grafo è ORIENTATO, usa semplicemente edges.add((u, v))
                    edge = tuple(sorted((u, v))) 
                    edges.add(edge)
                    
                    # 3. Registra i nodi univoci
                    nodes.add(u)
                    nodes.add(v)
                    
                except ValueError:
                    continue

        # Scrittura del file ripulito
        with open(output_path, 'w') as f:
            # La prima riga: numero_nodi numero_archi
            f.write(f"{len(nodes)} {len(edges)}\n")
            
            # Scrittura degli archi
            for u, v in sorted(edges):
                f.write(f"{u} {v}\n")
                
        print(f"Pulizia completata!")
        print(f"Nodi: {len(nodes)}, Archi: {len(edges)}")
        print(f"File salvato in: {output_path}")

    except FileNotFoundError:
        print("Errore: Il file di input non è stato trovato.")


clean_graph_dataset('../Datasets/linux.txt', '../Datasets/linux_clean.txt')

