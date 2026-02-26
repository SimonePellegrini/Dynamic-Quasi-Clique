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

# --- ESECUZIONE ---
# Assicurati che il file si trovi nella stessa cartella dello script
file_input = "./Datasets/Amazon0601.txt" 
file_output = "./Datasets/Amazon0601c.txt"

trasforma_grafo_con_header(file_input, file_output)