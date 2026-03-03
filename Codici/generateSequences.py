import os 
import random

def generate_sequences(file_name, number_of_sequences = 1, prob = 0.1):
    #open the file
    with open('../Datasets/'+file_name+".txt",'r') as file:
        lines = file.readlines()
    
    file_obj = open("./Sequences/StandardSequences/"+file_name+"_"+str(prob)+".txt", mode="w")

    n = int(lines[0].split(" ")[0])
    m = int(lines[0].split(" ")[1])

    file_obj.write(str(number_of_sequences)+"\n")
    file_obj.write(str(n)+ " "+str(m)+"\n")
    
    for i in range(number_of_sequences):
        insert_edges = [] 
        graph = lines[1:]
        random.shuffle(graph)
        #insert the firt half of the edges
        for i in range(1,m//2):
            u = graph[i].split(" ")[0]
            v = graph[i].split(" ")[1]
            file_obj.write("i "+str(u) +" "+ str(v))
            insert_edges.append((u,v))

        i = m//2+1
        while i<m :
            v = random.random()
            if v<prob:
                idx = random.randrange(len(insert_edges))
                edge_2_delete = insert_edges[idx]
                u = edge_2_delete[0]
                v = edge_2_delete[1]
                insert_edges[-1],insert_edges[idx] = insert_edges[idx],insert_edges[-1]
                insert_edges.pop()
                file_obj.write("d "+str(u) +" "+ str(v))
                
            else:
                u = graph[i].split(" ")[0]
                v = graph[i].split(" ")[1]
                file_obj.write("i "+str(u) +" "+ str(v))
                insert_edges.append((u,v))
                i+=1
        file_obj.write("e 0 0\n")
        
import random

def generate_mixed_sequences(file_name, number_of_sequences=1, prob=0.1, additional_operations=0.2, strategy="erdos"):
    with open('../Datasets/'+file_name+".txt", 'r') as file:
        lines = file.readlines()
    
    file_obj = open("./Sequences/MixedSequences/"+file_name+"_"+str(prob)+"_"+strategy+".txt", mode="w")

    n = int(lines[0].split(" ")[0])
    m = int(lines[0].split(" ")[1])

    file_obj.write(str(number_of_sequences)+"\n")
    file_obj.write(str(n)+ " "+str(m)+"\n")
    
    for seq_idx in range(number_of_sequences):
        insert_edges = []
        existing_edges = set() 
        
        graph = lines[1:]
        random.shuffle(graph)
        
        for edge_idx in range(1, m):
            parts = graph[edge_idx].strip().split(" ")
            
            u = parts[0]
            v = parts[1]
            file_obj.write(f"i {u} {v}\n")
            
            insert_edges.append((u, v))
            existing_edges.add((u, v))
            existing_edges.add((v, u))

        for op_idx in range(int(additional_operations * m)):
            v_rand = random.random()
            
            if v_rand < prob:         
                idx = random.randrange(len(insert_edges))
                edge_2_delete = insert_edges[idx]
                u = edge_2_delete[0]
                v = edge_2_delete[1]
                
                insert_edges[-1], insert_edges[idx] = insert_edges[idx], insert_edges[-1]
                insert_edges.pop()
 
                existing_edges.discard((u, v))
                existing_edges.discard((v, u))
                
                file_obj.write(f"d {u} {v}\n")
                
            else:
                
                if strategy=="erdos":
                    while True:
                        u_new = str(random.randint(0, n-1))
                        v_new = str(random.randint(0, n-1))
                        if u_new != v_new and (u_new, v_new) not in existing_edges:
                            break
                            
                elif strategy=="rich": 
                    while True:
                        rand_edge = random.choice(insert_edges)
                        u_new = rand_edge[random.randint(0, 1)]
                        v_new = str(random.randint(0, n-1))
                        
                        if u_new != v_new and (u_new, v_new) not in existing_edges:
                            break
                            
                elif strategy=="power":
                    while True:
                        rand_edge_1 = random.choice(insert_edges)
                        rand_edge_2 = random.choice(insert_edges)
                        u_new = rand_edge_1[random.randint(0, 1)]
                        v_new = rand_edge_2[random.randint(0, 1)]
                        
                        if u_new != v_new and (u_new, v_new) not in existing_edges:
                            break
                file_obj.write(f"i {u_new} {v_new}\n")
                insert_edges.append((u_new, v_new))
                existing_edges.add((u_new, v_new))
                existing_edges.add((v_new, u_new))
                
        file_obj.write("e 0 0\n")
        
    file_obj.close() 


datasets=["ego-facebook","email-enron","loc-gowalla","web-stanford","ca-hepph","ca-condmat"]
probabilities = [0.1,0.05]
types = ["erdos"]

for prob in probabilities:
    for name in datasets:
        for type in types:
            generate_mixed_sequences(name,5,prob,0.5,type)

for prob in probabilities:
    for name in datasets:
            generate_sequences(name,5,prob)


'''
import sys
from collections import defaultdict
import sys
from collections import defaultdict

def process_dynamic_graph(input_path, output_path):
    active_edges = set()   # archi attivi finali (non orientati)
    nodes = set()          # insieme dei nodi apparsi in operazioni valide
    output_lines = []      # memorizziamo le operazioni dinamiche valide
    num_op = 0
    with open(input_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) < 3:
                continue

            try:
                u = int(parts[0])
                v = int(parts[1])
                op = int(parts[2])
            except ValueError:
                continue  # salta header

            # rimuove self-loop
            if u == v:
                continue

            # rende non orientato
            a, b = min(u, v), max(u, v)
            edge = (a, b)

            if op == 1:
                if edge not in active_edges:
                    num_op+=1
                    active_edges.add(edge)
                    output_lines.append(f"i {a} {b}")
                    nodes.add(a)
                    nodes.add(b)

            elif op == -1:
                if edge in active_edges:
                    num_op+=1
                    active_edges.remove(edge)
                    output_lines.append(f"d {a} {b}")
                    nodes.add(a)
                    nodes.add(b)

    # Scrittura file finale
    with open(output_path, "w") as out:
        # prima riga: numero nodi e numero archi finali
        out.write(f"{len(nodes)} {num_op}\n")

        # poi tutte le operazioni dinamiche valide
        for line in output_lines:
            out.write(line + "\n")


if __name__ == "__main__":
    process_dynamic_graph("../Datasets/link-dynamic-plwiki.txt","./Sequences/StandardSequences/link-dynamic-plwiki_0.0.txt")
'''