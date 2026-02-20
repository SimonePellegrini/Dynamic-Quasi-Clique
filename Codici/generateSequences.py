import os 
import random

def generate_sequences(file_name, number_of_sequences = 1, prob = 0.1,additional_operations = 0.2):
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

        i = 0
        edges = m//2+1
        while i<m*additional_operations :
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
                u = graph[edges].split(" ")[0]
                v = graph[edges].split(" ")[1]
                file_obj.write("i "+str(u) +" "+ str(v))
                insert_edges.append((u,v))
                edges+=1
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


datasets=["ego-facebook","ca-CondMat","ca-HepPh","loc-gowalla","email-enron","web-stanford"]
probabilities = [0.0,0.1,0.2]
types = ["erdos","rich","power"]

for prob in probabilities:
    for name in datasets:
        for type in types:
            generate_mixed_sequences(name,10,prob,0.4,type)

for prob in probabilities:
    for name in datasets:
            generate_sequences(name,10,prob,0.4)
