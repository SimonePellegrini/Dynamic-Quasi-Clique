#ifndef CREDITS_ALGO
#define CREDITS_ALGO
#include <unordered_set>
#include <iostream>
#include <vector>
#include <queue>
#include <random>
#include <ctime>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <set>


using namespace std;

class CreditsAlgo{
    private: 
    int rounds = 0;
    int eta,delta,k,n,m,phi;
    double gamma,beta;
    vector<vector<int>> neighborhood;
    vector<vector<int>> signature;
    vector<int> a;
    vector<int> b;
    vector<int> credits;
    vector<set<int>> gamma_deg;
    vector<int> cur_sol;
    set<pair<int,int>> pq;
    vector<int> degree;
    int best_clique=0;
    int best_clique_size=0; 
    

    public:
        CreditsAlgo(int n,int m, int eta, int phi,int delta, double gamma,double beta,int k);
        vector<int> returnBestClique();
        double ct_score_minhash(int u, int v);
        void add_edge(int u, int v);
        int sift_num(int u);
        void update_signature(int u,int v); 
        int return_dim();
        double return_density();
};
int p = pow(2,31)-1;

CreditsAlgo::CreditsAlgo(int n,int m, int eta, int phi,int delta, double gamma,double beta,int k){
    /*
    input:
    n: il numero di nodi nel grafo
    m: il numero di archi nel grafo
    eta: un parametro che definisce il numero di rounds oltre i quali aggiornare la soluzione
    phi: un parametro che definisce il numero di nodi del vicinato a cui assegnare crediti
    delta: un parametro che definisce quanti elementi estrarre dalla coda con priorità
    gamma,beta,k: come l'algoritmo NBSim
    */

    this->n = n;
    this->m = m;
    this->eta = eta;
    this->phi = phi;
    this->delta = delta;
    this->gamma = gamma;
    this->beta = beta;
    this->k = k;

    //inizializzo i vari vettori
    neighborhood.resize(n);
    credits.resize(n);
    a.resize(k);
    b.resize(k);
    gamma_deg.resize(n);
    degree.resize(n);
    signature.resize(n);

    //inizializzo le k funzioni hash indipendenti
    for(int i = 0; i<k;i++){
        a[i] = rand()%(p-1)+1;
        b[i] = rand();
    }

    //inizializzo le firme min-hash
    for(int i=0; i<n; i++){
        signature[i].resize(k);
        for(int j=0; j<k; j++){
            signature[i][j] = int(((long long)(a[j]) * (long long)(i) + (long long)(b[j])) % p);
        }
    }

    //inizializzo la coda con priorità
    for(int i=0; i<n; i++){
        pq.insert({0,i});
    }
}

void CreditsAlgo::update_signature(int u,int v){
    //aggiorno la signature di N(u) 
    for(int i = 0; i<k; i++){
        signature[u][i] = min(signature[u][i],int(((long long)(a[i]) * (long long)(v) + (long long)(b[i])) % p));
    }
    //aggiorno la signature di N(v) 
    for(int i = 0; i<k; i++){
        signature[v][i] = min(signature[v][i],int(((long long)(a[i]) * (long long)(u) + (long long)(b[i])) % p));
    }
}

void CreditsAlgo::add_edge(int u, int v){
    //aggiungo l'arco al grafo
    neighborhood[v].emplace_back(u);
    neighborhood[u].emplace_back(v);
    degree[v]++;
    degree[u]++;
    
    
    //aggiorno la signature di N(u) e N(v), costo O(k)
    update_signature(u,v);

    //assegno i crediti a phi nodi
    int ind=0;
    for(int i = 0; i<min(phi,degree[u]); i++){     
        ind = neighborhood[u][rand()%degree[u]];
        if((degree[u]+1) >= gamma*(degree[ind]+1) && !gamma_deg[ind].count(u)){
            pq.erase({credits[ind],ind});
            credits[ind]+=1;
            pq.insert({credits[ind],ind});
            gamma_deg[ind].insert(u);
        }
        if((degree[ind]+1) >= gamma*(degree[u]+1) && !gamma_deg[u].count(ind)){
            pq.erase({credits[u],u});
            credits[u]+=1;
            pq.insert({credits[u],u});
            gamma_deg[u].insert(ind);
        }
    }
    for(int i = 0; i<min(phi,degree[v]); i++){ 
        ind = neighborhood[v][rand()%degree[v]];
        if((degree[v]+1) >= gamma*(degree[ind]+1) && !gamma_deg[ind].count(v)){
            pq.erase({credits[ind],ind});
            credits[ind]+=1;
            pq.insert({credits[ind],ind});
            gamma_deg[ind].insert(v);
        }
        if((degree[ind]+1) >= gamma*(degree[v]+1) && !gamma_deg[v].count(ind)){
            pq.erase({credits[v],v});
            credits[v]+=1;
            pq.insert({credits[v],v});
            gamma_deg[v].insert(ind);
        }
    }
    
    rounds++;
    //Ogni eta rounds io eseguo delta volte QCextract, costo O(delta*d_max)
    if(rounds >= eta){ 
        //estraggo i delta migliori candidati dalla coda con priorità
        auto it = pq.end();
        int count = 0;
        vector<pair<int,int>> candidates;
        candidates.resize(delta);
        while(it != pq.begin() && count<delta){
            it--;
            candidates[count] = *it;
            count++;
        }

        //valuto in batch se è possibile trovare quasi clique migliori
        for(auto best:candidates){
            int clique_len = sift_num(best.second);
            if(clique_len > best_clique_size){
                best_clique_size = clique_len;
                best_clique = best.second;
                if(best_clique_size != 0){
                    cur_sol.clear();
                    cur_sol.resize(best_clique_size);
                    cur_sol[0]=best_clique;
                    int j = 1;
                    for(int v : neighborhood[best_clique]) {
                        if (ct_score_minhash(best_clique, v) > gamma){
                            cur_sol[j]=v;
                            j++;
                        }
                    }
                }    
            }
        
        }
        rounds = 0;
    }
}

//QCExtract modificato per sostenere l'aggiornamento della chiave della coda con priorità
int CreditsAlgo::sift_num(int node) {
    int res = 0;
    for (int v : neighborhood[node]) {
       if (ct_score_minhash(node, v) > gamma) res++;
    }    
    double tmp = (double)res / (double)(degree[node]+1);
    //aggiorno con il numero di nodi sopra la soglia gamma
    pq.erase({credits[node], node});  
    credits[node] = res;
    pq.insert({credits[node], node}); 
    if (tmp < beta){
        return 0;
    }
    return res + 1;
}

//approssimazione del containment score utilizzando la firma min-hash
double CreditsAlgo::ct_score_minhash(int u, int v){
    int inter = 0;
    for(int i = 0; i<k; ++i){
        if(signature[u][i]==signature[v][i]){
            inter++;
        }
    }
    double js = (double)inter/(double)k;
    return (double)(degree[u] + degree[v] + 2) * js / (js + 1) / (double)(degree[u] + 1);
}

//funzione per ritornare la dimensione della migliore clique
int CreditsAlgo::return_dim(){
    return best_clique_size;
}

//funzione per ritornare la migliore clique
vector<int> CreditsAlgo::returnBestClique(){
    return cur_sol;
}

//funzione ausiliaria per valutare la densità della attuale soluzione
double CreditsAlgo::return_density(){
    if(best_clique_size==0) return 0;
    int cur_sol_edges = 0;
    for (int i : cur_sol) {
        for (int i2 : neighborhood[i]) {
            if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
        }
    }
    return double(cur_sol_edges) / double(best_clique_size* (best_clique_size-1));   
}

#endif