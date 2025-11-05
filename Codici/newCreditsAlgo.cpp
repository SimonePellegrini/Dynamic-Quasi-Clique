#ifndef CREDITS_ALGO_IMPROVED
#define CREDITS_ALGO_IMPROVED
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

int p = pow(2,31)-1;

class CredAlgoImproved{
    
    private:
    
        int k,n,m;
        float alpha,phi;
        double gamma,beta;
        int best_clique=0;
        int best_clique_size=0; 

        vector<vector<int>> neighborhood;
        vector<vector<int>> signature;
        vector<int> a;
        vector<int> b;
        vector<int> credits;
        vector<int> cur_sol;
        vector<int> degree;
        vector<int> past_credits; 

        void visit_node(int u);
    

    public:
        CredAlgoImproved(int n,int m, float phi,float alpha, double gamma,double beta,int k);
        vector<int> getCredits();
        vector<int> computeGammaDegree();
        vector<int> returnBestClique();
        double ct_score_minhash(int u, int v);
        void add_edge(int u, int v);
        int sift_num(int u);
        void update_signature(int u,int v); 
        int return_dim();
        double return_density();
};

CredAlgoImproved::CredAlgoImproved(int n,int m, float phi,float alpha, double gamma,double beta,int k){

    this->n = n;
    this->m = m;
    this->phi = phi;
    this->alpha = alpha;
    this->gamma = gamma;
    this->beta = beta;
    this->k = k;

    //inizializzo i vari vettori
    neighborhood.resize(n);
    credits.resize(n);
    a.resize(k);
    b.resize(k);
    degree.resize(n);
    signature.resize(n);
    past_credits.resize(n);

    for(int i = 0; i<n; i++){
        past_credits[i]++;
    }

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
}

void CredAlgoImproved::update_signature(int u,int v){
    
    //aggiorno la signature di N(u) 
    for(int i = 0; i<k; i++){
        signature[u][i] = min(signature[u][i],int(((long long)(a[i]) * (long long)(v) + (long long)(b[i])) % p));
    }

    //aggiorno la signature di N(v) 
    for(int i = 0; i<k; i++){
        signature[v][i] = min(signature[v][i],int(((long long)(a[i]) * (long long)(u) + (long long)(b[i])) % p));
    }
}

//visita il vicinato di un nodo e vede se è possibile estrarre qualcosa
void CredAlgoImproved::visit_node(int u){
    int qs_size = sift_num(u);
   
    if(qs_size>best_clique_size){
        best_clique = u;
        best_clique_size = qs_size;
        if(best_clique_size != 0){
                cur_sol.clear();
                cur_sol.resize(best_clique_size);
                cur_sol[0]=best_clique;
                int j = 1;
                for(int n : neighborhood[best_clique]) {
                    if (ct_score_minhash(best_clique, n) > gamma){
                        cur_sol[j]=n;
                        j++;
                    }
                }
            }  
    }
}

void CredAlgoImproved::add_edge(int u, int v){
    //aggiungo l'arco al grafo
    neighborhood[v].emplace_back(u);
    neighborhood[u].emplace_back(v);
    degree[v]++;
    degree[u]++;
    //aggiorno la signature di N(u) e N(v), costo O(k)
    update_signature(u,v);

    if((degree[u]+1) >= gamma*(degree[v]+1)){
        credits[v]+=1;
    }

    if((degree[v]+1) >= gamma*(degree[u]+1)){
        credits[u]+=1;
    }
    if((float)credits[u]/past_credits[u]>alpha && credits[u]>phi*best_clique_size){
        visit_node(u);
        past_credits[u] = credits[u];
    }

    if((float)credits[v]/past_credits[v]>alpha && credits[v]>phi*best_clique_size){ 
        visit_node(v);
        past_credits[v] = credits[v];
    }
}


//prende in input un nodo e restituisce la dimensione della quasi-clique estraibile da quel nodo

int CredAlgoImproved::sift_num(int node) {
    int res = 0;
    for (int v : neighborhood[node]) {
       if (ct_score_minhash(node, v) > gamma) res++;
    }    
    double tmp = (double)res / (double)(degree[node]+1);
    
    if (tmp < beta){
        return 0;
    }
    return res + 1;
}

//approssimazione del containment score utilizzando la firma min-hash

double CredAlgoImproved::ct_score_minhash(int u, int v){
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

int CredAlgoImproved::return_dim(){
    return best_clique_size;
}

//funzione per ritornare la migliore clique

vector<int> CredAlgoImproved::returnBestClique(){
    return cur_sol;
}

//funzione per ritornare tutti i crediti associati ai nodi del grafo

vector<int> CredAlgoImproved::getCredits(){
    return credits;
}

//calcola il gamma degree associato ai nodi del grafo

vector<int> CredAlgoImproved::computeGammaDegree(){
    vector<int> gammaDeg;
    gammaDeg.resize(n);
    for(int i = 0; i < n; i++){
        gammaDeg[i]++;
        for(auto x: neighborhood[i]){
            if(degree[x]+1 >= gamma*(degree[i]+1)){
                gammaDeg[i]++;
            }
        }
    }
    return gammaDeg;
}

//funzione ausiliaria per valutare la densità della attuale soluzione

double CredAlgoImproved::return_density(){
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