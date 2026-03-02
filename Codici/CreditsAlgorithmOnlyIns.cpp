#ifndef CREDITS_ALGO_IMPROVED2
#define CREDITS_ALGO_IMPROVED2
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
#include "DynamicMinHash.cpp"
#include <unordered_map>
#include <cmath>


using namespace std;

int p = pow(2,31)-1;

struct pair_hash {
    size_t operator()(const pair<int,int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

class creditsAlgorithm{
    
    private:
    
        int k,n,m;
        float alpha,phi;
        double gamma,beta;
        int best_clique_size=0; 
        int count_faults = 0;
        int count_removes = 0;
        vector<set<int>> neighborhood;
        vector<DynamicMinHash*> signature; 
        TabulationHash **hashes;
        vector<int> credits;
        vector<int> degree;
        vector<int> past_credits; 
        unordered_map<pair<int,int>, bool, pair_hash> received_credits;
        vector<int> sol_size;
        multiset<pair<int,int>> pq;
        vector<int> visit_node(int u);
        vector<int> best_clique;
        vector<int> obtain_solution(int node);
        vector<vector<int>> signatures;      // MinHash signatures
        vector<int> a,b;    

    public:
        creditsAlgorithm(int n,int m, float phi,float alpha, double gamma,double beta,int k);
        ~creditsAlgorithm();
        vector<int> getCredits();
        vector<int> getGammaDegree();
        vector<int> returnBestClique();
        double ct_score_minhash(int u, int v);
        void add_edge(int u, int v); 
        void incremental_minhash_update(int u,int v); 
        int return_dim();
        double return_density();
        void remove_edge(int u, int v);
        vector<int> count;
        double ct_scores(int v1, int v2);
};



creditsAlgorithm::~creditsAlgorithm(){

   for(int i=0; i<n; i++){ 
        if(signature[i] != nullptr) delete signature[i];
    }

    if(hashes != nullptr) {
        for(int i=0; i<k; i++) {
            delete hashes[i]; 
        }
        delete[] hashes;     
    }
};

creditsAlgorithm::creditsAlgorithm(int n,int m, float phi,float alpha, double gamma,double beta,int k){

    this->n = n;
    this->m = m;
    this->phi = phi;
    this->alpha = alpha;
    this->gamma = gamma;
    this->beta = beta;
    this->k = k;
    this->best_clique_size=0;
    a.resize(k);
    b.resize(k);
    //init the arrays
    neighborhood.resize(n);
    signatures.resize(n);
    credits.resize(n);
    degree.resize(n);
    signature.resize(n);
    past_credits.resize(n);
    count.resize(n);
    sol_size.resize(n);

    
    for(int i = 0; i<n; i++){
        past_credits[i]=1;
    }

        

    hashes = DynamicMinHash::createHashFunctions(k);
    //init the dynamic minhash signatures
    for(int i=0; i<n; i++){ 
        signature[i] = new DynamicMinHash(k,1, hashes);
        signature[i]->insert(i);
    }

     // inizializzo funzioni hash
    for(int i=0;i<k;i++){
        a[i] = rand()%(p-1)+1;
        b[i] = rand();
    }

    // inizializzo firme (self-loop)
    for(int i=0;i<n;i++){
        signatures[i].resize(k);
        for(int j=0;j<k;j++){
            signatures[i][j] =((long long)a[j]*i + b[j]) % p;
        }
        pq.insert({0,i});
    }

    for(int i = 0; i < n; ++i){
        pq.insert({0,i});
    }
}


//remove the edge and update credits and quasi-cliques
void creditsAlgorithm::remove_edge(int u, int v){
    
    //neighborhood[u].erase(v);
    //neighborhood[v].erase(u);
    degree[v]--;
    degree[u]--;
    
    //update the minhash signature
    if(signature[u]->remove(v)){
        signature[u]->insert(u);
        for(auto x:neighborhood[u]){
            signature[u]->insert(x);
        }
    }


    if(signature[v]->remove(u)){
        signature[v]->insert(v);
        for(auto x:neighborhood[v]){
            signature[v]->insert(x);
        }
    }

    //update the credits distribuited
    if(received_credits[{u,v}]){
        credits[v]--;
        received_credits[{u,v}]=false;
        count[v]++;
    }
    if(received_credits[{v,u}]){
        credits[u]--;
        received_credits[{v,u}]=false;
        count[u]++;
    }
    
    //check wheter it is necessary to update the quasi-clique

    if(past_credits[u]+count[u]>=(1+alpha)*past_credits[u] && sol_size[u]>=phi*(pq.rbegin()->first)){
        count[u] = 0;
        past_credits[u] = credits[u];
        pq.erase({sol_size[u],u});
        sol_size[u] = visit_node(u).size();
        pq.insert({sol_size[u],u});
    }

    if(past_credits[v]+count[v]>=(1+alpha)*past_credits[v] && sol_size[v]>=phi*(pq.rbegin()->first)){
        count[v] = 0;
        past_credits[v] = credits[v];
        pq.erase({sol_size[v],v});
        sol_size[v] = visit_node(v).size();
        pq.insert({sol_size[v],v});
    }
    
}

//add the edge to the graph and update credits and quasi-cliques, if necessary
void creditsAlgorithm::add_edge(int u, int v){
    
    //insert the edge in the graph
    neighborhood[v].insert(u);
    neighborhood[u].insert(v);
    degree[v] = neighborhood[v].size();
    degree[u] =  neighborhood[u].size();

    for(int i=0;i<k;i++){
        int hv = ((long long)a[i]*v + b[i]) % p;
        signatures[u][i] = min(signatures[u][i], hv);

        int hu = ((long long)a[i]*u + b[i]) % p;
        signatures[v][i] = min(signatures[v][i], hu);
    }
    
    //update credits
    if(degree[u] >= gamma*(degree[v])){
        credits[v]++;
        //received_credits[{u,v}]= true;
        //count[v]++;
    }

    if(degree[v] >= gamma*(degree[u])){
        credits[u]++;
        //received_credits[{v,u}]= true;
        //count[u]++;
    }
    
    //check wheter it is necessary to update the quasi-clique
    if(credits[u]>=(1+alpha)*past_credits[u] && credits[u]>=phi*best_clique_size){
        //count[u] = 0;
        past_credits[u] = credits[u];
        vector<int> sol = visit_node(u);
        if(sol_size[u]>best_clique_size){
            best_clique = sol;
            best_clique_size = sol.size();
        }
    }

    if(credits[v]>=(1+alpha)*past_credits[v] && credits[v]>=phi*best_clique_size){
        //count[v] = 0;
        past_credits[v] = credits[v];
        vector<int> sol = visit_node(v);
        if(sol.size()>best_clique_size){
            best_clique = sol;
            best_clique_size = sol.size();
        }
    }
}



//takes in input a node and try to extract a quasi-clique out of it

vector<int> creditsAlgorithm::visit_node(int node) {
    vector<int> sol = obtain_solution(node);
    int res = sol.size();
    double tmp = (double)res / (double)(degree[node]+1);
    
    if (tmp < beta){
        return {};
    }
    return sol;
}

//compute the set S associated to a node

vector<int> creditsAlgorithm::obtain_solution(int node){
    vector<int> sol;
    sol.push_back(node);
    for (int v : neighborhood[node]) {
       if (ct_score_minhash(node, v) >= gamma) sol.push_back(v);
    }    
    return sol;
}

//min-hash approximation of the containment score

double creditsAlgorithm::ct_score_minhash(int u, int v){
    int inter = 0;

    for(int i=0;i<k;i++){
        if(signatures[u][i] == signatures[v][i])
            inter++;
    }

    double js = (double)inter/(double)k;
    return (double)(degree[u] + degree[v] + 2) * js / (js + 1) / (double)(degree[u] + 1);
}

//function to compute the containment score

double creditsAlgorithm::ct_scores(int v1, int v2) {
    int res = 0;
    unordered_set<int> st;
    st.insert(v1);
    for (int num : neighborhood[v1]) st.insert(num);
    for (int num : neighborhood[v2]) if (st.count(num)) res++;
    if (st.count(v2)) res++;
    return (double)res / (double)(degree[v1] + 1);
}

//return the best clique mantained

int creditsAlgorithm::return_dim(){
    /*
    // take the node on the top of the priority que
    int node = pq.rbegin()->second;
    vector<int> solution = obtain_solution(node);
    best_clique = node;
    
    //return its size
    return solution.size();
    */
   return best_clique_size;
}

//return the number of credits associated to each node in the graph
vector<int> creditsAlgorithm::getCredits(){
    return credits;
}

//computes and return the gamma-degree associated to each node in the graph
vector<int> creditsAlgorithm::getGammaDegree(){
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

//return the density associated to the best clique found so far
#include <algorithm> // Necessario per utilizzare std::find

double creditsAlgorithm::return_density(){
    vector<int> cur_sol = best_clique;
    
    // Se la clique ha 0 o 1 nodo, la densità è 0
    if(cur_sol.size() <= 1){
        return 0.0;
    }
    
    int cur_sol_edges = 0;
    for (int i : cur_sol) {
        for (int i2 : cur_sol) {
            if(i != i2) {
                // Verifichiamo se i2 è presente nel vettore dei vicini di i
                bool i_has_i2 = std::find(neighborhood[i].begin(), neighborhood[i].end(), i2) != neighborhood[i].end();
                
                // Verifichiamo se i è presente nel vettore dei vicini di i2
                bool i2_has_i = std::find(neighborhood[i2].begin(), neighborhood[i2].end(), i) != neighborhood[i2].end();
                
                if(i_has_i2 && i2_has_i){
                    cur_sol_edges++;
                }
            }
        }
    }
        
    return double(cur_sol_edges) / double(cur_sol.size() * (cur_sol.size() - 1));  
}

#endif