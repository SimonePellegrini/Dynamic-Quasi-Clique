#ifndef NAIVE_ALGO
#define NAIVE_ALGO
#include <fstream>
#include <math.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>


using namespace std;

class DynamicNBSim
{
    private:
        double gamma,b;
        int n,m,i=0;
        vector<unordered_set<int>> graph;
        vector<unordered_map<int,double>> containmentScore;
        vector<set<int>> quasiClique;
        vector<int> NeighSize;
        set<pair<int,int>> pq;
        vector<double> rapporto;
        vector<int> quasiCliqueSize;
    
    public:
        DynamicNBSim(double gamma, double b,int n,int m);
        set<int> getBestClique();
        void Insert(int u, int v);
        double ct_score(int v1, int v2);
        int getBestCliqueSize();
        double sift_num(int node);
};

double DynamicNBSim::ct_score(int v1, int v2) {
    int res = 0;
    unordered_set<int> st;
    st.insert(v1);
    for (int num : graph[v1]) st.insert(num);
    for (int num : graph[v2]){
        if (st.count(num)) res++;
    }
    if (st.count(v2)) res++;
    return (double)res / (double)NeighSize[v1];
}

double DynamicNBSim::sift_num(int node) {
    int res = 0;
    
    for (int v : graph[node]) {
        if (ct_score(node,v) >= gamma) res++;
    }

    double tmp = (double)res / (double)(NeighSize[node]);
    if(tmp>=b){
        return res+1;
    }
    return 0;
 
}

int DynamicNBSim::getBestCliqueSize(){
    auto elem = *prev(pq.end());
    return elem.first;
}

DynamicNBSim::DynamicNBSim(double gamma,double b, int n, int m){
    /*
    input:
    gamma,b: parametri cinesi-like
    n: il numero di nodi nel grafo
    m: il numero di archi nel grafo
    */

    this->gamma = gamma;
    this->b = b;
    this->n = n;
    this->m = m;

    graph.resize(n);
    containmentScore.resize(n);
    NeighSize.resize(n);
    rapporto.resize(n);
    quasiCliqueSize.resize(n);

    for(int i = 0; i<n; i++){
        NeighSize[i] = 1;
    } 
    quasiClique.resize(n);
    for(int i = 0; i<n; i++){
        quasiClique[i].insert(i);
    }
}

void DynamicNBSim::Insert(int u,int v){
    int arr[2] = {u,v};
    for(int option: {0,1}){
    auto b_x = option^1;
    int x = arr[option];
    int y = arr[b_x];
    for(int w: graph[x])
    {
        if(graph[y].count(w)){
            //caso in cui w appartiene all'intersezione tra N(u) e N(v)

            //ricalcolo t(x,w)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x]+1)/(NeighSize[x]+1);
            //se il cont. score è maggiore di gamma e w non sta nella quasiClique di x, aggiungilo
            if(containmentScore[x][w]>=gamma && !(quasiClique[x].count(w))){ 
                quasiClique[x].insert(w);
            }

            //ricalcolo t(w,x) e aggiorno le quasi-clique di conseguenza
            containmentScore[w][x] = (containmentScore[w][x]*NeighSize[w]+1)/NeighSize[w];
            if(containmentScore[w][x]>=gamma && !(quasiClique[w].count(x))){
                //aggiorno la chiave nella coda con priorità
                pq.erase({quasiClique[w].size(), w});
                quasiClique[w].insert(x);
                quasiCliqueSize[w]++;
                double new_ratio = (double)(quasiClique[w].size() - 1) / NeighSize[w];
                if(new_ratio >= b) {
                    pq.insert({quasiClique[w].size(), w});
                }
            }
        }
        else{
            //caso in cui w non appartiene all'intersezione tra N(u) e N(v)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x])/(NeighSize[x]+1);
            if(containmentScore[x][w]<gamma && quasiClique[x].count(w)){   
                quasiClique[x].erase(w);
            }
        }   
    } 
    }
    NeighSize[u]++;
    NeighSize[v]++;
    graph[u].insert(v);
    graph[v].insert(u);
    containmentScore[u][v] = ct_score(u,v);
    containmentScore[v][u] = ct_score(v,u);
    //aggiorno le chiavi associate alle clique nella coda con priorità
    pq.erase({quasiCliqueSize[v], v});
    if(containmentScore[v][u]>=gamma){ 
        quasiClique[v].insert(u);
    } 
    quasiCliqueSize[v] = quasiClique[v].size();
    double ratiov = (double)(quasiClique[v].size() - 1) / NeighSize[v];
    if(ratiov >= b) {
        pq.insert({quasiClique[v].size(), v});    
    }
    
    pq.erase({quasiCliqueSize[u], u});
    if(containmentScore[u][v]>=gamma){
        quasiClique[u].insert(v);
    }
    quasiCliqueSize[u] = quasiClique[u].size();
    double ratiou = (double)(quasiClique[u].size() - 1) / NeighSize[u];
    if(ratiou >= b) {
        pq.insert({quasiClique[u].size(), u});    
    } 
}

#endif