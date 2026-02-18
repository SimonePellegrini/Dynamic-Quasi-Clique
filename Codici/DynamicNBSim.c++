#ifndef DYNAMICNBSIM
#define DYNAMICNBSIM
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
        void add_edge(int u, int v);
        double ct_score(int v1, int v2);
        int return_dim();
        void remove_edge(int u, int v);
        double return_density();
};

//compute the containment score
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

//return the density of the best clique found 
double DynamicNBSim::return_density(){
    int best_clique = pq.rbegin()->second;
    set<int> cur_sol = quasiClique[best_clique];
    if(cur_sol.size()==0){
        return 0;
    }
    int cur_sol_edges = 0;
    for (int i : cur_sol) {
        for (int i2 : graph[i]) {
            if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
        }
    }
    return double(cur_sol_edges) / double(cur_sol.size()* (cur_sol.size()-1));   
}

//return the size of the best quasi-clique found
int DynamicNBSim::return_dim(){
    auto elem = *prev(pq.end());
    return elem.first;
}


DynamicNBSim::DynamicNBSim(double gamma,double b, int n, int m){
    /*
    input:
    gamma,b: same parameters as FastNBSim
    n: number of nodes in the graph
    m: number of edges in the graph
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

//remove an edge from the graph and update dynamically the solution mantained
void DynamicNBSim::remove_edge(int u,int v){
    if(graph[u].count(v) && graph[v].count(u)){
    int arr[2] = {u,v};
    graph[u].erase(v);
    graph[v].erase(u);
    for(int option: {0,1}){
    auto b_x = option^1;
    int x = arr[option];
    int y = arr[b_x];
    for(int w: graph[x])
    { 
        if(graph[y].count(w)){
            //w belongs to the intersection of N(u) and N(v)

            //update t(x,w)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x]-1)/(NeighSize[x]-1);

            //update the quasi-clique associated to x, if necessary
            if(containmentScore[x][w]<gamma && quasiClique[x].count(w)){ 
                pq.erase({quasiClique[x].size(), x});
                quasiClique[x].erase(w);
                quasiCliqueSize[x]--;
                double new_ratio = (double)(quasiClique[x].size() - 1) / (NeighSize[x]-1);
                if(new_ratio >= b) {
                    pq.insert({quasiClique[x].size(), x});
                }
            }

            //update t(w,x) 
            containmentScore[w][x] = (containmentScore[w][x]*NeighSize[w]-1)/NeighSize[w];

            //update the quasi-clique associated to w, if necessary
            if(containmentScore[w][x]<gamma && quasiClique[w].count(x)){
                pq.erase({quasiClique[w].size(), w});
                quasiClique[w].erase(x);
                quasiCliqueSize[w]--;
                double new_ratio = (double)(quasiClique[w].size() - 1) / NeighSize[w];
                if(new_ratio >= b) {
                    pq.insert({quasiClique[w].size(), w});
                }
            }
        }
        else{
            //w doesn't belong to the interesction between N(u) and N(v)

            //update t(x,w)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x])/(NeighSize[x]-1);

            //update the quasi-clique associated to x, if necessary
            if(containmentScore[x][w]>=gamma && !(quasiClique[x].count(w))){   
                pq.erase({quasiClique[x].size(), x});
                quasiClique[x].insert(w);
                quasiCliqueSize[x]++;
                double new_ratio = (double)(quasiClique[x].size() - 1) / (NeighSize[x]-1);
                if(new_ratio >= b) {
                    pq.insert({quasiClique[x].size(), x});
                } 
            }
        }   
    } 
    }
    NeighSize[u]--;
    NeighSize[v]--; 
    containmentScore[u][v] = 0;
    containmentScore[v][u] = 0;
    
    //update the containment scores and the quasi-cliques associated to u and b
    pq.erase({quasiCliqueSize[v], v});
    quasiClique[v].erase(u);
    quasiCliqueSize[v] = quasiClique[v].size();
    double ratiov = (double)(quasiClique[v].size() - 1) / NeighSize[v];
    if(ratiov >= b) {
        pq.insert({quasiClique[v].size(), v});    
    }
    
    pq.erase({quasiCliqueSize[u], u});
    quasiClique[u].erase(v);
    quasiCliqueSize[u] = quasiClique[u].size();
    double ratiou = (double)(quasiClique[u].size() - 1) / NeighSize[u];
    if(ratiou >= b) {
        pq.insert({quasiClique[u].size(), u});    
    } 
    }
}

//add an edge to the graph and update dynamically the solution 
void DynamicNBSim::add_edge(int u,int v){
    if(!(graph[u].count(v))&&!(graph[v].count(u))){
    int arr[2] = {u,v};
    for(int option: {0,1}){
    auto b_x = option^1;
    int x = arr[option];
    int y = arr[b_x];
    for(int w: graph[x])
    {
        if(graph[y].count(w)){
            //w belongs to the intersection between N(u) and N(v)

            //update t(x,w)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x]+1)/(NeighSize[x]+1);

            //update the quasi-clique associated to x, if necessary
            if(containmentScore[x][w]>=gamma && !(quasiClique[x].count(w))){ 
                pq.erase({quasiClique[x].size(), x});
                quasiClique[x].insert(w);
                quasiCliqueSize[x]++;
                double new_ratio = (double)(quasiClique[x].size() - 1) / (NeighSize[x]+1);
                if(new_ratio >= b) {
                    pq.insert({quasiClique[x].size(), x});
                }
            }

            //update t(w,x)
            containmentScore[w][x] = (containmentScore[w][x]*NeighSize[w]+1)/NeighSize[w];

            //update the quasi-clique associated to w, if necessary
            if(containmentScore[w][x]>=gamma && !(quasiClique[w].count(x))){
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
            //update t(x,w)
            containmentScore[x][w] = (containmentScore[x][w]*NeighSize[x])/(NeighSize[x]+1);

            //update the quasi-clique associated to x, if necessary
            if(containmentScore[x][w]<gamma && quasiClique[x].count(w)){   
                pq.erase({quasiClique[x].size(), x});
                quasiClique[x].erase(w);
                quasiCliqueSize[x]--;
                double new_ratio = (double)(quasiClique[x].size() - 1) / (NeighSize[x]+1);
                if(new_ratio >= b) {
                    pq.insert({quasiClique[x].size(), x});
                } 
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
    
    //update the quasi-cliques and the containments scores of u and v
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
}

#endif