
#ifndef FastNBSim_h
#define FastNBSim_h
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
using namespace std;

class FastNBSim {
private:
    bool use_signature = false;
    int n, m;
    int cur_sol_vertex = 0;
    int cur_sol_size = 0;
    int cur_sol_edges = 0;
    int operations = 0;
    int k;
    int p = pow(2,31) - 1;
    double ts;
    double cur_sol_quasi;
    double f;
    double b;
    vector<vector<int>> neighbor;
    vector<int> degree;
    vector<int> kcore_size;
    vector<pair<int,int>> sorted_vertex;
    vector<int> cur_sol;
    vector<int> vec_a;
    vector<int> vec_b;
    vector<int> vis;
    vector<vector<int>> signature;
    
    void initial_variable();
    void get_k_signatures(int node);
    void sort_vertex();
    void iterate_vertex();
    int core_num(int node);
    int sift_num(int node);
    double ct_score(int v1, int v2);
    double ct_score_minhash(int v1, int v2);
    void get_res(bool density);
    
public:
    FastNBSim(int num_nodes, double threshold_ts, double threshold_b, int k_param = 0);
    void add_edge(int u, int v);
    void compute_result(bool density);
    void print_result(string filename = "batch");
    void clear_edges();
    void reset_computation();
    int getQuasiCliqueSize();
    double getQuasiCliqueDensity();
    vector<int> QuasiCliqueSizes();
};

// Implementazioni inline per evitare problemi di linking
FastNBSim::FastNBSim(int num_nodes, double threshold_ts, double threshold_b, int k_param) {
    n = num_nodes;
    m = 0;
    ts = threshold_ts;
    b = threshold_b;
    k = k_param;
    use_signature = (k > 0);
    initial_variable();
}

vector<int> FastNBSim::QuasiCliqueSizes(){
    vector<int> Qs_Size;
    Qs_Size.resize(n);
    for(int i=0; i<n; i++){
        Qs_Size[i] = sift_num(i);
    }
    return Qs_Size;

}

void FastNBSim::initial_variable() {
    degree.assign(n, 0);
    neighbor.resize(n);
    kcore_size.resize(n);
    vis.resize(n);
    if (use_signature) {
        vec_a.resize(k);
        vec_b.resize(k);
        signature.resize(n);
        for (int i = 0; i < k; i++) {
            vec_a[i] = rand() % (p - 1) + 1;
            vec_b[i] = rand();
        }
    }
}

void FastNBSim::add_edge(int u, int v) {
    if (u >= n || v >= n || u < 0 || v < 0) return;
    degree[u]++;
    degree[v]++;
    neighbor[u].push_back(v);
    neighbor[v].push_back(u);
    m++;
}

void FastNBSim::get_k_signatures(int node) {
    signature[node].clear();
    for (int index = 0; index < k; ++index) {
        int minhash = p;
        for (int j : neighbor[node]) {
            minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(j) + (long long)(vec_b[index])) % p));
        }
        minhash = min(minhash, int(((long long)(vec_a[index]) * (long long)(node) + (long long)(vec_b[index])) % p));
        signature[node].push_back(minhash);
    }
}

void FastNBSim::sort_vertex() {
    sorted_vertex.clear();
    for (int i = 0; i < n; ++i) {
        int num = core_num(i);
        if (num > 0) sorted_vertex.push_back({num, i});
    }
    sort(sorted_vertex.begin(), sorted_vertex.end(), greater<pair<int,int>>());
}

void FastNBSim::iterate_vertex() {
    for (int i = 0; i < sorted_vertex.size(); ++i) {
        auto p = sorted_vertex[i];
        if (p.first > cur_sol_size) {
            operations++;
            int num = sift_num(p.second);
            if (num > cur_sol_size) {
                cur_sol_size = num;
                cur_sol_vertex = p.second;
            }
        }
    }
}

int FastNBSim::core_num(int node) {
    int num = 0;
    for (int v : neighbor[node]) {
        if ((degree[v]+1) >= (degree[node]+1) * ts) num++;
    }
    return num;
}

double FastNBSim::ct_score(int v1, int v2) {
    int res = 0;
    unordered_set<int> st;
    st.insert(v1);
    for (int num : neighbor[v1]) st.insert(num);
    for (int num : neighbor[v2]) if (st.count(num)) res++;
    if (st.count(v2)) res++;
    return (double)res / (double)(degree[v1] + 1);
}

double FastNBSim::ct_score_minhash(int v1, int v2) {
    int res = 0;
    for (int i = 0; i < k; ++i) {
        if (signature[v1][i] == signature[v2][i]) res++;
    }
    double j_sim = (double)res / (double)k;
    return (double)(degree[v1] + degree[v2] + 2) * j_sim / (j_sim + 1) / (double)(degree[v1] + 1);
}

int FastNBSim::sift_num(int node) {
    int res = 0;
    if (use_signature){
       if (vis[node] == 0) {
           get_k_signatures(node);
           vis[node] = 1;
       }
       for (int v : neighbor[node]) {
           if (vis[v] == 0) {
               get_k_signatures(v);
               vis[v] = 1;
           }
           if (ct_score_minhash(node, v) >= ts) res++;
       }
    }
    else {
        for (int v : neighbor[node]) {
          if (ct_score(node, v) >= ts) res++;
        }
    }
    double tmp = (double)res / (double)(degree[node]+1);
    if (tmp < b) return 0;
    return res+1;
}

void FastNBSim::get_res(bool density) {
    sort_vertex();
    iterate_vertex();
    cur_sol.clear();
    cur_sol.push_back(cur_sol_vertex);
    cur_sol_edges = 0;
    int j = 0;
   
    for (int v : neighbor[cur_sol_vertex]) {
        if (use_signature) {
            if (ct_score_minhash(cur_sol_vertex, v) >= ts) cur_sol.push_back(v);
        }
        else {
            if (ct_score(cur_sol_vertex, v) >= ts) cur_sol.push_back(v);
        }
    }
        
    if(density){    
        for (int i : cur_sol) {
            for (int i2 : neighbor[i]) {
                if(find(cur_sol.begin(), cur_sol.end(), i2) != cur_sol.end()) cur_sol_edges++;
            }
        }
        cur_sol_quasi = double(cur_sol_edges) / double(cur_sol_size * (cur_sol_size -1));
    }
}

void FastNBSim::compute_result(bool density) {
    get_res(density);
}

int FastNBSim::getQuasiCliqueSize(){
    return cur_sol_size;
}

double FastNBSim::getQuasiCliqueDensity(){
    return cur_sol_quasi;
}

void FastNBSim::print_result(string filename) {
    cout << filename << " " << ts << " " << b << " " << k << " " << cur_sol_size << " " << cur_sol_quasi << endl;
}

void FastNBSim::clear_edges() {
    for (int i = 0; i < n; i++) {
        neighbor[i].clear();
        degree[i] = 0;
    }
    m = 0;
}

void FastNBSim::reset_computation() {
    cur_sol_vertex = 0;
    cur_sol_size = 0;
    cur_sol_edges = 0;
    operations = 0;
    cur_sol_quasi = 0;
    sorted_vertex.clear();
    cur_sol.clear();
    fill(vis.begin(), vis.end(), 0);
    if (use_signature) {
        for (int i = 0; i < n; i++) {
            signature[i].clear();
        }
    }
}

#endif 