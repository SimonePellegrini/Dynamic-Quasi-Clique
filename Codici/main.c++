#include <iostream>
#include <fstream>
#include "DynamicNBSim.c++"
#include <random>
#include <chrono>
#include <sstream>  
#include <cstdlib>  
#include <algorithm> 
#include "FastNBSim.c++"
#include "CreditsAlgorithm.cpp"
#include "math.h"
#include <format>

using namespace std;
using namespace std::chrono;
int n,m,number_of_experiments;
vector<pair<int,int>> graph;


void read_graph(string filename) {
    clock_t begin = clock();
    ifstream ifile;
    ifile.open(filename);
    int u,v;
    ifile >> n >> m;
    graph.resize(m);
    for (int i = 0; i < m; ++i) {
        ifile >> u >> v;
        graph[i]={u,v};
    }
    std::random_device rd;
    std::mt19937 g(rd());
    shuffle(graph.begin(), graph.end(), g);
    clock_t end = clock();
}
// struct to memorize the operations to be executed
struct GraphOp {
    char type;
    int u;
    int v;
};

vector<vector<GraphOp>> loadDataset(string fileName,float prob, string strategy, string type){
    
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    // open the sequences file
    if(strategy == "standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    } else {
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str +"_"+type+ ".txt");
    }

    if (!file.is_open()) {
        cout << "Error: the file doesn't exist!" << endl; 
    }

    file >> number_of_experiments >> n >> m;

    // save the operations to be executed
    cout << "Loading the dataset..." << endl;
    vector<vector<GraphOp>> all_experiments_ops(number_of_experiments);
    for (int i = 0; i < number_of_experiments; i++) {
        bool end = false;
        while (!end) {
            char type_of_op;
            int u, v;
            file >> type_of_op >> u >> v;
            
            if (type_of_op == 'i' || type_of_op == 'd') {
                all_experiments_ops[i].push_back({type_of_op, u, v});
            } else {
                end = true; 
            }
        }
    }
    file.close();
    return all_experiments_ops;
}

//computes the speed up obtained using the DynamicFastNBSim algorithm
//instead of the FastNBSim Algorithm

void confrontoAlgo1Baseline(string fileName,int num_nodes, int num_edges,int s_rate,int k){
    
    int n_iter = 5;
    double tot_time=0;
    int j = 0;
    double tempo_naive;

    ofstream speed_up("../Esperimenti/"+fileName+"/speedup.txt");
    ofstream tempoNBSim("../Esperimenti/"+fileName+"/tempoNBSim.txt");
    ifstream tempo_baseline_naive("../Esperimenti/"+fileName+"/naiveAvgTime.csv");

    tempo_baseline_naive>>tempo_naive;

    for(int i=0; i<n_iter; i++){
        j = 0;
        auto base =  FastNBSim(num_nodes,0.9,0.6,k);
        //cambia la permutazione del grafo
        read_graph("../Datasets/"+fileName+".txt");
        for(auto x: graph){
            auto start = chrono::high_resolution_clock::now();
            j++;
            //aggiungi l'arco al grafo
            base.add_edge(x.first,x.second);
            if(j%s_rate == 0){
                //usa l'algoritmo statico per trovare una buona quasi-clique
                base.compute_result(false);
            }
            auto end = chrono::high_resolution_clock::now();
            //resetta le variabili di FastNBSim (min-hash,ecc...)
            base.reset_computation();
            //somma il tempo di esecuzione dell'operazione
            tot_time+=std::chrono::duration<double, std::milli>(end - start).count();
        }   
    }
    tot_time = tot_time/n_iter;
    tempoNBSim<<tot_time/graph.size();
    speed_up<<((double)tot_time/graph.size())/tempo_naive;
    speed_up.close();
    tempo_baseline_naive.close();
    tempoNBSim.close();

}

//computes in batch the quasi-cliques obtain using the FastNBSim algorithm

void calcolaSoluzioniFastNBSIM(string fileName,int n,int m,int k){
    
    unordered_map<int,vector<double>> densità;
    unordered_map<int,vector<int>> dimensione;
    auto base =  FastNBSim(n,0.9,0.6,k);
    ofstream confronto_dimensioni_apx("../Esperimenti/"+fileName+"/confrontoDimensioniApx.csv");
    ofstream confronto_densità_apx("../Esperimenti/"+fileName+"/confrontoDensitàApx.csv");
    confronto_dimensioni_apx<<"numero di archi,dimensione"<<endl;
    confronto_densità_apx<<"numero di archi,densità"<<endl;
    int s_rate = m*0.1;
    int n_iter = 100;
    for(int i=0; i<n_iter; i++){
        int j = 0;
        for(auto x: graph){
            base.add_edge(x.first,x.second);
            if(++j%s_rate == 0 || j==m-1){
                base.compute_result(true);
                dimensione[j].push_back(base.getQuasiCliqueSize());
                densità[j].push_back(base.getQuasiCliqueDensity());
                base.reset_computation();
            }
            
        }
    }
    
    for(int i = 0; i < m; i+=s_rate) {
        double tot_densità=0;
        int tot_dimensione=0;
        int c_d = 0;
        int c_dim =0;
        for(auto x:densità[i]){
            if(x!=0){
                tot_densità+=x;
                c_d++;
            }
        }
        for(auto x:dimensione[i]){
            if(x!=0){
                tot_dimensione+=x;
                c_dim++;
            }
        }
        confronto_dimensioni_apx<<i<<","<<tot_dimensione/c_dim<<endl;
        confronto_densità_apx<<i<<","<<tot_densità/c_d<<endl;
    }
    double tot_densità=0;
    int tot_dimensione=0;
    int c_d = 0;
    int c_dim =0;
    for(auto x:densità[m-1]){
        if(x!=0){
            tot_densità+=x;
            c_d++;
        }
    }
    for(auto x:dimensione[m-1]){
        if(x!=0){
            tot_dimensione+=x;
            c_dim++;
        }
    }
    confronto_dimensioni_apx<<m-1<<","<<tot_dimensione/c_dim<<endl;
    confronto_densità_apx<<m-1<<","<<tot_densità/c_d<<endl;
    confronto_dimensioni_apx.close();
    confronto_densità_apx.close();
}

//computes the dimensions of the quasi-cliques obtained after inserting all the edges in the graph
//using the static algorithm

void getQsSize(string fileName,float gamma, float b){
    auto base =  FastNBSim(n,gamma,b);
    for(auto x: graph){
        base.add_edge(x.first,x.second);
    }
    ofstream qsSize("../Esperimenti/"+fileName+"/Qs_Size_b="+to_string(b)+".csv");
    vector<int> qs_size = base.QuasiCliqueSizes();
    qsSize<<qs_size[0];
    for(int i=1; i<qs_size.size(); i++){
        qsSize<<","<<qs_size[i];
    }
    qsSize.close();
}

//helps to 
template <typename T>
void saveToCSV(const std::string& filename, const std::vector<std::vector<T>>& data) {
    std::ofstream outFile(filename);

    if (!outFile.is_open() || data.empty()) return;
    size_t max_time_steps = 0;
    for (const auto& exp : data) {
        if (exp.size() > max_time_steps) max_time_steps = exp.size();
    }
    for (size_t t = 0; t < max_time_steps; ++t) {   
        for (size_t exp = 0; exp < data.size(); ++exp) {
            if (t < data[exp].size()) {
                outFile << data[exp][t];
            }
            if (exp < data.size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
    }
    outFile.close();
}



void countCreditsAndGammeDegree(string fileName, float phi, float alpha, float gamma, float b, int k, float prob) {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");

    if (!file.is_open()) {
        cout<<"error!"<<endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    vector<vector<int>> file_credits;
    vector<vector<int>> file_gamma;


    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
        bool end = false;
        char type_of_op;
        int u, v, edge_inserted = 0;

        while (!end) {
            file >> type_of_op>>u>>v;
            edge_inserted++;
            if (type_of_op == 'i') {
                algo_cred.add_edge(u, v);  
            } 
            else if (type_of_op == 'd') {
                algo_cred.remove_edge(u, v); 
            } 
            else {
                end = true;
            }
        }
        file_credits.push_back(algo_cred.getCredits());
        file_gamma.push_back(algo_cred.getGammaDegree());
    }
    saveToCSV("../Esperimenti/"+fileName+"/number_of_credits_"+p_str+".csv",file_credits);
    file.close();
}

void experimentDynamicOnly(string fileName, float gamma, float b, float prob, float perc = 0.1,string strategy="standard",string type="erdos") {
    string p_str = to_string(prob).substr(0, 3);

    //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);

    vector<vector<int>> dynamicSolution;
    vector<vector<float>> dynamicDensity;

    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_dyn = DynamicNBSim(gamma, b, n, m);
        int number_of_edges=0;
        dynamicSolution.push_back({});
        dynamicDensity.push_back({});

         for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_dyn.add_edge(op.u, op.v);  
            } else if (op.type == 'd') {
                algo_dyn.remove_edge(op.u, op.v); 
            }
            number_of_edges++;
            dynamicSolution.back().push_back(algo_dyn.return_dim());
            if(number_of_edges%(int)(m*perc)==0){
                dynamicDensity.back().push_back(algo_dyn.return_density());
            }
        }
    }
      if(strategy == "standard"){
        saveToCSV("../Esperimenti/"+fileName+"/standard/density_dynamic_"+p_str+".csv",dynamicDensity);
        saveToCSV("../Esperimenti/" + fileName + "/standard/results_dynamic_" + p_str + ".csv", dynamicSolution);
    }
    else{
        saveToCSV("../Esperimenti/"+fileName+"/mixed/density_dynamic_"+type+"_"+p_str+".csv",dynamicDensity);
        saveToCSV("../Esperimenti/" + fileName + "/mixed/results_dynamic_"+type+"_"+p_str + ".csv", dynamicSolution);
    }    
}

void experimentCreditsOnly(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,float perc = 0.1,string strategy="standard",string type="erdos") {
    string p_str = to_string(prob).substr(0, 3);
    //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    vector<vector<int>> creditsSolution;
    vector<vector<float>> creditsDensity;


    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
        creditsSolution.push_back({});
        creditsDensity.push_back({});
        bool end = false;
        int number_of_edges = 0;

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_cred.add_edge(op.u, op.v);
            } else if (op.type == 'd') {
                algo_cred.remove_edge(op.u, op.v); 
            }
            number_of_edges++;
            creditsSolution.back().push_back(algo_cred.return_dim());
            if(number_of_edges%(int)(m*perc)==0){
                creditsDensity.back().push_back(algo_cred.return_density());
            }
        
        }
    }
    if(strategy == "standard"){
        saveToCSV("../Esperimenti/"+fileName+"/standard/density_credits_"+p_str+"_"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/standard/results_credits_" + p_str +"_"+to_string(k)+ ".csv", creditsSolution);
    }
    else{
        saveToCSV("../Esperimenti/"+fileName+"/mixed/density_credits_"+type+"_"+p_str+"_"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/mixed/results_credits_"+type+"_"+p_str +"_alpha:"+to_string(alpha)+"_phi:"+to_string(phi)+"_k:"+to_string(k)+ ".csv", creditsSolution);
    }    
}

void experimentCreditsPerformance(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,string strategy="standard",string type="erdos",int q_freq=10) {
    string p_str = to_string(prob).substr(0, 3);
     //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    

    vector<vector<double>> executionTimes;

    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 
        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
        int n_op=0;
        auto start_time = std::chrono::high_resolution_clock::now();

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_cred.add_edge(op.u, op.v);  
            } else if (op.type == 'd') {
                algo_cred.remove_edge(op.u, op.v); 
            }
            n_op++;
            if(n_op==q_freq){
                algo_cred.return_dim();
                n_op=0;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;

        executionTimes.back().push_back(duration.count());
    }

    if(strategy == "standard"){
        saveToCSV("../Esperimenti/" + fileName + "/standard/times_credits_" + p_str +"_" +to_string(k) + ".csv", executionTimes);
    }
    else{
        saveToCSV("../Esperimenti/" + fileName + "/mixed/times_credits_"+type+"_" + p_str +"_alpha:"+to_string(alpha)+"_phi"+to_string(phi)+"_k:"+to_string(k)+ ".csv",executionTimes);
    } 
}

void experimentDynamicPerformance(string fileName, float gamma, float b, float prob,string strategy = "standard",string type="erdos",int q_freq=10) {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;

    //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    vector<vector<double>> executionTimes;

    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 

        auto algo_dyn = DynamicNBSim(gamma, b, n, m); 
        int n_op = 0;
        auto start_time = std::chrono::high_resolution_clock::now();

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_dyn.add_edge(op.u, op.v);  
            } else if (op.type == 'd') {
                algo_dyn.remove_edge(op.u, op.v); 
            }
            n_op++;
            if(n_op==q_freq){
                algo_dyn.return_dim();
                n_op=0;
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;
        executionTimes.back().push_back(duration.count());
    }
    if(strategy == "standard"){
        saveToCSV("../Esperimenti/" + fileName + "/standard/times_dynamic_" + p_str + ".csv", executionTimes);
    }
    else{
        saveToCSV("../Esperimenti/" + fileName + "/mixed/times_dynamic_"+type+"_" + p_str + ".csv", executionTimes);
    }
}


void gridSearchK(string fileName, float phi, float alpha, float gamma, float b, float prob, string strategy = "standard", string type = "erdos", int q_freq = 10) {
   
    cout << "Running dynamic baseline..." << endl;
    experimentDynamicOnly(fileName, gamma, b, prob, 0.1, strategy, type);
    experimentDynamicPerformance(fileName, gamma, b, prob, strategy, type, q_freq);
        
    vector<int> k_values = {32, 64, 80, 128, 150, 256};
    vector<double> alpha_values = {0.12};
    vector<double> phi_values = {0.6};

    for (int k : k_values) {
        for(double alpha : alpha_values){
            for(double phi : phi_values){
                experimentCreditsOnly(fileName, phi, alpha, gamma, b, k, prob, 0.1, strategy, type);
                experimentCreditsPerformance(fileName, phi, alpha, gamma, b, k, prob, strategy, type, q_freq);
            }
        }    
    }
    
    cout << "Grid search completed." << endl;
}




int main(int argc, const char * argv[]){
    /*prende dalla riga di comando i parametri in input*/
  
    float phi = atof(argv[1]);
    float alpha = atof(argv[2]);
    string fileName = argv[4];
    double gamma = 0.9;
    double b = 0.6;
    int k = atoi(argv[3]);
    double p = atof(argv[5]);
    /* Inizializza il dataset da leggere */
    read_graph("../Datasets/"+fileName+".txt");


    gridSearchK(fileName,phi,alpha,gamma,b,p,"standard");
    gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","power");
    gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","rich");
    gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","erdos");
    
    //gridSearchCreditsAlgorithm(fileName, gamma,b,p);

    cout << "Tutti gli esperimenti completati con successo!" << endl;
    return 0;
}

