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
int n,m;
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

void gridSearchCreditsAlgorithm(string fileName, float gamma, float b, float prob, string strategy = "standard") {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    
    // open the sequences file
    if(strategy == "standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    } else {
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str + ".txt");
    }

    if (!file.is_open()) {
        cout << "Errore: impossibile aprire il file di sequenza per la Grid Search!" << endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    // save the operations to be executed
    cout << "Caricamento operazioni in memoria..." << endl;
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

    // grid search parameters
    vector<int> k_values = {32,64,128,256};
    vector<float> phi_values = {0.1, 0.3, 0.5, 0.7,0.9};
    vector<float> alpha_values = {0.2,0.4,0.5,0.6,0.7, 0.8};

    ofstream summary("../Esperimenti/" + fileName + "/grid_search_summary_" + p_str + ".csv");
    summary << "k,phi,alpha,avg_size_over_time,avg_density_sampled,avg_time_ms\n";

    int total_combinations = k_values.size() * phi_values.size() * alpha_values.size();
    int current_comb = 0;
    
    // Intervall for the computation of the density
    int interval_10_percent = std::max(1, (int)(m * 0.1));

    cout << "Inizio Grid Search: " << total_combinations << " combinazioni totali." << endl;

    // Grid Search
    for (int k : k_values) {
        for (float phi : phi_values) {
            for (float alpha : alpha_values) {
                current_comb++;
                cout << "Test " << current_comb << "/" << total_combinations 
                     << " [k=" << k << ", phi=" << phi << ", alpha=" << alpha << "]" << endl;

                double total_time_ms = 0;
                
                double sum_of_sizes = 0;
                long long count_sizes = 0;
                
                double sum_of_densities = 0;
                long long count_densities = 0;

                for (int i = 0; i < number_of_experiments; i++) {
                    auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
                    int operations_done = 0;
                    
                    auto start_time = std::chrono::high_resolution_clock::now();
                    
                    for (const auto& op : all_experiments_ops[i]) {
                        operations_done++;

                        if (op.type == 'i') {
                            algo_cred.add_edge(op.u, op.v);  
                        } else if (op.type == 'd') {
                            algo_cred.remove_edge(op.u, op.v); 
                        }

                        sum_of_sizes += algo_cred.return_dim();
                        count_sizes++;

                        if (operations_done % interval_10_percent == 0) {
                            sum_of_densities += algo_cred.return_density();
                            count_densities++;
                        }
                    }

                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::milli> duration = end_time - start_time;
                    total_time_ms += duration.count();
                }

                // compute the average values
                double avg_time = total_time_ms / number_of_experiments;
                double avg_size_over_time = (count_sizes > 0) ? (sum_of_sizes / count_sizes) : 0;
                double avg_density_sampled = (count_densities > 0) ? (sum_of_densities / count_densities) : 0;

                summary << k << "," << phi << "," << alpha << "," 
                        << avg_size_over_time << "," << avg_density_sampled << "," << avg_time << "\n";
            }
        }
    }
    
    summary.close();
    cout << "Grid Search completata! File generato con le medie globali." << endl;
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

    //calcola gli speedup e inserisci le informazioni nei file
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


template <typename T>
void saveToCSV(const std::string& filename, const std::vector<std::vector<T>>& data) {
    std::ofstream outFile(filename);

    if (!outFile.is_open() || data.empty()) return;

    // 1. Trova l'esperimento più lungo (per sapere quante righe scrivere nel file)
    size_t max_time_steps = 0;
    for (const auto& exp : data) {
        if (exp.size() > max_time_steps) max_time_steps = exp.size();
    }

    // 2. Ciclo ESTERNO: Il Tempo (Le righe del tuo file Excel)
    for (size_t t = 0; t < max_time_steps; ++t) {
        
        // 3. Ciclo INTERNO: Gli Esperimenti (Le colonne del tuo file Excel)
        for (size_t exp = 0; exp < data.size(); ++exp) {
            
            // Controlla se l'esperimento 'exp' ha ancora dati al tempo 't'
            if (t < data[exp].size()) {
                outFile << data[exp][t];
            }
            // Altrimenti lascia vuoto (o metti 0 se preferisci)
            
            // Metti la virgola se non è l'ultima colonna
            if (exp < data.size() - 1) {
                outFile << ",";
            }
        }
        // Fine della riga temporale, andiamo a capo
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

void experimentDynamicOnly(string fileName, float gamma, float b, float prob, float perc = 0,string strategy="standard") {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    if(strategy=="standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    }
    else{
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str + ".txt");
    }
    
    if (!file.is_open()) {
        cout<<"error!"<<endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    vector<vector<int>> dynamicSolution;
    vector<vector<float>> dynamicDensity;

    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_dyn = DynamicNBSim(gamma, b, n, m);
        dynamicSolution.push_back({});
        dynamicDensity.push_back({});

        bool end = false;
        char type_of_op;
        int u, v, number_of_edges = 0;

        while (!end) {
            file >> type_of_op>>u>>v;
            if (type_of_op == 'i') {
                number_of_edges++;
                algo_dyn.add_edge(u, v); 
                dynamicSolution.back().push_back(algo_dyn.return_dim());
                if(number_of_edges%(int)(m*perc)==0){
                    dynamicDensity.back().push_back(algo_dyn.return_density());
                }
            } 
            else if (type_of_op == 'd') {
                algo_dyn.remove_edge(u, v);
                if(strategy=="mixed"){
                    dynamicSolution.back().push_back(algo_dyn.return_dim());
                    if(number_of_edges%(int)(m*perc)==0){
                        dynamicDensity.back().push_back(algo_dyn.return_density());
                    } 
                }
            } 
            else {
                end = true;
            }
        }
    }
    saveToCSV("../Esperimenti/"+fileName+"/density_credits_"+p_str+".csv",dynamicDensity);
    saveToCSV("../Esperimenti/" + fileName + "/results_dynamic_" + p_str + ".csv", dynamicSolution);
    file.close();
}

void experimentCreditsOnly(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,float perc = 0.1,string strategy="standard") {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    if(strategy=="standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    }
    else{
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str + ".txt");
    }

    if (!file.is_open()) {
        cout<<"error!"<<endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    vector<vector<int>> creditsSolution;
    vector<vector<float>> creditsDensity;


    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
        creditsSolution.push_back({});
        creditsDensity.push_back({});
        bool end = false;
        char type_of_op;
        int u, v, edge_inserted = 0;

        while (!end) {
            file >> type_of_op>>u>>v;
            edge_inserted++;
            if (type_of_op == 'i') {
                algo_cred.add_edge(u, v);  
                creditsSolution.back().push_back(algo_cred.return_dim());
                if(edge_inserted%(int)(m*perc)==0){
                    creditsDensity.back().push_back(algo_cred.return_density());
                }
            } 
            else if (type_of_op == 'd') {
                algo_cred.remove_edge(u, v); 
                if(strategy=="mixed"){
                    creditsSolution.back().push_back(algo_cred.return_dim());
                    if(edge_inserted%(int)(m*perc)==0){
                        creditsDensity.back().push_back(algo_cred.return_density());
                    }
                }
            } 
            else {
                end = true;
            }
        }
    }
    saveToCSV("../Esperimenti/"+fileName+"/density_credits_"+p_str+".csv",creditsDensity);
    saveToCSV("../Esperimenti/" + fileName + "/results_credits_" + p_str + ".csv", creditsSolution);
    file.close();
}

void experimentCreditsPerformance(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,string strategy="standard") {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    if(strategy=="standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    }
    else{
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str + ".txt");
    }
    if (!file.is_open()) {
        cout << "Errore: impossibile aprire il file di sequenza!" << endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    vector<vector<double>> executionTimes;

    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 
        bool end = false;
        char type_of_op;
        int u, v;

        

        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);
        
        auto start_time = std::chrono::high_resolution_clock::now();

        while (!end) {
            file >> type_of_op >> u >> v;
            
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

        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;

        executionTimes.back().push_back(duration.count());
    }

    saveToCSV("../Esperimenti/" + fileName + "/times_credits_" + p_str + ".csv", executionTimes);
    file.close();
}

void experimentDynamicPerformance(string fileName, float gamma, float b, float prob,string strategy = "standard") {
    string p_str = to_string(prob).substr(0, 3);
    std::ifstream file;
    if(strategy=="standard"){
        file.open("./Sequences/StandardSequences/" + fileName + "_" + p_str + ".txt");
    }
    else{
        file.open("./Sequences/MixedSequences/" + fileName + "_" + p_str + ".txt");
    }
    if (!file.is_open()) {
        cout << "Errore: impossibile aprire il file di sequenza!" << endl;
        return; 
    }

    int number_of_experiments, n, m;
    file >> number_of_experiments >> n >> m;

    vector<vector<double>> executionTimes;

    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 
        bool end = false;
        char type_of_op;
        int u, v;

        auto algo_dyn = DynamicNBSim(gamma, b, n, m); 
        
        auto start_time = std::chrono::high_resolution_clock::now();

        while (!end) {
            file >> type_of_op >> u >> v;
            
            if (type_of_op == 'i') {
                algo_dyn.add_edge(u, v);  
            } 
            else if (type_of_op == 'd') {
                algo_dyn.remove_edge(u, v);
            } 
            else {
                end = true;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;

        executionTimes.back().push_back(duration.count());
    }

    saveToCSV("../Esperimenti/" + fileName + "/times_dynamic_" + p_str + ".csv", executionTimes);
    file.close();
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

    /*
    experimentDynamicOnly(fileName,gamma,b,p,0.1,"mixed");
    experimentCreditsOnly(fileName,phi,alpha,gamma,b,k,p,0.1,"mixed");
    experimentCreditsPerformance(fileName,phi,alpha,gamma,b,k,p,"mixed");
    experimentDynamicPerformance(fileName,gamma,b,p,"mixed");
    */
    gridSearchCreditsAlgorithm(fileName, gamma,b,p);

    cout << "Tutti gli esperimenti completati con successo!" << endl;
    return 0;
}

