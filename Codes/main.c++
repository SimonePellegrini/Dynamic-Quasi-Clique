#include <iostream>
#include <fstream>
#include "DynamicNBSim.c++"
#include <random>
#include <chrono>
#include <sstream>  
#include <cstdlib>  
#include <algorithm> 
#include "CreditsAlgorithm.cpp"
#include "CreditsAlgorithmOnlyIns.cpp"
#include "math.h"
#include <unordered_map>
#include <unordered_map>

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
        return {};
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
    cout<<"file read!"<<endl;
    return all_experiments_ops;
}

template <typename T>
void saveToCSV(const std::string& filename, const std::vector<std::vector<T>>& data) {
    std::ofstream outFile(filename);

    if (!outFile.is_open() || data.empty()){
        cout<<"Error:the file has not been saved correctly!"<<endl;
        return;
    }
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

void countCreditsAndGammeDegree(string fileName, float phi, float alpha, float gamma, float b, int k, float prob, string strategy="standard", string type="erdos") {
    // Calculate the probability string for saving
    string p_str = to_string((int)(prob*100));

    // Load the dataset using the helper function
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName, prob, strategy, type);

    // If loading fails, all_experiments_ops will be empty
    if (all_experiments_ops.empty()) {
        cout << "Error: impossible to load the dataset for credits!" << endl;
        return;
    }

    vector<vector<int>> file_credits;
    vector<vector<int>> file_gamma;

    // number_of_experiments, n, m are global variables initialized by loadDataset
    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_cred = creditsAlgorithm(n, m, phi, alpha, gamma, b, k);

        // Iterate directly over the pre-loaded operations in memory
        for (const auto& op : all_experiments_ops[i]) {
            if (op.type == 'i') {
                algo_cred.add_edge(op.u, op.v);  
            } 
            else if (op.type == 'd') {
                algo_cred.remove_edge(op.u, op.v); 
            }
        }
        
        file_credits.push_back(algo_cred.getCredits());
        file_gamma.push_back(algo_cred.getGammaDegree());
    }
    
    saveToCSV("../Esperimenti/"+fileName+"/number_of_credits_"+p_str+".csv", file_credits);
    saveToCSV("../Esperimenti/"+fileName+"/gamma_degree_"+p_str+".csv", file_gamma);
}

void experimentDynamicOnly(string fileName, float gamma, int b, float prob, float perc = 0.1,string strategy="standard",string type="erdos") {
    string p_str = to_string((int)(prob*100));

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
    string p_str = to_string((int)(prob*100));
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
        saveToCSV("../Esperimenti/"+fileName+"/standard/density_credits_"+p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/standard/results_credits_" + p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv" , creditsSolution);
    }
    else{
        saveToCSV("../Esperimenti/"+fileName+"/mixed/density_credits_"+type+"_"+p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/mixed/results_credits_"+type+"_"+p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv", creditsSolution);
    }    
}


void experimentCreditsOnlyOnlyIns(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,float perc = 0.1,string strategy="standard",string type="erdos") {
    string p_str = to_string((int)(prob*100));
    //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    vector<vector<int>> creditsSolution;
    vector<vector<float>> creditsDensity;


    for (int i = 0; i < number_of_experiments; i++) {
        auto algo_cred = creditsAlgorithmOnlyIns(n, m, phi, alpha, gamma, b, k);
        creditsSolution.push_back({});
        creditsDensity.push_back({});
        bool end = false;
        int number_of_edges = 0;

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_cred.add_edge(op.u, op.v);
            } 
            number_of_edges++;
            creditsSolution.back().push_back(algo_cred.return_dim());
            if(number_of_edges%(int)(m*perc)==0){
                creditsDensity.back().push_back(algo_cred.return_density());
            }
        }
    }

    
      
    if(strategy == "standard"){
        saveToCSV("../Esperimenti/"+fileName+"/standard/density_credits_"+p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/standard/results_credits_" + p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv" , creditsSolution);
    }
    else{
        saveToCSV("../Esperimenti/"+fileName+"/mixed/density_credits_"+type+"_"+p_str+"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+".csv",creditsDensity);
        saveToCSV("../Esperimenti/" + fileName + "/mixed/results_credits_"+type+"_"+p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv", creditsSolution);
    }    
}

void experimentCreditsPerformance(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,string strategy="standard",string type="erdos",int q_freq=10) {
    string p_str = to_string((int)(prob*100));
    
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
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;

        executionTimes.back().push_back(duration.count());
    }

    if(strategy == "standard"){
        saveToCSV("../Esperimenti/" + fileName + "/standard/times_credits_" + p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k" +to_string(k) + ".csv", executionTimes);
    }
    else{
        saveToCSV("../Esperimenti/" + fileName + "/mixed/times_credits_"+type+"_" + p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv",executionTimes);
    } 
}

void experimentCreditsPerformanceOnlyIns(string fileName, float phi, float alpha, float gamma, float b, int k, float prob,string strategy="standard",string type="erdos",int q_freq=10) {
    string p_str = to_string((int)(prob*100));
    
     //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    

    vector<vector<double>> executionTimes;

    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 
        auto algo_cred = creditsAlgorithmOnlyIns(n, m, phi, alpha, gamma, b, k);
        int n_op=0;
        auto start_time = std::chrono::high_resolution_clock::now();

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_cred.add_edge(op.u, op.v);  
            }
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> duration = end_time - start_time;

        executionTimes.back().push_back(duration.count());
    }

    if(strategy == "standard"){
        saveToCSV("../Esperimenti/" + fileName + "/standard/times_credits_" + p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k" +to_string(k) + ".csv", executionTimes);
    }
    else{
        saveToCSV("../Esperimenti/" + fileName + "/mixed/times_credits_"+type+"_" + p_str +"_alpha"+to_string(alpha)+"_phi"+to_string(phi)+"_k"+to_string(k)+ ".csv",executionTimes);
    } 
}

void experimentDynamicPerformance(string fileName, float gamma, float b, float prob,string strategy = "standard",string type="erdos",int q_freq=10, bool only_insert = false) {
    string p_str = to_string((int)(prob*100));
    std::ifstream file;

    
    //load the dataset
    vector<vector<GraphOp>> all_experiments_ops = loadDataset(fileName,prob,strategy,type);
    vector<vector<double>> executionTimes;
    
    for (int i = 0; i < number_of_experiments; i++) {
        executionTimes.push_back({}); 
        auto algo_dyn = DynamicNBSim(gamma, b, n, m); 

        auto start_time = std::chrono::high_resolution_clock::now();

        for (const auto& op : all_experiments_ops[i]) {    
            if (op.type == 'i') {
                algo_dyn.add_edge(op.u, op.v);  
            } else if (op.type == 'd') {
                algo_dyn.remove_edge(op.u, op.v); 
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


void gridSearchK(string fileName, float phi, float alpha, float gamma, float b, float prob, string strategy = "standard", string type = "erdos", int q_freq = 10,bool run_dyn=true) {
    if(run_dyn){
        cout << "Running dynamic baseline..." << endl;
        experimentDynamicOnly(fileName, gamma, b, prob, 0.1, strategy, type);
        experimentDynamicPerformance(fileName, gamma, b, prob, strategy, type, q_freq);
    }
    vector<int> k_values = {16,32,64,128};
    vector<double> alpha_values = {0.1,0.3,0.6};
    vector<double> phi_values = {0.8};

    for (int k : k_values) {
        for(double alpha : alpha_values){
            for(double phi : phi_values){
                if(prob!=0){
                    experimentCreditsOnly(fileName, phi, alpha, gamma, b, k, prob, 0.1, strategy, type);
                    experimentCreditsPerformance(fileName, phi, alpha, gamma, b, k, prob, strategy, type, q_freq);
                }
                else{
                    experimentCreditsOnlyOnlyIns(fileName, phi, alpha, gamma, b, k, prob, 0.1, strategy, type);
                    experimentCreditsPerformanceOnlyIns(fileName, phi, alpha, gamma, b, k, prob, strategy, type, q_freq);
                }
            }
        }    
    }
    
    cout << "Grid search completed." << endl;
}




int main(int argc, const char * argv[]){

  
    if (argc < 8) {
        cerr << "Use: " << argv[0] << " <phi> <alpha> <k> <fileName> <p> <gamma> <b>" << endl;
        return 1;
    }
    float phi = atof(argv[1]);
    float alpha = atof(argv[2]);
    string fileName = argv[4];
    double gamma = atof(argv[6]);
    double b = atof(argv[7]);
    int k = atoi(argv[3]);
    double p = atof(argv[5]);
    
    read_graph("../Datasets/"+fileName+".txt");

    

    cout<<fileName<<endl;
    //cout<<"standard"<<endl;
    //countCreditsAndGammeDegree(fileName,phi,alpha,gamma,b,k,p);
    gridSearchK(fileName,phi,alpha,gamma,b,p,"standard","erdos",10,false);
    cout<<"erdos"<<endl;
    //gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","erdos");
    cout<<"rich"<<endl;
    gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","rich",10,false);
    cout<<"power"<<endl;
    //gridSearchK(fileName,phi,alpha,gamma,b,p,"mixed","power");
    
    //gridSearchCreditsAlgorithm(fileName, gamma,b,p);

    cout << "Experiments completed!" << endl;
    return 0;
}

