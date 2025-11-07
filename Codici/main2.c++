#include <iostream>
#include <fstream>
#include "DynamicNBSim.c++"
#include <random>
#include <chrono>
#include <sstream>  
#include <cstdlib>  
#include <algorithm> 
#include "FastNBSim.c++"
#include "newCreditsAlgo.cpp"

using namespace std;
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



//esegue la grid search per trovare i migliori parametri

void gridSearch(string fileName){
    //parametri utilizzati nella grid search
    float alpha_param[] = {0,1,1.03,1.06,1.09,1.12};
    float phi_param[] = {0,0.3,0.5,0.8,1,1.2,1.5,2};
    
    ifstream tempo_baseline_naive("../Esperimenti/"+fileName+"/naiveAvgTime.csv");
    ifstream qc_baseline_naive("../Esperimenti/"+fileName+"/naiveSolution.csv");
    string riga;
    vector<int> quasi_clique_sizes;

    getline(qc_baseline_naive, riga);

    while (getline(qc_baseline_naive, riga)) {
        stringstream ss(riga);
        string archi, size;
    
        getline(ss, archi, ',');
        getline(ss, size, ',');
    
        quasi_clique_sizes.push_back(stoi(size));
    }
    double tempo_naive;
    int n_iter = 100;
    tempo_baseline_naive>>tempo_naive;
    ofstream confronto_tempi("../Esperimenti/"+fileName+"/confrontoTempi.txt");
    ofstream confronto_dimensioni("../Esperimenti/"+fileName+"/confrontoDimensioni.txt");
        for(auto phi: phi_param){
            for(auto alpha: alpha_param){
                confronto_tempi << "(phi = "<<phi<<" / alpha = " << alpha <<" )"<<endl;
                confronto_dimensioni << "(phi = "<<phi<<" / alpha = " << alpha<<" )"<<endl;
                double media_tempi = 0;
                double media_dimensioni = 0;
                for(int j = 0; j<n_iter; j++){
                    int count = 0;
                    int i = 0;
                    read_graph("../Datasets/"+fileName+".txt");
                    auto c_a = CredAlgoImproved(n,m,phi,alpha,0.9,0.6,512);
                    auto start = chrono::high_resolution_clock::now();
                    for(auto x: graph){
                        c_a.add_edge(x.first,x.second);
                        if(c_a.return_dim()>(float)quasi_clique_sizes[i]*0.8){
                            count++;
                        }
                        i++;
                    }
                    auto end = chrono::high_resolution_clock::now();
                    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
                    media_tempi+=(double)tempo_naive/((double)duration.count()/graph.size());   
                    media_dimensioni+=(double)count/graph.size();
                }
                confronto_tempi<<(double)media_tempi/n_iter<<endl;  
                confronto_dimensioni<<(double)media_dimensioni/100<<endl;
            }
        }
    confronto_tempi.close();
    confronto_dimensioni.close();

    
    tempo_baseline_naive.close();
    tempo_baseline_naive.close();
    qc_baseline_naive.close();
}

//calcola le dimensioni delle quasi clique estratte dall'algoritmo dinamico Naive e i tempi di inserimento medi
void getSoluzioneDynamicNBSim(string fileName){
    std::ofstream solution("../Esperimenti/"+fileName+"/naiveSolution.csv");
    std::ofstream tempi_naive("../Esperimenti/"+fileName+"/naiveAvgTime.csv");
    auto n_a = DynamicNBSim(0.9,0.6,n,m);
    int i = 0;
    solution<<"numero di archi,size quasi-clique-naive"<<endl;
    auto start = chrono::high_resolution_clock::now();
    for(auto x: graph){
        n_a.Insert(x.first,x.second);
        solution<<i<<","<<n_a.getBestCliqueSize()<<endl;  
        i++;
    }
    auto end = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    tempi_naive<<(double)duration.count()/graph.size();
    tempi_naive.close();
    solution.close();
}

//confronta i tempi di inserimento medi ottenuti utilizzando
//l'algoritmo dei cinesi "in batch" e l'algoritmo dinamico naive

void confrontoAlgo1Baseline(string fileName,int num_nodes, int num_edges,int s_rate,int k){
    
    int i=0;
    double tot_time=0;
    int n_iter = 10;
    ofstream speed_up("../Esperimenti/"+fileName+"/speedup.txt");
    ofstream tempoNBSim("../Esperimenti/"+fileName+"/tempoNBSim.txt");
    ifstream tempo_baseline_naive("../Esperimenti/"+fileName+"/naiveAvgTime.csv");
    double tempo_naive;
    tempo_baseline_naive>>tempo_naive;
    for(int i=0; i<n_iter; i++){
        auto base =  FastNBSim(num_nodes,0.9,0.6,k);
        for(auto x: graph){
            auto start = chrono::high_resolution_clock::now();
            
            base.add_edge(x.first,x.second);
            if(i++%s_rate == 0){
                base.compute_result(false);
            }
            auto end = chrono::high_resolution_clock::now();
            base.reset_computation();
            tot_time+=chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
    }
    tot_time = tot_time/n_iter;
    tempoNBSim<<tot_time/graph.size();
    speed_up<<((double)tot_time/graph.size())/tempo_naive;
    speed_up.close();
    tempo_baseline_naive.close();
    tempoNBSim.close();
}

//Calcola le soluzioni (dimensioni dell quasi clique e densità) ottenute 
//utilizzando l'algoritmo dei cinesi in "batch"

void calcolaSoluzioniFastNBSIM(string fileName,int n,int m,int k){
    
    unordered_map<int,vector<double>> densità;
    unordered_map<int,vector<int>> dimensione;
    ofstream confronto_dimensioni_apx("../Esperimenti/"+fileName+"/confrontoDimensioniApx.csv");
    ofstream confronto_densità_apx("../Esperimenti/"+fileName+"/confrontoDensitàApx.csv");
    confronto_dimensioni_apx<<"numero di archi,dimensione"<<endl;
    confronto_densità_apx<<"numero di archi,densità"<<endl;
    int s_rate = m*0.1;
    int n_iter = 5;
    for(int i=0; i<n_iter; i++){
        int j = 0;
        auto base =  FastNBSim(n,0.9,0.6,k);
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

//Calcola le dimensioni delle quasi-clique estratte da ogni nodo del grafo utilizzando l'algoritmo dei cinesi 
//sul grafo finale

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


//calcola il numero di crediti associati a ogni nodo alla fine dell'inserimento di ogni arco

void credsAlgoSolutions(string fileName, int n, int m,float phi_default,float alpha_default,double gamma,double b, int k){
    std::ofstream file_varianza("../Esperimenti/"+fileName+"/solVarianza"+to_string(k)+".csv"); 
    std::ofstream file_densità("../Esperimenti/"+fileName+"/solDensità"+to_string(k)+".csv"); 
    int s_rate = 0.1*m;
    unordered_map<int,vector<double>> densità;
    unordered_map<int,vector<double>> soluzioni;
    int num_iter = 5;
    for(int j=0; j<num_iter; j++){
        read_graph("../Datasets/"+fileName+".txt");
        auto c_a = CredAlgoImproved(n,m,phi_default,alpha_default,gamma,b,k);
        int i = 0;
        for(auto x: graph){
            c_a.add_edge(x.first,x.second);
            soluzioni[i].push_back(c_a.return_dim());
            if(i%s_rate==0 || i==m-1){
                densità[i].push_back(c_a.return_density());
            }
            i++;
        }
    }

    file_densità<<"numero archi,densità media"<<endl;
    for(int i = 0; i<m; i+=s_rate){
        file_densità<<i<<",";
        double densità_media = 0;
        int n_c = 0;
        for(auto x: densità[i]){
            if(x!=0){
                densità_media+=x;
                n_c++;
            }
        }
        file_densità<<(double)densità_media/n_c<<endl;
    }
    double densità_media = 0;
    int n_c = 0;
    for(auto x: densità[m-1]){
        if(x!=0){
            densità_media+=x;
            n_c++;
        }
    }
    file_densità<<m-1<<","<<(double)densità_media/n_c<<endl;

    file_varianza<<"numero di archi";
    for(int i=0; i<num_iter; i++){
        file_varianza<<",t"<<to_string(i);
    }
    file_varianza<<endl;

    for(int i=0; i<m;i++){
        file_varianza<<to_string(i);
        for(auto x:soluzioni[i]){
            file_varianza<<","<<x;
        }
        file_varianza<<endl;
    }
    file_varianza.close();
    file_densità.close();
}

//calcola lo speed up medio nei tempi di inserimento rispetto all'algoritmo dinamico Naive
void calcolaSpeedUpMedio(string fileName,int n, int m, float phi_default,float alpha_default,double gamma,double b, int k){
    
    int i=0;
    double tempo_medio_totale=0;
    ofstream speed_up("../Esperimenti/"+fileName+"/Credits_speedup.txt");
    ofstream creditsAvgTime("../Esperimenti/"+fileName+"/CreditsAvgTime.txt");
    ifstream tempo_baseline_naive("../Esperimenti/"+fileName+"/naiveAvgTime.csv");
    double tempo_naive;
    tempo_baseline_naive>>tempo_naive;
    for(int j=0; j<5; j++){
        read_graph("../Datasets/"+fileName+".txt");
        auto credAlgo =  CredAlgoImproved(n,m,phi_default,alpha_default,gamma,b,k);
        auto start = chrono::high_resolution_clock::now();
        for(auto x: graph){   
            credAlgo.add_edge(x.first,x.second);
            credAlgo.return_dim();
        }
        auto end = chrono::high_resolution_clock::now();
        tempo_medio_totale+=((double)chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/graph.size());
    }
    double tempo_medio = tempo_medio_totale/5;
    creditsAvgTime<<tempo_medio;
    speed_up<<(double)tempo_naive/tempo_medio;
    tempo_baseline_naive.close();
    creditsAvgTime.close();
    speed_up.close();
}

int main(int argc, const char * argv[]){
    /*prende dalla riga di comando i parametri in input*/
  
    float phi_default = atof(argv[1]);
    float alpha_default = atof(argv[2]);
    string fileName = argv[4];
    double gamma = 0.9;
    double b = 0.6;
    int k = atoi(argv[3]);
    /* Inizializza il dataset da leggere */
    read_graph("../Datasets/"+fileName+".txt");

    /* 
     * Calcola le dimensioni delle quasi clique estratte 
     * dall'algoritmo dinamico Naive e i tempi di inserimento medi
     */
    getSoluzioneDynamicNBSim(fileName);

    /* 
     * Esegue la grid search per trovare i migliori parametri
     */
    //gridSearch(fileName);

    /* 
     * Confronta i tempi di inserimento medi ottenuti utilizzando 
     * l'algoritmo dei cinesi "in batch" e l'algoritmo dinamico naive
     */
    confrontoAlgo1Baseline(fileName, n, m, 10, 8);

    /* 
     * Calcola le soluzioni (dimensioni e densità) ottenute 
     * utilizzando l'algoritmo dei cinesi in "batch"
     */
    //calcolaSoluzioniFastNBSIM(fileName, n, m, 0);

    /* 
     * Calcola le dimensioni delle quasi-clique estratte da ogni nodo del grafo 
     * utilizzando l'algoritmo dei cinesi sul grafo finale
     */
    //getQsSize(fileName, gamma, b);

    /* 
     * Calcola il numero di crediti associati a ogni nodo alla fine 
     * dell'inserimento di ogni arco
     */
    //credsAlgoSolutions(fileName, n, m, phi_default, alpha_default, gamma, b, k);

    /* 
     * Calcola lo speed up medio nei tempi di inserimento 
     * rispetto all'algoritmo dinamico Naive
     */
    calcolaSpeedUpMedio(fileName, n, m, phi_default, alpha_default, gamma, b, k);

    cout << "Tutti gli esperimenti completati con successo!" << endl;
    return 0;
}

