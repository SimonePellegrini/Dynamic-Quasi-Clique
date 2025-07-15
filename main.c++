#include <iostream>
#include <fstream>
#include "DynamicNBSim.c++"
#include <random>
#include <chrono>
#include "creditsAlgo.cpp"
#include <sstream>  
#include <cstdlib>  
#include <algorithm> 
#include "FastNBSim.c++"

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
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
}

void gridSearch(string fileName){
    int eta_param[] = {500,1250,2000};
    int phi_param[] = {3,9,27};
    int delta_param[] = {75,150,225,300};
    ifstream tempo_baseline_naive("../Plotting/data/"+fileName+"/naiveAvgTime.csv");
    ifstream qc_baseline_naive("../Plotting/data/"+fileName+"/naiveSolution.csv");
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
    
    tempo_baseline_naive>>tempo_naive;
    ofstream confronto_tempi("../Plotting/data/"+fileName+"/confrontoTempi.txt");
    ofstream confronto_dimensioni("../Plotting/data/"+fileName+"/confrontoDimensioni.txt");
    for(auto eta: eta_param){
        for(auto phi: phi_param){
            for(auto delta: delta_param){
                confronto_tempi << "( eta = "<<eta<<" / phi = "<<phi<<" / delta = " << delta<<" )"<<endl;
                confronto_dimensioni << "( eta = "<<eta<<" / phi = "<<phi<<" / delta = " << delta<<" )"<<endl;
                double media_tempi = 0;
                double media_dimensioni = 0;
                for(int j = 0; j<5; j++){
                    int count = 0;
                    int i = 0;
                    auto c_a = creditsAlgo(n,m,eta,phi,delta,0.9,0.6,512);
                    auto start = chrono::high_resolution_clock::now();
                    for(auto x: graph){
                        c_a.add_edge(x.first,x.second);
                        if(c_a.return_dim()>quasi_clique_sizes[i]*0.8){
                            count++;
                        }
                        i++;
                    }
                    auto end = chrono::high_resolution_clock::now();
                    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
                    media_tempi+=(double)tempo_naive/((double)duration.count()/graph.size());   
                    media_dimensioni+=(double)count/graph.size();
                }
                confronto_tempi<<(double)media_tempi/5<<endl;  
                confronto_dimensioni<<(double)media_dimensioni/5<<endl;
            }
        }
    }
    confronto_tempi.close();
    confronto_dimensioni.close();

    ofstream confronto_densità("../Plotting/data/"+fileName+"/confrontoDensità.txt");
    for(auto eta: eta_param){
        for(auto phi: phi_param){
            for(auto delta: delta_param){
                confronto_densità << "( eta = "<<eta<<" / phi = "<<phi<<" / delta = " << delta<<" )"<<endl;
                auto c_a = creditsAlgo(n,m,eta,phi,delta,0.9,0.6,512);
                int count = 0;
                int s_rate = 0.1*m;
                int i = 0;
                for(auto x: graph){
                    c_a.add_edge(x.first,x.second);
                    if(i%s_rate == 0){
                        double d = c_a.return_density();
                        if(d!=0){
                            i++;
                        }
                        if(d>0.83){
                            count++;
                        }
                    }
                }
                confronto_densità<<(double)count/i<<endl;
            }
        }
    }
    confronto_densità.close();
    tempo_baseline_naive.close();
    tempo_baseline_naive.close();
    qc_baseline_naive.close();
}

void getSoluzioneDynamicNBSim(string fileName){
    std::ofstream solution("../Plotting/data/"+fileName+"/naiveSolution.csv");
    std::ofstream tempi_naive("../Plotting/data/"+fileName+"/naiveAvgTime.csv");
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

void getSoluzioneCreditsAlgo(string fileName, int n, int m,int eta_default,int phi_default,int delta_default,double gamma,double b, int k){
    std::ofstream file_varianza("../Plotting/data/"+fileName+"/solVarianza"+to_string(k)+".csv"); 
    std::ofstream file_densità("../Plotting/data/"+fileName+"/solDensità"+to_string(k)+".csv"); 
    int s_rate = 0.1*m;
    unordered_map<int,vector<double>> densità;
    unordered_map<int,vector<double>> soluzioni;
    for(int j=0; j<5; j++){
        auto c_a = creditsAlgo(n,m,eta_default,phi_default,delta_default,gamma,b,k);
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
    for(int i=0; i<5; i++){
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

void confrontoAlgo1Baseline(string fileName,int num_nodes, int num_edges,int s_rate,int k){
    auto base =  FastNBSim(num_nodes,0.6,0.9,k);
    int i=0;
    double tot_time=0;
    ofstream speed_up("../Plotting/data/"+fileName+"/speedup.txt");
    ifstream tempo_baseline_naive("../Plotting/data/"+fileName+"/naiveAvgTime.csv");
    double tempo_naive;
    tempo_baseline_naive>>tempo_naive;
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
    cout<<((double)tot_time/graph.size())/tempo_naive;
    speed_up<<((double)tot_time/graph.size())/tempo_naive;
    speed_up.close();
    tempo_baseline_naive.close();
}

void effettoK(string fileName,int n, int m,int eta_default,int phi_default,int delta_default,double gamma,double b){
    int k_list[] = {4,8,16,32};
    ofstream k_time("../Plotting/data/"+fileName+"/kSpeedUp.csv");
    ifstream tempo_baseline_naive("../Plotting/data/"+fileName+"/naiveAvgTime.csv");
    double tempo_naive;
    tempo_baseline_naive>>tempo_naive;
    k_time<<"k,tempo medio"<<endl;

    for(auto k: k_list){
        double tempo_medio_totale=0;
        for(int j=0; j<5; j++){
            auto credAlgo =  creditsAlgo(n,m,eta_default,phi_default,delta_default,gamma,b,k);
            auto start = chrono::high_resolution_clock::now();
            for(auto x: graph){   
                credAlgo.add_edge(x.first,x.second);
                credAlgo.return_dim();
            }
            auto end = chrono::high_resolution_clock::now();
            tempo_medio_totale+=((double)chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/graph.size());
        }
        double tempo_medio = tempo_medio_totale/5;
        k_time<<k<<","<<tempo_naive/tempo_medio<<endl;
    }
    
    int k_list_grafici[] = {4,8,16,32};    
    for(auto k:k_list_grafici){
        getSoluzioneCreditsAlgo(fileName,n,m,eta_default,phi_default,delta_default,gamma,b,k);
    }

    tempo_baseline_naive.close();
    k_time.close();
}


void calcolaSpeedUpMedio(string fileName,int n, int m,int eta_default,int phi_default,int delta_default,double gamma,double b, int k){
    
    int i=0;
    double tempo_medio_totale=0;
    ofstream speed_up("../Plotting/data/"+fileName+"/Credits_speedup.txt");
    ofstream creditsAvgTime("../Plotting/data/"+fileName+"/CreditsAvgTime.txt");
    ifstream tempo_baseline_naive("../Plotting/data/"+fileName+"/naiveAvgTime.csv");
    double tempo_naive;
    tempo_baseline_naive>>tempo_naive;
    for(int j=0; j<5; j++){
        auto credAlgo =  creditsAlgo(n,m,eta_default,phi_default,delta_default,gamma,b,k);
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

void calcolaSoluzioniFastNBSIM(string fileName,int n,int m,int k){
    
    unordered_map<int,vector<double>> densità;
    unordered_map<int,vector<int>> dimensione;
    ofstream confronto_dimensioni_apx("../Plotting/data/"+fileName+"/confrontoDimensioniApx.csv");
    ofstream confronto_densità_apx("../Plotting/data/"+fileName+"/confrontoDensitàApx.csv");
    confronto_dimensioni_apx<<"numero di archi,dimensione"<<endl;
    confronto_densità_apx<<"numero di archi,densità"<<endl;
    int s_rate = m*0.1;
    for(int i=0; i<5; i++){
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

int main(int argc, const char * argv[]){
    /*prende dalla riga di comando i parametri in input*/
    int eta_default = atoi(argv[1]);
    int phi_default = atoi(argv[2]);
    int delta_default = atoi(argv[3]);
    string fileName = argv[5];
    double gamma = 0.9;
    double b = 0.6;
    int k = atoi(argv[4]);

    /*inizializza il dataset da leggere*/
    read_graph("../Datasets/"+fileName+".txt");
    /* calcola la soluzione ottenuta utilizzando DynamicNBSim e il tempo medio di inseritmento*/
    getSoluzioneDynamicNBSim(fileName);
    /*grid search utilizzata per la selezione dei parametri*/
    gridSearch(fileName);
    /*Valuta la densità media e le dimensioni medie delle clique estratte usando CreditsAlgo*/
    getSoluzioneCreditsAlgo(fileName,n,m,eta_default,phi_default,delta_default,gamma,b,k);
    /*Permette di confrontare il running time di DynamicNBSim e NBSim (se k=0) o FastNBSim (se k>0)*/
    int frequenza_query = 100;
    confrontoAlgo1Baseline(fileName,n,m,100,k);
    /*permette di calcolare lo speed up ottenuto utilizzando CreditsAlgo invece di DynamicNBSim*/
    calcolaSpeedUpMedio(fileName,n,m,eta_default,phi_default,delta_default,gamma,b,k);
    /*Permette di valutare gli effetti sulle prestazioni di CreditsAlgo del parametro k*/
    effettoK(fileName,n,m,eta_default,phi_default,delta_default,gamma,b);
    /*Calcola le soluzioni utilizzando l'algoritmo FastNBSim in batch ad intervalli regolari.*/
    calcolaSoluzioniFastNBSIM(fileName,n,m,k); 
}

