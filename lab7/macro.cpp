/*
KADD lab7_v1
Kacper Skelnik
*/
#include <TMath.h>

using namespace std;
TRandom *gen = new TRandom();

double kolo(double x, double y){
    return x*x + y*y;
}

bool w_kole(double x, double y, double r){
    if (kolo(x,y) <= r){ return true; }
    else { return false; }    
}

double pi(int N, double r){
    double x[N];
    double y[N];
    for(int i = 0; i < N; i++){
        x[i] = gen->Rndm();
        y[i] = gen->Rndm();
    }

    vector<double> x_wkole, x_poza, y_wkole, y_poza;
    for(int i = 0; i < N; i++){
        if(w_kole(x[i], y[i], r) == true){
            x_wkole.push_back(x[i]);
            y_wkole.push_back(y[i]);
        }
        else{
            x_poza.push_back(x[i]);
            y_poza.push_back(y[i]);
        }
    }

    //cout << "pi = " << 4*(double)x_wkole.size()/N << endl;
    cout << "Wzgledna dokladnosc calki = " << 1/TMath::Sqrt(N) << endl;

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1000, 600);
	c->SetWindowSize(1100, 700);

    TMultiGraph *g_graph = new TMultiGraph();

    TGraph *g_kolo = new TGraph (x_wkole.size(), &x_wkole[0], &y_wkole[0]);
    g_kolo->SetMarkerColor(2);
    g_kolo->SetMarkerStyle(9);
    g_kolo->SetMarkerSize(1);

    TGraph *g_poza = new TGraph (x_poza.size(), &x_poza[0], &y_poza[0]);
    g_poza->SetMarkerColor(1);
    g_poza->SetMarkerStyle(9);
    g_poza->SetMarkerSize(1);

    g_graph->Add(g_kolo);
    g_graph->Add(g_poza);

    g_graph->Draw("AP");


    return 4*(double)x_wkole.size()/N;
}

int macro(){ 
    int N = 500000;
    double r = 1;
    
    double Pi = pi(N,r);
    cout << "pi = " << Pi << endl;

return 0;
}