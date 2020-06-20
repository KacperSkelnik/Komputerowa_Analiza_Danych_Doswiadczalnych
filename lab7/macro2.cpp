/*
KADD lab7_v2
Kacper Skelnik
*/
#include <TMath.h>

using namespace std;
TRandom *gen = new TRandom();

Double_t myFunc(double *x, double *par){ 
    return (1/(par[1]*TMath::Sqrt(2* TMath::Pi())))*TMath::Exp(-((x[0]-par[0])/par[1])*((x[0]-par[0])/par[1])/(2));
}   

double losujVonNeumann(TF1* g, double min, double max){
    double x;
    double y;
    
    bool run = true;
    while(run == true){
        x = gen->Uniform(min, max);
        y = gen->Uniform(0, g->GetMaximum(min, max));

        if( y <= g->Eval(x) ){
            run = false;
            return x;
        }
    }
    return x;
}

double wydajnoscVonNeumann(TF1* g, double min, double max, int N){
    double x[N];
    double y[N];

    for(int i = 0; i < N; i++){
        x[i] = gen->Uniform(min, max);
        y[i] = gen->Uniform(0, g->GetMaximum(min, max));
    }

    vector<double> x_w, x_poza, y_w, y_poza;
    for(int i = 0; i < N; i++){

        if( y[i] <= g->Eval(x[i]) ){    
            x_w.push_back(x[i]);
            y_w.push_back(y[i]);
        }
        else{
            x_poza.push_back(x[i]);
            y_poza.push_back(y[i]);
        }
    }
    return (double)x_w.size()/N;
}

double calkaVonNeumann(TF1* g, double min, double max, int N){
    double wydajnosc = wydajnoscVonNeumann(g, min, max, N);
    
    return wydajnosc * (max - min) * (g->GetMaximum(min, max));
}


int macro2(){ 
    int N = 20000;

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1000, 600);
	c->SetWindowSize(1100, 700);
	
	TF1 *fun = new TF1("f1",myFunc, 0,10,2);
	fun->SetParameter(0,5);
	fun->SetParameter(1,2);

    double x = losujVonNeumann(fun, 0, 10);

    cout << "wylosowana jedna liczba pseudolosowa = " << x << endl;
    cout << "wydajnosc metody = " << wydajnoscVonNeumann(fun, 0, 10, N) << endl;
    cout << "\n";
    cout << "całka monte carlo = " << calkaVonNeumann(fun, 2, 4, N) << endl;
    cout << "Wzgledna dokladnosc calki = " << 1/TMath::Sqrt(N) << endl;
    cout << "wydajnosc metody na przedziale (2,4) = " << wydajnoscVonNeumann(fun, 2, 4, N) << endl;
    cout << "\n";
    cout << "całka wbudowana = " << fun->Integral(2,4) << endl;


    c->cd(1);

	TH1D* neumannHist = new TH1D("neumannHist", "density", 70, 0, 10);
	for (int i = 0; i < N; i++) {
		neumannHist->Fill(losujVonNeumann(fun, 0, 10));
	}

    neumannHist->Scale(1/neumannHist->Integral(), "width");
    neumannHist->SetLineColor(2);
	neumannHist->SetStats(kFALSE);
	neumannHist->GetXaxis()->SetTitle("x");
	neumannHist->GetYaxis()->SetTitle("y");
	neumannHist->Draw("HIST");

    TH1D* random = new TH1D("random", "density", 70, 0, 10);
	for (int i = 0; i < N; i++) {
		random->Fill(fun->GetRandom(0, 10));
	}

    random->Scale(1/random->Integral(), "width");
    random->SetLineColor(4);
	random->SetStats(kFALSE);
	random->Draw("HIST,SAME");

    fun->SetLineColor(8);
    fun->Draw("SAME");

return 0;
}