/*
KADD lab4
Kacper Skelnik
*/

#include <iostream>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TF12.h>

using namespace std;

TF2 *fun;
Double_t c1;
double x,y;

Double_t myFunc(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return TMath::Cos(x[0]/2)*TMath::Sin(x[1]);
    }
    else return 0;
}

Double_t my_normed_Func(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return par[2]*TMath::Cos(x[0]/2)*TMath::Sin(x[1]);
    }
    else return 0;
}

void fill_hist(int N, TH2D *histogram){
    for(int i=0; i<N; i++){
	fun->GetRandom2(x,y);
	histogram->Fill(x,y);
	}
}

void normalize_hist(TH2D *histogram){
    Double_t scale = 1/histogram->Integral(0,TMath::Pi(), 0,TMath::Pi());
    histogram->Scale(scale);
}

int macro(){    

    auto canvas = new TCanvas("c","canvas");
	canvas->SetCanvasSize(600, 600);
	canvas->SetWindowSize(700, 700);
	canvas->Divide(1,1);
	
	TF2 *raw_fun = new TF2("f1",myFunc, 0,TMath::Pi(), 0,TMath::Pi(), 2);
	raw_fun->SetParameter(0,0);
	raw_fun->SetParameter(1,TMath::Pi());

    c1 = 1/raw_fun->TF2::Integral(0,TMath::Pi(), 0,TMath::Pi());
    cout << "stała A = " << c1 << endl;
    cout << endl;

    fun = new TF2("Gestosc prawdopodobienstwa",my_normed_Func, 0,TMath::Pi(), 0,TMath::Pi(), 3);
	fun->SetParameter(0,0);
	fun->SetParameter(1,TMath::Pi());
    fun->SetParameter(2,c1);
    fun->GetXaxis()->SetTitle("x");
    fun->GetYaxis()->SetTitle("y");
    fun->GetZaxis()->SetTitle("f(x,y)");

    canvas->cd();
    fun->Draw("surf1");
    
    auto canvas2 = new TCanvas("c1","canvas2");
	canvas2->SetCanvasSize(1000, 1000);
	canvas2->SetWindowSize(1100, 1100);
	canvas2->Divide(2,2);

    //histogram1
    TH2D *histogram = new TH2D("histogram","tytul histogramu",10, 0,TMath::Pi(), 10, 0, TMath::Pi());
    fill_hist(100, histogram);
    normalize_hist(histogram);
    histogram->SetTitle("histogram N = 100");
	histogram->GetXaxis()->SetTitle("x");			
	histogram->GetYaxis()->SetTitle("y");
	histogram->GetZaxis()->SetTitle("f(x,y)");
    canvas2->cd(1);
	histogram->Draw("lego2");

    cout << "histogram N = 100" << endl;
    cout << "E(X) = " << histogram->TH1::GetMean(1) << endl;
    cout << "E(Y) = " << histogram->TH1::GetMean(2) << endl;
    cout << "Std. Dev X = " << histogram->TH1::GetRMS(1) << endl;
    cout << "Std. Dev Y = " << histogram->TH1::GetRMS(2) << endl;
    cout << "cov(X,Y) = " << histogram->TH2::GetCovariance(1,2) << endl;
    cout << "rho(X,Y) = " << histogram->TH2::GetCorrelationFactor(1,2) << endl;
    cout << endl;

    //histogram2
    TH2D *histogram2 = new TH2D("histogram2","tytul histogramu2",10, 0,TMath::Pi(), 10, 0, TMath::Pi());
    fill_hist(1000, histogram2);
    normalize_hist(histogram2);
    histogram2->SetTitle("histogram N = 1000");
	histogram2->GetXaxis()->SetTitle("x");			
	histogram2->GetYaxis()->SetTitle("y");
	histogram2->GetZaxis()->SetTitle("f(x,y)");
	canvas2->cd(2);
    histogram2->Draw("lego2");

    cout << "histogram N = 1000" << endl;
    cout << "E(X) = " << histogram2->TH1::GetMean(1) << endl;
    cout << "E(Y) = " << histogram2->TH1::GetMean(2) << endl;
    cout << "Std. Dev X = " << histogram2->TH1::GetRMS(1) << endl;
    cout << "Std. Dev Y = " << histogram2->TH1::GetRMS(2) << endl;
    cout << "cov(X,Y) = " << histogram2->TH2::GetCovariance(1,2) << endl;
    cout << "rho(X,Y) = " << histogram2->TH2::GetCorrelationFactor(1,2) << endl;
    cout << endl;

    //histogram3
    TH2D *histogram3 = new TH2D("histogram3","tytul histogramu3",10, 0,TMath::Pi(), 10, 0, TMath::Pi());
    fill_hist(10000, histogram3);
    normalize_hist(histogram3);
    histogram3->SetTitle("histogram N = 10000");
	histogram3->GetXaxis()->SetTitle("x");			
	histogram3->GetYaxis()->SetTitle("y");
	histogram3->GetZaxis()->SetTitle("f(x,y)");
	canvas2->cd(3);
	histogram3->Draw("lego2");

    cout << "histogram N = 10000" << endl;
    cout << "E(X) = " << histogram3->TH1::GetMean(1) << endl;
    cout << "E(Y) = " << histogram3->TH1::GetMean(2) << endl;
    cout << "Std. Dev X = " << histogram3->TH1::GetRMS(1) << endl;
    cout << "Std. Dev Y = " << histogram3->TH1::GetRMS(2) << endl;
    cout << "cov(X,Y) = " << histogram3->TH2::GetCovariance(1,2) << endl;
    cout << "rho(X,Y) = " << histogram3->TH2::GetCorrelationFactor(1,2) << endl;
    cout << endl;

    //histogram4
    TH2D *histogram4 = new TH2D("histogram4","tytul histogramu4",10, 0,TMath::Pi(), 10, 0, TMath::Pi());
    fill_hist(100000, histogram4);
    normalize_hist(histogram4);
    histogram4->SetTitle("histogram N = 100000");
	histogram4->GetXaxis()->SetTitle("x");			
	histogram4->GetYaxis()->SetTitle("y");
	histogram4->GetZaxis()->SetTitle("f(x,y)");
	canvas2->cd(4);
	histogram4->Draw("lego2");

    cout << "histogram N = 100000" << endl;
    cout << "E(X) = " << histogram4->TH1::GetMean(1) << endl;
    cout << "E(Y) = " << histogram4->TH1::GetMean(2) << endl;
    cout << "Std. Dev X = " << histogram4->TH1::GetRMS(1) << endl;
    cout << "Std. Dev Y = " << histogram4->TH1::GetRMS(2) << endl;
    cout << "cov(X,Y) = " << histogram4->TH2::GetCovariance(1,2) << endl;
    cout << "rho(X,Y) = " << histogram4->TH2::GetCorrelationFactor(1,2) << endl;
    cout << endl;

    cout<<"E(X) pierwotnej funkcji: "<< fun->TF1::Mean(0,TMath::Pi()) <<endl;

return 0;
}