/*
KADD lab3
Kacper Skelnik
*/

#include <iostream>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TF12.h>

using namespace std;

Double_t sigma_x = 0.5;
Double_t sigma_y = 1;
Double_t c1;
TF2 *fun;
double dx1 = 0.1;
double dx2 = 0.2;

Double_t myFunc(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return TMath::Exp(-((x[0]*x[0])/(2*sigma_x*sigma_x))-((x[1]*x[1])/(2*sigma_y*sigma_y)));
    }
    else return 0;
}
Double_t my_normed_Func(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return par[2]*TMath::Exp(-((x[0]*x[0])/(2*sigma_x*sigma_x))-((x[1]*x[1])/(2*sigma_y*sigma_y)));
    }
    else return 0;
}

double cumulative(double *x, double *par){
	fun->SetParameter(0,-2);
	fun->SetParameter(1,2);
	fun->SetParameter(2,c1);
	return fun->TF2::Integral(-2,x[0],-2,x[1]);
}

int number_of_i(double xp, double xk, double dx){
    return (xk-xp)/dx;
}

Double_t to_calc_g(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return par[2]*TMath::Exp(-(x[0]*x[0])/(2*sigma_y*sigma_y));
    }
    else return 0;
}

Double_t to_calc_h(double *x, double *par){ 
    if ((x[0] >= par[0] && x[0] <= par[1]) && (x[1] >= par[0] && x[1] <= par[1])){
        return par[2]*TMath::Exp(-(x[0]*x[0])/(2*sigma_x*sigma_x));
    }
    else return 0;
}

double *x(double dx, int N){
    int k, m;

    double x2[N/2];
    double x[N];
    for(int q = 0; q <= N/2; q++){
        x[q] = -q*dx;
    }
    for(int q = 0; q <= N/2; q++){
        x2[q] = q*dx;
    }
    for(m=0, k=N/2; k<=N && m<=N/2; m++, k++){
		x[k]=x2[m];
    }

    return x;
}

double *g_x(double dx){
    int N = number_of_i(-2, 2, dx);

    double x1[N];
    double *x_tymczasowe1 = x(dx, N);
    for(int i= 0; i < N; i ++){ x1[i] = x_tymczasowe1[i];    }

    double y[N];

    double y1 = 0;
    fun = new TF2("tymczasowa",to_calc_g, -2,2,-2,2, 3);
	fun->SetParameter(0,-2);
	fun->SetParameter(1,2);
    fun->SetParameter(2,c1);
    y1 = fun->TF2::Integral(-2,2);

    for(int i = 0; i < N;  i++){
        y[i] = y1*TMath::Exp(-((x1[i]*x1[i])/(2*sigma_x*sigma_x)));
    }

    return y;
}

double *h_y(double dx){
    int N = number_of_i(-2, 2, dx);

    double x1[N];
    double *x_tymczasowe1 = x(dx, N);
    for(int i= 0; i < N; i ++){ x1[i] = x_tymczasowe1[i];    }

    double y[N];
    
    double y1 = 0;
    fun = new TF2("tymczasowa",to_calc_h, -2,2,-2,2, 3);
	fun->SetParameter(0,-2);
	fun->SetParameter(1,2);
    fun->SetParameter(2,c1);
    y1 = fun->TF2::Integral(-2,2);

    for(int i = 0; i < N;  i++){
        y[i] = y1*TMath::Exp(-((x1[i]*x1[i])/(2*sigma_y*sigma_y)));
    }

    return y;
}

int macro(){

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1000, 600);
	c->SetWindowSize(1100, 700);
	c->Divide(2,2);
	
	TF2 *raw_fun = new TF2("f1",myFunc, -2,2,-2,2, 2);
	raw_fun->SetParameter(0,-2);
	raw_fun->SetParameter(1,2);

    c1 = 1/raw_fun->TF2::Integral(-2,2, -2,2);
    cout << "stała: " << c1 << endl;
        
    fun = new TF2("Gestosc prawdopodobienstwa",my_normed_Func, -2,2,-2,2, 3);
	fun->SetParameter(0,-2);
	fun->SetParameter(1,2);
    fun->SetParameter(2,c1);
    fun->GetXaxis()->SetTitle("x^2");
    fun->GetYaxis()->SetTitle("y^2");
    c->cd(1);
    fun->Draw("surf1");

    TF2 *cdf = new TF2("Dystrybuanta", cumulative, -2,2 , -2,2, 3);
    cdf->GetXaxis()->SetTitle("x^2");
    cdf->GetYaxis()->SetTitle("y^2");
	c->cd(2);
	cdf->Draw("surf1");

    /*c->cd(3);
    TF12 *gx = new TF12("f12",fun,0,"x");
    gx->SetTitle("g(x)");
    gx->Draw("C");
    
    c->cd(4);
    TF12 *hy = new TF12("f12",fun,0,"y");
    hy->SetTitle("h(y)");
    hy->Draw("C");*/

    TMultiGraph *g_graph = new TMultiGraph();
    g_graph->SetTitle("g(x)");
    g_graph->GetXaxis()->SetTitle("x");
    g_graph->GetYaxis()->SetTitle("g(x)");

    TMultiGraph *h_graph = new TMultiGraph();
    h_graph->SetTitle("h(y)");
    h_graph->GetXaxis()->SetTitle("y");
    h_graph->GetYaxis()->SetTitle("h(y)");



    int N = number_of_i(-2, 2, dx1);
    double x1[N];
    double *x_tymczasowe1 = x(dx1, N);
    for(int i= 0; i < N; i ++){ x1[i] = x_tymczasowe1[i];    }

    double g1[N];
    double *y_tymczasowe1 = g_x(dx1);
    for(int i= 0; i < N; i ++){ g1[i] = y_tymczasowe1[i];    }
    TGraph *gx1 = new TGraph (N, x1, g1);
    gx1->SetMarkerColor(632);
    gx1->SetMarkerSize(1);
    gx1->SetMarkerStyle(22);

    double h1[N];
    y_tymczasowe1 = h_y(dx1);
    for(int i= 0; i < N; i ++){ h1[i] = y_tymczasowe1[i];    }
    TGraph *hy1 = new TGraph (N, x1, h1);
    hy1->SetMarkerColor(632);
    hy1->SetMarkerSize(1);
    hy1->SetMarkerStyle(22);
    
    
    
    int N2 = number_of_i(-2, 2, dx2);
    double x2[N2];
    double *x_tymczasowe2 = x(dx2, N2);
    for(int i= 0; i < N2; i ++){ x2[i] = x_tymczasowe2[i];    }

    double g2[N2];
    double *y_tymczasowe2 = g_x(dx2);
    for(int i= 0; i < N; i ++){ g2[i] = y_tymczasowe2[i];    }
    TGraph *gx2 = new TGraph (N2, x2, g2);
    gx2->SetMarkerColor(4);
    gx2->SetMarkerSize(1);
    gx2->SetMarkerStyle(21);

    double h2[N2];
    y_tymczasowe2 = h_y(dx2);
    for(int i= 0; i < N2; i ++){ h2[i] = y_tymczasowe2[i];    }
    TGraph *hy2 = new TGraph (N2, x2, h2);
    hy2->SetMarkerColor(4);
    hy2->SetMarkerSize(1);
    hy2->SetMarkerStyle(21);



    c->cd(3);
    g_graph->Add(gx1);
    g_graph->Add(gx2);
    g_graph->Draw("AP");
    auto legend1 = new TLegend(0.1,0.9,0.3,0.7);
    legend1->AddEntry(hy1,"krok = 0.1","p");
    legend1->AddEntry(hy2,"krok = 0.2","p");
    legend1->Draw();

    c->cd(4);
    h_graph->Add(hy1);
    h_graph->Add(hy2);
    h_graph->Draw("AP");
    auto legend2 = new TLegend(0.1,0.9,0.3,0.7);
    legend2->AddEntry(hy1,"krok = 0.1","p");
    legend2->AddEntry(hy2,"krok = 0.2","p");
    legend2->Draw();

    return 0;
}