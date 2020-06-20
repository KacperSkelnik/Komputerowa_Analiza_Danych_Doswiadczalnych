/*
KADD lab8
Kacper Skelnik
*/
#include <TMath.h>

using namespace std;
TRandom *gen = new TRandom();
std::default_random_engine generator;

Double_t Gauss(double *x, double *par){ 
    if (x[0] > -10 && x[0] < 10){
        return (1/(par[1]*TMath::Sqrt(2* TMath::Pi())))*TMath::Exp(-((x[0]-par[0])/par[1])*((x[0]-par[0])/par[1])/(2));
    }
    else{   return 0;    }
}

Double_t Exp(double *x, double *par){ 
    if (x[0] > 0 && x[0] < 20){
        return TMath::Exp(-x[0]/5);
    }
    else{   return 0;    }
}

Double_t Uniform(double* x, double* par){
    if(x[0] >= 0 && x[0] <= 2){
        return 0.5;
    }
    else{
        return 0;
    }
}

void splot(TF1 **tab, int k, TH1D* h, int n){
    TF1 **tab2 = tab;

    for(int i = 0; i < k-1; i++) {
        TF1Convolution *f_conv = new TF1Convolution(tab2[i], tab2[i+1], -10, 20, true);
        f_conv->SetRange(-10,20);
        f_conv->SetNofPointsFFT(n);

        TF1* func_conv = new TF1("f", *f_conv, -10, 20, f_conv->GetNpar() ); 
        func_conv->SetNpx(n);
        func_conv->SetParameters(tab2[i]->GetParameter(0), tab2[i]->GetParameter(1), tab2[i+1]->GetParameter(0), tab2[i+1]->GetParameter(1));

        tab2[i+1] = func_conv;
    }

    for (int i = 0;  i < n; i++) {
        h->Fill(tab2[k-1]->GetRandom(-10, 20));
    }
}

void deskaGaltona_pascal(TH1D* hist, int n, int l){
    l = l+1;
    long long **trojkatPascala;
    trojkatPascala= new long long *[l];
    for (int j=0;j<l;j++){
    trojkatPascala[j]=new long long [j+1];
    trojkatPascala[j][0]=1;
    trojkatPascala[j][j]=1;

    for (int i=0; i<j-1; i++){
        trojkatPascala[j][i+1]=trojkatPascala[j-1][i]+trojkatPascala[j-1][i+1];
        //cout << trojkatPascala[j][i+1] << " ";
        //if(i == j-2) {cout << endl;}
        }
    }

    int sum = 0;
    vector<double> prob;

    for (int r = 0; r < l; r++){    sum += trojkatPascala[l-1][r];    }

    for (int j = 0; j < l; j++){    prob.push_back((double)trojkatPascala[l-1][j]/sum);    }

    std::discrete_distribution<> d(prob.begin(), prob.end());

    for (int i = 0; i <n; i++){
        hist->Fill(d(generator));
    }
}

void deskaGaltona(TH1D * hist, int n, int l, float p){
	for(int i = 0; i<n; i++){
		int bin = 0;
		for(int j = 0; j<l; j++){
			if(gen->Uniform(0, 1) > p){
				bin+=1;
			}
			else{
				bin-=1;
			}
		}
		hist->Fill(bin);
	}	
}

int macro(){ 
    int n = 10000;
    
    TF1 *gauss1 = new TF1("gauss", Gauss, -10, 20, 2);
    gauss1->SetParameter(0,0);
	gauss1->SetParameter(1,2);

    TF1 *gauss2 = new TF1("gauss2", Gauss, -10, 20, 2);
    gauss2->SetParameter(0,3);
	gauss2->SetParameter(1,2);

    TF1 *exp = new TF1("exp", Exp, -10, 20, 0);

    TF1 *U = new TF1("U", Uniform, 0, 20, 0);


    TH1D *gauss_gauss = new TH1D("h1", "Gauss-Gauss", 100, -10, 20);
    TF1 *gaus_gaus_tab[3] = {gauss1, gauss2};
    TF1 **tab1 = gaus_gaus_tab;
    splot(tab1, 2   , gauss_gauss, n);

    TH1D *gauss_exp = new TH1D("h2", "Gauss-Exp", 100, -10, 20);
    TF1 *gaus_exp_tab[2] = {gauss1, exp};
    TF1 **tab2 = gaus_exp_tab;
    splot(tab2, 2, gauss_exp, n);

    TH1D *uniform_uniform = new TH1D("h3", "uniform_x2", 100, 0, 10);
    TF1 *uniform_uniform_tab[2] = {U, U};
    TF1 **tab3 = uniform_uniform_tab;
    splot(tab3, 2, uniform_uniform, n);

    TH1D *uniform_x3 = new TH1D("h4", "uniform_x3", 100, 0, 10);
    TF1 *uniform_x3_tab[3] = {U, U, U};
    TF1 **tab4 = uniform_x3_tab;
    splot(tab4, 3, uniform_x3, n);

    TH1D *uniform_x4 = new TH1D("h5", "uniform_x4", 100, 0, 10);
    TF1 *uniform_x4_tab[4] = {U, U, U, U};
    TF1 **tab5 = uniform_x4_tab;
    splot(tab5, 4, uniform_x4, n);
   

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1350, 900);
	c->SetWindowSize(1400, 950);
    c->Divide(4,2);

    c->cd(1);
    gauss_gauss->Draw();

    c->cd(2);
    gauss_exp->Draw();

    c->cd(3);
    uniform_uniform->Draw();

    c->cd(4);
    uniform_x3->Draw();

    c->cd(5);
    uniform_x4->Draw();


    int l = 10;
    float p = 0.5;

    TH1D *galton_pascal = new TH1D("h6", "deska Galtona - trojkat pascala", 100, 0, 10);
    deskaGaltona_pascal(galton_pascal, n, l);

    TH1D *galton_sym = new TH1D("h7", "deska Galtona - symulacja", 100, -10, 10);
    deskaGaltona(galton_sym, n, l, p);

    c->cd(6);
    galton_pascal->Draw();

    c->cd(7);
    galton_sym->Draw();

return 0;
}