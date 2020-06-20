/*
KADD lab10
Kacper Skelnik
*/
#include <TMath.h>

using namespace std;

Double_t Teor(double* x, double* par){
    //return TMath::Sqrt(TMath::Pi()) * par[0] * TMath::Exp(par[1]*par[1]/(4*par[2])) * TMath::Erf( (par[1] + 2*par[2]* x[0]) / (2 * TMath::Pi()) ) / ( 2* TMath::Pi() );

    return par[0]*TMath::Exp(-par[1]*x[0]*x[0] - par[2]*x[0]);
}

double chi2(TGraphErrors* Tge, TF1* Tf1){
    double T;

    for (int i=0; i < Tge->GetN(); i++){
        T += TMath::Power((Tge->GetY()[i] - Tf1->Eval(Tge->GetX()[i])),2)/(Tf1->Eval(Tge->GetX()[i]));
    }

    return T;
}

bool test_chi2(double T, double alfa, int stopnie_swobody) {
		double lambda = double(stopnie_swobody) / 2;
		TF1* chi_sq_cdf = new TF1("chi_sq_cdf", "TMath::Gamma([0], x / 2)", 0, 30);
		chi_sq_cdf->SetParameter(0, lambda);
		double chi_val = chi_sq_cdf->GetX(1 - alfa);

		if (T < chi_val){
            cout << "Poziom istotnosci: " << alfa << endl;
            cout << "Wartość statystyki testowej: " << T << endl;
            cout << "Liczba stopni swobody: " << stopnie_swobody << endl;
            cout << "Wartość krytyczna to: " << chi_val << endl;
            cout << "Nie ma powodów do odrzucenia hipotezy" << endl;

            return true;
		}
		else{
            cout << "Poziom istotnosci: " << alfa << endl;
            cout << "Wartość statystyki testowej: " << T << endl;
            cout << "Liczba stopni swobody: " << stopnie_swobody << endl;
            cout << "Wartość krytyczna to: " << chi_val << endl;
            cout << "Nalezy odrzucic hipoteze" << endl;

			return false;
		}
}

int getNDF(TGraphErrors *Tge, TF1 *Tf1){
    return Tge->GetN() - Tf1->GetNpar() - 1;
}


int macro(){

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1350, 900);
	c->SetWindowSize(1400, 950);
    c->SetLogy();

    auto gr = new TGraphErrors("dane.dat", "%lg %lg %lg %lg", "\t");
    
    TF1 *fun_teor = new TF1("f1",Teor, 0, 3, 3);
    fun_teor->SetParameters(1, 1, 1);

    gr->Fit("f1", "L");
    fun_teor = gr->GetFunction("f1");

    gr->Draw("AP");

    test_chi2(chi2(gr, fun_teor), 0.01, getNDF(gr, fun_teor));

return 0;
}