/*
KADD lab9
Kacper Skelnik
*/
#include <TMath.h>

using namespace std;

Double_t chi_2(double *x, double *par){ 
    return (1/(TMath::Power(2,par[0]/2)*TMath::Gamma(par[0]/2)))*TMath::Power(x[0],(par[0]/2)-1)*TMath::Exp(-x[0]/2);
}

Double_t Uniform(double* x, double* par){
    if(x[0] >= 0 && x[0] <= 1){
        return 1;
    }
    else{
        return 0;
    }
}

void splot(TF1 **tab, int k, TH1D* h, int n){
    for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < k; j++) {
			double y = tab[j]->GetRandom();
			sum += y;
		}
		h->Fill(sum);
	}
}

int macro(){

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1350, 900);
	c->SetWindowSize(1400, 950);
    c->Divide(2,2);

    c->cd(1);

    TF1 *chi2[46];
    chi2[0] = new TF1("chi2", chi_2, 0, 100, 1);
    chi2[0]->SetParameter(0,5);
    chi2[0]->SetLineColor(1);
    chi2[0]->SetLineWidth(4);
    chi2[0]->Draw();
    
    for (int i = 1; i <= 45; i++){
        chi2[i] = new TF1("chi2", chi_2, 0, 100, 1);
        chi2[i]->SetParameter(0,i+5);
        chi2[i]->SetLineColor(i+1);
        chi2[i]->SetLineWidth(4);
        chi2[i]->Draw("SAME");
    }

    c->cd(2);

    TGraph *g = (TGraph*)chi2[0]->DrawIntegral("P");
    g->Draw();

    for (int i = 1; i <= 45; i++){
        TGraph *g = (TGraph*)chi2[i]->DrawIntegral("P");
        g->Draw("SAME");
    }

    c->cd(3);
    
    int n = 10000;

    TF1 *U = new TF1("U", Uniform, 0, 1, 0);

    TH1D *uniform[12];
    TF1 *fit[12];

    Double_t chi22, ndf;

    for(int i = 2; i <= 11; i ++){
        uniform[i] = new TH1D("", "zad2", 100, 0, 10);
        TF1 *uniform_uniform_tab[i];
        for(int j = 0; j < i; j++){
            uniform_uniform_tab[j] = U;
        }
        TF1 **tab = uniform_uniform_tab;
        splot(tab, i, uniform[i], n);

        uniform[i]->SetMarkerStyle(kFullCircle);
        uniform[i]->SetMarkerSize(0.8);
        if(i != 11) {   uniform[i]->SetMarkerColor(i-1);  }
        else{   uniform[i]->SetMarkerColor(12);  }

        uniform[i]->Fit("gaus", "Q", "SAME");

        fit[i] = uniform[i]->GetFunction("gaus");
        if(i != 11) {   fit[i]->SetLineColor(i-1);  }
        else{   fit[i]->SetLineColor(12);  }

        uniform[i]->Draw("SAME P E");
    }

    for(int i = 2; i <= 11; i ++){
        cout << "For n = " << i << "\t" << "chi2/ndf: " << fit[i]->GetChisquare() << " / " << fit[i]->GetNDF() << " = " << fit[i]->GetChisquare()/fit[i]->GetNDF() << endl;
    }

return 0;
}