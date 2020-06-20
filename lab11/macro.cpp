/*
KADD lab11
Kacper Skelnik
*/
#include <TMath.h>
#include <TMatrixD.h>

using namespace std;

void wypisz(TMatrixD *mac){
    for (int i=0; i<mac->GetNrows(); i++){
        for(int j=0; j<mac->GetNcols(); j++){
            cout << (*mac)(i,j) << "\t";
            if(j == mac->GetNcols()-1){ cout << endl; }
        }
    }
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
            cout << endl;

            return true;
		}
		else{
            cout << "Poziom istotnosci: " << alfa << endl;
            cout << "Wartość statystyki testowej: " << T << endl;
            cout << "Liczba stopni swobody: " << stopnie_swobody << endl;
            cout << "Wartość krytyczna to: " << chi_val << endl;
            cout << "Nalezy odrzucic hipoteze" << endl;
            cout << endl;

			return false;
		}
}

int getNDF(TGraphErrors *Tge, TF1 *Tf1){
    return Tge->GetN() - Tf1->GetNpar() - 1;
}

double dopasuj(int r, int n, vector<double> X, vector<double> Y, vector<double> sigma, double *wsp){
    TMatrixD *A = new TMatrixD(n,r+1);
    TMatrixD *H = new TMatrixD(n,n); 
    TMatrixD *c = new TMatrixD(n,1);

    for (int i=0; i<n; i++){
        for (int j=0; j <= r; j++){
            (*A)(i,j) = pow(X[i],j);
        }

        for (int j=0; j<n; j++){
            if(i==j){
                (*H)(i,j) = 1/sigma[i];
            }
            else{
                (*H)(i,j) = 0;
            }
        }

        (*c)(i,0) = Y[i];
    }

    TMatrixD *Yprim = new TMatrixD(*H, TMatrixD::kMult, *c);
    TMatrixD *Aprim = new TMatrixD(*H, TMatrixD::kMult, *A);

    TMatrixD *AprimT = new TMatrixD(TMatrixD::kTransposed, *Aprim);
    TMatrixD *AprimTAprim = new TMatrixD(*AprimT, TMatrixD::kMult, *Aprim);
    TMatrixD *AprimTAprimInv = new TMatrixD(TMatrixD::kInverted, *AprimTAprim);
    TMatrixD *AprimTYprim = new TMatrixD(*AprimT, TMatrixD::kMult, *Yprim);

    TMatrixD *Xprim = new TMatrixD(*AprimTAprimInv, TMatrixD::kMult, *AprimTYprim);
    TMatrixD *ni = new TMatrixD(*A, TMatrixD::kMult, *Xprim);

    //wypisz(ni);
    
    for (int i=0; i<r+1; i++){
        wsp[i] = (*Xprim)(i,0);
    }

    double M = 0;
    for (int i=0; i<n; i++){
        M += (*c)(i, 0) - (*ni)(i, 0);
    }

	return M;
}

int macro(){
    double val1;
    double val2;
    vector<double> X;
    vector<double> Y;
    vector<double> sigma;

    ifstream ifile;
	ifile.open("dane2.dat");
	while (ifile >> val1 >> val2) 
	{
		X.push_back(val1);
		Y.push_back(val2);
		sigma.push_back(TMath::Sqrt(val2));
	}
	ifile.close();

    int r = 5;
    int n = X.size();
    double wsp[r];

    TF1* f0 = new TF1("st. 0", "[0]",-1,1);
	TF1* f1 = new TF1("st. 1", "[0]+x*[1]",-1,1);
	TF1* f2 = new TF1("st. 2", "[0]+x*[1]+TMath::Power(x,2)*[2]",-1,1);
	TF1* f3 = new TF1("st. 3", "[0]+x*[1]+TMath::Power(x,2)*[2]+TMath::Power(x,3)*[3]",-1,1);
	TF1* f4 = new TF1("st. 4", "[0]+x*[1]+TMath::Power(x,2)*[2]+TMath::Power(x,3)*[3]+TMath::Power(x,4)*[4]",-1,1);
	TF1* f5 = new TF1("st. 5", "[0]+x*[1]+TMath::Power(x,2)*[2]+TMath::Power(x,3)*[3]+TMath::Power(x,4)*[4]+TMath::Power(x,5)*[5]",-1,1);

    TF1 **ff = new TF1*[r]; 
    ff[0] = f0;
	ff[1] = f1;
	ff[2] = f2;
	ff[3] = f3;
	ff[4] = f4;
	ff[5] = f5;

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1350, 900);
	c->SetWindowSize(1400, 950);

    c->cd(1);
    double err[n];
    for (int i=0; i<n; i++){ err[i] = 0; }
    TGraphErrors *gr = new TGraphErrors(X.size(), &(X[0]), &(Y[0]), err, &(sigma[0]));
    gr->SetMarkerStyle(20);
	gr->Draw("AP");

    for(int i =0; i<=r; i++)
	{
		double M = dopasuj(i, n, X, Y, sigma, wsp);

        for (int j=0; j<r; j++) {ff[i]->SetParameter(j, wsp[j]);}
        ff[i]->SetLineColor(i+1);
        ff[i]->Draw("SAME");
        
        cout<<"Wspolczynniki dopasowania: " << endl;
        for (int j=0; j<i+1; j++){
            cout << wsp[j] << ", ";
        }
        cout<<endl;
        cout<<"Dopasowanie wielomianem stopnia "<<i<<endl;
        test_chi2(M, 0.05, getNDF(gr, ff[i]));
	}

return 0;
}