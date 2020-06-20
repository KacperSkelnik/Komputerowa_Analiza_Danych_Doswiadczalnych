/*
KADD lab5
Kacper Skelnik
*/

#include <TMatrixT.h>
#include <TMath.h>

using namespace std;

int macro(){    

    ifstream plik1, plik2, plik3;
	plik1.open("dane1.dat");
	plik2.open("dane2.dat");
	plik3.open("dane3.dat");

    double val1, val2, val3;

    std::vector < double > tab1;
    std::vector < double > tab2;
    std::vector < double > tab3;
    
    while (plik1 >> val1) {
        tab1.push_back(val1);
    }
    while (plik2 >> val2) {
        tab2.push_back(val2);
    }
    while (plik3 >> val3) {
        tab3.push_back(val3);
    }

    double tab1_max = *max_element(tab1.begin(), tab1.end()); 
    double tab1_min = *min_element(tab1.begin(), tab1.end()); 

    double tab2_max = *max_element(tab2.begin(), tab2.end()); 
    double tab2_min = *min_element(tab2.begin(), tab2.end()); 

    double tab3_max = *max_element(tab3.begin(), tab3.end()); 
    double tab3_min = *min_element(tab3.begin(), tab3.end());

    plik1.close();
	plik2.close();
	plik3.close();

    plik1.open("dane1.dat");
	plik2.open("dane2.dat");
	plik3.open("dane3.dat");


    TH1D* h1 = new TH1D("h1", "Hist1", 100, tab1_min, tab1_max);
    TH1D* h2 = new TH1D("h2", "Hist2", 100, tab2_min, tab2_max);
    TH1D* h3 = new TH1D("h3", "Hist3", 100, tab3_min, tab3_max);
    

    TH2D* h12 = new TH2D("h12", "Hist12", 100, tab1_min, tab1_max, 100, tab2_min, tab2_max);
	TH2D* h13 = new TH2D("h13", "Hist13", 100, tab1_min, tab1_max, 100, tab3_min, tab3_max);
	TH2D* h23 = new TH2D("h23", "Hist23", 100, tab2_min, tab2_max, 100, tab3_min, tab3_max);

    while (plik1 >> val1 && plik2 >> val2 && plik3 >> val3) {
        h1->Fill(val1);
        h2->Fill(val2);
        h3->Fill(val3);

		h12->Fill(val1, val2);
		h13->Fill(val1, val3);
		h23->Fill(val2, val3);
	}

    plik1.close();
	plik2.close();
	plik3.close();
    

    auto canvas = new TCanvas("c","canvas");
	canvas->SetCanvasSize(1000, 1000);
	canvas->SetWindowSize(1100, 1100);
	canvas->Divide(2,3);

    canvas->cd(1);
    h1->Draw();
    canvas->cd(3);
    h2->Draw();
    canvas->cd(5);
    h3->Draw();

    canvas->cd(2);
    h12->Draw("COLZ");
    canvas->cd(4);
    h13->Draw("COLZ");
    canvas->cd(6);
    h23->Draw("COLZ");


	cout << "Srednia dane1: " << h1->GetMean(1) << endl;
	cout << "Srednia dane2: " << h2->GetMean(1) << endl;
	cout << "Srednia dane3: " << h3->GetMean(1) << endl;

    cout << "Odchylenie standardowe dane1: " << h1->GetStdDev(1) << endl;
    cout << "Odchylenie standardowe dane2: " << h2->GetStdDev(1) << endl;
    cout << "Odchylenie standardowe dane3: " << h3->GetStdDev(1) << endl;

    Double_t cov11 = h12->GetCovariance(1, 1);
	Double_t cov22 = h23->GetCovariance(1, 1);
	Double_t cov33 = h23->GetCovariance(2, 2);
	Double_t cov12 = h12->GetCovariance(1, 2);
	Double_t cov13 = h13->GetCovariance(1, 2);
	Double_t cov23 = h23->GetCovariance(1, 2);

    TMatrixD MCX(3, 3);
	MCX(0, 0) = h12->GetCovariance(1, 1);
	MCX(1, 1) = h23->GetCovariance(1, 1);
	MCX(2, 2) = h23->GetCovariance(2, 2);
	MCX(0, 1) = MCX(1, 0) = h12->GetCovariance(1, 2);
	MCX(0, 2) = MCX(2, 0) = h13->GetCovariance(1, 2);
	MCX(1, 2) = MCX(2, 1) = h23->GetCovariance(1, 2);

    cout << "Macierz kowiariancji zmiennych X1, X2 i X3" << endl;
    MCX.Print();

    std::vector < int > sizes;
    sizes.push_back(tab1.size());
    sizes.push_back(tab2.size());
    sizes.push_back(tab3.size());

    int NY = *min_element(sizes.begin(), sizes.end());
    double Y1[NY];
    double Y2[NY];
    for(int i = 0; i < NY ; i++){
        Y1[i] = 2*tab1[i] + 5*tab2[i] + tab3[i];
        Y2[i] = 3 + 0.5*tab1[i] + 4*tab2[i];
    }

    cout << "Średnia wartoś Y1: " << TMath::Mean(NY, Y1) << endl;
    cout << "Średnia wartoś Y2: " << TMath::Mean(NY, Y2) << endl << endl;

    cout << "Odchylenie standardowe wartości Y1: " << TMath::StdDev(NY, Y1) << endl;
    cout << "Odchylenie standardowe wartości Y2: " << TMath::StdDev(NY, Y2) << endl << endl;

    double Y1_max = TMath::MaxElement(NY, Y1);
    double Y1_min = TMath::MinElement(NY, Y1);

    double Y2_max = TMath::MaxElement(NY, Y2);
    double Y2_min = TMath::MinElement(NY, Y2);

    TH2D* hY12 = new TH2D("hY12", "HistY12", 100, Y1_min, Y1_max, 100, Y2_min, Y2_max);
    TH1D* hY1 = new TH1D("Y1", "HistY1", 100, Y1_min, Y1_max);
    TH1D* hY2 = new TH1D("Y2", "HistY2", 100, Y2_min, Y2_max);


    for(int i = 0; i < NY ; i++){
        hY12->Fill(Y1[i], Y2[i]);
        
        hY1->Fill(Y1[i]);
        hY2->Fill(Y2[i]);
    }

    TMatrixD MCY(2, 2);
	MCY(0, 0) = hY12->GetCovariance(1, 1);
	MCY(1, 1) = hY12->GetCovariance(2, 2);
	MCY(0, 1) = MCY(1, 0) = hY12->GetCovariance(1, 2);

    cout << "Macierz kowiariancji zmiennych Y1 i Y2" << endl;
    MCY.Print();

    cout << "Wspóczunnik korelacji pomiędzy zmiennymi Y1 i Y2: " << hY12->GetCorrelationFactor(1, 2) << endl;

    auto canvas2 = new TCanvas("c2","canvas2");
	canvas2->SetCanvasSize(1000, 1000);
	canvas2->SetWindowSize(1100, 1100);
	canvas2->Divide(2,2);

    canvas2->cd(1);
    hY1->Draw();
    canvas2->cd(3);
    hY2->Draw();

    canvas2->cd(2);
    hY12->Draw("COLZ");

return 0;
}