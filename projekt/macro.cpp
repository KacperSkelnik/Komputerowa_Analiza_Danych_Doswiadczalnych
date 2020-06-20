/*
KADD projekt
Kacper Skelnik
*/

using namespace std;

void read_file(string path, vector<double>& X, vector<vector<double>>& rozklady_norm, vector<vector<double>>& rozklady){
    /*
    funkcja odpowiedzialna za pobieranie danych z plikow
    wczutuje dane:
    X - znormalizowany rozklad
    rozklady - tablica na wszystkie rozklady
    rozklady_norm - tablica na wszystkie znormalizowane
    */
    int val1;
    int val2;

    ifstream ifile;
	ifile.open(path);
	while (ifile >> val1 >> val2){
        if(val2 != 0){  X.push_back(val2);  }
	}

    rozklady.push_back(X);

    double X_max = *max_element(X.begin(), X.end());
    for(int i=0; i<X.size(); i++){
        X[i] = X[i]/X_max;
    }
    rozklady_norm.push_back(X);

	ifile.close();
}

TGraph * def_graph(vector<double>& Y){
    /*
    funkcja odpowiedzialna za definicje TGraph
    */
    vector<double> year;

    int min = 2019 - Y.size() +1; // trzeba tak z racji na rozna ilosc lat w danych
    int max = 2019;

    for(int i=min; i<=max; i++){
        year.push_back(i);
    }

    TGraph *gr = new TGraph(Y.size(), &year[0], &Y[0]);

    return gr;
}

void draw(TGraph *gr, string title, int color){
    /*
    funkcja odpowiedzialna za rysowanie wykresow
    */
    const char *t = title.c_str();
    gr->SetTitle(t);
    gr->GetXaxis()->SetTitle("rok");
    gr->GetXaxis()->CenterTitle();
    gr->GetYaxis()->SetTitle("unormowana liczba zliczen");
    gr->GetYaxis()->CenterTitle();

    gr->GetXaxis()->SetRangeUser(2003,2019);

    gr->SetMarkerColor(color);
    gr->SetMarkerSize(1);
    gr->SetMarkerStyle(20);
    gr->Draw("AP");
}

double std_dev(vector<double>& X){
    /*
    funkcja odpowiedzialna za liczenie odchylenia standardowego
    */
    double x;
    double x2;
    for(int i=0; i<X.size(); i++){
        x2 += X[i]*X[i]/X.size();
        x += X[i]/X.size();
    }

    return TMath::Sqrt(x2 - x*x);
}

double statystyka_chi2(vector<double>& O, vector<double>& E){
    /*
    funkcja odpowiedzialna za liczenie statystkyki chi2
    */
    vector<double> O_tmp;
    vector<double> E_tmp;

    // ponizej jest sposob na to zeby obslugiwac sytuacje gdy mamy rozne dlugosci tablic
    if(O.size() > E.size()){
        for(int i=0; i<E.size(); i++){
            E_tmp.push_back(E[i]);
        }
        for(int i=O.size(); i>int(O.size()-E.size()); i--){
            O_tmp.insert(O_tmp.begin(), O[i]); //dodaj na poczatek (bo for leci od konca)
        }
    }
    else if (O.size() < E.size()){
        for(int i=0; i<O.size(); i++){
            O_tmp.push_back(O[i]);
        }
        for(int i=E.size(); i>int(E.size()-O.size()); i--){
            E_tmp.insert(E_tmp.begin(), E[i]);
        }
    }
    else{
        for(int i=0; i<O.size(); i++){
            O_tmp.push_back(O[i]);
            E_tmp.push_back(E[i]);
        }
    }

    if(E_tmp.size() != O_tmp.size()){cout << "cos jest zle " << endl;}
    
    double chi2;
    double std = std_dev(E_tmp);

    if(O.size() >= E_tmp.size()){
        for(int i=E_tmp.size(); i>0; i--){
            chi2 += TMath::Power((O_tmp[i-1] - E_tmp[i-1])/std,2);
        }
    }
    else if(O.size() < E.size()){
        for(int i=O.size(); i>0; i--){
            chi2 += TMath::Power((O_tmp[i-1] - E_tmp[i-1])/std,2);
        }
    }

    return chi2;
}

bool test_chi2(double T, double alfa, int stopnie_swobody) {
    /*
    funkcja odpowiedzialna za zwracanie wyniku testu chi2
    */
    double lambda = double(stopnie_swobody) / 2;
    TF1* chi_sq_cdf = new TF1("chi_sq_cdf", "TMath::Gamma([0], x / 2)", 0, 30);
    chi_sq_cdf->SetParameter(0, lambda);
    double chi_val = chi_sq_cdf->GetX(1 - alfa);

    if (T < chi_val){
        //cout << "Poziom istotnosci: " << alfa << endl;
        //cout << "Wartość statystyki testowej: " << T << endl;
        //cout << "Liczba stopni swobody: " << stopnie_swobody << endl;
        //cout << "Wartość krytyczna to: " << chi_val << endl;
        //cout << "Nie ma powodów do odrzucenia hipotezy" << endl;

        return true;
    }
    else{
        //cout << "Poziom istotnosci: " << alfa << endl;
        //cout << "Wartość statystyki testowej: " << T << endl;
        //cout << "Liczba stopni swobody: " << stopnie_swobody << endl;
        //cout << "Wartość krytyczna to: " << chi_val << endl;
        //cout << "Nalezy odrzucic hipoteze" << endl;

        return false;
    }
}

double R2(vector<double>& Y, vector<double>& Y_teor){
    /*
    funkcja odpowiedzialna za liczenie statystyki R2 punktów do punktów
    */
    double mean;
    for(int i=0; i<Y.size(); i++){
        mean += Y[i]/Y.size();
    }

    double R2_licznik;
    double R2_mianownik;
    for(int i=0; i<Y.size(); i++){
        R2_mianownik += (Y[i] - mean)*(Y[i] - mean);
    }
    for(int i=Y_teor.size()-1; i>=int(Y_teor.size()-Y.size()); i--){ //zawsze od konca i jezeli jedna tablica jest wieksza, to nie dojdzie do 0
        R2_licznik += (Y_teor[i] - mean)*(Y_teor[i] - mean);
    }
    
    return R2_licznik/R2_mianownik;
}

double R2_fit(vector<double>& Y, TF1 *Y_teor){
    /*
    funkcja odpowiedzialna za liczenie statystyki R2 punktów do funkcji
    */
    double mean;
    for(int i=0; i<Y.size(); i++){
        mean += Y[i]/Y.size();
    }

    double R2_licznik;
    double R2_mianownik;
    for(int i=0; i<Y.size(); i++){
        R2_mianownik += (Y[i] - mean)*(Y[i] - mean);
    }
    for(int i=2019-Y.size()+1; i<=2019; i++){
        R2_licznik += (Y_teor->Eval(i) - mean)*(Y_teor->Eval(i) - mean); 
    }

    return R2_licznik/R2_mianownik;
}

double test_Kolmogorowa_Smirnowa(vector<double>& X, vector<double>& Y){
    /*
    funkcja odpowiedzialna za liczenie statystyki Kolmogorowa-Smirnowa
    opcja z roznymi dlugosciami tablic rozwiazana tak jak w chi2
    */
    vector<double> X_tmp;
    vector<double> Y_tmp;

    if(Y.size() > X.size()){
        for(int i=0; i<X.size(); i++){
            X_tmp.push_back(X[i]);
        }
        for(int i=Y.size(); i>int(Y.size()-X.size()); i--){
            Y_tmp.insert(Y_tmp.begin(), Y[i]);
        }
    }
    else if (Y.size() < X.size()){
        for(int i=0; i<Y.size(); i++){
            Y_tmp.push_back(Y[i]);
        }
        for(int i=X.size(); i>int(X.size()-Y.size()); i--){
            X_tmp.insert(X_tmp.begin(), X[i]);
        }
    }
    else{
        for(int i=0; i<Y.size(); i++){
            Y_tmp.push_back(Y[i]);
            X_tmp.push_back(X[i]);
        }
    }

    if(X_tmp.size() != Y_tmp.size()){cout << "cos jest zle " << endl;}

    int sumX = 0;
    int sumY = 0;
    sumX = accumulate(X_tmp.begin(), X_tmp.end(), 0);
    sumY = accumulate(Y_tmp.begin(), Y_tmp.end(), 0);

    vector<double> dystrybuantaX;
    vector<double> dystrybuantaY;
    for(int i=0; i<X_tmp.size(); i++){
        dystrybuantaX.push_back(X_tmp[i]/sumX);
        dystrybuantaY.push_back(Y_tmp[i]/sumY);
    }
    
    double test = 0;
    for(int i=0; i<X_tmp.size(); i++){
            test += abs(dystrybuantaX[i] - dystrybuantaY[i]);
    }

    double statystyka = TMath::Sqrt(X_tmp.size()*Y_tmp.size()/(X_tmp.size()+Y_tmp.size()))*test;

    return statystyka;
}

int macro(){
    vector<double> niedobor_seksualny;
    vector<double> dluzsza_nieobecnosc;
    vector<double> naganny_stosunek_do_czlonkow_rodziny;
    vector<double> niedochowanie_wiernosci_malzenskiej;
    vector<double> roznice_swiatopogladowe;
    vector<double> hazard;
    vector<double> narkotyki;
    vector<double> nieporozumienia_na_tle_finansowym;
    vector<double> trudnosci_mieszkaniowe;
    vector<double> naduzywanie_alkoholu;
    vector<double> niezgodnosc_charakterow;

    vector<vector<double>> rozklady_norm;
    vector<vector<double>> rozklady;

    read_file("dane/niedobor_seksualny.txt", niedobor_seksualny, rozklady_norm, rozklady);
    read_file("dane/dluzsza_nieobecnosc.txt", dluzsza_nieobecnosc, rozklady_norm, rozklady);
    read_file("dane/naganny_stosunek_do_czlonkow_rodziny.txt", naganny_stosunek_do_czlonkow_rodziny, rozklady_norm, rozklady);
    read_file("dane/niedochowanie_wiernosci_malzenskiej.txt", niedochowanie_wiernosci_malzenskiej, rozklady_norm, rozklady);
    read_file("dane/roznice_swiatopogladowe.txt", roznice_swiatopogladowe, rozklady_norm, rozklady);
    read_file("dane/hazard.txt", hazard, rozklady_norm, rozklady);
    read_file("dane/narkotyki.txt", narkotyki, rozklady_norm, rozklady);
    read_file("dane/nieporozumienia_na_tle_finansowym.txt", nieporozumienia_na_tle_finansowym, rozklady_norm, rozklady);
    read_file("dane/trudnosci_mieszkaniowe.txt", trudnosci_mieszkaniowe, rozklady_norm, rozklady);
    read_file("dane/naduzywanie_alkoholu.txt", naduzywanie_alkoholu, rozklady_norm, rozklady);
    read_file("dane/niezgodnosc_charakterow.txt", niezgodnosc_charakterow, rozklady_norm, rozklady);

    for(int i=0; i<rozklady.size(); i++){ //sprawdzanie czy sie wczytalo
        if (rozklady[i].empty() == 1){
            cout << "vector " << i << " pusty" <<  endl;
        }
        else if (rozklady_norm[i].empty() == 1){
            cout << "vector znormalizowany " << i << " pusty" <<  endl;
        }    
    }

    vector<double> zbiorczy;
    for(int i=0; i<17; i++){
            zbiorczy.push_back(0);
    }

    for(int i=0; i<rozklady.size(); i++){ // tworzenie rozkladu zbiorczego
        for(int j=0; j<rozklady[i].size(); j++){
            zbiorczy[j] = zbiorczy[j] + rozklady[i][j];
        }
    }

    double zbiorczy_max = *max_element(zbiorczy.begin(), zbiorczy.end()); // tworzenie rozkladu zbiorczego unormowanego
    for(int i=0; i<zbiorczy.size(); i++){
            zbiorczy[i] = zbiorczy[i]/zbiorczy_max;
    }

    TGraph *gr_niedobor_seksualny = def_graph(niedobor_seksualny);
    TGraph *gr_dluzsza_nieobecnosc = def_graph(dluzsza_nieobecnosc);
    TGraph *gr_naganny_stosunek_do_czlonkow_rodziny = def_graph(naganny_stosunek_do_czlonkow_rodziny);
    TGraph *gr_niedochowanie_wiernosci_malzenskiej = def_graph(niedochowanie_wiernosci_malzenskiej);
    TGraph *gr_roznice_swiatopogladowe = def_graph(roznice_swiatopogladowe);
    TGraph *gr_hazard = def_graph(hazard);
    TGraph *gr_narkotyki = def_graph(narkotyki);
    TGraph *gr_nieporozumienia_na_tle_finansowym = def_graph(nieporozumienia_na_tle_finansowym);
    TGraph *gr_trudnosci_mieszkaniowe = def_graph(trudnosci_mieszkaniowe);
    TGraph *gr_naduzywanie_alkoholu = def_graph(naduzywanie_alkoholu);
    TGraph *gr_niezgodnosc_charakterow = def_graph(niezgodnosc_charakterow);

    vector<string> names = { "niedobor_seksualny", "dluzsza_nieobecnosc", "naganny_stosunek_do_czlonkow_rodziny", "niedochowanie_wiernosci_malzenskiej","roznice_swiatopogladowe",
     "hazard", "narkotyki", "nieporozumienia_na_tle_finansowym", "trudnosci_mieszkaniowe", "naduzywanie_alkoholu", "niezgodnosc_charakterow"  }; //tablica przechowywujaca nazwy rozkladow

    TGraph *gr_zbiorczy = def_graph(zbiorczy);

    auto c = new TCanvas("c","canvas");
	c->SetCanvasSize(1350, 900);
	c->SetWindowSize(1400, 950);
    c->Divide(4,3);
    

	TF1* f3 = new TF1("f3", "[0]+x*[1]+TMath::Power(x,2)*[2]+TMath::Power(x,3)*[3]",2003,2019); // def wielomianu 3 stopnia do dopasowywania
    f3->SetParameters(-0.01, -0.01, -0.01, -0.01);


    c->cd(1);
    draw(gr_niedobor_seksualny, "niedobor seksualny", 922);                     //narysuj
    gr_niedobor_seksualny->Fit("pol1");                                         //dopasuj    
    TF1 *fit_niedobor_seksualny = gr_niedobor_seksualny->GetFunction("pol1");   //wez dopasowana funkcje
    c->cd(2);
    draw(gr_dluzsza_nieobecnosc, "dluzsza nieobecnosc", 1);
    gr_dluzsza_nieobecnosc->Fit("pol1");
    TF1 *fit_dluzsza_nieobecnosc = gr_dluzsza_nieobecnosc->GetFunction("pol1");
    c->cd(3);
    draw(gr_naganny_stosunek_do_czlonkow_rodziny, "naganny stosunek do czlonkow rodziny", 1);
    gr_naganny_stosunek_do_czlonkow_rodziny->Fit("pol1");
    TF1 *fit_naganny_stosunek_do_czlonkow_rodziny = gr_naganny_stosunek_do_czlonkow_rodziny->GetFunction("pol1");
    c->cd(4);
    draw(gr_niedochowanie_wiernosci_malzenskiej, "niedochowanie wiernosci malzenskiej", 1);
    gr_niedochowanie_wiernosci_malzenskiej->Fit("pol1");
    TF1 *fit_niedochowanie_wiernosci_malzenskiej = gr_niedochowanie_wiernosci_malzenskiej->GetFunction("pol1");
    c->cd(5);
    draw(gr_roznice_swiatopogladowe, "roznice swiatopogladowe", 4);
    gr_roznice_swiatopogladowe->Fit(f3);
    TF1 *fit_roznice_swiatopogladowe = gr_roznice_swiatopogladowe->GetFunction("f3");
    c->cd(6);
    draw(gr_hazard, "hazard", 2);
    gr_hazard->Fit("pol1");
    TF1 *fit_hazard = gr_hazard->GetFunction("pol1");
    c->cd(7);
    draw(gr_narkotyki, "narkotyki", 2);
    gr_narkotyki->Fit("pol1");
    TF1 *fit_narkotyki = gr_narkotyki->GetFunction("pol1");
    c->cd(8);
    draw(gr_nieporozumienia_na_tle_finansowym, "nieporozumienia na tle finansowym", 1);
    gr_nieporozumienia_na_tle_finansowym->Fit("pol1");
    TF1 *fit_nieporozumienia_na_tle_finansowym = gr_nieporozumienia_na_tle_finansowym->GetFunction("pol1");
    c->cd(9);
    draw(gr_trudnosci_mieszkaniowe, "trudnosci mieszkaniowe", 4);
    gr_trudnosci_mieszkaniowe->Fit(f3);
    TF1 *fit_trudnosci_mieszkaniowe = gr_trudnosci_mieszkaniowe->GetFunction("f3");
    c->cd(10);
    draw(gr_naduzywanie_alkoholu, "naduzywanie alkoholu", 1);
    gr_naduzywanie_alkoholu->Fit("pol1");
    TF1 *fit_naduzywanie_alkoholu = gr_naduzywanie_alkoholu->GetFunction("pol1");
    c->cd(11);
    draw(gr_niezgodnosc_charakterow, "niezgodnosc charakterow", 4);
    gr_niezgodnosc_charakterow->Fit(f3);
    TF1 *fit_niezgodnosc_charakterow = gr_niezgodnosc_charakterow->GetFunction("f3");
    c->cd(12);
    draw(gr_zbiorczy, "zbiorczy", 1);
    gr_zbiorczy->Fit("pol1");
    TF1 *fit_zbiorczy = gr_zbiorczy->GetFunction("pol1");


    int macierz[rozklady_norm.size()][rozklady_norm.size()];
    cout << endl;
    cout << "macierz wynikow testu Chi2" << endl;
    for(int i=0; i<rozklady_norm.size(); i++){  //tworzenie macierzy z wynikami chi2
        for(int j=0; j<rozklady_norm.size(); j++){
            if(i != j){
                if(test_chi2(statystyka_chi2(rozklady_norm[i], rozklady_norm[j]), 0.01, 1) == 0){
                    macierz[i][j] = 0;
                }
                else{
                    macierz[i][j] = 1;
                }
            }
            else{   macierz[i][j] = 1; }
        }
    }

    for(int i=0; i<rozklady_norm.size(); i++){ //wypisywanie macierzy
        for(int j=0; j<rozklady_norm.size(); j++){
            cout << macierz[i][j] << " ";
            if(j==rozklady_norm.size()-1){  cout << " " << names[i] << endl; }
        }
    }
    cout << endl;
    for(int i=0; i<20; i++){ //petla leci linijka po linijce jezeli jest co wypisywac to to robi napisy pod X
        for(int j=0; j<11; j++){
            char ch;
            if( names[j].length() > i ){
                ch = names[j][i];
                cout << ch << " "; 
            }
            else{ cout << "  "; }
        }
        cout << endl;
    }

    cout << endl;
    cout << "wyniki R2 z wykresem zbiorczym (jako funkcja teoretyczna)" << endl;
    for(int i=0; i<rozklady_norm.size(); i++){ //test R2 dla wszystkich rozkladow ze zbiorczym
        cout << names[i] << ": " << R2(rozklady_norm[i], zbiorczy) << endl;
    }

    cout << endl;
    cout << "wyniki R2 z rozkladow z ich dopasowanymi krzywymi" << endl; //test R2 dla wszystkich rozkladow z funkcja do nich dopasowana
    cout << names[0] << ": " << R2_fit(niedobor_seksualny, fit_niedobor_seksualny) << endl;
    cout << names[1] << ": " << R2_fit(dluzsza_nieobecnosc, fit_dluzsza_nieobecnosc) << endl;
    cout << names[2] << ": " << R2_fit(naganny_stosunek_do_czlonkow_rodziny, fit_naganny_stosunek_do_czlonkow_rodziny) << endl;
    cout << names[3] << ": " << R2_fit(niedochowanie_wiernosci_malzenskiej, fit_niedochowanie_wiernosci_malzenskiej) << endl;
    cout << names[4] << ": " << R2_fit(roznice_swiatopogladowe, fit_roznice_swiatopogladowe) << endl;
    cout << names[5] << ": " << R2_fit(hazard, fit_hazard) << endl;
    cout << names[6] << ": " << R2_fit(narkotyki, fit_narkotyki) << endl;
    cout << names[7] << ": " << R2_fit(nieporozumienia_na_tle_finansowym, fit_nieporozumienia_na_tle_finansowym) << endl;
    cout << names[8] << ": " << R2_fit(trudnosci_mieszkaniowe, fit_trudnosci_mieszkaniowe) << endl;
    cout << names[9] << ": " << R2_fit(naduzywanie_alkoholu, fit_naduzywanie_alkoholu) << endl;
    cout << names[10] << ": " << R2_fit(niezgodnosc_charakterow, fit_niezgodnosc_charakterow) << endl;
    cout << "zbiorczy" << ": " << R2_fit(zbiorczy, fit_zbiorczy) << endl;

    cout << endl;
    cout << "wyniki R2 z  krzywa dopasowana do wykresu zbiorczego" << endl;
    for(int i=0; i<rozklady_norm.size(); i++){ //test R2 dla wszystkich rozkladow z funckja dopasowana do zbiorczego
        cout << names[i] << ": " << R2_fit(rozklady_norm[i], fit_zbiorczy) << endl;
    }


    // to samo co dla chi2 ale dla Kolmogorowa Smirnowa
    double macierz_K_S[rozklady_norm.size()][rozklady_norm.size()];
    cout << endl;
    cout << "macierz wynikow testu Kolmogorowa Smirnowa" << endl;
    for(int i=0; i<rozklady.size(); i++){
        for(int j=0; j<rozklady.size(); j++){
            if(i != j){
                if(test_Kolmogorowa_Smirnowa(rozklady[i], rozklady[j]) < 0.673){
                    macierz_K_S[i][j] = 1;
                }
                else{
                    macierz_K_S[i][j] = 0;
                }
            }
            else{ macierz_K_S[i][j] = 1; }
        }
    }
    for(int i=0; i<rozklady.size(); i++){
        for(int j=0; j<rozklady.size(); j++){
            cout << macierz_K_S[i][j] << " ";
            if(j==rozklady.size()-1){  cout << " " << names[i] << endl; }
        }
    }
    cout << endl;
    for(int i=0; i<20; i++){
        for(int j=0; j<11; j++){
            char ch;
            if( names[j].length() > i ){
                ch = names[j][i];
                cout << ch << " ";
            }
            else{ cout << "  "; }
        }
        cout << endl;
    }

return 0;
}