// C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <math.h>

// ROOT libraries
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>

using namespace std;
void Ajuste_Tr(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    string Direccion_ajustes = "Parametros_Traza_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
        A_fix.push_back(A);
        B_fix.push_back(B);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    char Traza_name[Num_N][50]; 
    char Traza_Name_histo[Num_N][80];
    TH1F* Traza_hist[Num_N];
    for (int i = 0; i < Num_N; i++){
        sprintf(Traza_name[i], "Traza_%i", Particulas[i]);
        sprintf(Traza_Name_histo[i], "Traza_%i ; Residuo ; Counts", Particulas[i]);
        Traza_hist[i]= new TH1F(Traza_name[i], Traza_Name_histo[i], Particulas[i]/2, -0.25, 0.25);
    }
    // Open the file
      for (int i = 0; i < Num_N; i++){
            string filename = "Scaling/ETH_Hopping_";
            filename += to_string(Particulas[i]) + ".txt"; 
            ifstream infile(filename);
            if (!infile.is_open()) {
                std::cerr << "Unable to open file '" << filename << "'" << std::endl;
            return;
            }
            std::string line;
            while (std::getline(infile, line)) {
                std::stringstream iss(line);
                double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
                iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
                if (eid>-3.3 && eid<-2.3) {
                    try {
                        //cout << "eid: " << eid << "traza" << traza2<< endl;
                        double Residuo = eid*A_fix[i]+B_fix[i]-traza3; 
                        Traza_hist[i]->Fill(Residuo);
                    } catch (const std::exception& e) {
                        std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                    }
                }
            }
    }
}

void Ajuste_Autoval_Corr(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    char Autoval_name[Num_N][50]; 
    char Autoval_Name_histo[Num_N][80];
    TH1F* Autoval_hist[Num_N];
    for (int i = 0; i < Num_N; i++){
        sprintf(Autoval_name[i], "Autovalor_%i", Particulas[i]);
        sprintf(Autoval_Name_histo[i], "Autovalor_%i ; Residuo ; Counts", Particulas[i]);
        Autoval_hist[i]= new TH1F(Autoval_name[i], Autoval_Name_histo[i], Particulas[i]/2, -0.25, 0.25);
    }
// Open the file
      for (int i = 0; i < Num_N; i++){
            string filename = "Scaling/ETH_Corriente_";
            filename += to_string(Particulas[i]) + ".txt"; 
            ifstream infile(filename);
            if (!infile.is_open()) {
                std::cerr << "Unable to open file '" << filename << "'" << std::endl;
            return;
            }
            std::string line;
            while (std::getline(infile, line)) {
                std::stringstream iss(line);
                double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
                iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
                if (eid>-3.3 && eid<-2.3) {
                    try {
                        //cout << "eid: " << eid << "traza" << traza2<< endl;
                        Autoval_hist[i]->Fill(eig1_2);
                    } catch (const std::exception& e) {
                        std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                    }
                }
            }
    }
}

void Ajuste_Autoval_Rot(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    char Autoval_name[Num_N][3][50]; 
    char Autoval_Name_histo[Num_N][3][80];
    TH1F* Autoval_hist[Num_N][3];
    for (int i = 0; i < Num_N; i++){
        for (int j = 0; j < 3; j++){
            sprintf(Autoval_name[i][j], "Autovalor_%i_N_%i",j, Particulas[i]);
            sprintf(Autoval_Name_histo[i][j], "Autovalor_%i_%i ; Residuo ; Counts",j, Particulas[i]);
            Autoval_hist[i][j] = new TH1F(Autoval_name[i][j], Autoval_Name_histo[i][j], Particulas[i], -0.15, 0.15);
        }
    }
     // Open the file
      for (int i = 0; i < Num_N; i++){
            string filename = "Scaling/ETH_Rotacion_";
            filename += to_string(Particulas[i]) + ".txt"; 
            ifstream infile(filename);
            if (!infile.is_open()) {
                std::cerr << "Unable to open file '" << filename << "'" << std::endl;
            return;
            }
            std::string line;
            while (std::getline(infile, line)) {
                std::stringstream iss(line);
                double eid,eR,traza2,re_eig1_2,re_eig2_2, im_eig1_2, im_eig2_2, traza3, re_eig1_3,re_eig2_3,re_eig3_3, im_eig1_3, im_eig2_3, im_eig3_3, num;
                iss >> eid >> eR >> traza2 >> re_eig1_2 >> re_eig2_2 >>  im_eig1_2 >> im_eig2_2   >> traza3 >> re_eig1_3 >> re_eig2_3 >> re_eig3_3 >>  im_eig1_3 >> im_eig2_3 >> im_eig3_3 >> num;            
                if (eid>-3.2 && eid<-2.6) {
                    try {
                        //cout << "eid: " << eid << "traza" << traza2<< endl;
                        Autoval_hist[i][0]->Fill(re_eig1_3);
                        Autoval_hist[i][1]->Fill(re_eig2_3);
                        Autoval_hist[i][2]->Fill(re_eig3_3);
                    } catch (const std::exception& e) {
                        std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                    }
                }
            }
    }
}

void Ordenar(string Operador){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Ajustes_Parametros_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    int Num_N = Particulas.size();
    for (int i=0; i < Num_N; i++){
        string filename = "Scaling/ETH_" + Operador+ "_";
        string newfilename = filename+ "new";
        filename += to_string(Particulas[i]) + ".txt"; 
        newfilename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        ofstream outfile(newfilename);
        std::string line;
        if (Operador != "Rotacion"){
            while (std::getline(infile, line)) {
                std::stringstream iss(line);
                double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
                iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
                if (num < 2.5){
                    outfile << eid << eR << traza2 << eig1_2 << eig2_2 << traza3 << eig3_3 << eig1_3 << eig2_3  << std::endl;
                }
                else{
                    //Movidas que comprobar;
                    if (abs(eig1_2) >= abs(eig2_2)){
                        double maximo = eig1_2; 
                        if (maximo != eig3_3 && maximo != eig2_3 && maximo != eig1_3){
                        cout <<"N:" << Particulas[i] << " eid:" << eid << "maximo en 2 algo raro" << endl;
                    }
                    }
                    else{
                        double maximo = eig2_2;
                        if (maximo != eig3_3 && maximo != eig2_3 && maximo != eig1_3){
                        cout <<"N:" << Particulas[i] << " eid:" << eid << "maximo en 2 algo raro" << endl;

                    }
                        //cout <<"N:" << Particulas[i] << " eid:" << eid << "maximo en 2 algo raro" << endl;
                    }
                    
                    outfile << eid << eR << traza2 << eig1_2 << eig2_2 << traza3 << eig3_3 << eig1_3 << eig2_3  << std::endl;
                }
            }
        }
    }

}


void Mean_Rot(TFile* outputFile, double mean_re[9], double mean_im[9]){
    //Meto N A y B
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        string filename = "Scaling/ETH_Rotacion_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid,eR,traza2,re_eig1_2,re_eig2_2, im_eig1_2, im_eig2_2, traza3, re_eig1_3,re_eig2_3,re_eig3_3, im_eig1_3, im_eig2_3, im_eig3_3, num;
            iss >> eid >> eR >> traza2 >> re_eig1_2 >> re_eig2_2 >>  im_eig1_2 >> im_eig2_2   >> traza3 >> re_eig1_3 >> re_eig2_3 >> re_eig3_3 >>  im_eig1_3 >> im_eig2_3 >> im_eig3_3 >> num;            
            if (eid>-3 && eid<-2.6) {
                Eventos_zona += 3;
                try {
                    mean_re[i] += re_eig1_3 + re_eig2_3 + re_eig3_3;
                    mean_im[i] += im_eig1_3 + im_eig2_3 + im_eig3_3;
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                   
                } catch (const std::exception& e) {
                    std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        mean_re[i] = mean_re[i]/Eventos_zona;
        mean_im[i] = mean_im[i]/Eventos_zona;
        cout << "N: " << Particulas[i] << " Media real " << mean_re[i] << " Media imag " << mean_im[i] << endl;
    }
}

void Sigma_Rot(TFile*outputFile, double mean_re[9], double mean_im[9]){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    double sigma_re[Num_N];
    double sigma_im[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        string filename = "Scaling/ETH_Rotacion_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid,eR,traza2,re_eig1_2,re_eig2_2, im_eig1_2, im_eig2_2, traza3, re_eig1_3,re_eig2_3,re_eig3_3, im_eig1_3, im_eig2_3, im_eig3_3, num;
            iss >> eid >> eR >> traza2 >> re_eig1_2 >> re_eig2_2 >>  im_eig1_2 >> im_eig2_2   >> traza3 >> re_eig1_3 >> re_eig2_3 >> re_eig3_3 >>  im_eig1_3 >> im_eig2_3 >> im_eig3_3 >> num;            
            if (eid>-3.1 && eid<-2.8) {
                Eventos_zona += 3;
                try {
                    sigma_re[i] += TMath::Power(re_eig1_3-mean_re[i],2)+TMath::Power(re_eig2_3-mean_re[i],2)+TMath::Power(re_eig3_3-mean_re[i],2);
                    sigma_im[i] += TMath::Power(im_eig1_3-mean_im[i],2)+TMath::Power(im_eig2_3-mean_im[i],2)+TMath::Power(im_eig3_3-mean_im[i],2);
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                    //cout <<eid <<" " <<Particulas[i] << " " << TMath::Power(re_eig1_3-mean_re[i],2) << " " <<TMath::Power(im_eig1_3-mean_im[i],2) << " " << TMath::Power(re_eig2_3-mean_re[i],2) << " " << TMath::Power(im_eig2_3-mean_im[i],2)<< " "<< TMath::Power(re_eig3_3-mean_re[i],2)<< " "<< TMath::Power(im_eig3_3-mean_im[i],2) <<" "<< sigma_re[i] << " " << sigma_im[i] << endl;
                } catch (const std::exception& e) {
                    std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        sigma_re[i] = sqrt(sigma_re[i]/Eventos_zona);
        sigma_im[i] = sqrt(sigma_im[i]/Eventos_zona);
        cout << "N: " << Particulas[i] << " Sigma real " << sigma_re[i] << " Sigma imag " << sigma_im[i] << endl;
    }
}


void Radio_Rot(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        string filename = "Scaling/ETH_Rotacion_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        string outpurfilename = "Orden2_"+to_string(Particulas[i])+".dat";
        std::ofstream outputRadio(outpurfilename);
        if (!outputRadio.is_open()) {
            std::cerr << "Unable to open file for calibration output: " << outpurfilename << std::endl;
            return;
        }
        double Autovalor[2];
        double DiffAutoval;
        int Paso = 0;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid,eR,traza2,re_eig1_2,re_eig2_2, im_eig1_2, im_eig2_2, traza3, re_eig1_3,re_eig2_3,re_eig3_3, im_eig1_3, im_eig2_3, im_eig3_3, num;
            iss >> eid >> eR >> traza2 >> re_eig1_2 >> re_eig2_2 >>  im_eig1_2 >> im_eig2_2   >> traza3 >> re_eig1_3 >> re_eig2_3 >> re_eig3_3 >>  im_eig1_3 >> im_eig2_3 >> im_eig3_3 >> num;            
            double Valores[7];
            Valores[0] = eid;
            
            if (abs(im_eig1_3) < 0.0001){
                Valores[1] = re_eig1_3;
                Valores[2] = im_eig1_3;
                Valores[3] = re_eig2_3;
                Valores[4] = abs(im_eig2_3);
                Valores[5] = re_eig3_3;
                Valores[6] = -abs(im_eig3_3);
                //outputRadio << eid << " " << re_eig1_3<< endl;
            } else if (abs(im_eig2_3) < 0.0001){
                Valores[1] = re_eig2_3;
                Valores[2] = im_eig2_3;
                Valores[3] = re_eig1_3;
                Valores[4] = abs(im_eig1_3);
                Valores[5] = re_eig3_3;
                Valores[6] = -abs(im_eig3_3);
                //outputRadio << eid << " " << re_eig2_3<<endl;
            } else if(abs(im_eig3_3) < 0.0001){
                Valores[1] = re_eig3_3;
                Valores[2] = im_eig3_3;
                Valores[3] = re_eig2_3;
                Valores[4] = abs(im_eig2_3);
                Valores[5] = re_eig1_3;
                Valores[6] = -abs(im_eig1_3);
                //outputRadio << eid << " " << re_eig3_3<< endl;
                }
            
            outputRadio << Valores[0] << " " << Valores[1]  << " " << Valores[2]  << " " << Valores[3]  << " " << Valores[4] << " " << Valores[5] << " " << Valores[6] << endl; 
        }
    }
}

void Positivos_Corriente(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        string filename = "Scaling/ETH_Corriente_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        string outpurfilename = "Corriente_Positivos_"+to_string(Particulas[i])+".dat";
        std::ofstream outputRadio(outpurfilename);
        if (!outputRadio.is_open()) {
            std::cerr << "Unable to open file for calibration output: " << outpurfilename << std::endl;
            return;
        }
        double Autovalor[2];
        double DiffAutoval;
        int Paso = 0;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid, eR, traza2, eig1_2, eig2_2, traza3, eig1_3, eig2_3, eig3_3, num;
            iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2  >> traza3 >> eig1_3 >> eig2_3 >>eig3_3  >> num;            
            outputRadio << eid << " " << abs(eig1_2) << " " << -abs(eig1_2) << "\n";
        }
    }
}

void Sigma_Corriente(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    string outpurfilename = "Sigma_Corriente_Centro.dat";
    std::ofstream outputRadio(outpurfilename);
    if (!outputRadio.is_open()) {
         std::cerr << "Unable to open file for calibration output: " << outpurfilename << std::endl;
         return;
    }

    string MaximoFilename = "Maximo_Corriente_Centro.dat";
    std::ofstream outputMaximo(MaximoFilename);
    if (!outputMaximo.is_open()) {
         std::cerr << "Unable to open file for calibration output: " << MaximoFilename << std::endl;
         return;
    }
     //Defino los histogramas
    int Num_N = Particulas.size();
    double sigma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        Maximo[i] = 0;
        string filename = "Corriente_Positivos_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        double Energia;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid, eig1_2, eig2_2;
            iss >> eid >>  eig1_2 >> eig2_2; 
            if (-3.1 < eid && eid < -2.3){
                 sigma[i] += eig1_2*eig1_2;
                 Eventos_zona += 1;
                 if (eig1_2 > Maximo[i]){
                     Maximo[i] = eig1_2; 
                     Energia = eid;
                 }
                 //cout << Particulas[i] << " " << eid << " " <<   eig1_2*eig1_2 << "\n";
            } 
        }
        sigma[i] = sqrt(sigma[i]/Eventos_zona);
        outputRadio << Particulas[i] << " " << sigma[i] << "\n";
        outputMaximo << Particulas[i] << " " << Energia << " " << Maximo[i] << "\n";
        cout << Particulas[i] << " Sigma " << sigma[i] << " E_Max " << Energia << " Max " << Maximo[i] << " Numero de eventos " << Eventos_zona << "\n";
    }
}


void Sigma_Corriente_Histo(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    string outpurfilename = "Sigma_Corriente_Centro_Histo.dat";
    std::ofstream outputRadio(outpurfilename);
    if (!outputRadio.is_open()) {
         std::cerr << "Unable to open file for calibration output: " << outpurfilename << std::endl;
         return;
    }

     //Defino los histogramas
    int Num_N = Particulas.size();
    int Bineado[9] = {4,5,6,9,11,14,18,23,29};
    char Corriente_name[Num_N][50]; 
    char Corriente_Name_histo[Num_N][80];
    TH1F* Corriente_hist[Num_N];
    for (int i = 0; i < Num_N; i++){
        sprintf(Corriente_name[i], "Autovelores_%i", Particulas[i]);
        sprintf(Corriente_Name_histo[i], "Autovalores_%i ; Residuo ; Counts", Particulas[i]);
        Corriente_hist[i]= new TH1F(Corriente_name[i], Corriente_Name_histo[i], Particulas[i]/2, -0.2, 0.2);
    }
    // Open the file
    for (int i = 0; i < Num_N; i++){
        cout << Particulas[i] << "\n";
        string filename = "Corriente_Positivos_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        double Energia;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid, eig1_2, eig2_2;
            iss >> eid >>  eig1_2 >> eig2_2; 
            if (-3.1 < eid && eid < -2.8){
                 Corriente_hist[i]-> Fill(eig1_2);
                 cout << eig1_2 << "\n";
                 cout << eig2_2 << "\n";
                 Corriente_hist[i]-> Fill(eig2_2);
            } 
        }
        //sigma[i] = sqrt(sigma[i]/Eventos_zona);
        //outputRadio << Particulas[i] << " " << sigma[i] << "\n";
        //cout << Particulas[i] << " Sigma " << sigma[i]  << " Numero de eventos " << Eventos_zona << "\n";
    }
}

void Ajuste_Traza_Hopping(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    string Direccion_ajustes = "Parametros_Traza_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
        A_fix.push_back(A);
        B_fix.push_back(B);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Num_eventos = 0;
        double media =0;
        string filename = "Scaling/ETH_Hopping_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)){
            std::stringstream iss(line);
            double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
            iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
            if (eid>-3.1 && eid<-2.8) {
                try {
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                     //cout << "Hola" << "\n";
                     Num_eventos += 1;
                     double Residuo = eid*A_fix[i]+B_fix[i]-traza3; 
                     media += Residuo;
                } catch (const std::exception& e) {
                      std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        media = media/Num_eventos;
        cout << Particulas[i] << " " << media << "\n";
    }
}

void Ajuste_Traza_Hopping_sigma(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    std::vector<double> Mean;
    string Direccion_ajustes = "Parametros_Traza_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
        A_fix.push_back(A);
        B_fix.push_back(B);
    }
    infile_ajuste.close();

    string Direccion_media = "Media_Traza_Hopping.dat";   
    ifstream infile_media(Direccion_media);
    if (!infile_media.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_media << "'" << std::endl;
        return;
    }
    cout << Direccion_media << endl;
    std::string line2;
    while (std::getline(infile_media, line2)) {
        std::stringstream iss(line2);
        double N, media;
        iss >> N >> media;        
        Mean.push_back(media);
        //cout << N <<" " << media << endl;
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    double sigma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        sigma[i] = 0;
        Maximo[i] = 0;
        double  Energia;
        int Eventos_zona = 0;
        string filename = "Scaling/ETH_Hopping_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)){
            std::stringstream iss(line);
            double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
            iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
            if (eid>-3.1 && eid<-2.8) {
                try {
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                     //cout << "Hola" << "\n";
                     Eventos_zona += 1;
                     double Residuo = eid*A_fix[i]+B_fix[i]-traza3; 
                     sigma[i] += (Residuo-Mean[i])*(Residuo-Mean[i]);
                     if (abs(Residuo)>Maximo[i]){
                         Maximo[i] = abs(Residuo);
                         Energia = eid;
                     }
                } catch (const std::exception& e) {
                      std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        sigma[i] = sqrt(sigma[i]/Eventos_zona);
       // cout << Particulas[i] << " " << sigma[i] << "\n";
        cout << Particulas[i] << "  " << Energia<< " " << Maximo[i] << "\n";
    }
}

void Ordenar_Autoval_Hopping(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    std::vector<double> Mean;
    string Direccion_ajustes = "Parametros_Traza_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    infile_ajuste.close();
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0;
        string filename = "Scaling/ETH_Hopping_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        string filenameOrdenado = "ETH_Hopping_Autoval_Full_";
        filenameOrdenado += to_string(Particulas[i]) + ".dat";
        std::ofstream outputAutoval(filenameOrdenado);
        if (!outputAutoval.is_open()) {
            std::cerr << "Unable to open file for calibration output: " << filenameOrdenado << std::endl;
        return;
        }

        std::string line;
        while (std::getline(infile, line)){
            std::stringstream iss(line);
            double eid,eR,traza2,eig1_2,eig2_2,traza3,eig1_3,eig2_3,eig3_3, num;
            iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3>> num;            
            //if (eid>-3.1 && eid<-2.8) {
                try {
                    outputAutoval << eid << " " << eig1_3 << "\n";
                    outputAutoval << eid << " " << eig2_3 << "\n";
                    outputAutoval << eid << " " << eig3_3 << "\n";
                } catch (const std::exception& e) {
                      std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            //}
        }
    }
}

void Ajuste_Autoval_Hopping(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    string Direccion_ajustes = "Parametros_Autovalores_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
        A_fix.push_back(A);
        B_fix.push_back(B);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    // Open the file
    for (int i = 0; i < Num_N; i++){
        double media =0;
        int  Eventos_zona = 0;
        string filename = "ETH_Hopping_Autoval_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)){
            std::stringstream iss(line);
            double eid, eig;
            iss >> eid >> eig;            
            if (eid>-3.1 && eid<-2.3) {
                try {
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                     //cout << "Hola" << "\n";
                     Eventos_zona += 1;
                     double Residuo = eid*A_fix[i]+B_fix[i]-eig; 
                     media += Residuo;
                } catch (const std::exception& e) {
                      std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        cout << Eventos_zona << "\n";
        cout << Particulas[i] << " " << media << "\n";
    }
}

void Ajuste_Autoval_Hopping_sigma(TFile* outputFile){
    //Meto N A y B
    std::vector<int> Particulas;
    std::vector<double> A_fix;
    std::vector<double> B_fix;
    std::vector<double> Mean;
    string Direccion_ajustes = "Parametros_Autovalores_Hopping.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
        A_fix.push_back(A);
        B_fix.push_back(B);
    }
    infile_ajuste.close();

    string Direccion_media = "Media_Autovalores_Hopping.dat";   
    ifstream infile_media(Direccion_media);
    if (!infile_media.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_media << "'" << std::endl;
        return;
    }
    cout << Direccion_media << endl;
    std::string line2;
    while (std::getline(infile_media, line2)) {
        std::stringstream iss(line2);
        double N, media;
        iss >> N >> media;        
        Mean.push_back(media);
        cout << N <<" " << media << endl;
    }

    //Defino los histogramas
    int Num_N = Particulas.size();
    double sigma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        sigma[i] = 0;
        Maximo[i] = 0;
        double Energia = 0;
        int Eventos_zona = 0;
        string filename = "ETH_Hopping_Autoval_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)){
            std::stringstream iss(line);
            double eid, eig;
            iss >> eid >> eig;           
            if (eid>-3.1 && eid<-2.8) {
                try {
                    //cout << "eid: " << eid << "traza" << traza2<< endl;
                     //cout << "Hola" << "\n";
                     Eventos_zona += 1;
                     double Residuo = eid*A_fix[i]+B_fix[i]-eig; 
                     sigma[i] += (Residuo-Mean[i])*(Residuo-Mean[i]);
                     cout << Residuo << endl;
                     if (abs(Residuo) > Maximo[i]){
                         Maximo[i] = abs(Residuo); 
                         //cout << Maximo[i] << "\n";
                         Energia = eid;
                     }
                } catch (const std::exception& e) {
                      std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        //cout << "Eventos: " << Eventos_zona << "\n";
        //sigma[i] = sqrt(sigma[i]/Eventos_zona);
        //cout << Particulas[i] << " " << sigma[i] << "\n";
        cout << Particulas[i] << "  " << Energia<< " " << Maximo[i] << "\n";
    }
}

void Sigma_Rot_Radio(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    //Defino los histogramas
    int Num_N = Particulas.size();
    double sigma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Eventos_zona = 0; 
        Maximo[i] = 0;
        sigma[i] = 0;

        double Energia;
        string filename = "Orden_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid, radio, im_radio, re_eig_1, im_eig_1, re_eig_2, im_eig_2;
            iss >> eid >> radio >> im_radio >> re_eig_1 >> im_eig_1 >> re_eig_2 >> im_eig_2;
            if (eid>-3.1 && eid<-2.3) {
                Eventos_zona += 1;
                try {
                    sigma[i] += radio*radio;
                    if (abs(radio)> Maximo[i]){
                        Maximo[i] = abs(radio);
                        Energia =eid;
                   }
                } catch (const std::exception& e) {
                    std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        sigma[i] = sqrt(sigma[i]/Eventos_zona);
        cout << " " << Particulas[i] << " "<< sigma[i] << endl;
        //cout << " " << Particulas[i] << " "<< Energia <<" " <<Maximo[i] << endl;
    }
}


void Traza_Densidad(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    infile_ajuste.close();

    string Direccion_densidad = "micro_hopping.dat";   
    ifstream infile_densidad(Direccion_densidad);
    if (!infile_densidad.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_densidad << "'" << std::endl;
        return;
    }
    std::vector<double> Micro;
    std::vector<double> Densidad;
    std::vector<double> Energia_Micro;
    std::string line_densidad;
    while (std::getline(infile_densidad, line_densidad)) {
        std::stringstream iss(line_densidad);
        double e, rho, Micro_hopping;
        iss >> e >> rho >> Micro_hopping;
        Energia_Micro.push_back(e);
        Densidad.push_back(rho);
        Micro.push_back(Micro_hopping);
    }

    int Num_N = Particulas.size();
    double sigma[Num_N];
    double media[Num_N];
    double suma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Paso_Micro = 0;
        int Eventos_zona = 0; 
        Maximo[i] = 0;
        sigma[i] = 0;
        media[i] = 0;
        suma[i] = 0;
        double Energia;
        string filename = "Scaling/ETH_Hopping_";
        filename += to_string(Particulas[i]) + ".txt"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid, eR, traza2, eig1_2, eig2_2, traza3, eig1_3, eig2_3, eig3_3, num;
            iss >> eid >> eR >> traza2 >> eig1_2 >> eig2_2 >> traza3 >> eig1_3 >> eig2_3 >> eig3_3 >> num;
            if (eid>-3.1 && eid<-2.3) {
                Eventos_zona += 1;
                try {
                    while (Energia_Micro[Paso_Micro] < eid ){
                        Paso_Micro +=1;
                    }
                    double Porcentaje_1,Porcentaje_2;
                    Porcentaje_1 = 1-abs(eid-Energia_Micro[Paso_Micro])*1000;
                    Porcentaje_2 = 1-abs(eid-Energia_Micro[Paso_Micro-1])*1000;
                    double Interpolado = Porcentaje_1*Micro[Paso_Micro]+Porcentaje_2*Micro[Paso_Micro-1];
                    double Residuo = 3*Interpolado - traza3;
                    suma[i] += Residuo;
                    //cout << "Particula: " << Particulas[i] << " Energia-id " << eid << " Energia_Micro " << Energia_Micro[Paso_Micro] << " Intervalo 1 " << Porcentaje_1 <<" Energia anterior" <<  Energia_Micro[Paso_Micro-1] << " Intervalo 2 " << Porcentaje_2 << endl; 
                    sigma[i] += Residuo*Residuo;
                    if (abs(Residuo)> Maximo[i]){
                        Maximo[i] = abs(Residuo);
                        Energia =eid;
                   }
                } catch (const std::exception& e) {
                    std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        media[i] = suma[i]/Eventos_zona;
        sigma[i] = sqrt(sigma[i]/Eventos_zona-media[i]*media[i]);
        cout << " " << Particulas[i]  << " "<<sigma[i] << endl;
        //cout << " " << Particulas[i] << " "<< Energia <<" " <<Maximo[i] << endl;
    }
}

void Autovalores_Densidad(TFile*outputFile){
    std::vector<int> Particulas;
    string Direccion_ajustes = "Parametros_Traza_Hopping_3.txt";   
    ifstream infile_ajuste(Direccion_ajustes);
    if (!infile_ajuste.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_ajustes << "'" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile_ajuste, line)) {
        std::stringstream iss(line);
        double N, A, B;
        iss >> N >> A >> B;        
        Particulas.push_back(N);
    }
    infile_ajuste.close();

    string Direccion_densidad = "micro_hopping.dat";   
    ifstream infile_densidad(Direccion_densidad);
    if (!infile_densidad.is_open()) {
            std::cerr << "Unable to open file '" << Direccion_densidad << "'" << std::endl;
        return;
    }
    std::vector<double> Micro;
    std::vector<double> Densidad;
    std::vector<double> Energia_Micro;
    std::string line_densidad;
    while (std::getline(infile_densidad, line_densidad)) {
        std::stringstream iss(line_densidad);
        double e, rho, Micro_hopping;
        iss >> e >> rho >> Micro_hopping;
        Energia_Micro.push_back(e);
        Densidad.push_back(rho);
        Micro.push_back(Micro_hopping);
    }

    int Num_N = Particulas.size();
    double sigma[Num_N];
    double media[Num_N];
    double suma[Num_N];
    double Maximo[Num_N];
    // Open the file
    for (int i = 0; i < Num_N; i++){
        int Paso_Micro = 0;
        int Eventos_zona = 0; 
        Maximo[i] = 0;
        sigma[i] = 0;
        media[i] = 0;
        suma[i] = 0;
        double Energia;
        string filename = "ETH_Hopping_Autoval_Full_";
        filename += to_string(Particulas[i]) + ".dat"; 
        ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Unable to open file '" << filename << "'" << std::endl;
        return;
        }
        std::string line;
        while (std::getline(infile, line)) {
            std::stringstream iss(line);
            double eid,eig;
            iss >> eid >> eig;
            if (eid>-3.1 && eid<-2.3) {
                Eventos_zona += 1;
                try {
                    while (Energia_Micro[Paso_Micro] < eid ){
                        Paso_Micro +=1;
                    }
                    double Porcentaje_1,Porcentaje_2;
                    Porcentaje_1 = 1-abs(eid-Energia_Micro[Paso_Micro])*1000;
                    Porcentaje_2 = 1-abs(eid-Energia_Micro[Paso_Micro-1])*1000;
                    double Interpolado = Porcentaje_1*Micro[Paso_Micro]+Porcentaje_2*Micro[Paso_Micro-1];
                    double Residuo = Interpolado - eig;
                    suma[i] += Residuo;
                    //cout << "Particula: " << Particulas[i] << " Energia-id " << eid << " Energia_Micro " << Energia_Micro[Paso_Micro] << " Intervalo 1 " << Porcentaje_1 <<" Energia anterior" <<  Energia_Micro[Paso_Micro-1] << " Intervalo 2 " << Porcentaje_2 << endl; 
                    sigma[i] += Residuo*Residuo;
                    if (abs(Residuo)> Maximo[i]){
                        Maximo[i] = abs(Residuo);
                        Energia = eid;
                   }
                } catch (const std::exception& e) {
                    std::cerr << "Error extracting channel or strip information: " << e.what() << std::endl;
                }
            }
        }
        media[i] = suma[i]/Eventos_zona;
        sigma[i] = sqrt(sigma[i]/Eventos_zona-media[i]*media[i]);
        cout << " " << Particulas[i]  <<" " << sigma[i] << endl;
        //cout << " " << Particulas[i] << " "<< Energia <<" " <<Maximo[i] << endl;
    }
}


int Ajuste_Traza() {
    //string Operador = "Hopping";
    //Ordenar(Operador);
    string Name = "Scaling";
    std::string outputFileName = "./Analized/"+Name+"-Analized.root";
    TFile* outputFile = new TFile(outputFileName.c_str(), "recreate");
    //Autovalores_Densidad(outputFile);
    //Traza_Densidad(outputFile);
    //Sigma_Rot_Radio(outputFile);
    //Ajuste_Autoval_Hopping_sigma(outputFile);
    //Ajuste_Traza_Hopping_sigma(outputFile);
    //Ordenar_Autoval_Hopping(outputFile);
   //Sigma_Corriente(outputFile);
    //Sigma_Corriente_Histo(outputFile);
    //Positivos_Corriente(outputFile);
   //Ajuste_Tr(outputFile);
    //Ajuste_Autoval_Corr(outputFile);
    //Ajuste_Autoval_Rot(outputFile);
    //new TBrowser;
    /*
    double mean_1[9] = {0};
    double mean_2[9] = {0};
    Mean_Rot(outputFile, mean_1, mean_2);
    Sigma_Rot(outputFile, mean_1, mean_2);*/
    Radio_Rot(outputFile);
    return 0;
}

