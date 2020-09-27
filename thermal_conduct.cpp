// Created by Siddarth Achar - PhD student : University of Pittsburgh
// Project collaborated with Meiirbek Islamov (Phd student : University of Pittsburgh)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <stdlib.h>
// #include <stdio.h>
using namespace std;

ifstream data_in("Heat_flux-100K.txt");
ofstream data_out1("hcacf.txt");
ofstream data_out2("k_data.txt");


int main()
{
    double dt = 0.002*5;
    double L = 7.42;
    int N = 200001;
    int M = 2000; //correlation length
    double T_d = 100/121.0;
    double V = L*L*L;
    double heatf[N][3];
    double heatx, heaty, heatz;
    double hcacfx[M], hcacfy[M], hcacfz[M]; // Heat current autocorrelation function
    double t[M]; //correlation time
    double hcacf; // average heat current autocorrelation function
    double k_b = 1; //Boltzmann's constant
    double kt;
    double hh = 0;


    for(int i=0; i<N; ++i){
        for(int j=0; j<3;++j)
        {
            data_in >> heatf[i][j];
        }
    }
 
    for (int m=0; m<M; m++){
        
        heatx = 0;
        heaty = 0;
        heatz = 0;
        for (int n=0; n<(N-m); n++){
            
            heatx = heatx + heatf[n][0] * heatf[m+n][0];
            heaty = heaty + heatf[n][1] * heatf[m+n][1];
            heatz = heatz + heatf[n][2] * heatf[m+n][2];
            
        }
        hh+=heatx+heaty+heatz;
        
        hcacfx[m] = heatx/(N-m);
        hcacfy[m] = heaty/(N-m);
        hcacfz[m] = heatz/(N-m);
        
        data_out1<<m<<"    "<<hcacfx[m]<<"    "<<hcacfy[m]<<"    "<<hcacfz[m]<<endl;

    }
    
    // Thermal conductivity calculation
    
    
    for (int i=0; i<M; ++i){ // correlation time list
        t[i] = i*dt*2.14;
    }
    
    for (int i=0; i<M; ++i){
        
        hcacf = (hcacfx[i]+hcacfy[i]+hcacfz[i])/3;
        
        if (i == 0) {
            
            kt =  hcacf / 2 / V * dt / (k_b * T_d * T_d)*0.0190;
            data_out2<<t[i]<<"    "<<hcacf<<"    "<<kt<<endl;
        }
        
        else{
            
            kt =  kt + hcacf/V * dt / (k_b * T_d * T_d)*0.0190;
            data_out2<<t[i]<<"    "<<hcacf<<"    "<<kt<<endl;
        }
    }
    
    
    

    // LJ thermal conductivity scale - 0.0189, 20K-1.2W/mk
    
  


    
}

