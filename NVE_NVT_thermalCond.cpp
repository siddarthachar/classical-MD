// Created by Siddarth Achar - PhD student : University of Pittsburgh
// Project collaborated with Meiirbek Islamov (Phd student : University of Pittsburgh)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// ofstream f_ke("Kinetic.txt");
// ofstream f_pe("Potential.txt");
// ofstream f_temp("Temperature.txt");
// ofstream f_pres("Pressure.txt");
// ofstream f_tot("TotalE.txt");
// ofstream f_mom("Mom.txt");
// ofstream f_tempd("Dimless_Temp.txt");
ofstream f_heatf("Heat_flux.txt");
// ofstream f_pos("Last_Position.txt");
ifstream data_in("liquid256.txt");
int N = 256;
double Tset = 90;

double Rc = 2.5;
double Rc2 = Rc * Rc;
double Rc2_1 = 1/Rc2;
double Rc6_1 = Rc2_1 * Rc2_1 * Rc2_1; //1/rij^6
double Rc12_1 = Rc6_1 * Rc6_1; //1/rij^12

double URc = 4*(Rc12_1 - Rc6_1);
double dURc = (24/Rc)*(-2*Rc12_1 + Rc6_1);
double preFijRc = -(24/Rc)*(2*Rc12_1 - Rc6_1);

double L = 7.14; //box length

double V = L*L*L;
double NearestImage(double rtemp){
    if(rtemp >= L/2){
        rtemp -= L;
    }else if(rtemp < -L/2){
        rtemp += L;
    }
    return rtemp;
}

double InsideCell(double rtemp){
    if (rtemp <= 0){
        rtemp += L;
    } else if (rtemp > L){
        rtemp -= L;
    }
    return rtemp;
}


// MAIN -------------- MAIN ------------- MAIN ------------- MAIN ----------- MAIN
int main(){
    double Fcurr[N][3];
    
    int m, n;
    for (m = 0; m<N; m++)
    {
        for (n = 0; n<3; n++)
            Fcurr[m][n] = 0.0;
    }
    double rcurr[N][3], vcurr[N][3], Ucurr = 0, Fij[3], Ui_arr[N],Ui, avgTemp = 0, Temp, Vir, Pres, avgPres = 0, Par = 4.24*pow(10.0,7.0), Tempd, jterm2[3],jterm1[3],totei, jterm2_temp[3],heat_flux[3],temp_j2, KEi;
    double vix, viy, viz, vjx, vjy, vjz, rij, rijx, rijy, rijz, rij2, r2_1, r6_1, r12_1, Uij, preFij, Utot = 0; // For the force loop
    double vprev, rprev; // For vel verlet
    double eta_prev = 0, eta_curr = 0;

    for(int i=0; i<N; ++i){
        for(int j=0; j<3;++j)
            data_in >> rcurr[i][j];
    }
    // velocity initialization - change to random initialization
    
    double a = -1, b = 1;
      double step = (b-a)/RAND_MAX;
      double sum[3] = {0};

      for(int i = 0; i < N; ++i){
    for(int j = 0; j<3; ++j){
      // printf("%g\n", a + rand()*step);
      vcurr[i][j] = a + rand()*step;
      sum[j] += vcurr[i][j];
  }}

      double sumN[3] ;
      sumN[0] = sum[0]/N;
      sumN[1] = sum[1]/N;
      sumN[2] = sum[2]/N;
      double sum_sq = 0;

      for(int i = 0; i < N; ++i)
        for(int j = 0; j<3; ++j){
              vcurr[i][j] = vcurr[i][j] - sumN[j];
             sum_sq += vcurr[i][j] * vcurr[i][j];}

     double KE_scale = sqrt(2*320/sum_sq), ke_now = 0;
     // printf("%g\n", KE_scale);

     for(int i = 0; i < N; ++i){
    for(int j = 0; j<3; ++j){
      vcurr[i][j] = KE_scale*vcurr[i][j];
          ke_now += vcurr[i][j] * vcurr[i][j];}}

printf("-------------------------------------lslsls------%g\n", ke_now/2);

// For the first time step ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    Utot = 0;
    Vir = 0;
    // Calculate initial Force for position rcurr
    for(int i=0; i<N; ++i){
        Fij[0] = 0;
        Fij[1] = 0;
        Fij[2] = 0;
        vix = vcurr[i][0]; 
        viy = vcurr[i][1];/
        viz = vcurr[i][2];/
        
        Ui = 0; 
        
        for(int j = i+1; j<N && j !=i; ++j){
            rijx = rcurr[i][0] - rcurr[j][0];
            rijy = rcurr[i][1] - rcurr[j][1];
            rijz = rcurr[i][2] - rcurr[j][2];
            
            vjx = vcurr[j][0];
            vjy = vcurr[j][1];
            vjz = vcurr[j][2];

            rijx = NearestImage(rijx);
            rijy = NearestImage(rijy);
            rijz = NearestImage(rijz);


            rij2 = rijx*rijx + rijy*rijy + rijz*rijz; // rij^2
            rij = pow(rij2,0.5);
            
            if (rij <= Rc){

                r2_1 = 1/(rij2);                      // 1/rij^2
                r6_1 = r2_1 * r2_1 * r2_1;            //1/rij^6
                r12_1 = r6_1 * r6_1;                  //1/rij^12

                Uij = 4*(r12_1 - r6_1) - URc - (rij - Rc)*dURc;
                preFij = (24*r2_1)*(2*r12_1 - r6_1); // Fij

                Fij[0] = (preFij+preFijRc/rij)*rijx; // Fijx
                Fij[1] = (preFij+preFijRc/rij)*rijy; // Fijy
                Fij[2] = (preFij+preFijRc/rij)*rijz; // Fijz

                
                temp_j2 = Fij[0]*(vix+vjx) + Fij[1]*(viy+vjy) +Fij[2]*(viz+vjz); 
                jterm2_temp[0] += rijx*temp_j2;
                jterm2_temp[1] += rijy*temp_j2;
                jterm2_temp[2] += rijz*temp_j2;

                
                Fcurr[i][0] += Fij[0];
                Fcurr[j][0] -= Fij[0];

                Fcurr[i][1] += Fij[1];
                Fcurr[j][1] -= Fij[1];

                Fcurr[i][2] += Fij[2];
                Fcurr[j][2] -= Fij[2];

            }else{
                Uij = 0;

                Fij[0] = 0;
                Fij[1] = 0;
                Fij[2] = 0;
            }
        }
        Ui_arr[i] = Ui; //pot en for each atom stored
        
        jterm2[0] += jterm2_temp[0]; // sum of jterm2 for all atoms
        jterm2[1] += jterm2_temp[1];
        jterm2[2] += jterm2_temp[2];
            
            jterm2_temp[0] = 0;
            jterm2_temp[1] = 0;
            jterm2_temp[2] = 0;
    }// Calculate initial Force for position rcurr

    // Calculate Kinetic energy

    float KE = 0;
    for(int i = 0; i < N; ++i){
        KEi = (vcurr[i][0]*vcurr[i][0] + vcurr[i][1]*vcurr[i][1] + vcurr[i][2]*vcurr[i][2])/2;
        KE += KEi;
        totei = KEi+0.5*Ui_arr[i];                //1/2 mv2 + 1/2 Ui
        jterm1[0] += vcurr[i][0]*totei;
        jterm1[1] += vcurr[i][1]*totei;
        jterm1[2] += vcurr[i][2]*totei;}
    
    
    heat_flux[0] = jterm1[0] + 0.5*jterm2[0];
    heat_flux[1] = jterm1[1] + 0.5*jterm2[1];
    heat_flux[2] = jterm1[2] + 0.5*jterm2[2];
        
    
    // Calculate Kinetic energy
    Temp = (2*KE*121/(3*(N-1)));
    Tempd = (2*KE/(3*(N-1)));
    Pres = Par*(N*Temp/V + Vir/3/V);
    cout<<"Kinetic Energy-----------------------: "<<KE<<"\n";
    cout<<"Potential Energy: "<<Utot/2<<"\n";
    cout<<"Total Energy: "<<(KE+Utot/2)<<"\n";
    cout<<"Temperature: "<<Temp<<"\n";
    cout<<"Pressure: "<<Pres<<"\n";
    cout<<"Heat flux x:"<<heat_flux[0]<<"\n";
    cout<<"Heat flux y:"<<heat_flux[1]<<"\n";
    cout<<"Heat flux z:"<<heat_flux[2]<<"\n";
    cout<<"Eta : "<<eta_curr<<"\n";
    // f_ke << KE<<"\n";
    // f_pe << Utot/2<<"\n";
    // f_tot << (KE+Utot/2)<<"\n";
    // f_temp << Temp<<"\n";
    // f_tempd << Tempd<<"\n";
    // f_pres << Pres<<"\n";
    // f_mom<< 0.0 << " "<< 0.0 << " "<< 0.0 << "\n";
    // f_heatf<< heat_flux[0] << " "<< heat_flux[1] << " "<< heat_flux[2] << "\n";

    // cout<<"\n";

// ---------------------------------------------------------------------------------------------

    
    
    
    
    
    
    
    

// For every time step after the first ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    double tau2 = 0.05*0.05, num, den;
    float dt = 0.002;
    float dt2 = dt * dt;
    float t_last = 3000, tot_step = t_last/dt;
    double Fprev[N][3], rcurr_temp, vhalf[N][3];
    for(int i=0; i < tot_step; ++i){ // Time loop-----------
        eta_prev = eta_curr;
        for(int j = 0; j < N; ++j){

            for(int k =0; k<3; ++k){
                Fprev[j][k] = Fcurr[j][k];
                rprev = rcurr[j][k];
                vprev = vcurr[j][k];
                vhalf[j][k] = vprev + (Fprev[j][k] - eta_prev*vprev)*dt/2;
                //update eta
                rcurr_temp = rprev + vhalf[j][k]*dt;
                rcurr[j][k] = InsideCell(rcurr_temp);
                Fcurr[j][k] *= 0; // Initialize all the forces to zero for the next time step
                //Overwrite velocities for the next step
            }
        }
        if (i >= 300000){ eta_curr = 0;} // For NVE
        if (i < 300000) {eta_curr = eta_prev + dt*(Temp/Tset - 1)/tau2;} // For NVT

    //------------------------------/ Calculate Force: insert loop / ------------------
        Utot = 0;
        Vir = 0;
        for(int j=0; j<N; ++j){
//             Fij[3] = {0};
            Fij[0] = 0;
            Fij[1] = 0;
            Fij[2] = 0;

            
            vix = vcurr[j][0];
            viy = vcurr[j][1];
            viz = vcurr[j][2];
            
            Ui = 0; 
            for(int k = 0; k<N && k !=j; ++k){
                rijx = rcurr[j][0] - rcurr[k][0];
                rijy = rcurr[j][1] - rcurr[k][1];
                rijz = rcurr[j][2] - rcurr[k][2];
                
                vjx = vcurr[k][0];
                vjy = vcurr[k][1];
                vjz = vcurr[k][2];


                rijx = NearestImage(rijx);
                rijy = NearestImage(rijy);
                rijz = NearestImage(rijz);

                rij2 = rijx*rijx + rijy*rijy + rijz*rijz; // rij^2
                rij = sqrt(rij2);
                
                if (rij <= Rc){


                    r2_1 = 1/(rij2);                      // 1/rij^2
                    r6_1 = r2_1 * r2_1 * r2_1;            //1/rij^6
                    r12_1 = r6_1 * r6_1;                  //1/rij^12

                    Uij = 4*(r12_1 - r6_1) - URc - (rij - Rc)*dURc;
                    preFij = (24*r2_1)*(2*r12_1 - r6_1); // Fij

                    Fij[0] = (preFij+preFijRc/rij)*rijx; // Fijx
                    Fij[1] = (preFij+preFijRc/rij)*rijy; // Fijy
                    Fij[2] = (preFij+preFijRc/rij)*rijz; // Fijz

                    Utot += Uij;
                    Ui += Uij;
                    
                    Vir += Fij[0]*rijx + Fij[1]*rijy + Fij[2]*rijz;

                    Fcurr[j][0] += Fij[0];
                    Fcurr[k][0] -= Fij[0];

                    Fcurr[j][1] += Fij[1];
                    Fcurr[k][1] -= Fij[1];

                    Fcurr[j][2] += Fij[2];
                    Fcurr[k][2] -= Fij[2];
             
                    
                    

                }else{
                    Uij = 0;

                    Fij[0] = 0;
                    Fij[1] = 0;
                    Fij[2] = 0;
                    }
            }
            Ui_arr[j] = Ui;
            
                
        }// Calculate Force: insert loop


        double mom[3];
        mom[0] = 0.0;
        mom[1] = 0.0;
        mom[2] = 0.0;
        for (int j = 0; j<N; ++j){
            for(int k=0; k<3; ++k){

                num = vhalf[j][k] + Fcurr[j][k]*dt/2;
                den = 1 + eta_curr*dt/2;
                vcurr[j][k] = num/den;
                mom[k] += vcurr[j][k];
            }
        }
        
        //------------------------------/ Calculate Force: insert loop / ------------------
    
                for(int j=0; j<N; ++j){
       
                    Fij[0] = 0;
                    Fij[1] = 0;
                    Fij[2] = 0;

                    vix = vcurr[j][0];
                    viy = vcurr[j][1];
                    viz = vcurr[j][2];
                    
                    Ui = 0; 
                    for(int k = 0; k<N && k !=j; ++k){
                        rijx = rcurr[j][0] - rcurr[k][0];
                        rijy = rcurr[j][1] - rcurr[k][1];
                        rijz = rcurr[j][2] - rcurr[k][2];
                        
                        vjx = vcurr[k][0];
                        vjy = vcurr[k][1];
                        vjz = vcurr[k][2];

                        rijx = NearestImage(rijx);
                        rijy = NearestImage(rijy);
                        rijz = NearestImage(rijz);

                        rij2 = rijx*rijx + rijy*rijy + rijz*rijz; // rij^2
                        rij = sqrt(rij2);

                        if (rij <= Rc){


                            r2_1 = 1/(rij2);                      // 1/rij^2
                            r6_1 = r2_1 * r2_1 * r2_1;            //1/rij^6
                            r12_1 = r6_1 * r6_1;                  //1/rij^12
                            Uij = 4*(r12_1 - r6_1) - URc - (rij - Rc)*dURc;
                            preFij = (24*r2_1)*(2*r12_1 - r6_1); // Fij

                            Fij[0] = (preFij+preFijRc/rij)*rijx; // Fijx
                            Fij[1] = (preFij+preFijRc/rij)*rijy; // Fijy
                            Fij[2] = (preFij+preFijRc/rij)*rijz; // Fijz
                            
                            temp_j2 = Fij[0]*(vix+vjx) + Fij[1]*(viy+vjy) +Fij[2]*(viz+vjz); 

                            jterm2_temp[0] += rijx*temp_j2;
                            jterm2_temp[1] += rijy*temp_j2;
                            jterm2_temp[2] += rijz*temp_j2;
                            
                        }else{


                            Fij[0] = 0;
                            Fij[1] = 0;
                            Fij[2] = 0;
                            }
                    }
                    
                    jterm2[0] += jterm2_temp[0];
                    jterm2[1] += jterm2_temp[1];
                    jterm2[2] += jterm2_temp[2];
                    
                    jterm2_temp[0] = 0;
                    jterm2_temp[1] = 0;
                    jterm2_temp[2] = 0;
                }// Calculate Force: insert loop

        float KE = 0;
        for(int j = 0; j < N; ++j){
        KEi = (vcurr[j][0]*vcurr[j][0] + vcurr[j][1]*vcurr[j][1] + vcurr[j][2]*vcurr[j][2])/2;
        KE += KEi;
        totei = KEi+0.5*Ui_arr[j];                //1/2 mv2 + 1/2 Ui
        jterm1[0] += vcurr[j][0]*totei;
        jterm1[1] += vcurr[j][1]*totei;
        jterm1[2] += vcurr[j][2]*totei;
        
        }
        
        Temp = (2*KE*121/(3*(N-1)));
        Pres = Par*(N*Temp/V/121 + Vir/3/V);
        Tempd = (2*KE/(3*(N-1)));
        cout<<"\n";
        heat_flux[0] = jterm1[0] + 0.5*jterm2[0];
        heat_flux[1] = jterm1[1] + 0.5*jterm2[1];
        heat_flux[2] = jterm1[2] + 0.5*jterm2[2];
        
        avgTemp += Temp;
        avgPres += Pres;
        
        cout<<i<<" "<<tot_step<<"\n";
        cout<<"Momentum x: "<<mom[0]<<"\n";
        cout<<"Momentum y: "<<mom[1]<<"\n";
        cout<<"Momentum z: "<<mom[2]<<"\n";
        cout<<"Kinetic Energy: "<<KE<<"\n";
        cout<<"Potential Energy: "<<Utot<<"\n";
        cout<<"Total Energy: "<<(KE+Utot)<<"\n";
        cout<<"Temperature: "<<Temp<<"\n";
        cout<<"Pressure: "<<Pres<<"\n";
        cout<<"Heat flux x:"<<heat_flux[0]<<"\n";
        cout<<"Heat flux y:"<<heat_flux[1]<<"\n";
        cout<<"Heat flux z:"<<heat_flux[2]<<"\n";
        cout<<"Eta : "<<eta_curr<<"\n";
        if (i%5==0 && i>500000){ // Collect data only for the last 1mil steps - 0.5 mil to 1.5mil
	// f_ke << KE<<"\n";
        // f_pe << Utot/2<<"\n";
        // f_tot << (KE+Utot/2)<<"\n";
        // f_temp << Temp<<"\n";
        // f_tempd << Tempd<<"\n";
        // f_pres << Pres<<"\n";
        // f_mom<< mom[0]<< " "<< mom[1] << " "<< mom[2] << "\n";
        //if(i>500000){
        f_heatf<< heat_flux[0] << " "<< heat_flux[1] << " "<< heat_flux[2] << "\n";
        }

        cout<<"\n";
//         jterm1[3] = {0};
//         jterm2[3] = {0};
        
        jterm1[0] = 0.0;
        jterm1[1] = 0.0;
        jterm1[2] = 0.0;
        
        jterm2[0] = 0.0;
        jterm2[1] = 0.0;
        jterm2[2] = 0.0;
        
        cout<<"zeroing:"<<jterm2[0]<<"\n";
       // print KE, PE, TE and Temp
        
    }
    cout<<"Average Temperature: "<<avgTemp/tot_step;
    cout<<"Average Pressure: "<<avgPres/tot_step;
        // for (int ii = 0;ii < N; ii++){
        //     f_pos<< rcurr[ii][0]<<" "<< rcurr[ii][1] << " "<< rcurr[ii][2] << "\n";}
}

