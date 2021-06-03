// mclusters.C  Nuevo codigo basado en dat6_0 para construir clusters con todo tipo de particulas y con informacion mas precisa: tiempos y energias separados para electrones y muones.
// !!!!!
// Cuando puedas cambiar px, py y pz por pmod y pt!!!!
// !!!!!
// JAG.Ene.21. Intento de redefinir los clusters con todas las particulas y separando componente electronica y muonica.
//Programa muy robusto. ANALISIS DE MULTIPLICIDAD. Identificacion de particulas a menos de 1 metro de radio entre ellas. Con este programa queremos conocer los efectos sitematicos existente en la cascada.
//Metodo recursivo para valores n muy grandes, para un unico valor no funciona..
// Posible metodo utilizar los array y si son excesivamente grandes pasarse al metodo recursivo

#define mclusters_cxx
#include "mclusters.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#define mele 0.000511 // GeV
#define mmu 0.10566 // GeV
#define mgam 0.0 // GeV
#define mp 0.9383  // GeV
#define mn 0.9396  // GeV
#define mother 0.14 // GeV. mainly pions
#define rclust 1.0 //m. Distancia de separación entre partículas
#define limite 5.0 // m. limite de separacion entre interacciones

void mclusters::Loop()
{
    
// Analisis de cascadas atmosfericas de rayos cosmicos de distinta energia con incidencia vertical.
/*
 *  Intervalos de Energia.c
 *
 *  Created by Juan A. Garzon on 23/01/13
 *  Copyright LabCAF. IGFAE/USC. All rights reserved.
 *  Modified by G.Kornakov on  11/03/2013
 *  Improved by JAGarzon several times
 *  Modified by Yanis Fontenla Barba on 23/05/2017 for Corsika
*/
    
//   In a ROOT session, you can do:
//      root> .L dat6.C
//      root> dat6 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    ishow is the global entry number in the chain
//    itree is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    ishow for TChain::GetEntry
//    itree for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(ishow);       //read all branches
//by  b_branchname->GetEntry(itree); //read only this branch

    if (fChain == 0) return;

    Long64_t nshows = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    Long64_t fbytes = 0, fb = 0;
    Long64_t pbytes = 0, pb = 0;
    Long64_t hbytes = 0, hb = 0;
    Long64_t tbytes = 0, tb = 0;
    // Identificadores de particulas en Corsika:
    //  1        gammas
    //  2 3      e- e+
    //  5 6      mu- mu+
    //  7 8 9    pi0 pi+ pi-
    //  11 12    K+ K-
    //  13 14    n p
    // =======================================

    Int_t   ncont=0, n=0, nsecp, icont=0, ntsep=0, iclust=0;
    Int_t   label=-1, cmult=0, cuentmas=0;
    Int_t   pid=0, pid1=0, pid2=0, pidmin=1000000000, pidmax=0;
    Int_t   idp;  // I introduce my own particle identificator
    Int_t   ngam=0, nele=0, nmu=0, nn=0, np=0, noth=0;
    Int_t   ntgam=0, ntele=0, ntmu=0, ntn=0, ntp=0, ntoth=0;
    Long64_t id1=.0, id2=.0, idclus=.0;
    Float_t id, x, y, x1, x2, y1, y2, e, r, r1, r2, t ; //
    Float_t theta=0., theta1=0., theta2=0., phi=0, phi1=0, phi2=0,
            theta_tan=0, theta_tan1=0,  thetamean=0., tttheta=0., sigmatheta=0.;
    Float_t time =0., t1=0., t2=0., tmin = 1000000000, tmax=0.,
            rmin = 1000000000, rmax=0,
            thetamin = 1000000000,thetamax = 0, phimin = 1000000000,  phimax = 0;
    Float_t shwh = 0, shwth = 0, shwph = 0, pt = 0;
    Float_t px, py, pz, px1, py1, pz1, px2, py2, pz2, pmod,
            psq, psq1, psq2; // momenta particles
    Float_t dsq, rsq2, dxsq, dysq, sigmar;
    Float_t ene, ene1, ene2,
            kene, kene1, kene2, kenesq, kene2sq,
            kmax, kmin=1000000000, kkmax, kkmin=1000000000;
    Float_t xmax=0, xmin=0, ymax=0, ymin=0;
    Float_t xcentral=0, sigmaX=0, ycentral=0, sigmaY=0;
    Float_t xmean = 0, xmean_old=0, xmeansq = 0, xmean_oldsq=0, xsigmasq=0;
    Float_t ymean = 0, ymean_old=0, ymeansq = 0, ymean_oldsq=0, ysigmasq=0;
    Float_t tmean = 0, tmean_old=0, tmeansq = 0, tmean_oldsq=0, tsigmasq=0;
    Float_t kmean = 0, kmean_old=0, kmeansq = 0, kmean_oldsq=0, ksigmasq=0;
    Float_t rad2g;
    rad2g = 180/TMath::Pi();

    // Definimos las energias en intervalos logaritmicos de 0.25

    fstream archivo;
    //        archivo.open("C_P_Electrons_1E15_10.txt", fstream::out);
    archivo.open("dat6_out.txt", fstream::out);
        archivo <<
        "#NShow"      << "\t" <<
        "NClust"      << "\t" <<
        "SzClust"     << "\t" <<
        "IdClust"     << "\t" <<
        "DistOr/m"    << "\t" <<
        "Id_{Ti}"     << "\t" <<
        "x_{Ti}/m"    << "\t" <<
        "y_{Ti}/m"    << "\t" <<
        "time_{i}"    << "\t" <<
        "Zen_{Ti}/º"  << "\t" <<
        "Azim_{Ti}/º" << "\t" <<
        "KEne_{Ti}/GeV" << "\t" <<
        "Id_{Tf}"     << "\t" <<
        "x_{Tf}/m"    << "\t" <<
        "y_{Tf}/m"    << "\t" <<
        "time_{f}"    << "\t" <<
        "Zen_{Tf}/º"  << "\t" <<
        "Azim_{Tf}/º" << "\t" <<
        "KEne_{Tf}/GeV" << "\t" <<
        "XClust/m"   << "\t" <<
        "sXClust/m"  << "\t" <<
        "YClust/m"   << "\t" <<
        "sYClust/m"  << "\t" <<
        "Tmean/ns"   << "\t" <<
        "sTmean/ns"  << "\t" <<
        "KMax/GeV"   << "\t" <<
        "KMin/GeV"   << "\t" <<
        "Kmean/GeV"  << "\t" <<
        "sKmean/GeV" << endl;

nshows = 1;

   for (Long64_t ishow=0; ishow<nshows; ishow++) {     // Loop in showers
       
       icont = -1;
       Long64_t itree = LoadTree(ishow);

       if (itree < 0) break;

       nb = fChain->GetEntry(itree);   nbytes += nb;
       // if (Cut(itree) < 0) continue;

       fb  = b_shower_FirstHeight->GetEntry(itree);   fbytes += fb;
       hb  = b_shower_Theta->GetEntry(itree);         hbytes += hb;
       tb  = b_shower_Phi->GetEntry(itree);           tbytes += tb;

       shwh  = shower_FirstHeight/100; // h en metros
       shwth = shower_Theta * rad2g;;
       shwph = shower_Phi * rad2g;;

       Float_t arr0[100000] = { [0 ... 99999] = -10. };
       Float_t arr1[100000] = { [0 ... 99999] = -10. };
       Float_t arr2[100000] = { [0 ... 99999] = -10. };
       Float_t arr3[100000] = { [0 ... 99999] = -10. };
       Float_t arr4[100000] = { [0 ... 99999] = -10. };
       Float_t arr5[100000] = { [0 ... 99999] = -10. };
       Float_t arr6[100000] = { [0 ... 99999] = -10. };
       Float_t arr7[100000] = { [0 ... 99999] = -10. };
       Float_t arr8[100000] = { [0 ... 99999] = -10. };
       Float_t arr9[100000] = { [0 ... 99999] = -10. };
       Float_t arr10[100000] = { [0 ... 99999] = -10. };
       Float_t ilab[100000] = { [0 ... 99999] =  -1. };
       
       // First loop in particles
       nsecp = particle__;  // secondary part. in the shower
       for(Int_t ip=0; ip < particle__; ip++){
           
           x         = particle__x[ip];
           y         = particle__y[ip];
           r         = sqrt(x*x + y*y);
           time      = particle__Time[ip];
           px        = particle__Px[ip];
           py        = particle__Py[ip];
           pz        = particle__Pz[ip];
           psq       = px*px + py*py + pz*pz;
           pmod      = sqrt(psq);
           //pt        = ((px*y)-(py*x))/(r*psq); ???
           pt        = sqrt(px*px + py*py);
           theta     = acos(pz/pmod) * rad2g;
           phi       = atan2(py,px) * rad2g;
           //          theta_tan = atan2(y/x)*(180/TMath::Pi());
           
           id = particle__ParticleID[ip];
   //      Select type of particle
           if(id==1){ // GAMMAS
               idp = 1;
               kene = pmod;
               ene  = pmod;
               ngam++;
               ntgam++;
           }
           else if( (id==2) || (id==3) ){ // ELECTRONS
               idp = 100;
               ene  = sqrt(psq + mele*mele);
               kene = ene - mele;
               // if(kene1<0.1) continue;
               nele++;
               ntele++;
           }
           else if( (id==5) || (id==6) ){ // MUONS
               idp = 10000;
               ene  = sqrt(psq + mmu*mmu);
               kene = ene - mmu;
               // if(kene1<0.1) continue;
               nmu++;
               ntmu++;
           }
           else if( id==13 ){ // NEUTRONS
               idp = 100000;
               ene  = sqrt(psq + mn*mn);
               kene = ene - mn;
               // if(kene1<0.1) continue;
               nn++;
               ntn++;
           }
           else if( id==14){ // PROTONS
               idp = 1000000;
               ene  = sqrt(psq + mp*mp);
               kene = ene - mp;
               // if(kene1<0.1) continue;
               np++;
               ntp++;
           }
           else{
               idp = 10000000;
               ene  = sqrt(psq + mother*mother);
               kene = ene - mother;
               // if(kene1<0.1) continue;
               noth++;
               ntoth++;
           }
           icont++;     // particle

           arr0[icont] = idp;
           arr1[icont] = x/100; //m
           arr2[icont] = y/100; //m
           arr3[icont] = time;
           arr4[icont] = px;
           arr5[icont] = py;
           arr6[icont] = pz;
           arr7[icont] = theta;
           arr8[icont] = phi;
           arr9[icont] = ene;
           arr10[icont] = kene;
       //}
       }    // end loop in particles
       
       // cout << nsecp << ": " << ngam << " " << nele << " " << nmu <<  " " << nn <<  " " << np <<  " " << noth << endl;
       
       ngam = 0 ; nele=0; nmu=0; nn=0; np=0; noth=0;
       ntsep = ntsep + nsecp;   // total nb. secondary part.
       
       // Multiplicidad
       for(Int_t i=0; i < icont; i++){

           kene1 = 0;
           ene1  = 0;

           if(ilab[i]!=-1) continue;

           pid     = arr0[i];
           x1      = arr1[i];
           y1      = arr2[i];
           t1      = arr3[i];
           px1     = arr4[i];
           py1     = arr5[i];
           pz1     = arr6[i];
           theta1  = arr7[i];
           phi1    = arr8[i];
           ene1    = arr9[i];
           kene1   = arr10[i];
           r1      = sqrt(x1*x1 + y1*y1);
           psq1    = px1*px1 + py1*py1 + pz1*pz1;
           
           //cout << " pid " <<  pid  << endl;

           // ELECTRONS
           /*
            id1      = 1000; // Codigo de contaje para electrones. Hasta 999 electrones
            ene1  = sqrt(psq1 + mele*mele);
            kene1 = ene1 - mele;
            if(kene1<0.1) continue;
            */

           // MUONS
           id1      = 1000000; // Codigo de contaje para muones.
           ene1  = sqrt(psq1 + mmu*mmu);
           kene1 = ene1 - mmu;
           idclus   = id1;
           xmean    = x1;
           ymean    = y1;
           tmean    = t1;
           kmean    = kene1;
           ilab[i]++;

           for(Int_t j=i+1; j < icont; j++){      // busca en el resto de la lista

               kene2 = 0;
               ene2  = 0;

               if(ilab[j]!=-1) continue;

               iclust++;

               dxsq   = (arr1[j]-xmean) * (arr1[j]-xmean);
               dysq   = (arr2[j]-ymean) * (arr2[j]-ymean);
               rsq2   = dxsq + dysq;
               //sigmar = (sqrt(dxsq)*xsigmasq + sqrt(dysq)*ysigmasq)/sqrt(rsq2);
               sigmar  = sqrt(dxsq*xsigmasq + dysq*ysigmasq) / sqrt(rsq2);  // JAG
           
               if((sqrt(rsq2)<=rclust) && (sigmar<=limite)){   // es necesario sigmar?

                   label++;
                   cmult++;

                   pid2    = arr0[j];
                   t2      = arr3[j];
                   x2      = arr1[j];
                   y2      = arr2[j];
                   r2      = sqrt(x2*x2 + y2*y2);
                   px2     = arr4[j];
                   py2     = arr5[j];
                   pz2     = arr6[j];
                   psq2   = px2*px2 + py2*py2 + pz2*pz2;
                   theta2  = arr7[j];
                   phi2    = arr8[j];

                   // Calculo del nuevo centro del cluster y su anchura
                   // Cuidado con la forma recursiva Hay que sumarle 1 al conteo de
                   // multiplicidad ya que falta contar la 1º interaccion

                   xmean_old = xmean;
                   ymean_old = ymean;
                   tmean_old = tmean;
                   xmean_oldsq = xmean_old * xmean_old;
                   ymean_oldsq = ymean_old * ymean_old;
                   tmean_oldsq = tmean_old * tmean_old;

                   // sustituyo las formulas de Yanis por otras equivalentes
                   xmean  = xmean + (x2 - xmean)/(cmult+1);
                   ymean  = ymean + (y2 - ymean)/(cmult+1);
                   tmean  = tmean + (t2 - tmean)/(cmult+1);
                   xmeansq = xmean * xmean;
                   ymeansq = ymean * ymean;
                   tmeansq = tmean * tmean;
                   

                   xsigmasq = xsigmasq + xmean_oldsq - xmeansq +
                           ((x2*x2 - xsigmasq - xmean_oldsq)/(cmult+1));
                   ysigmasq = ysigmasq + ymean_oldsq - ymeansq +
                           ((y2*y2 - ysigmasq - ymean_oldsq)/(cmult+1));
               
                   // Reescala de los valores temporales. Si no, aparecen incoherencias
                   // el resultado debido a la precision        <--- ? JAG
               
                   tmean_old = tmean_old - tmean;
                   t2        = t2 - tmean;
                   tsigmasq = tsigmasq + tmean_oldsq - tmeansq +
                            ((t2*t2 - tsigmasq - tmean_oldsq)/(cmult+1));
   /*
                   // ELECTRONS
                   id2     = 1000; // Codigo de contaje para electrones
                   ene2  = sqrt(psq2 + mele*mele);
                   kene2 = ene2 - mele;
   */
                   // MUONS
               
                   id2       = 1000000; // Codigo de contaje para muones
                   ene2   = sqrt(psq2 + mmu*mmu);
                   kene2  = ene2 - mmu;
                   kene2sq = kene2 * kene2;
                   kmean_old = kmean;
                   kmean_oldsq = kmean_old * kmean_old;
                   kmean     = kmean + (kene2 - kmean)/(cmult+1);
                   kmeansq  = kmean * kmean;
                   ksigmasq = ksigmasq + kmean_oldsq - kmeansq +
                           ((kene2sq-ksigmasq - kmean_oldsq)/(cmult+1));

                   // search the slowest particle
                   if(t2 > tmax){
                       tmax     = t2;
                       xmax     = x2;
                       ymax     = y2;
                       pidmax   = id2;
                       kmax     = kene2;
                       thetamax = theta2;
                       phimax   = phi2;
                   }
                   // search the fastest particle
                   if(t2 < tmin){
                       tmin     = t2;
                       xmin     = x2;
                       ymin     = y2;
                       pidmin   = id2;
                       kmin     = kene2;
                       thetamin = theta2;
                       phimin   = phi2;
                   }
                   /*
                    if(kene2 > kkmax){ kkmax = kene2;}
                    if(kene2 < kkmin){ kkmin = kene2;}
                    */
               
                   idclus = idclus + id2;
               }
           
               ilab[j] = label;
               label     = -1;
           }  //  end for j

       /*
        if(kene1 > kkmax){kkmax  = kene1; }
        if(kene1 < kkmin){kkmin  = kene1; }
        */
           
       if(t1 > tmax){
           tmax     = t1;
           xmax     = x1;
           ymax     = y1;
           pidmax   = id1;
           kmax     = kene1;
           thetamax = theta1;
           phimax   = phi1;
       }

       if(t1 < tmin){
           tmin     = t1;
           xmin     = x1;
           ymin     = y1;
           pidmin   = id1;
           kmin     = kene1;
           thetamin = theta1;
           phimin   = phi1;
       }

       /*
        xsigmasq = sqrt(xsigmasq);
        ysigmasq = sqrt(ysigmasq);
        tsigmasq = sqrt(tsigmasq);
        ksigmasq = sqrt(ksigmasq);
        */
       // << " hasta aqui ha llegado " << endl;
       archivo << ishow+1 << "\t" << iclust+1 << "\t" << cmult+1 << "\t" <<
         idclus  << "\t" << r1  << "\t" << pidmin << "\t" <<
         xmin << "\t" << ymin << "\t" << tmin  << "\t" <<
         thetamin << "\t" << phimin << "\t" << kmin << "\t" << pidmax << "\t" <<
         xmax << "\t" << ymax << "\t" << tmax  << "\t" <<
         thetamax << "\t" << phimax << "\t" << kmax << "\t" <<
         xmean  << "\t" << sqrt(xsigmasq)  << "\t" <<
         ymean  << "\t" << sqrt(ysigmasq)  << "\t" <<
         tmean  << "\t" << sqrt(tsigmasq)  << "\t" <<
         kkmax << "\t"  << kkmin << "\t" << "\t" <<
         kmean << "\t"  << sqrt(ksigmasq) << endl;

       iclust    = 0.;
       idclus    = 0.;
       x1        = 0.;
       y1        = 0.;
       x2        = 0.;
       y2        = 0.;
       id1       = 0.;
       id2       = 0.;
       pid1      = 0.;
       pid2      = 0.;
       label     = -1;
       cmult     = 0.;
       tmax      = 0.;
       tmin      = 1000000000.;
       t2        = 0.;
       t1        = 0.;
       xmax      = 0.;
       xmin      = 1000000000.;
       xmean     = 0.;
       xcentral  = 0.;
       sigmaX    = 0.;
       xsigmasq = 0.;
       xmean_old = 0.;
       ymax      = 0.;
       ymin      = 1000000000.;
       ymean     = 0.;
       ycentral  = 0.;
       sigmaY    = 0.;
       ysigmasq = 0.;
       ymean_old = 0.;
       pidmax    = 0.;
       pidmin    = 1000000000.;
       kmin      = 1000000000.;
       phimin    = 1000000000.;
       thetamin  = 1000000000.;
       thetamax  = 0.;
       phimax    = 0.;
       kmax      = 0.;
       tmean     = 0.;
       tsigmasq = 0.;
       tmean_old = 0.;
       kmean     = 0.;
       ksigmasq = 0.;
       kmean_old = 0.;
       kkmax     = 0.;
       kkmin     = 1000000000.;
       dxsq     = 0.;
       dysq     = 0.;
       rsq2     = 0.;
       sigmar    = 0.;
       psq2     = 0.;
       r1        = 0.;
       }

       for(Int_t i=0; i < 100000; i++){
           ilab[i]=-1;
           arr0[i]=-10.;
           arr1[i]=-10.;
           arr2[i]=-10.;
           arr3[i]=-10.;
           arr4[i]=-10.;
           arr5[i]=-10.;
           arr6[i]=-10.;
           arr7[i]=-10.;
           arr8[i]=-10.;
       }
   }  // end for in showers
    
printf("*** Proceso completado");
archivo.close();
cout << endl;
cout << "* Numero de cascadas : " << nshows << endl;
//cout << endl;
cout << "* Numero total de secundarios :" << ntsep << endl;
cout << "* Distribucion (gemnpo): " << ntgam << " " << ntele << " " << ntmu <<  " " << ntn <<  " " << ntp <<  " " << ntoth << endl;
    cout << "* Distribucion relativa * 1000: " << 1000*ntgam/ntsep << " " << 1000*ntele/ntsep << " " << 1000*ntmu/ntsep <<  " " << 1000*ntn/ntsep <<  " " << 1000*ntp/ntsep <<  " " << 1000*ntoth/ntsep << endl;
cout << endl;

// Estimamos densidades de particulas y sigmas correspondientes

}
