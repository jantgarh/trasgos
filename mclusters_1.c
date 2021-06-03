// mclusters_1.C  Nuevo codigo basado en dat6_0 para construir clusters con todo tipo de particulas y con informacion mas precisa: tiempos y energias separados para electrones y muones.
//   Version operativa 23.abr.21
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
#define rclust 2.0 //m. Distancia de separación entre partículas
#define limite 1.0 // m. limite de separacion entre interacciones

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

    Int_t   ncont=0, n=0, nsecp, icont=0, ntsep=0;
    Int_t   label=-1, icshow, nclust=0, cmult, cmultp1;
    Int_t   pid=0, pid1=0, pid2=0, pidmin=1000000000, pidmax=0;
    Int_t   idp;  // I introduce my own particle identificator
    Int_t   ngam=0, nele=0, nmu=0, nn=0, np=0, noth=0;
    Int_t   ntgam=0, ntele=0, ntmu=0, ntn=0, ntp=0, ntoth=0;
    Long64_t id1=.0, id2=.0, idclus=0, tidclus=0;
    Float_t id, x, y, dx, dy, x1, x2, y1, y2, e, r, r1, r2, t;
    Float_t theta=0., theta1=0., theta2=0., phi=0, phi1=0, phi2=0,
            theta_tan=0, theta_tan1=0,  thetamean=0., tttheta=0., sigmatheta=0.;
    Float_t time =0., t1=0., t2=0., tmin = 1000000000, tmax=0.,
            rmin = 1000000000, rmax=0,
            thetamin = 1000000000,thetamax = 0, phimin = 1000000000,  phimax = 0;
    Float_t shwh = 0, shwth = 0, shwph = 0, pt = 0;
    Float_t px, py, pz, px1, py1, pz1, px2, py2, pz2, pmod,
            psq, psq1, psq2; // momenta particles
    Float_t dxsq, dysq, drsq, ds, sigr=0;
    Float_t ene, ene1, ene2,
            kene, kene1, kene2, kenesq, kene2sq,
            kmax, kmin=1000000000, kkmax, kkmin=1000000000;
    Float_t xmax=0, xmin=0, ymax=0, ymin=0;
    Float_t xcentral=0, sigmaX=0, ycentral=0, sigmaY=0;
    Float_t xmean = 0, xmeansq = 0, xmold=0, xmoldsq=0, sigxsq=0;
    Float_t ymean = 0, ymeansq = 0, ymold=0, ymoldsq=0, sigysq=0;
    Float_t tmean = 0, tmeansq = 0, tmold=0, tmoldsq=0, sigtsq=0;
    Float_t kmean = 0, kmeansq = 0, kmold=0, kmoldsq=0, sigksq=0;
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

nshows = 2;

   for (Long64_t ishow=0; ishow<nshows; ishow++) {     // Loop in showers
       
       cout << " *** new shower *** " << endl;
       
       icont = -1;   // index for total nb. of secondaries
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
       cmult = 1;
       icshow = 0;
       particle__ = 100;
       nsecp = particle__;  // nb. secondary part. in shower
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
               idp = 1000;
               ene  = sqrt(psq + mele*mele);
               kene = ene - mele;
               // if(kene1<0.1) continue;
               nele++;
               ntele++;
           }
           else if( (id==5) || (id==6) ){ // MUONS
               idp = 1000000;
               ene  = sqrt(psq + mmu*mmu);
               kene = ene - mmu;
               // if(kene1<0.1) continue;
               nmu++;
               ntmu++;
           }
           else if( id==13 ){ // NEUTRONS
               idp = 1000000000;   // nucleon
               ene  = sqrt(psq + mn*mn);
               kene = ene - mn;
               // if(kene1<0.1) continue;
               nn++;
               ntn++;
           }
           else if( id==14){ // PROTONS
               idp = 1000000000;   // nucleon
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
       
       // Looking for clusters. First loop in the shower
       for(Int_t iss=0; iss<nsecp; iss++){

           // kene1 = 0;
           // ene1  = 0;

           if(ilab[iss]!=-1) continue;

           pid     = arr0[iss];
           x1      = arr1[iss];
           y1      = arr2[iss];
           t1      = arr3[iss];
           px1     = arr4[iss];
           py1     = arr5[iss];
           pz1     = arr6[iss];
           theta1  = arr7[iss];
           phi1    = arr8[iss];
           ene1    = arr9[iss];
           kene1   = arr10[iss];
           r1      = sqrt(x1*x1 + y1*y1);
           psq1    = px1*px1 + py1*py1 + pz1*pz1;
           
           // cout << "pid, cmult: " << pid << " " << cmult << endl;
           
           idclus   = pid;
           xmean    = x1;
           ymean    = y1;
           tmean    = t1;
           kmean    = kene1;
           
           ilab[iss]++;
           
           // cout << "* "<< ishow << " " << iss << " "<< ilab[iss] << endl;
           
           // = 0;
           // Looking for clusters. Look for more particles
           for(Int_t jss=iss+1; jss < nsecp; jss++){

               if(ilab[jss]!=-1) continue;

               dx   = arr1[jss] - xmean;
               dxsq = dx * dx ;
               dy   = arr2[jss] - ymean;
               dysq = dy * dy ;
               drsq = dxsq + dysq;

               //sigr = (sqrt(dxsq)*sigxsq + sqrt(dysq)*sigysq)/sqrt(drsq);
               // no entiendo muy bien este corte.
               //sigr = sqrt(sigxsq * dxsq/drsq + sigysq * dysq/drsq) ;  // JAG
           
               // new particle in the cluster
               if(sqrt(drsq)<=rclust && sigr<=limite){
                   
                   //cout << " dx, dy, sr: " << dx << " " << dy << " " << sigr << endl;
 
                   label++;   // ?
                   cmult++;
                   
                   pid2   = arr0[jss];
                   t2     = arr3[jss];
                   x2     = arr1[jss];
                   y2     = arr2[jss];
                   r2     = sqrt(x2*x2 + y2*y2);
                   px2    = arr4[jss];
                   py2    = arr5[jss];
                   pz2    = arr6[jss];
                   psq2   = px2*px2 + py2*py2 + pz2*pz2;
                   theta2 = arr7[jss];
                   phi2   = arr8[jss];
                   ene2   = arr9[jss];
                   kene2  = arr10[jss];
                   
                   // cout << "pid2, cmult " << pid2 << " " << cmult << endl;
                   
                   idclus = idclus + pid2;

                   // Calculo del nuevo centro del cluster y su anchura
                   // Cuidado con la forma recursiva Hay que sumarle 1 al conteo de multiplicidad ya que falta contar la 1º interaccion

                   xmold   = xmean;
                   ymold   = ymean;
                   tmold   = tmean;
                   xmoldsq = xmold * xmold;
                   ymoldsq = ymold * ymold;
                   tmoldsq = tmold * tmold;
                   // Formulas recursivas del valor medio y la dispersion
                   // sustituyo las formulas de Yanis por otras equivalentes
                   xmean   = xmean + (x2 - xmean)/cmult;
                   ymean   = ymean + (y2 - ymean)/cmult;
                   tmean   = tmean + (t2 - tmean)/cmult;
                   xmeansq = xmean * xmean;
                   ymeansq = ymean * ymean;
                   tmeansq = tmean * tmean;
                
                   sigxsq = sigxsq + xmoldsq - xmeansq +
                           ((x2*x2 - sigxsq - xmoldsq)/cmult);
                   sigysq = sigysq + ymoldsq - ymeansq +
                           ((y2*y2 - sigysq - ymoldsq)/cmult);
                   sigr   = sqrt(sigxsq + sigysq) ;
               
                   // Reescala de los valores temporales. Si no, aparecen incoherencias
                   // el resultado debido a la precision        <--- ? JAG
               
                   tmold  = tmold - tmean;
                   t2     = t2 - tmean;
                   sigtsq = sigtsq + tmoldsq - tmeansq +
                            ((t2*t2 - sigtsq - tmoldsq)/cmult);
               
                   kene2sq = kene2 * kene2;
                   kmold = kmean;
                   kmoldsq = kmold * kmold;
                   kmean     = kmean + (kene2 - kmean)/cmult;
                   kmeansq  = kmean * kmean;
                   sigksq = sigksq + kmoldsq - kmeansq +
                           ((kene2sq-sigksq - kmoldsq)/cmult);

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
               
               // cout << " * lab cmult " << label << " " << cmult << endl;
               
               ilab[jss] = label;
               label     = -1;
           }  //  end for jss
           
           if (cmult > 1){
               nclust++;
               icshow ++;
               tidclus = tidclus + idclus;  // suma logica de los clusters en el shower
               
               //*
               cout << "*** icshow " << icshow <<  endl;
               cout << "**** nclust " << nclust << endl;
               cout << "* idclus " << idclus << endl;
               cout << "* cmult "  << cmult <<  endl;
               cout << "(xmean, ymean), sigr " << xmean << " " << ymean << " " << sigr << endl;
               cout << endl;
               
               archivo << ishow+1 << "\t" << nclust << "\t" << cmult << "\t" <<
                 idclus  << "\t" << r1  << "\t" << pidmin << "\t" <<
                 xmin << "\t" << ymin << "\t" << tmin  << "\t" <<
                 thetamin << "\t" << phimin << "\t" << kmin << "\t" << pidmax << "\t" <<
                 xmax << "\t" << ymax << "\t" << tmax  << "\t" <<
                 thetamax << "\t" << phimax << "\t" << kmax << "\t" <<
                 xmean  << "\t" << sqrt(sigxsq)  << "\t" <<
                 ymean  << "\t" << sqrt(sigysq)  << "\t" <<
                 tmean  << "\t" << sqrt(sigtsq)  << "\t" <<
                 kkmax << "\t"  << kkmin << "\t" << "\t" <<
                 kmean << "\t"  << sqrt(sigksq) << endl;
               
               //*/
           }  // endif cmult > 1:  new cluster

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
       }  // endif

       if(t1 < tmin){
           tmin     = t1;
           xmin     = x1;
           ymin     = y1;
           pidmin   = id1;
           kmin     = kene1;
           thetamin = theta1;
           phimin   = phi1;
       }  // endif

       /*
        sigxsq = sqrt(sigxsq);
        sigysq = sqrt(sigysq);
        sigtsq = sqrt(sigtsq);
        sigksq = sqrt(sigksq);
        */

           
       //nclust    = 0.;
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
       cmult     = 1.;
       tmax      = 0.;
       tmin      = 1000000000.;
       t2        = 0.;
       t1        = 0.;
       xmax      = 0.;
       xmin      = 1000000000.;
       xmean     = 0.;
       xcentral  = 0.;
       sigmaX    = 0.;
       sigxsq = 0.;
       xmold = 0.;
       ymax      = 0.;
       ymin      = 1000000000.;
       ymean     = 0.;
       ycentral  = 0.;
       sigmaY    = 0.;
       sigysq = 0.;
       ymold = 0.;
       pidmax    = 0.;
       pidmin    = 1000000000.;
       kmin      = 1000000000.;
       phimin    = 1000000000.;
       thetamin  = 1000000000.;
       thetamax  = 0.;
       phimax    = 0.;
       kmax      = 0.;
       tmean     = 0.;
       sigtsq    = 0.;
       tmold = 0.;
       kmean     = 0.;
       sigksq    = 0.;
       kmold = 0.;
       kkmax     = 0.;
       kkmin     = 1000000000.;
       dxsq      = 0.;
       dysq      = 0.;
       drsq      = 0.;
       sigr      = 0.;
       psq2      = 0.;
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
       
       cout << endl;
       cout << "* icshow " << icshow <<  endl;
       cout << "* Suma de clusters ids: " << tidclus << endl;
       cout << endl;
       tidclus   = 0;
       
   }  // end for in showers
    
printf("*** Proceso completado");
archivo.close();
cout << endl;
cout << "* Numero de cascadas : " << nshows << endl;
//cout << endl;
cout << "* Numero total de secundarios: " << ntsep << endl;
cout << "* Numero total de clusters: " << nclust << endl;
cout << "* Distribucion (gemnpo): " << ntgam << " " << ntele << " " << ntmu <<  " " << ntn <<  " " << ntp <<  " " << ntoth << endl;
    cout << "* Distribucion relativa*1000: " << 1000*ntgam/ntsep << " " << 1000*ntele/ntsep << " " << 1000*ntmu/ntsep <<  " " << 1000*ntn/ntsep <<  " " << 1000*ntp/ntsep <<  " " << 1000*ntoth/ntsep << endl;
cout << endl;

// Estimamos densidades de particulas y sigmas correspondientes

}
