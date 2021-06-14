// mclusters.C  Nuevo codigo basado en dat6_0 para construir clusters con todo tipo de particulas y con informacion mas precisa: tiempos y energias separados para electrones y muones.
// !!!!! Valido para simulaciones sin fijar la altura
// JAG.Ene.21. Intento de redefinir los clusters con todas las particulas y separando componente electronica y muonica.
//Programa muy robusto. ANALISIS DE MULTIPLICIDAD. Identificacion de particulas a menos de una cierta distancia entre ellas. Con este programa queremos conocer los efectos sitematicos existente en la cascada.

#define mclusters_cxx
#include "mclusters.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#define mele 0.000511 // GeV
#define mmu 0.10566 // GeV
#define mgam 0.0  // GeV
#define mp 0.9383 // GeV
#define mn 0.9396 // GeV
#define mother 0.14 // GeV. mainly pions
#define distmx 2.0 //m. Distancia maxima de separación entre partículas
#define sigrmx 1.0 //m. Dispersion radial maxima del cluster
#define mpcr 1   // A mass of primary cosmic ray
#define epcr 5  // Log(Energy/GeV) primary cosmic ray
#define zhpcr 0 // Zenith angle of primary cosmic ray
//#define hghtfi 0 // Height/km of first interaction
#define fout1 "xxx.txt" //
#define fout2 "xxy.txt" //

// ----------------------------------------------------------------
void mclusters::Loop()
{
    
// Analisis de cascadas atmosfericas de rayos cosmicos
/*
 *
 *  Copyright LabCAF. IGFAE/USC. All rights reserved.
 *  Created by Juan A. Garzon on 23/01/13
 *  Modified by G.Kornakov on  11/03/2013
 *  Modified by Yanis Fontenla Barba on 23/05/2017 for Corsika
 *  Modified by JAGarzon on Apr.2021 for cluster analysis
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

    Long64_t nshow = fChain->GetEntriesFast();
    cout << " - nShowers in file: " << nshow << endl;
    nshow = 1;   //-
    cout << " - nShowers analized: " << nshow << endl;
    cout << endl;
 
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
    //  13 14    n p
// ================================================================================
    Int_t   ncont=0, n=0, nsecp, icont=0, ntsep=0;
    Int_t   tag=-1, icshow, iclust=0, cmult, cmultp1,
            iclems=0, iclmus=0, iclmxs=0, iclots=0,    // cluster count in showers
            iclemt=0, iclmut=0, iclmxt=0, iclott=0;    // total cluster count
    Int_t   pid, idp1, idp2, pidfst, pidlst;
    Int_t   ngam=0, nele=0, nmu=0, nn=0, np=0, noth=0, ngclst =0, neclst=0, nmclst=0;
    Int_t   nel50=0, nel70=0, nel100=0, nel150=0, nel250=0;
    Int_t   ntgam=0, ntele=0, ntmu=0, ntn=0, ntp=0, ntoth=0;
    Int_t   ngams=0, neles=0, nmus=0, nns=0, nps=0, nots=0;
    Int_t   ifg, ife=0, ifm=0;
    Long64_t idgam=1, idele=1000, idmu=1000000, idn=100000000, ido=100000000;
    Long64_t idp=0, id1=0, id2= 0, idclus=0, sidclst=0, tidclus=0;
    Float_t rx, ry, dx, dy, x1, x2, y1, y2, r, rsq,
            xf=0, xl=0, yf=0, yl=0, dfl=0, afl=0,
            xfe=0, xle=0, yfe=0, yle=0,
            xfm=0, xlm=0, yfm=0, ylm=0;
    Float_t theta=0., zen1=0., zen2=0., phi=0, azh1=0, azh2=0;
    Float_t time =0., t0, t1=0., t2=0.,
            tfst, tlst, tfstg=0, tlstg=0, tfste=0, tlste=0, tfstm=0, tlstm=0,
            zhf=0, zhl=0, azf=0, azl=0,
            zhfste=0, zhlste=0, zhfstm=0, zhlstm=0,
            azfste=0, azlste=0, azfstm=0, azlstm=0;
    Float_t px, py, pz, pt, px1, py1, pz1, px2, py2, pz2, pmod,
            pm1, pm2, pmfste, pmlste, pmfstm, pmlstm, psq; // particle momenta
    Float_t dxsq, dysq, drsq, dr, sigr=0;
    Float_t xolde=0, xolesq=0, xmeane=0, xmnesq=0, sgxesq=0, sigre =0;
    Float_t yolde=0, yolesq=0, ymeane=0, ymnesq=0, sgyesq=0;
    Float_t xoldm=0, xolmsq=0, xmeanm=0, xmnmsq=0, sgxmsq=0, sigrm =0;
    Float_t yoldm=0, yolmsq=0, ymeanm=0, ymnmsq=0, sgymsq=0;
    Float_t xmean = 0, xmeansq = 0, xmold=0, xmoldsq=0, sigxsq=0;
    Float_t ymean = 0, ymeansq = 0, ymold=0, ymoldsq=0, sigysq=0;
    Float_t tmean = 0, tmeansq = 0, tmold=0, tmoldsq=0, sigt=0, sigtsq=0;
    Float_t nxf, nyf, nzf, nxl, nyl, nzl;
    Float_t hghtfi;
    Float_t rad2g;
    rad2g = 180/TMath::Pi();

    fstream file1;
    fstream file2;
    file1.open(fout1, fstream::out); //-    Main output
    file2.open(fout2, fstream::out); //-    Shower summary

        file1 <<
        "#PrimCR"    << "\t" <<
        "EnePCR"    << "\t" <<
        "ZthPCR"    << "\t" <<
        "HghPCR"    << "\t" <<
        "IShow"     << "\t" <<
        "IClust"    << "\t" <<
        "IdClust"   << "\t" <<
        "SidClst"   << "\t" <<    // Short ID of cluster
        "NparClt"   << "\t" <<
        "NelClst"   << "\t" <<
        "NmuClst"   << "\t" <<
        "PidfCst"   << "\t" <<
        "PidlCst"   << "\t" <<
        "XmClust"   << "\t" <<
        "YmClust"   << "\t" <<
        "RWdtClt"   << "\t" <<
        "TfClust"   << "\t" <<
        "TlClust"   << "\t" <<
        "TmClust"   << "\t" <<
        "sTClust"   << "\t" <<
        "TFirstE"   << "\t" <<
        "TLastE"    << "\t" <<
        "ZhFrstE"   << "\t" <<
        "AzFrstE"   << "\t" <<
        "ZhLastE"   << "\t" <<
        "AzlastE"   << "\t" <<
        "PmFrstE"   << "\t" <<
        "PmLastE"   << "\t" <<
        "TFrstM"    << "\t" <<
        "TLastM"    << "\t" <<
        "ZhfrstM"   << "\t" <<
        "AzfrstM"   << "\t" <<
        "ZhlastM"   << "\t" <<
        "AzlastM"   << "\t" <<
        "PmFrstM"   << "\t" <<
        "PmLastM"   << "\t" <<
        "DisFtLt"   << "\t" <<
        "AngFtLt"   << endl;
    
        file2 <<
        "#IShow"  << "\t" <<
        "NSecPS"  << "\t" <<
        "NClstS"  << "\t" <<
        "NGamS"   << "\t" <<
        "NeleS"   << "\t" <<
        "%50MV"   << "\t" <<
        "%70MV"   << "\t" <<
        "%100MV"  << "\t" <<
        "%150MV"  << "\t" <<
        "%250MV"  << "\t" <<
        "SigRES"  << "\t" <<
        "NMuS"    << "\t" <<
        "SigRMS"  << "\t" <<
        "NNeutS"  << "\t" <<
        "NProtS"  << "\t" <<
        "NOthrS"  << "\t" <<    // Short ID of cluster
        "NClElS"  << "\t" <<
        "NClMuS"  << "\t" <<
        "NClMxS"  << "\t" <<
        "NClOtS"  << endl;
// ********************************************************************************
    
   for (Long64_t ishow=0; ishow<nshows; ishow++) {     // Loop in showers
       
       //cout << endl;
       //cout << " ***** new shower " << endl;
       //cout << endl;
       
       icont = -1;   // index for total nb. of secondaries
       Long64_t itree = LoadTree(ishow);
       
       if (itree < 0) break;

       nb = fChain->GetEntry(itree);   nbytes += nb;
       // if (Cut(itree) < 0) continue;

       fb  = b_shower_FirstHeight->GetEntry(itree);   fbytes += fb;
       hb  = b_shower_Theta->GetEntry(itree);         hbytes += hb;
       tb  = b_shower_Phi->GetEntry(itree);           tbytes += tb;

       hghtfi  = shower_FirstHeight/100000; // h en kilometros

       Float_t arr0[100000] = { [0 ... 99999] = -10. };
       Float_t arr1[100000] = { [0 ... 99999] = -10. };
       Float_t arr2[100000] = { [0 ... 99999] = -10. };
       Float_t arr3[100000] = { [0 ... 99999] = -10. };
       Float_t arr4[100000] = { [0 ... 99999] = -10. };
       Float_t arr5[100000] = { [0 ... 99999] = -10. };
       Float_t arr6[100000] = { [0 ... 99999] = -10. };
       Float_t itag[100000] = { [0 ... 99999] =  -1. };
       //Float_t mshow[nshows][14];
       
       // First loop in particles
       cmult = 1;
       icshow = 0;
       // particle__ = 10;  //-
       nsecp = particle__;  // nb. secondary part. in shower
       
       // cout << " nsecp " << nsecp << endl;
       
       // Look for the fastest particle in the shower
       t0 = 1000000;
       for(Int_t ip=0; ip < particle__; ip++){
           time      = particle__Time[ip];
           if(time < t0){
               t0 = time;
           }
           // cout << t0 << endl;
       }
       
       // scans all particles
       for(Int_t ip=0; ip < particle__; ip++){
           rx        = particle__x[ip]/100;   // meters
           ry        = particle__y[ip]/100;   // meters
           r         = sqrt(rx*rx + ry*ry);
           time      = particle__Time[ip] - t0;
           px        = particle__Px[ip];
           py        = particle__Py[ip];
           pz        = particle__Pz[ip];
           psq       = px*px + py*py + pz*pz;
           pmod      = sqrt(psq);
           //pt        = ((px*y)-(py*x))/(r*psq); ???
           pt        = sqrt(px*px + py*py);
           theta     = acos(pz/pmod); // * rad2g;
           phi       = atan2(py,px); // * rad2g;
           
           pid = particle__ParticleID[ip];

   //      Select type of particle
           if(pid==1){ // GAMMAS
               idp = idgam;
               ngam++;
               ntgam++;
               ngams++;
           }
           else if( (pid==2) || (pid==3) ){ // ELECTRONS
               idp = idele;
               nele++;
               ntele++;
               neles++;
               if (pmod < 0.07){nel50++;}
               else if (pmod>0.07 && pmod<0.1){nel70++;}
               else if (pmod>0.1  && pmod<0.15){nel100++;}
               else if (pmod>0.14 && pmod<0.25){nel150++;}
               else {nel250++;}
               
               //   Electron radial distribution
               xolde  = xmeane;
               yolde  = ymeane;
               xolesq = xolde * xolde;
               yolesq = yolde * yolde;
               xmeane = xmeane + (rx - xmeane)/neles;
               ymeane = ymeane + (ry - ymeane)/neles;
               xmnesq = xmeane * xmeane;
               ymnesq = ymeane * ymeane;
               sgxesq = sgxesq + xolesq - xmnesq +
                       ((rx*rx - sgxesq - xolesq)/neles);
               sgyesq = sgyesq + yolesq - ymnesq +
                       ((ry*ry - sgyesq - yolesq)/neles);
               sigre   = sqrt((sgxesq + sgyesq/2)) ;
           }
           else if( (pid==5) || (pid==6) ){ // MUONS
               idp = idmu;
               nmu++;
               ntmu++;
               nmus++;
               //   Muon radial distribution
               xoldm  = xmeanm;
               yoldm  = ymeanm;
               xolmsq = xoldm * xoldm;
               yolmsq = yoldm * yoldm;
               xmeanm = xmeanm + (rx - xmeanm)/nmus;
               ymeanm = ymeanm + (ry - ymeanm)/nmus;
               xmnmsq = xmeanm * xmeanm;
               ymnmsq = ymeanm * ymeanm;
               sgxmsq = sgxmsq + xolmsq - xmnmsq +
                       ((rx*rx - sgxmsq - xolmsq)/nmus);
               sgymsq = sgymsq + yolmsq - ymnmsq +
                       ((ry*ry - sgymsq - yolmsq)/nmus);
               sigrm   = sqrt((sgxmsq + sgymsq)/2) ;
               
           }
           else if( pid==13 ){ // NEUTRONS
               idp = idn;   // nucleon
               nn++;
               ntn++;
               nns++;
           }
           else if( pid==14){ // PROTONS
               idp = idn;   // nucleon
               np++;
               ntp++;
               nps++;
           }
           else{
               idp = ido;
               noth++;
               ntoth++;
               nots++;
           }   // endif in pid
           
           icont++;     // nb of particles

           arr0[icont] = idp;
           arr1[icont] = rx;
           arr2[icont] = ry;
           arr3[icont] = time;
           arr4[icont] = pmod;
           arr5[icont] = theta;
           arr6[icont] = phi;
           
       //}
       }    // end of first scan of particles

       ntsep  = ntsep + nsecp;   // total nb. secondary part.
              
       // Loop in the first particle
       for(Int_t ifp=0; ifp<nsecp; ifp++){

           //if(itag[ifp]!=-1) continue;  // Creo que no hace falta

           idp1  = arr0[ifp];
           x1    = arr1[ifp];
           y1    = arr2[ifp];
           t1    = arr3[ifp];
           pm1   = arr4[ifp];
           zen1  = arr5[ifp];
           azh1  = arr6[ifp];
           
           xf = x1; xl = x1; yf = y1; yl = y1; xmean = x1; ymean = y1;
           zhf = zen1; zhl = zen1; azf = azh1; azl = azh1;
           tfst = t1; tlst = t1; tmean = t1; sigt=0; sigr = 0;
           pidfst = idp1; pidlst = idp1;
           
           if(idp1 == idgam){        // gamma
               tfstg = t1;
               tlstg = t1;
               ngclst++;
               ifg = 1;
           }
           
           if(idp1 == idele){        // electron
               xfe     = x1;
               xle     = x1;
               yfe     = y1;
               yle     = y1;
               tfste   = t1;
               tlste   = t1;
               pmfste  = pm1;
               pmlste  = pm1;
               zhfste  = zen1;
               zhlste  = zen1;
               azfste  = azh1;
               azlste  = azh1;
               neclst++;
               ife = 1;
           }
           
           if(idp1 == idmu){      // muon
               xfm     = x1;
               xlm     = x1;
               yfm     = y1;
               ylm     = y1;
               tfstm   = t1;
               tlstm   = t1;
               pmfstm  = pm1;
               pmlstm  = pm1;
               zhfstm  = zen1;
               zhlstm  = zen1;
               azfstm  = azh1;
               azlstm  = azh1;
               nmclst++;
               ifm = 1;
           }
                      
           //itag[ifp]++;
           
           // Look for clusters in the remaining particles
           for(Int_t isp=ifp+1; isp < nsecp; isp++){
               if(itag[isp]!=-1) continue;

               dx   = arr1[isp] - xmean; dxsq = dx * dx ;
               dy   = arr2[isp] - ymean; dysq = dy * dy ;
               drsq = dxsq + dysq; dr = sqrt(drsq);
           
               // new cluster and new particle
               if( dr <= distmx && sigr <= sigrmx){
                   
                   if(cmult==1){idclus = idp1;}
                   
                   tag++;   // ?
                   cmult++;
                   
                   idp2   = arr0[isp]; idclus = idclus + idp2;
                   x2   = arr1[isp]; y2 = arr2[isp]; t2 = arr3[isp];
                   pm2  = arr4[isp];
                   zen2 = arr5[isp]; azh2 = arr6[isp];
                   
                   // cout << " --- cmult - idp1, idp2 " << cmult << " " << idp1 << " " << idp2 << endl;
                   
                   if(t2<tfst){
                       xf = x2; yf = y2; tfst = t2;
                       zhf = zen2; azf = azh2;
                       pidfst = idp2;

                   }
                   
                   if(t2>tlst){
                       xl = x2; yf = y2; tlst = t2;
                       zhl = zen2; azl = azh2;
                       pidlst = idp2;
                   }
                
                   if(idp2 == idgam){         // gamma
                       ngclst++;
                       if(ifg==0){
                           tfstg   = t2;
                           tlstg   = t2;
                           ifg = 1;
                       }
                       if(t2<tfstg){ tfstg  = t2;}
                       if(t2>tlstg){ tlstg  = t2;}
                   }
                   if(idp2 == idele){        // electron
                       neclst++;
                       if(ife==0){
                           xfe     = x2;
                           xle     = x2;
                           yfe     = y2;
                           yle     = y2;
                           tfste   = t2;
                           tlste   = t2;
                           pmfste  = pm2;
                           pmlste  = pm2;
                           zhfste  = zen2;
                           zhlste  = zen2;
                           azfste  = azh2;
                           azlste  = azh2;
                           ife = 1;
                       }
                       if(t2<tfste){
                           xfe   = x2; yfe = y2; tfste = t2;
                           pmfste = pm2; zhfste = zen2; azfste = azh2;
                       }
                       if(t2>tlste){
                           xle    = x2; yle = y2; tlste = t2;
                           pmlste = pm2; zhlste = zen2; azlste = azh2;
                       }
                   }
                   if(idp2 == idmu){      // muon
                       nmclst++;
                       if(ifm == 0){
                           xfm     = x2;
                           xlm     = x2;
                           yfm     = y2;
                           ylm     = y2;
                           tfstm   = t2;
                           tlstm   = t2;
                           pmfstm  = pm2;
                           pmlstm  = pm2;
                           zhfstm  = zen2;
                           zhlstm  = zen2;
                           azfstm  = azh2;
                           azlstm  = azh2;
                           ifm = 1;
                       }
                       if(t2<tfstm){
                           xfm    = x2; yfm = y2; tfstm = t2;
                           pmfstm = pm2, zhfstm = zen2; azfstm = azh2;
                       }
                       
                       if(t2>tlstm){
                           xlm    = x2; ylm = y2; tlstm = t2;
                           pmlstm = pm2, zhlstm = zen2; azlstm = azh2;
                       }
                   }
                   
                   // Calculo del nuevo centro del cluster y su anchura
                   // Formulas recursivas del valor medio y la dispersion
                   
                   xmold   = xmean;
                   ymold   = ymean;
                   tmold   = tmean;
                   xmoldsq = xmold * xmold;
                   ymoldsq = ymold * ymold;
                   tmoldsq = tmold * tmold;
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
                   sigtsq = sigtsq + tmoldsq - tmeansq +
                            ((t2*t2 - sigtsq - tmoldsq)/cmult);
                   sigt   = sqrt(sigtsq);
                   sigr   = sqrt((sigxsq + sigysq)/2) ;
                   
               }  // ends adding particle to the cluster
               
               // cout << " * lab cmult " << tag << " " << cmult << endl;
               
               itag[isp] = tag;
               tag       = -1;
           }  //  End of isp for;
           
           if (cmult > 1){      //  New cluster found
               iclust++;
               icshow ++;
               tidclus = tidclus + idclus;  // suma logica de los clusters en el shower
               
               // Clasificamos el cluster . Ignoramos las particulas pesadas
               
               if (idclus < idmu && neclst > 1){   // clean EM cluster with >1 electron
                   sidclst = 1;
                   iclems ++;
                   iclemt ++;
                   xf = xfe; yf = yfe; xl = xle; yl = yle;
                   zhf = zhfste; zhl = zhlste; azf = azfste; azl = azlste;
               }
               else if (nmclst > 1 && neclst == 0 && ngclst == 0){   //  cluster with >1 muons
                   sidclst = 2;
                   iclmus ++;
                   iclmut ++;
                   xf = xfm; yf = yfm; xl = xlm; yl = ylm;
                   zhf = zhfstm; zhl = zhlstm; azf = azfstm; azl = azlstm;
               }
               else if (neclst>0 && nmclst>0) {          // mixt cluster with e's & mu's
                   sidclst = 3;
                   iclmxs ++;
                   iclmxt ++;
                   if(tfste<tfstm){
                      xf  = xfe; yf = yfe;
                      zhf = zhfste ; azf = azfste ;
                   }
                   else{
                      xf  = xfm; yf = yfm;
                      zhf = zhfstm ; azf = azfstm ;
                   }
                   if(tlste>tlstm){
                      xl  = xle; yl = yle;
                      zhl = zhlste ; azl = azlste ;
                   }
                   else{
                      xl  = xlm; yl  = ylm;
                      zhl = zhlstm ; azl = azlstm ;
                   }
               }
               else{iclots ++; iclott ++;}
               
               nxf = sin(zhf) * cos(azf) ;
               nyf = sin(zhf) * sin(azf) ;
               nzf = cos(zhf);
               nxl = sin(zhl) * cos(azl) ;
               nyl = sin(zhl) * sin(azl) ;
               nzl = cos(zhl);
               afl = acos(nxf*nxl + nyf*nyl + nzf*nzl); // * rad2g;
               dfl = sqrt((xl-xf)*(xl-xf) + (yl-yf)*(yl-yf));
 
               /*   ---  New cluster
               cout << "*** ICluster " << icshow << " - itclust " << iclust << endl;
               cout << "* ClMult "  << cmult <<  " - ClustID " << idclus << endl;
               cout << "--- ShClustID (0-4):   " << sidclst << endl;
               cout << "-- iclems, iclmus, iclmxs, iclots: " << iclems << " " << iclmus << " " << iclmxs << " " << iclots << endl;
               //cout << "---dfl, afl:   " << dfl << " " << afl << endl;
               //cout<< "--- tfst tlst: " << tfst << "   "<< tlst << endl;
               //cout << "(xmean, ymean), sigr " << xmean << " " << ymean << " " << sigr << endl;
               cout << endl;
               */
               
               //*
               file1 << mpcr << "\t"
                  << epcr   << "\t"
                  << zhpcr  << "\t"
                  << hghtfi << "\t"
                  << ishow+1<< "\t"
                  << iclust << "\t"
                  << idclus << "\t"
                  << sidclst<< "\t"
                  << cmult  << "\t"
                  << neclst << "\t"
                  << nmclst << "\t"
                  << pidfst << "\t"
                  << pidlst << "\t"
                  << xmean  << "\t"
                  << ymean  << "\t"
                  << sigr   << "\t"
                  << tfst   << "\t"
                  << tlst   << "\t"
                  << tmean  << "\t"
                  << sigt   << "\t"
                  << tfste  << "\t"
                  << tlste  << "\t"
                  << zhfste << "\t"
                  << azfste << "\t"
                  << zhlste << "\t"
                  << azlste << "\t"
                  << pmfste << "\t"
                  << pmlste << "\t"
                  << tfstm  << "\t"
                  << tlstm  << "\t"
                  << zhfstm << "\t"
                  << azfstm << "\t"
                  << zhlstm << "\t"
                  << azlstm << "\t"
                  << pmfstm << "\t"
                  << pmlstm << "\t"
                  << dfl    << "\t"
                  << afl    << endl;
                    //*/
               
           }  // endif cmult > 1:  new cluster
                   
           tag     = -1; idclus  = 0.; sidclst=0; tmean=0, cmult = 1.;
           sigxsq  = 0.; sigysq  = 0.; sigtsq  = 0;
           ifg     = 0; ife     = 0; ifm     = 0;
           ngclst  = 0; neclst  = 0; nmclst  = 0;
           tfstg   = 0; tlstg   = 0; tfste   = 0; tlste   = 0;
           zhfste  = 0; zhlste  = 0; azfste  = 0; azlste  = 0;
           pmfste  = 0; pmlste  = 0;
           tfstm   = 0; tlstm   = 0;
           zhfstm  = 0; zhlstm  = 0; azfstm  = 0; azlstm  = 0;
           nxf = 0; nxl= 0; nyf = 0; nyl= 0; nzf = 0; nzl= 0;
           xf = 0; xl = 0; yf = 0; yl = 0;
           
       }  // ends ifp loop

       for(Int_t i=0; i < 100000; i++){
           itag[i]=-1;
           arr0[i]=-10.;
           arr1[i]=-10.;
           arr2[i]=-10.;
           arr3[i]=-10.;
           arr4[i]=-10.;
           arr5[i]=-10.;
           arr6[i]=-10.;
       }
       
       /*
       //cout << endl;
       cout << "*** Shower summary " << endl;
       cout << "* IShower - NClust - NSecP:  " << ishow+1 << "  -  " << icshow <<  "  -  " << nsecp << endl;
       //cout << "* Suma de clusters IDs (onmmeeeggg): " << tidclus << endl;
       cout << "* ShIdClusters: 2Ele, 2MU, 2Ch, Oth: " << iclems << " " << iclmus << " " << iclmxs << " " << iclots << endl;
       cout << "* Distribucion (gemnpo): " << ngams << " " << neles << " " << nmus <<  " " << nns <<  " " << nps <<  " " << nots << endl;
           cout << "* Distribucion relativa/%: " << 100*ngams/nsecp << " " << 100*neles/nsecp << " " << 100*nmus/nsecp <<  " " << 100*nns/nsecp <<  " " << 100*nps/nsecp <<  " " << 100*nots/nsecp << endl;
       cout << endl;
       */
       
       if (neles==0){neles =-1;}
       //*
       file2 <<ishow+1<< "\t"
          << nsecp   << "\t"
          << icshow  << "\t"
          << ngams   << "\t"
          << neles   << "\t"
          << 100*nel50/neles  << "\t"
          << 100*nel70/neles  << "\t"
          << 100*nel100/neles << "\t"
          << 100*nel150/neles << "\t"
          << 100*nel250/neles << "\t"
          << sigre   << "\t"
          << nmus    << "\t"
          << sigrm   << "\t"
          << nns     << "\t"
          << nps     << "\t"
          << nots    << "\t"
          << iclems  << "\t"
          << iclmus  << "\t"
          << iclmxs  << "\t"
          << iclots  << endl;
        //*/
       
       // cout << "SigREls - SigRMus" << sigre << "  " << sigrm << endl;
       
       ngam = 0; nele=0; nmu=0; nn=0; np=0; noth=0;
       nel50=0; nel70=0; nel100=0; nel150=0; nel250=0;
       ngams=0; neles=0; nmus=0; nns=0; nps=0; nots=0;
       xolde=0; xolesq=0; xmeane=0; xmnesq=0; sgxesq=0; sigre =0;
       yolde=0; yolesq=0; ymeane=0; ymnesq=0; sgyesq=0;
       
       tidclus = 0;
       iclems  = 0;  // count of EM clusters in the shower
       iclmus  = 0;  // count of muon clusters in the shower
       iclmxs  = 0;  // count of mixed clusters in the shower
       iclots  = 0;  // count of other type of clusters
       
   }  // end for in showers

    /*
    for (Long64_t ishow=0; ishow<nshows; ishow++) {     // Loop in showers
        mshow[ishow][6]= sigrm;
    }
     */
    
cout << endl;
printf("******  Proceso completado  ******");
file1.close();
file2.close();
    
    if (ntsep==0){ntsep=-1;}
    
cout << " - Parametros iniciales. DistMx - SigrMx: " << distmx << "  " << sigrmx << endl;
cout << "* NShowers: " << nshows << endl;
//cout << endl;
cout << "* Numero total de clusters: " << iclust << endl;
cout << "* Numero total de secundarios: " << ntsep << endl;
cout << "* Contador (gemnpo): " << ntgam << " " << ntele << " " << ntmu <<  " " << ntn <<  " " << ntp <<  " " << ntoth << endl;
cout << "* Distribucion/% (gemnpo): " << 100*ntgam/ntsep << " " << 100*ntele/ntsep << " " << 100*ntmu/ntsep <<  " " << 100*ntn/ntsep <<  " " << 100*ntp/ntsep <<  " " << 100*ntoth/ntsep << endl;
cout << "* NParticulas/shower (gemnpo): " << ntgam/nshows << " " << ntele/nshows << " " << ntmu/nshows <<  " " << ntn/nshows <<  " " << ntp/nshows <<  " " << ntoth/nshows << endl;
cout << "* ShIdClust: 2Ele, 2MU, 2Ch, Oth: " << iclemt << " " << iclmut << " " << iclmxt << " " << iclott << endl;
cout << "* ShIdClust/shower: 2Ele, 2MU, 2Ch, Oth: " << iclemt/nshows << " " << iclmut/nshows << " " << iclmxt/nshows << " " << iclott/nshows << endl;
cout << endl;

//TH2D *h1 = new TH2D("h1", "cluster distribution", 100, 0., 1000., 100, 0., 1000.);
    
// Estimamos densidades de particulas y sigmas correspondientes

}
/*
 void bigpm() {
    TFile *f = new TFile("cernstaff.root");
    TTree *T = (TTree*)f->Get("T");
    TCanvas *c1 = new TCanvas("c1");
    T->Draw("Cost:Age>>hist","","goff");
    TH2F *h = new TH2F("h","Cost vs Age",60,10,70,20,0,20000);
    h->Draw();
    c1->Update();
    T->SetMarkerStyle(3);
    T->SetMarkerColor(3); T.Draw("Cost:Age","Grade==3","same");
    T->SetMarkerColor(4); T.Draw("Cost:Age","Grade==4","same");
    T->SetMarkerColor(5); T.Draw("Cost:Age","Grade==5","same");
    T->SetMarkerColor(6); T.Draw("Cost:Age","Grade==6","same");
    T->SetMarkerColor(1); T.Draw("Cost:Age","Grade==10","same");
    T->SetMarkerColor(2); T.Draw("Cost:Age","Grade==12","same");
 }
 */

