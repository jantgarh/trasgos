//  mclusth.C  Analisis de clusters con interaccion en diferentes alturas
//Programa muy robusto. ANALISIS DE CLUSTERS. Identificacion de particulas a menos de una cierta distancia entre ellas.
// Programa para datos con alturas fijas
// Elige los clusters mediante el algoritmo iterativo puesto a punto por Yanis

#define mclusth_cxx
#include "mclusth.h"
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
#define distmx 2.0  //m. Distancia maxima de separación entre partículas
#define sigrmx 1.0  //m. Dispersion radial maxima del cluster
#define mpcr 1  // A mass of primary cosmic ray
//#define epcr 2  // Log(Energy/GeV) primary cosmic ray
//#define dene 0.25 // Energy interval for the total flux estimation
// #define hghtfi 15 // Height/km of first interaction
#define zhpcr 0 // Zenith angle of primary cosmic ray
#define spindex -2.7  // Spectral index
#define f0 1.8E4  // Intercept hidrogen
//#define f0he 5.0E3 // Intercept helium
//#define f0mZ 1.0E3 // Intercept medium Z
//#define f0hZ 1.5E3 // Intercept hight Z

//#define fout1 "soclst_p_02_15km.txt" //   cluster output summary
//#define fout2 "soshws_p_02_15km.txt" //   shower output summary
//#define fout3 "sopart_p_02_15km.txt" //   particle output summary
#define fout1 "clstsum.txt" //
#define fout2 "showsum.txt" //
#define fout3 "partsum.txt" //

void mclusth::Loop()
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
//      root> .L mclusth.C
//      root> mclusth t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//    Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    
   if (fChain == 0) return;
    
   Long64_t nshows = fChain->GetEntriesFast();
    
   cout << " nShower in file: " << nshows << endl;
   nshows = 10;   //-
   cout << " nShowers analyzed " << nshows << endl;

   Long64_t nbytes = 0, nb = 0;
   Long64_t fbytes = 0, fb = 0;
   Long64_t ebytes = 0, eb = 0;
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
    Int_t   nel50=0, nel100=0, nel150=0, nel200=0, nel300=0,
    nel500, nel1000, nel2000;
    Int_t   ntel50=0, ntel100=0, ntel150=0, ntel200=0, ntel300=0,
    ntel500, ntel1000, ntel2000;
    Int_t   ntgam=0, ntele=0, ntmu=0, ntn=0, ntp=0, ntoth=0;
    Int_t   ngams=0, neles=0, nmus=0, nns=0, nps=0, nots=0;
    Int_t   ifg, ife=0, ifm=0;
    Long64_t idgam=1, idele=1000, idmu=1000000, idn=100000000, ido=100000000;
    Long64_t idp=0, id1=0, id2= 0, idclus=0, sidclst=0, tidclus=0;
    Float_t rx, ry, dx, dy, x1, x2, y1, y2, r, rsq,
            xf=0, xl=0, yf=0, yl=0, distfl=0, angfl=0,
            xfe=0, xle=0, yfe=0, yle=0,
            xfm=0, xlm=0, yfm=0, ylm=0;
    Float_t theta=0., zen1=0., zen2=0., phi=0, azh1=0, azh2=0;
    Float_t time =0., t0, t1=0., t2=0.,
            tfst, tlst, tfstg=0, tlstg=0, tfste=0, tlste=0, tfstm=0, tlstm=0,
            zhf=0, zhl=0, azf=0, azl=0,
            zhfste=0, zhlste=0, zhfstm=0, zhlstm=0,
            azfste=0, azlste=0, azfstm=0, azlstm=0;
    Float_t px=0, py=0, pz=0, pt=0, px1=0, py1=0, pz1=0, px2=0, py2=0, pz2=0,
            pmod, pm1, pm2, psq,
            pmfste, pmlste, pmfstm, pmlstm; // particle momenta
    Float_t dxsq, dysq, drsq, dr, sigr=0;
    Float_t xolde=0, xolesq=0, xmeane=0, xmnesq=0, sgxesq=0, sigre =0;
    Float_t yolde=0, yolesq=0, ymeane=0, ymnesq=0, sgyesq=0;
    Float_t xoldm=0, xolmsq=0, xmeanm=0, xmnmsq=0, sgxmsq=0, sigrm =0;
    Float_t yoldm=0, yolmsq=0, ymeanm=0, ymnmsq=0, sgymsq=0;
    Float_t xmean = 0, xmeansq = 0, xmold=0, xmoldsq=0, sigxsq=0;
    Float_t ymean = 0, ymeansq = 0, ymold=0, ymoldsq=0, sigysq=0;
    Float_t tmean = 0, tmeansq = 0, tmold=0, tmoldsq=0, sigt=0, sigtsq=0;
    Float_t nxf, nyf, nzf, nxl, nyl, nzl;
    Float_t hghtfi, mnhgt=0;  // first height, mean hght
    Float_t emin, emax, gm1, flint;   // energy limits, spindex-1, flux integral
    Float_t rad2g;
    rad2g = 180/TMath::Pi();
    
    Float_t epcr, dene=0.25, mepcr=0;
    

    fstream file1;
    fstream file2;
    fstream file3;
    file1.open(fout1, fstream::out); //-    Clusters summary
    file2.open(fout2, fstream::out); //-    Shower summary
    file3.open(fout3, fstream::out); //-    Particle summary

    file1 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    file2 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    file3 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;

    file1 <<
    "#PrCR"    << "\t" <<
    "EnePCR"    << "\t" <<
    "ZthPCR"    << "\t" <<
    "HghPCR"    << "\t" <<
    "IShow"     << "\t" <<
    "IClust"    << "\t" <<
    "IdClst"   << "\t" <<
    "SIdClt"   << "\t" <<    // Short ID of cluster
    "NparCt"   << "\t" <<
    "NelClt"   << "\t" <<
    "NmuClt"   << "\t" <<
    "PidFPr"   << "\t" <<
    "PidLPr"   << "\t" <<
    "XmClst"   << "\t" <<
    "YmClst"   << "\t" <<
    "RadWCt"   << "\t" <<
    "TfClst"   << "\t" <<
    "TlClst"   << "\t" <<
    "TmClst"   << "\t" <<
    "sTClst"   << "\t" <<
    "TFstEl"   << "\t" <<
    "TLstEl"    << "\t" <<
    "ZhFstE"   << "\t" <<
    "AzFstE"   << "\t" <<
    "ZhLstE"   << "\t" <<
    "AzlstE"   << "\t" <<
    "PmFstE"   << "\t" <<
    "PmstE"   << "\t" <<
    "TFrstM"    << "\t" <<
    "TLastM"    << "\t" <<
    "ZhfstM"   << "\t" <<
    "AzfstM"   << "\t" <<
    "ZhlstM"   << "\t" <<
    "AzlstM"   << "\t" <<
    "PmFstM"   << "\t" <<
    "PmLstM"   << "\t" <<
    "DistFL"   << "\t" <<
    "AnglFL"   << endl;

    file2 <<
    "#PrimCR" << "\t" <<
    "EnePCR"  << "\t" <<
    "ZthPCR"  << "\t" <<
    "HghPCR"  << "\t" <<
    "IShow"   << "\t" <<
    "NSPar"   << "\t" <<
    "NClst"   << "\t" <<
    "NGam"    << "\t" <<
    "Nele"    << "\t" <<
    "RSigEl"  << "\t" <<
    "NMuS"    << "\t" <<
    "RSigMu"  << "\t" <<
    "NNeut"   << "\t" <<
    "NProt"   << "\t" <<
    "NOthr"   << "\t" <<
    "NClEl"   << "\t" <<
    "NClMu"   << "\t" <<
    "NClMx"   << "\t" <<
    "NClRm"   << "\t" <<
    "Nel50S"  << "\t" <<
    "Ne100S"  << "\t" <<
    "Ne150S"  << "\t" <<
    "Ne200S"  << "\t" <<
    "Ne300S"  << "\t" <<
    "Ne500S"  << "\t" <<
    "Nel1kS"  << "\t" <<
    "Nel2kS"  <<endl ;

    file3 <<
    "#PrimCR"  << "\t" <<
    "EnePCR"   << "\t" <<
    "ZthPCR"   << "\t" <<
    "HghPCR"   << "\t" <<
    "NShow"    << "\t" <<
    "NtSePs"   << "\t" <<
    "NtClts"   << "\t" <<
    "NtGam"    << "\t" <<
    "Ntele"    << "\t" <<
    "NMu"      << "\t" <<
    "NNeut"    << "\t" <<
    "NProt"    << "\t" <<
    "NOthr"    << "\t" <<
    "NClEl"    << "\t" <<
    "NClMu"    << "\t" <<
    "NClMx"    << "\t" <<
    "NClRm"    << "\t" <<
    "Ntel50"   << "\t" <<
    "Nte100"   << "\t" <<
    "Nte150"   << "\t" <<
    "Nte200"   << "\t" <<
    "Nte300"   << "\t" <<
    "Nte500"   << "\t" <<
    "Ntel1K"   << "\t" <<
    "Ntel2K"   <<endl ;

// *******************************************************************************
//nshows = 1;   //-
// cout << "nshows " << nshows << endl;
       
   for (Long64_t ishow=0; ishow<nshows; ishow++) {

       // cout << " *** ishow " << ishow << endl;
       
       Long64_t itree = LoadTree(ishow);
       
       if (itree < 0) break;
       
       nb = fChain->GetEntry(itree);   nbytes += nb;
       
       //if (Cut(ientry) < 0) continue;
       //eb  = b_shower_Energy->GetEntry(itree);        ebytes += eb;
       fb  = b_shower_FirstHeight->GetEntry(itree);   fbytes += fb;
       hb  = b_shower_Theta->GetEntry(itree);         hbytes += hb;
       tb  = b_shower_Phi->GetEntry(itree);           tbytes += tb;
       
       icont = -1; // index for total nb. of secondaries

       hghtfi = shower_FirstHeight/100000; // h en kilometros
       mnhgt = mnhgt + (hghtfi-mnhgt)/(ishow+1);
       epcr = log10(shower_Energy); // h en kilometros
       mepcr = mepcr + (epcr-mepcr)/(ishow+1);
       
       
       // cout << "LogEPCR, MeanEne " << epcr << " " << mlge << endl;

       Float_t arr0[100000] = { [0 ... 99999] = -10. };
       Float_t arr1[100000] = { [0 ... 99999] = -10. };
       Float_t arr2[100000] = { [0 ... 99999] = -10. };
       Float_t arr3[100000] = { [0 ... 99999] = -10. };
       Float_t arr4[100000] = { [0 ... 99999] = -10. };
       Float_t arr5[100000] = { [0 ... 99999] = -10. };
       Float_t arr6[100000] = { [0 ... 99999] = -10. };
       Float_t itag[100000] = { [0 ... 99999] =  -1. };

       // First loop in particles
       cmult = 1;
       icshow = 0;
       // particle__ = 100;  //-
       nsecp = particle__;  // nb. secondary part. in shower
    
       // Look for the fastest particle in the shower
       t0 = 1000000;
       for(Int_t ip=0; ip < particle__; ip++){
           time      = particle__Time[ip];
           if(time < t0){
               t0 = time;
           }
       }

   //cout << " ishow, nsecp, t0: " << ishow << "   " << nsecp << "   " << t0 << endl;
   
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
       
       //cout << "*** ishow, ip, pid " << ishow << "   " << ip << "   " << pid << endl;
       
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
           
           // cout << "*** ishow, ip, pmod " << ishow << "   " << ip << "   " << pmod << endl;
           
           if (pmod < 0.10){nel50++; ntel50++;}
           else if (pmod>0.10 && pmod<0.15){nel100++; ntel100++;}
           else if (pmod>0.15 && pmod<0.20){nel150++; ntel150++;}
           else if (pmod>0.20 && pmod<0.30){nel200++; ntel200++;}
           else if (pmod>0.30 && pmod<0.50){nel300++; ntel300++;}
           else if (pmod>0.50 && pmod<1.00){nel500++; ntel500++;}
           else if (pmod>1.00 && pmod<2.00){nel1000++;ntel1000++;}
           else {nel2000++; ntel2000++;
           }
           
           //cout << " *** ishow, ip " << ishow << "   " << ip  << endl;
           //cout << "     *** pmod_ele: " << pmod << endl;

           
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
           
           //cout << "*** noth, pid: " << noth << " " << pid << endl;
           
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

       //cout << endl;
       //cout << "********* 1: itag  " << itag[ifp] << endl;
       
       if(itag[ifp]!=-1) continue;

       idp1  = arr0[ifp];
       x1    = arr1[ifp];
       y1    = arr2[ifp];
       t1    = arr3[ifp];
       pm1   = arr4[ifp];
       zen1  = arr5[ifp];
       azh1  = arr6[ifp];
       
       //cout << "*** ifp, idp1 " << ifp << "   " << idp1 << endl;
       
       xf = x1; xl = x1; yf = y1; yl = y1; xmean = x1; ymean = y1;
       zhf = zen1; zhl = zen1; azf = azh1; azl = azh1;
       tfst = t1; tlst = t1; tmean = t1; sigt=0; sigr = 0;
       pidfst = idp1; pidlst = idp1;
       
       if(idp1 == idgam){        // gamma
           ngclst++;
           ifg = 1;
           tfstg = t1;
           tlstg = t1;
       }
       
       if(idp1 == idele){        // electron
           neclst++;
           ife = 1;
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
       }
       
       if(idp1 == idmu){      // muon
           nmclst++;
           ifm = 1;
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
       }
                  
       itag[ifp]++;

       // Look for clusters in the remaining particles
       for(Int_t isp=ifp+1; isp < nsecp; isp++){
           
           //cout << "  ** look for clusters: itag  " << itag[ifp] << endl;
           
           if(itag[isp]!=-1) continue;

           dx   = arr1[isp] - xmean; dxsq = dx * dx ;
           dy   = arr2[isp] - ymean; dysq = dy * dy ;
           drsq = dxsq + dysq; dr = sqrt(drsq);
       
           // new cluster and new particle
           if( dr <= distmx && sigr <= sigrmx){
               
               if(cmult==1){idclus = idp1;}
               
                // cout << "   * new cluster: tag  " << tag << endl;
               
               tag++;   // cluster found
               cmult++;
               
               idp2   = arr0[isp];
               idclus = idclus + idp2;
               x2   = arr1[isp]; y2 = arr2[isp]; t2 = arr3[isp];
               pm2  = arr4[isp];
               zen2 = arr5[isp]; azh2 = arr6[isp];
               
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
           angfl = acos(nxf*nxl + nyf*nyl + nzf*nzl); // * rad2g;
           distfl = sqrt((xl-xf)*(xl-xf) + (yl-yf)*(yl-yf));

           /*   ---  New cluster
           cout << "*** ICluster " << icshow << " - itclust " << iclust << endl;
           cout << "* ClMult "  << cmult <<  " - ClustID " << idclus << endl;
           cout << "--- ShClustID (0-4):   " << sidclst << endl;
           cout << "-- iclems, iclmus, iclmxs, iclots: " << iclems << " " << iclmus << " " << iclmxs << " " << iclots << endl;
           //cout << "---distfl, angfl:   " << distfl << " " << angfl << endl;
           //cout<< "--- tfst tlst: " << tfst << "   "<< tlst << endl;
           //cout << "(xmean, ymean), sigr " << xmean << " " << ymean << " " << sigr << endl;
           cout << endl;
           */
           
           //*
           file1 << mpcr << "\t"
              << mepcr   << "\t"
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
//             << fixed << setprecision(3)
              << setprecision(4)
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
              << distfl << "\t"
              << angfl  << endl;
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
   
       if (neles==0){neles=-1;}
       
   //*
   file2
      << mpcr << "\t"
      << epcr   << "\t"
      << zhpcr  << "\t"
      << hghtfi << "\t"
      << ishow+1<< "\t"
      << nsecp   << "\t"
      << icshow  << "\t"
      << ngams   << "\t"
      << neles   << "\t"
      << setprecision(4)
      << sigre   << "\t"
      << setprecision(0)
      << nmus    << "\t"
      << setprecision(4)
      << sigrm   << "\t"
      << setprecision(0)
      << nns     << "\t"
      << nps     << "\t"
      << nots    << "\t"
      << iclems  << "\t"
      << iclmus  << "\t"
      << iclmxs  << "\t"
      << iclots  << "\t"
      << nel50   << "\t"
      << nel100  << "\t"
      << nel150  << "\t"
      << nel200  << "\t"
      << nel300  << "\t"
      << nel500  << "\t"
      << nel1000 << "\t"
      << nel2000 << endl;
      //*/

     // cout << "SigREls - SigRMus" << sigre << "  " << sigrm << endl;
   
     ngam=0; nele=0; nmu=0; nn=0; np=0; noth=0;
     nel50=0; nel100=0; nel150=0; nel200=0; nel300=0;
     nel500=0;nel1000=0; nel2000=0;
     ngams=0; neles=0; nmus=0; nns=0; nps=0; nots=0;
     xolde=0; xolesq=0; xmeane=0; xmnesq=0; sgxesq=0; sigre =0;
     yolde=0; yolesq=0; ymeane=0; ymnesq=0; sgyesq=0;
   
     tidclus = 0;
     iclems  = 0;  // count of EM clusters in the shower
     iclmus  = 0;  // count of muon clusters in the shower
     iclmxs  = 0;  // count of mixed clusters in the shower
     iclots  = 0;  // count of other type of clusters
       
     // cout << " ishow, nsecp " << ishow << "   " << nsecp  << endl;
       
   } // end of for-loop in showers
    
    emin = mepcr - dene/2;
    emax = mepcr + dene/2;
    gm1  = spindex - 1;
    flint = (f0/gm1) * (pow(emax,-gm1) - pow(emin,-gm1));

file3<< mpcr  << "\t"
    << epcr    << "\t"
    << zhpcr   << "\t"
    << mnhgt   << "\t"
    << nshows  << "\t"
    << ntsep   << "\t"
    << iclust  << "\t"
    << ntgam   << "\t"
    << ntele   << "\t"
    << ntmu    << "\t"
    << ntn     << "\t"
    << ntp     << "\t"
    << ntoth   << "\t"
    << iclemt  << "\t"
    << iclmut  << "\t"
    << iclmxt  << "\t"
    << iclott  << "\t"
    << ntel50   << "\t"
    << ntel100  << "\t"
    << ntel150  << "\t"
    << ntel200  << "\t"
    << ntel300  << "\t"
    << ntel500  << "\t"
    << ntel1000 << "\t"
    << ntel2000 << endl;
    
cout << endl;
printf("******  Proceso completado  ******");
file1.close();
cout << endl;
file2.close();
cout << endl;
file3.close();
cout << endl;
    
cout << " - Parametros iniciales. MxDist - MxSigr: " << distmx << "  " << sigrmx << endl;
cout << "* NShowers: " << nshows << endl;
//cout << endl;
cout << "* Numero total de secundarios: " << ntsep << endl;
cout << "* Contador (gemnpo): " << ntgam << " " << ntele << " " << ntmu <<  " " << ntn <<  " " << ntp <<  " " << ntoth << endl;
cout << "* Distribucion/% (gemnpo): " << 100*ntgam/ntsep << " " << 100*ntele/ntsep << " " << 100*ntmu/ntsep <<  " " << 100*ntn/ntsep <<  " " << 100*ntp/ntsep <<  " " << 100*ntoth/ntsep << endl;
cout << "* NParticulas/shower (gemnpo): " << ntgam/nshows << " " << ntele/nshows << " " << ntmu/nshows <<  " " << ntn/nshows <<  " " << ntp/nshows <<  " " << ntoth/nshows << endl;
    cout << " --- " << endl;
    cout << "* Numero total de clusters: " << iclust << endl;
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

