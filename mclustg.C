//Programa muy robusto. ANALISIS DE CLUSTERS. Identificacion de particulas a menos de una cierta distancia entre ellas.
// Programa para datos con alturas fijas
// Elige los clusters a traves de una cuadricula de tamaño variable


#define mclustg_cxx
#include "mclustg.h"
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
#define epcr 2  // Log(Energy/GeV) primary cosmic ray
#define dene 0.25 // Energy interval for the total flux estimation
// #define hghtfi 15 // Height/km of first interaction
#define zhpcr 0 // Zenith angle of primary cosmic ray
#define spindex 2.7  // Spectral index
#define f0 1.8E4  // Flux Intercept hidrogen /m^2/s/sr
#define f0he 5.0E3 // Flux Intercept helium  /m^2/s/sr
#define f0mZ 1.0E3 // Flux Intercept medium Z /m^2/s/sr
#define f0hZ 1.5E3 // Flux Intercept hight Z  /m^2/s/sr

#define fout1 "xout1.txt" // p_e20_clsts.txt" //
#define fout2 "xout2.txt" //"p_e20_shows.txt" //
#define fout3 "xout3.txt" //"p_e20_parts.txt" //
#define fout4 "xout4.txt" //"p_e20_pgrid.txt" //
#define fout5 "xout5.txt" //"p_e20_cgrid.txt" //
//#define fout6 "xout1.txt" //"p_e20_cgride.txt" //
//#define fout7 "xout1.txt" //"p_e20_cgridm.txt" //
//#define fout8 "xout1.txt" //"p_e20_cgridx.txt" //
#define fout9 "xout9.txt" //"p_e20_rpart.txt" //
#define fouta "xouta.txt" //"p_e20_rclst.txt" //
#define foutb "xoutb.txt" //"p_e20_rclte.txt" //
#define foutc "xoutc.txt" //"p_e20_rcltm.txt" //
#define foutd "xoutd.txt" //"p_e20_rcltx.txt" //

//#define ngrid 1000
#define nbox 200
#define len 1000. //m. Size of the box

// using namespace std;

void mclustg::Loop()
{
// Analisis de cascadas atmosfericas de rayos cosmicos
/*
 *  Copyright LabCAF. IGFAE/USC. All rights reserved.
 *  Created by Juan A. Garzon on 23/01/13
 *  Modified by G.Kornakov on  11/03/2013
 *  Modified by Yanis Fontenla Barba on 23/05/2017 for Corsika and cluster analysis
 *  Improved by JAGarzon on Apr.2021 for cluster analysis
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

   Long64_t nshow = fChain->GetEntriesFast();
    cout << " nshow: " << nshow << endl;
   nshow = 1;   //-
   cout << "nShowers " << nshow << endl;

    // Estimated flux factor integral
    Float_t emin, emax, gm1, flint;   // energy limits, spindex-1, flux integral
    emin = epcr - dene/2;
    emax = epcr + dene/2;
    gm1  = spindex - 1;
    flint= (f0/gm1) * (pow(pow(10,emin),-gm1) - pow(pow(10,emax),-gm1));  // flux integral factor
    
    // cout << " flint  " << flint << endl;
    
    // Identificadores de particulas en Corsika:
    //  1        gammas
    //  2 3      e- e+
    //  5 6      mu- mu+
    //  13 14    n p
// ================================================================================
    Int_t   ncont=0, n=0, nsecp=0, nasep=0, icont=0, ntsep=0;
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
    Float_t rad2g;
    rad2g = 180/TMath::Pi();
    
    fstream ofclst;
    fstream ofshws;
    fstream ofpart;
    fstream ofpgrd;
    fstream ofcgrd;
    //fstream ofcgrde;
    //fstream ofcgrdm;
    //fstream ofcgrdx;
    fstream ofrpart;
    fstream ofrclst;
    fstream ofrclte;
    fstream ofrcltm;
    fstream ofrcltx;
    
    ofclst.open(fout1, fstream::out);  //    Clusters summary
    ofshws.open(fout2, fstream::out);  //    Shower summary
    ofpart.open(fout3, fstream::out);  //    Particle summary
    ofpgrd.open(fout4, fstream::out);  //    Particle grid
    ofcgrd.open(fout5, fstream::out);  //    Cluster grid
    ofrpart.open(fout9, fstream::out); //    Cluster grid
    ofrclst.open(fouta, fstream::out); //    Cluster grid
    ofrclte.open(foutb, fstream::out); //    Cluster grid
    ofrcltm.open(foutc, fstream::out); //    Cluster grid
    ofrcltx.open(foutd, fstream::out); //    Cluster grid

    ofclst << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    ofshws << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    ofpart << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    ofpgrd << "# Particle grid. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofcgrd << "# Cluster grid. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofrpart<< "# Particle radial dist. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofrclst<< "# Clusters radial dist. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofrclte<< "# EM Clusters radial dist. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofrcltm<< "# Muon Clusters radial dist. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    ofrcltx<< "# Mixed Clusters radial dist. Mass, logE/GeV " << "\t" << mpcr << "\t" << epcr << endl;
    
    ofclst <<
    "#PrCR"    << "\t" <<
    "EnePCR"   << "\t" <<
    "ZthPCR"   << "\t" <<
    "HghPCR"   << "\t" <<
    "IShow"    << "\t" <<
    "IClust"   << "\t" <<
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
    "TLstEl"   << "\t" <<
    "ZhFstE"   << "\t" <<
    "AzFstE"   << "\t" <<
    "ZhLstE"   << "\t" <<
    "AzlstE"   << "\t" <<
    "PmFstE"   << "\t" <<
    "PmstE"    << "\t" <<
    "TFrstM"   << "\t" <<
    "TLastM"   << "\t" <<
    "ZhfstM"   << "\t" <<
    "AzfstM"   << "\t" <<
    "ZhlstM"   << "\t" <<
    "AzlstM"   << "\t" <<
    "PmFstM"   << "\t" <<
    "PmLstM"   << "\t" <<
    "DistFL"   << "\t" <<
    "AnglFL"   << endl;

    ofshws <<
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
    "NCmss"   << "\t" <<
    "Nel50S"  << "\t" <<
    "Ne100S"  << "\t" <<
    "Ne150S"  << "\t" <<
    "Ne200S"  << "\t" <<
    "Ne300S"  << "\t" <<
    "Ne500S"  << "\t" <<
    "Nel1kS"  << "\t" <<
    "Nel2kS"  <<endl ;

    ofpart <<
    "#PrimCR"  << "\t" <<
    "EnePCR"   << "\t" <<
    "ZthPCR"   << "\t" <<
    "HghPCR"   << "\t" <<
    "Flint"    << "\t" <<
    "NTShow"   << "\t" <<
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
    
    ofpgrd <<
    "#Ishow"  << "\t" <<
    "nBox" << "\t" <<
    "...... nBox*nBox particle contents......" <<endl ;
    
    ofcgrd <<
    "#Ishow"  << "\t" <<
    "nBox" << "\t" <<
    "...... nBox*nBox cluster contents......" <<endl ;
    
    ofrpart <<
    "#Nshws"  << "\t" <<
    "Flint" << "\t" <<
    "LenGrd" << "\t" <<
    "nBox" << "\t" <<
    "...... Particle radial distribution......" <<endl;
    
    ofrclst <<
    "#Nshws"  << "\t" <<
    "Flint" << "\t" <<
    "LenGrd" << "\t" <<
    "nBox" << "\t" <<
    "...... Cluster radial distribution......" <<endl ;
    
    ofrclte <<
    "#Nshws"  << "\t" <<
    "Flint" << "\t" <<
    "LenGrd" << "\t" <<
    "nBox" << "\t" <<
    "...... EM cluster radial distribution......" <<endl ;
    
    ofrcltm <<
    "#Nshws"  << "\t" <<
    "Flint" << "\t" <<
    "LenGrd" << "\t" <<
    "nBox" << "\t" <<
    "...... Muon cluster radial distribution......" <<endl ;
    
    ofrcltx <<
    "#Nshws"  << "\t" <<
    "Flint" << "\t" <<
    "LenGrd" << "\t" <<
    "nBox" << "\t" <<
    "...... Mxd cluster radial distribution......" <<endl ;
    
// *******************************************************************************
// *******************************************************************************
    
    Int_t pgrid[nshow][nbox][nbox],cgrid[nshow][nbox][nbox];
    //Int_t cgride[nshow][nbox][nbox],cgridm[nshow][nbox][nbox],cgridx[nshow][nbox][nbox];
    Float_t rgrid[nbox][nbox], rsign,
    rpart[nbox], rclst[nbox], rclte[nbox], rcltm[nbox], rcltx[nbox];
    Float_t cellw, celdw, lminx=-len/2., lminy=-len/2., lmind= -len/sqrt(2.), xcell, ycell ;
    Int_t irow, icol, irows, icols, idia;

    cellw = len/nbox;
    celdw = cellw * sqrt(2);   // cell diagonal width
 
    for(Int_t i=0; i<nshow; i++ )
        for(Int_t j=0; j<nbox; j++ )
    for(Int_t k=0; k<nbox; k++  ){
                pgrid[i][j][k] = 0.;
                cgrid[i][j][k] = 0.;
                //cgride[i][j][k] = 0.;
                //cgridm[i][j][k] = 0.;
                //cgridx[i][j][k] = 0.;

                //cout << "a[" << j << "][" << k << "]: ";
                //std::cout<<"matriz="<< pgrid[i][j][k]<<endl;
                }
    
    for(Int_t i=0; i<nshow; i++ )
        for(Int_t j=0; j<nbox; j++ )
            for(Int_t k=0; k<nbox; k++  ){
                pgrid[i][j][k] = 0.;
                cgrid[i][j][k] = 0.;
                //cgride[i][j][k] = 0.;
                //cgridm[i][j][k] = 0.;
                //cgridx[i][j][k] = 0.;

                //cout << "a[" << j << "][" << k << "]: ";
                //std::cout<<"matriz="<< pgrid[i][j][k]<<endl;
                }
    
    for(Int_t i=0; i<nbox; i++ )
        for(Int_t j=0; j<nbox; j++  ){
            xcell = (i+0.5)*cellw - len/2 - 0.000001;   // avoiding zero value
            ycell = (j+0.5)*cellw - len/2 - 0.000001;
            rgrid[i][j] = sqrt(pow(xcell,2) + pow(ycell,2));
            rsign = ((ycell-nbox*xcell)/abs(ycell-nbox*xcell)) ;  // we add a sign depending on the x position
            rgrid[i][j] = rsign * rgrid[i][j];
                // cout << "--- rgrid:  " << rgrid[i][j] << endl;
            }
    
    for(Int_t i=0; i<nbox; i++  ){
        rpart[i]=0;  // radial distribution of all particles
        rclst[i]=0;  // radial distribution of all particles
        rclte[i]=0;  // radial distribution of all particles
        rcltm[i]=0;  // radial distribution of all particles
        rcltx[i]=0;  // radial distribution of all particles
    }
        
    for (Long64_t ishow=0; ishow<nshow; ishow++) {
       // cout << " *** new shower " << ishow+1 << endl;
       Long64_t itree = LoadTree(ishow);
       
       if (itree < 0) break;
        
        Long64_t nbytes = 0, nb = 0;
        Long64_t fbytes = 0, fb = 0;
        Long64_t ebytes = 0, eb = 0;
        Long64_t hbytes = 0, hb = 0;
        Long64_t tbytes = 0, tb = 0;
       
       nb = fChain->GetEntry(itree);   nbytes += nb;
       // if (Cut(ientry) < 0) continue;
       //eb  = b_shower_Energy->GetEntry(itree);        ebytes += eb;
       fb  = b_shower_FirstHeight->GetEntry(itree);   fbytes += fb;
       hb  = b_shower_Theta->GetEntry(itree);         hbytes += hb;
       tb  = b_shower_Phi->GetEntry(itree);           tbytes += tb;
       
       icont = -1; // index for total nb. of secondaries
    
       hghtfi = shower_FirstHeight/100000; // h en kilometros
       mnhgt   = mnhgt + (hghtfi-mnhgt)/(ishow+1);
       cout << "ene, hghtfi: " << shower_Energy << " " << hghtfi << endl;
       
       // cout << "LogEPCR, MeanEne " << epcr << " " << mlge << endl;

       Float_t arr0[100000] = { [0 ... 99999] = -10. };
       Float_t arr1[100000] = { [0 ... 99999] = -10. };
       Float_t arr2[100000] = { [0 ... 99999] = -10. };
       Float_t arr3[100000] = { [0 ... 99999] = -10. };
       Float_t arr4[100000] = { [0 ... 99999] = -10. };
       Float_t arr5[100000] = { [0 ... 99999] = -10. };
       Float_t arr6[100000] = { [0 ... 99999] = -10. };
       Float_t itag[100000] = { [0 ... 99999] =  -1. };
        //cout << " Punto 2" << endl;
// =========================================================================
       cmult = 1;
       icshow = 0;
       particle__ = 20;  //-
       nsecp = particle__;  // nb. secondary part. in shower
    
       // Look for the fastest particle in the shower
       t0 = 1000000;
       for(Int_t ip=0; ip < particle__; ip++){
           time      = particle__Time[ip];
           if(time < t0){
               t0 = time;
           }
       }
   
   // First loop in  all particles ===================
   for(Int_t ip=0; ip < particle__; ip++){
       
       pid = particle__ParticleID[ip];
       
       //cout << " Punto 3a" << endl;
       
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
       
       irow = int((rx - lminx)/cellw);
       icol = int((ry - lminy)/cellw);
              
       if (irow<0 | icol<0 | irow+1>nbox | icol+1>nbox){
           itag[ip]++;
           continue;
       }
       
       icont++;     // nb of particles
       nasep++;  // nb. of accepted secondaries
       pgrid[ishow][icol][irow]++;
       idia = int((rgrid[irow][icol] - lmind)/celdw);
       rpart[idia]++;
       
       //ofpgrd << setprecision(4)
       //<< rx << "\t"
       //<< ry <<endl;

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
           
           if (pmod < 0.10){nel50++; ntel50++;}
           else if (pmod>0.10 && pmod<0.15){nel100++; ntel100++;}
           else if (pmod>0.15 && pmod<0.20){nel150++; ntel150++;}
           else if (pmod>0.20 && pmod<0.30){nel200++; ntel200++;}
           else if (pmod>0.30 && pmod<0.50){nel300++; ntel300++;}
           else if (pmod>0.50 && pmod<1.00){nel500++; ntel500++;}
           else if (pmod>1.00 && pmod<2.00){nel1000++;ntel1000++;}
           else {nel2000++; ntel2000++;
           }
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
           // cout << " ***** neutron" << endl;
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
       }   // endif in pid
       
       // cout << "****** idp " << idp << endl;
       
       arr0[ip] = idp;
       arr1[ip] = rx;
       arr2[ip] = ry;
       arr3[ip] = time;
       arr4[ip] = pmod;
       arr5[ip] = theta;
       arr6[ip] = phi;
    
       
   }    // end of first scan of particles
        
   ntsep  = ntsep + nsecp;   // total nb. secondary part.
          
   // Second loop in the accepted particles ==========
   for(Int_t ifp=0; ifp<nasep; ifp++){
       
       if(itag[ifp]!=-1) continue;
       itag[ifp]++;

       idp1  = arr0[ifp];
       
       x1    = arr1[ifp];
       y1    = arr2[ifp];
       t1    = arr3[ifp];
       pm1   = arr4[ifp];
       zen1  = arr5[ifp];
       azh1  = arr6[ifp];
       
       irow = int((x1 - lminx)/cellw);
       icol = int((y1 - lminy)/cellw);
       //ofpgrd << setprecision(4)
       //<< rx << "\t"
       //<< ry <<endl;
       
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

       // Loop in the remaining particles
       for(Int_t isp=ifp+1; isp < nsecp; isp++){
           
           if(itag[isp]!=-1) continue;
           
           idp2 = arr0[isp];
           x2   = arr1[isp] ;
           y2   = arr2[isp] ;
           irows = int((x2 - lminx)/cellw);
           icols = int((y2 - lminy)/cellw);
       
           if( irows == irow && icols == icol){   // new particle in the cell
               
               if (cmult==1) {
                   itag[isp]++;
                   idclus = idp1;
                   iclust ++;
                   icshow ++;
                   //cgrid[ishow][irow][icol]= idclus;
                   cgrid[ishow][irow][icol]++;                   
                   idia = int((rgrid[irow][icol] - lmind)/celdw);
                   rclst[idia]++;
               }
               
               itag[isp]++;
               cmult++;
               idclus = idclus + idp2;
               //rclst[idia]++;
               
               //cout << "  ifp, isp, idp1, idp2: " << ifp << " " << isp << " " << idp1<< " " << idp2 << endl;
               
               t2 = arr3[isp];    pm2  = arr4[isp];
               zen2 = arr5[isp];  azh2 = arr6[isp];
               
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
                              
           }  // ends adding particle to the cluster
           
           //itag[isp] = tag;
           //tag = -1;
           
       }  //  End of isp for;
       
       if (cmult > 1){      //  New cluster found
           
           // cout << " ****** iclust, irow, icol, cmult: " << iclust << " " << irow  << " " << icol << " " << cmult << endl;
           //cout << "---irow, icol, idclus: " << irow << " " << icol << " " << idclus << endl;
           //cgrid[ishow][irow][icol] = idclus;
           tidclus = tidclus + idclus; // suma logica de todos los clusters en el shower
           
           // Clasificamos el cluster . Ignoramos las particulas pesadas
       
           if (neclst>1 && nmclst ==0 ){   // clean EM cluster with >1 electron
               sidclst = 1;
               //cgride[ishow][irow][icol]++;
               rclte[idia]++;
               iclems ++;
               iclemt ++;
               xf = xfe; yf = yfe; xl = xle; yl = yle;
               zhf = zhfste; zhl = zhlste; azf = azfste; azl = azlste;
           }
           else if (nmclst > 1 && neclst == 0){   //  cluster with >1 muons
               sidclst = 2;
               //cgridm[ishow][irow][icol]++;
               rcltm[idia]++;
               iclmus ++;
               iclmut ++;
               xf = xfm; yf = yfm; xl = xlm; yl = ylm;
               zhf = zhfstm; zhl = zhlstm; azf = azfstm; azl = azlstm;
           }
           else if (neclst>0 && nmclst>0) {          // mixt cluster with e's & mu's
               sidclst = 3;
               //cgridx[ishow][irow][icol]++;
               rcltx[idia]++;
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
           ofclst << mpcr << "\t"
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
              //  << fixed << setprecision(3)
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
       
   ofshws
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
       
     //Int_t inte=0, res[nbox][nbox];

     ofpgrd << ishow+1 << "\t" << nbox << "\t" ;
     ofcgrd << ishow+1 << "\t" << nbox << "\t" ;
     for (Int_t row=0; row<nbox; row++){
        for (Int_t col=0; col<nbox; col++){
            //cout<<"cgrid["<< ishow<<"]["<<row<<"]["<<col<<"]="<<cgrid[ishow][row][col] <<endl;
            ofpgrd << pgrid[ishow][row][col] << "\t";
            ofcgrd << cgrid[ishow][row][col] << "\t";
            //inte = inte + pgrid[ishow][row][col];
        }
     }
    ofpgrd << endl;
    ofcgrd << endl;
    
        /*
        for (Int_t col=0; col<nbox; col++){
           for (Int_t row=0; row<nbox; row++){
               cout<<"cgrid["<< ishow<<"]["<<
               row<<"]["<<col<<"]="<< cgrid[ishow][row][col] <<endl;
           }
        }
        
        cout << endl;
        
        for (Int_t col=0; col<nbox; col++){
           for (Int_t row=0; row<nbox; row++){
               cout<<"cgride["<< ishow<<"]["<<
               row<<"]["<<col<<"]="<< cgride[ishow][row][col] <<endl;
           }
        }
        
        cout << endl;
        
        for (Int_t col=0; col<nbox; col++){
           for (Int_t row=0; row<nbox; row++){
               cout<<"cgridm["<< ishow<<"]["<<
               row<<"]["<<col<<"]="<< cgridm[ishow][row][col] <<endl;
           }
        }
        
        cout << endl;
        
        for (Int_t col=0; col<nbox; col++){
           for (Int_t row=0; row<nbox; row++){
               cout<<"cgridx["<< ishow<<"]["<<
               row<<"]["<<col<<"]="<< cgridx[ishow][row][col] <<endl;
           }
        }
        
        cout << endl;
         */
        
        //ofpgrd << setprecision(4)
        //<< rx << "\t"
        //<< ry <<endl;
       
    } // end of for-loop in showers
    
    ofrpart << nshow << "\t" << flint << "\t" << len << "\t" << nbox << "\t";
    ofrclst << nshow << "\t" << flint << "\t" << len << "\t" << nbox << "\t";
    ofrclte << nshow << "\t" << flint << "\t" << len << "\t" << nbox << "\t";
    ofrcltm << nshow << "\t" << flint << "\t" << len << "\t" << nbox << "\t";
    ofrcltx << nshow << "\t" << flint << "\t" << len << "\t" << nbox << "\t";
    for (Int_t i=0; i<nbox; i++){
        ofrpart << rpart[i] << "\t";
        ofrclst << rclst[i] << "\t";
        ofrclte << rclte[i] << "\t";
        ofrcltm << rcltm[i] << "\t";
        ofrcltx << rcltx[i] << "\t";
        // cout << "*** rpart   " << rpart[i] << endl;
    }
    ofrpart << endl;
    ofrclst << endl;
    ofrclte << endl;
    ofrcltm << endl;
    ofrcltx << endl;

ofpart<< mpcr  << "\t"
    << epcr    << "\t"
    << zhpcr   << "\t"
    << mnhgt   << "\t"
    << flint   << "\t"
    << nshow   << "\t"
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

    ofclst.close();
    ofshws.close();
    ofpart.close();
    
    printf("******  Proceso completado  ******");
    cout << endl;
    
    cout << " - Parametros iniciales. MxDist - MxSigr: " << distmx << "  " << sigrmx << endl;
    cout << "* NShowers: " << nshow << endl;
    //cout << endl;
    cout << "* Numero total de secundarios: " << ntsep << endl;
    cout << "* Numero total aceptados: " << nasep << " / "<< 100 * nasep/ntsep << " %" <<endl;
    cout << "* Contador (gemnpo): " << ntgam << " " << ntele << " " << ntmu <<  " " << ntn <<  " " << ntp <<  " " << ntoth << endl;
    cout << "* Distribucion/% (gemnpo): " << 100*ntgam/ntsep << " " << 100*ntele/ntsep << " " << 100*ntmu/ntsep <<  " " << 100*ntn/ntsep <<  " " << 100*ntp/ntsep <<  " " << 100*ntoth/ntsep << endl;
    cout << "* NParticulas/shower (gemnpo): " << ntgam/nshow << " " << ntele/nshow << " " << ntmu/nshow <<  " " << ntn/nshow <<  " " << ntp/nshow <<  " " << ntoth/nshow << endl;
        cout << " --- " << endl;
        cout << "* Numero total de clusters: " << iclust << endl;
    cout << "* ShIdClust: 2Ele, 2MU, 2Ch, Oth: " << iclemt << " " << iclmut << " " << iclmxt << " " << iclott << endl;
    cout << "* ShIdClust/shower: 2Ele, 2MU, 2Ch, Oth: " << iclemt/nshow << " " << iclmut/nshow << " " << iclmxt/nshow << " " << iclott/nshow << endl;
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

