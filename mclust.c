// mclust.C
// Nueva version basada en mclusters.C para estudio general de clusters a distinta altura y diferentes masas. Incorporadas algunas de las opciones de mclusth y mclustg
//
// !!!!! Valido para simulaciones sin fijar la altura
// JAG.Ene.21. Intento de redefinir los clusters con todas las particulas y separando componente electronica y muonica.
//Programa muy robusto. ANALISIS DE MULTIPLICIDAD. Identificacion de particulas a menos de una cierta distancia entre ellas. Con este programa queremos conocer los efectos sitematicos existente en la cascada.

#define mclust_cxx
#include "mclust.h"
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
#define fout1 "of_clustc.txt" // output files
#define fout2 "of_clusts.txt" //
#define fout3 "of_clustp.txt" //
#define fout4 "of_clustm.txt" //

void mclust::Loop()
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

   if (fChain == 0) return;
    
   Long64_t nshows = fChain->GetEntriesFast();
   cout << "- nShower in file:   " << nshows << endl;
   //nshows = 10;   //-
   cout << "- nShowers analyzed: " << nshows << endl;

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
    Int_t   pid, idp1, idp2, idpf, idpl;
    Int_t   ngam=0, nele=0, nmu=0, nn=0, np=0, noth=0, ngclst =0, neclst=0, nmclst=0;
    Int_t   nel50=0, nel100=0, nel150=0, nel200=0, nel300=0,
    nel500, nel1000, nel2000;
    Int_t   ntel50=0, ntel100=0, ntel150=0, ntel200=0, ntel300=0,
    ntel500, ntel1000, ntel2000;
    Int_t   ntgam=0, ntele=0, ntmu=0, ntn=0, ntp=0, ntoth=0;
    Int_t   ngams=0, neles=0, nmus=0, nns=0, nps=0, nots=0;
    Int_t   ifg, ife=0, ifm=0;
    Long64_t idgam=1, idele=1000, idmu=1000000, idn=100000000, ido=100000000;
    Long64_t idp=0, idclus=0, sidclst=0, tidclus=0;
    Float_t rx, ry, dx, dy, x1, x2, y1, y2, r, rsq,
            xf=0, xl=0, yf=0, yl=0, dcc=0, dtc=0, dac=0,
            xfe=0, xle=0, yfe=0, yle=0,
            xfm=0, xlm=0, yfm=0, ylm=0;
    Float_t azh=0., zen1=0., zen2=0., aza=0, azh1=0, azh2=0;
    Float_t time =0., t0, t1=0., t2=0.,
            tf, tl, tfg=0, tlg=0, tfe=0, tle=0, tfm=0, tlm=0,
            zhf=0, zhl=0, azf=0, azl=0,
            zhfe=0, zhle=0, zhfm=0, zhlm=0,
            azfe=0, azle=0, azfm=0, azlm=0;
    //Float_t tfc=0, tlc=0, zhfc=0, zhlc=0, azfc=0, azlc=0;
    Float_t px=0, py=0, pz=0, pt=0, px1=0, py1=0, pz1=0, px2=0, py2=0, pz2=0,
            pmod, pm1, pm2, psq,
            pmfe, pmle, pmfm, pmlm; // particle momenta
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
    Int_t ndist=8, ndelt=6, ndang=6, iner=0, idelt=0, idang=0;
    Int_t vdccem[ndist], vdccmu[ndist], vdccmx[ndist],
            vdtcem[ndelt], vdtcmu[ndelt], vdtcmx[ndelt],
            vdacem[ndang], vdacmu[ndang], vdacmx[ndang];
    Int_t mclem[ndist][ndelt][ndang], mclmu[ndist][ndelt][ndang], mclmx[ndist][ndelt][ndang];
    Float_t iee1=0.1, iee2=0.15, iee3=0.20, iee4=0.30, iee5=0.50, iee6=1.0, iee7=2.0;
    Int_t id1=50, id2=100, id3=150, id4=200, id5=300, id6=500, id7=1000;
    Int_t it1=1, it2= 3, it3 = 6, it4=20, it5=40;
    Int_t ia1=2, ia2= 5, ia3 = 9, ia4=15, ia5=25;
    Int_t jd, jt, ja; // indexes for dist, time and angle

    Float_t rad2g;
    rad2g = 180/TMath::Pi();
    Float_t epcr, dene=0.25, mepcr=0;
    
    fstream file1;
    fstream file2;
    fstream file3;
    fstream file4;
    
    file1.open(fout1, fstream::out); //- Clusters summary
    file2.open(fout2, fstream::out); //- Shower summary
    file3.open(fout3, fstream::out); //- Particle summary
    file4.open(fout4, fstream::out); //- Micro summary

    file1 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    file2 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    file3 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;
    file4 << "# DistMx, SigRMx: " << "\t" <<distmx  << "\t" << sigrmx << endl;

    file1 <<
    "#PrCR"    << "\t" <<
    "EnePCR"   << "\t" <<
    "HghPCR"   << "\t" <<
    "IShow"    << "\t" <<
    "IClust"   << "\t" <<
    "IdClst"   << "\t" <<
    "SIdClt"   << "\t" <<    // Short ID of cluster
    "NparCt"   << "\t" <<
    "NelClt"   << "\t" <<
    "NmuClt"   << "\t" <<
    "idpfPr"   << "\t" <<
    "idplPr"   << "\t" <<
    "XmClst"   << "\t" <<
    "YmClst"   << "\t" <<
    "RadWCt"   << "\t" <<
    "TfClst"   << "\t" <<
    "TlClst"   << "\t" <<
    "TmClst"   << "\t" <<
    "sTClst"   << "\t" <<
    "tfel"   << "\t" <<
    "tlel"   << "\t" <<
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
    "dcc"   << "\t" <<
    "dac/º"  << endl;

    file2 <<
    "#PrimCR" << "\t" <<
    "EnePCR"  << "\t" <<
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
    "MnHPCR"   << "\t" <<
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
    
    file4 <<
    "#PrimCR"  << "\t" <<
    "EnePCR"   << "\t" <<
    "NShow"    << "\t" <<
    "Ntele"    << "\t" <<
    "NtMu"     << "\t" <<
    "NClEl"    << "\t" <<
    "NClMu"    << "\t" <<
    "NClMx"    << "\t" <<
    "Nte100"   << "\t" <<
    "Nte200"   << "\t" <<
    "Nte500"   << "\t" <<
    "Ntel1K"   << "\t" <<
    "NCElR1"   << "\t" <<
    "NCElR2"   << "\t" <<
    "NCElR3"   << "\t" <<
    "NCElR4"   << "\t" <<
    "NCMuR1"   << "\t" <<
    "NCMuR2"   << "\t" <<
    "NCMuR3"   << "\t" <<
    "NCMuR4"   << "\t" <<
    "NCMxR1"   << "\t" <<
    "NCMxR2"   << "\t" <<
    "NCMxR3"   << "\t" <<
    "NCMxR4"   << "\t" <<
    "NCElT1"   << "\t" <<
    "NCElT2"   << "\t" <<
    "NCElT3"   << "\t" <<
    "NCMuT1"   << "\t" <<
    "NCMuT2"   << "\t" <<
    "NCMuT3"   << "\t" <<
    "NCMxT1"   << "\t" <<
    "NCMxT2"   << "\t" <<
    "NCMxT3"   << "\t" <<
    "NCElA1"   << "\t" <<
    "NCElA2"   << "\t" <<
    "NCElA3"   << "\t" <<
    "NCMuA1"   << "\t" <<
    "NCMuA2"   << "\t" <<
    "NCMuA3"   << "\t" <<
    "NCMxA1"   << "\t" <<
    "NCMxA2"   << "\t" <<
    "NCMxA3"   << endl ;
    
    for (Int_t ie=0; ie<ndist; ie++)
        for (Int_t it=0; it<ndelt; it++)
            for (Int_t ia=0; ia< ndang; ia++)
    {
        mclem[ie][it][ia]=0;
        mclmu[ie][it][ia]=0;
        mclmx[ie][it][ia]=0;
    }
    
    for (Int_t ie=0; ie<ndist; ie++){
        vdccem[ie]=0; vdccmu[ie]=0; vdccmx[ie]=0;}
    for (Int_t it=0; it<ndelt; it++){
        vdtcem[it]=0; vdtcmu[it]=0; vdtcmx[it]=0;}
    for (Int_t ia=0; ia< ndang; ia++){
        vdacem[ia]=0; vdacmu[ia]=0; vdacmx[ia]=0;}
// *******************************************************************************
    
   for (Long64_t ishow=0; ishow<nshows; ishow++) {

       // cout << " *** new shower " << endl;
       
       Long64_t itree = LoadTree(ishow);
       
       if (itree < 0) break;
       
       nb = fChain->GetEntry(itree);   nbytes += nb;
       
       // if (Cut(ientry) < 0) continue;
       
       //eb  = b_shower_Energy->GetEntry(itree);        ebytes += eb;
       
       //cout << " Ene " << shower_Energy << endl;
       
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
       azh     = acos(pz/pmod); // * rad2g;
       aza       = atan2(py,px); // * rad2g;
       
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
           
           if (pmod < iee1){nel50++; ntel50++;}
           else if (pmod>iee1 && pmod<iee2){nel100++; ntel100++;}
           else if (pmod>iee2 && pmod<iee3){nel150++; ntel150++;}
           else if (pmod>iee3 && pmod<iee4){nel200++; ntel200++;}
           else if (pmod>iee4 && pmod<iee5){nel300++; ntel300++;}
           else if (pmod>iee5 && pmod<iee6){nel500++; ntel500++;}
           else if (pmod>iee6 && pmod<iee7){nel1000++;ntel1000++;}
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
           //cout << "-- pid  " << pid << endl;
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
       arr5[icont] = azh;
       arr6[icont] = aza;
       
   //}
   }    // end of first scan of particles

   ntsep  = ntsep + nsecp;   // total nb. secondary part.
          
   // ===================    Second loop in the first particle
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
       tf = t1; tl = t1; tmean = t1; sigt=0; sigr = 0;
       idpf = idp1; idpl = idp1;
       
       if(idp1 == idgam){        // gamma
           ngclst++;
           ifg = 1;
           tfg = t1;
           tlg = t1;
       }
       
       if(idp1 == idele){        // electron
           neclst++;
           ife = 1;
           xfe   = x1;
           xle   = x1;
           yfe   = y1;
           yle   = y1;
           tfe   = t1;
           tle   = t1;
           pmfe  = pm1;
           pmle  = pm1;
           zhfe  = zen1;
           zhle  = zen1;
           azfe  = azh1;
           azle  = azh1;
       }
       
       if(idp1 == idmu){      // muon
           nmclst++;
           ifm = 1;
           xfm   = x1;
           xlm   = x1;
           yfm   = y1;
           ylm   = y1;
           tfm   = t1;
           tlm   = t1;
           pmfm  = pm1;
           pmlm  = pm1;
           zhfm  = zen1;
           zhlm  = zen1;
           azfm  = azh1;
           azlm  = azh1;
       }
                  
       itag[ifp]++;

       // Look for clusters in the remaining particles
       for(Int_t isp=ifp+1; isp < nsecp; isp++){
           
           //cout << "  ** look for clusters: itag  " << itag[ifp] << endl;
           
           if(itag[isp]!=-1) continue;

           dx   = arr1[isp] - xmean; dxsq = dx * dx ;
           dy   = arr2[isp] - ymean; dysq = dy * dy ;
           drsq = dxsq + dysq; dr = sqrt(drsq);
       
           // ==================================== fills cluster
           
           if(dr <= distmx && sigr <= sigrmx){
               
               if(cmult==1){idclus = idp1;}
               // cout << "   * new cluster: tag  " << tag << endl;
               tag++;   // cluster found
               cmult++;
               
               idp2   = arr0[isp];
               idclus = idclus + idp2;
               x2   = arr1[isp]; y2 = arr2[isp]; t2 = arr3[isp];
               pm2  = arr4[isp];
               zen2 = arr5[isp]; azh2 = arr6[isp];
               
               if(t2<tf){
                   tf = t2; xf = x2; yf = y2;
                   zhf = zen2; azf = azh2;
                   idpf = idp2;
               }
               
               if(t2>tl){
                   tl  = t2; xl = x2; yf = y2;
                   zhl = zen2; azl = azh2;
                   idpl = idp2;
               }
                   
               nxf = sin(zhf) * cos(azf) ;
               nyf = sin(zhf) * sin(azf) ;
               nzf = cos(zhf);
               nxl = sin(zhl) * cos(azl) ;
               nyl = sin(zhl) * sin(azl) ;
               nzl = cos(zhl);
               dac = acos(nxf*nxl + nyf*nyl + nzf*nzl) * rad2g;
               dcc = sqrt(((xl-xf)*(xl-xf) + (yl-yf)*(yl-yf))/2);
               dtc = tl - tf;

               if(idp2 == idgam){         // gamma
                   ngclst++;
                   if(ifg==0){
                       tfg = t2;
                       tlg = t2;
                       ifg = 1;
                   }
                   if(t2<tfg){ tfg = t2;}
                   if(t2>tlg){ tlg = t2;}
               }
               if(idp2 == idele){        // electron
                   neclst++;
                   if(ife==0){
                       xfe = x2; xle = x2; yfe = y2; yle = y2;
                       tfe = t2; tle = t2; pmfe = pm2; pmle = pm2;
                       zhfe= zen2; zhle = zen2; azfe = azh2; azle  = azh2;
                       ife = 1;
                   }
                   if(t2<tfe){
                       xfe = x2; yfe = y2; tfe = t2;
                       pmfe = pm2; zhfe = zen2; azfe = azh2;
                   }
                   if(t2>tle){
                       xle    = x2; yle = y2; tle = t2;
                       pmle = pm2; zhle = zen2; azle = azh2;
                   }
               }
               if(idp2 == idmu){      // muon
                   nmclst++;
                   if(ifm == 0){
                       xfm = x2; xlm = x2; yfm = y2; ylm = y2;
                       tfm = t2; tlm = t2; pmfm = pm2; pmlm = pm2;
                       zhfm = zen2; zhlm = zen2; azfm = azh2; azlm = azh2;
                       ifm = 1;
                   }
                   if(t2<tfm){
                       tfm = t2; xfm = x2; yfm = y2;
                       pmfm = pm2, zhfm = zen2; azfm = azh2;
                   }
                   
                   if(t2>tlm){
                       tlm = t2; xlm = x2; ylm = y2;
                       pmlm = pm2, zhlm = zen2; azlm = azh2;
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
               dcc   = sqrt((xmean*xmean + ymean*ymean)/2);
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
       
           if (dcc < id1){jd=0;}
           else if (dcc>id1 && dcc<id2){jd=1;}
           else if (dcc>id2 && dcc<id3){jd=2;}
           else if (dcc>id3 && dcc<id4){jd=3;}
           else if (dcc>id4 && dcc<id5){jd=4;}
           else if (dcc>id5 && dcc<id6){jd=5;}
           else if (dcc>id6 && dcc<id7){jd=6;}
           else {jd=7;}
           
           if (dtc < it1){jt=0;}
           else if (dtc> it1 && dtc<it2){jt=1;}
           else if (dtc> it2 && dtc<it3){jt=2;}
           else if (dtc> it3 && dtc<it4){jt=3;}
           else if (dtc> it4 && dtc<it5){jt=4;}
           else {jt=5;}
           
           if (dac < ia1){ja=0;}
           else if (dac> ia1 && dac<ia2){ja=1;}
           else if (dac> ia2 && dac<ia3){ja=2;}
           else if (dac> ia3 && dac<ia4){ja=3;}
           else if (dac> ia4 && dac<ia5){ja=4;}
           else {ja=5;}
           
       } //  End of isp for;
           
       //*
       if (cmult > 1){      //  ======== Stores new cluster
           iclust++;
           icshow ++;
           tidclus = tidclus + idclus;  // suma logica de los clusters en el shower
           
           // Clasificamos el cluster . Ignoramos las particulas pesadas
           
           //cout << "dcc, dtc, dac " << dcc << " " << dtc << " " << dac << endl;
       

           if (idclus < idmu && neclst > 1){   // clean EM cluster with >1 electron and no heavy particles
               sidclst = 1;
               iclems ++;
               iclemt ++;
               xf = xfe; yf = yfe; xl = xle; yl = yle;
               tf = tfe; tl = tle;
               zhf = zhfe; zhl = zhle; azf = azfe; azl = azle;
               mclem[jd][jt][ja]++;
               vdccem[jd] ++;
               vdtcem[jt] ++;
               vdacem[ja] ++;
           }
           else if (nmclst > 1 && neclst == 0){   //  cluster with >1 muons and perhaps a few gammas
               sidclst = 2;
               iclmus ++;
               iclmut ++;
               xf = xfm; yf = yfm; xl = xlm; yl = ylm;
               tf = tfm; tl = tlm;
               zhf = zhfm; zhl = zhlm; azf = azfm; azl = azlm;
               mclmu[jd][jt][ja]++;
               vdccmu[jd] ++;
               vdtcmu[jt] ++;
               vdacmu[ja] ++;
           }
           else if (neclst>0 && nmclst>0) {          // mixt cluster with e's & mu's
               sidclst = 3;
               iclmxs ++;
               iclmxt ++;
               if(tfe<tfm){
                  xf  = xfe; yf = yfe;
                  zhf = zhfe ; azf = azfe ;
               }
               else{
                  xf  = xfm; yf = yfm;
                  zhf = zhfm ; azf = azfm ;
               }
               if(tle>tlm){
                  xl  = xle; yl = yle;
                  zhl = zhle ; azl = azle ;
               }
               else{
                  xl  = xlm; yl  = ylm;
                  zhl = zhlm ; azl = azlm ;
               }
               
               mclmx[jd][jt][ja]++;
               vdccmx[jd] ++;
               vdtcmx[jt] ++;
               vdacmx[ja] ++;
           }
           else{iclots ++; iclott ++;}
       
           /*
           //   ---  New cluster
           cout << "*** ICluster " << icshow << " - itclust " << iclust << endl;
           cout << "* ClMult "  << cmult <<  " - ClustID " << idclus << endl;
           cout << "--- ShClustID (0-4):   " << sidclst << endl;
           cout << "-- iclems, iclmus, iclmxs, iclots: " << iclems << " " << iclmus << " " << iclmxs << " " << iclots << endl;
           //cout << "---dcc, dac:   " << dcc << " " << dac << endl;
           //cout<< "--- tf tl: " << tf << "   "<< tl << endl;
           //cout << "(xmean, ymean), sigr " << xmean << " " << ymean << " " << sigr << endl;
           cout << endl;
           */
           
           
           file1 << mpcr << "\t"
              << mepcr   << "\t"
              << hghtfi << "\t"
              << ishow+1<< "\t"
              << iclust << "\t"
              << idclus << "\t"
              << sidclst<< "\t"
              << cmult  << "\t"
              << neclst << "\t"
              << nmclst << "\t"
              << idpf << "\t"
              << idpl << "\t"
//             << fixed << setprecision(3)
              << setprecision(4)
              << xmean  << "\t"
              << ymean  << "\t"
              << sigr   << "\t"
              << tf   << "\t"
              << tl   << "\t"
              << tmean  << "\t"
              << sigt   << "\t"
              << tfe  << "\t"
              << tle  << "\t"
              << zhfe << "\t"
              << azfe << "\t"
              << zhle << "\t"
              << azle << "\t"
              << pmfe << "\t"
              << pmle << "\t"
              << tfm  << "\t"
              << tlm  << "\t"
              << zhfm << "\t"
              << azfm << "\t"
              << zhlm << "\t"
              << azlm << "\t"
              << pmfm << "\t"
              << pmlm << "\t"
              << dcc  << "\t"
              << dac  << endl;
           
       }  // endif cmult > 1:  new cluster
       
       //*/
       
       tag     = -1; idclus  = 0.; sidclst=0; tmean=0, cmult = 1.;
       sigxsq  = 0.; sigysq  = 0.; sigtsq  = 0;
       ifg     = 0; ife     = 0; ifm     = 0;
       ngclst  = 0; neclst  = 0; nmclst  = 0;
       tfg   = 0; tlg   = 0; tfe   = 0; tle   = 0;
       zhfe  = 0; zhle  = 0; azfe  = 0; azle  = 0;
       pmfe  = 0; pmle  = 0;
       tfm   = 0; tlm   = 0;
       zhfm  = 0; zhlm  = 0; azfm  = 0; azlm  = 0;
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
   //} // ***
       
    emin = mepcr - dene/2;
    emax = mepcr + dene/2;
    gm1  = spindex - 1;
    flint = (f0/gm1) * (pow(emax,-gm1) - pow(emin,-gm1));

file3<< mpcr  << "\t"
    << epcr    << "\t"
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
    
file4 << mpcr  << "\t"
    << epcr    << "\t"
    << nshows  << "\t"
    << ntele   << "\t"
    << ntmu    << "\t"
    << iclemt  << "\t"
    << iclmut  << "\t"
    << iclmxt  << "\t"
    << ntel50 + ntel100  << "\t"
    << ntel150 + ntel200  << "\t"
    << ntel300 + ntel500  << "\t"
    << ntel1000 + ntel2000  << "\t"
    << vdccem[0] << "\t"
    << vdccem[1] + vdccem[2] << "\t"
    << vdccem[3] + vdccem[4] << "\t"
    << vdccem[5] + vdccem[6] + vdccem[7] << "\t"
    << vdccmu[0] << "\t"
    << vdccmu[1] + vdccmu[2] << "\t"
    << vdccmu[3] + vdccmu[4] << "\t"
    << vdccmu[5] + vdccmu[6] + vdccmu[7] << "\t"
    << vdccmu[0] << "\t"
    << vdccmx[1] + vdccmx[2] << "\t"
    << vdccmx[3] + vdccmx[4] << "\t"
    << vdccmx[5] + vdccmx[6] + vdccmx[7] << "\t"
    << vdtcem[0] << "\t"
    << vdtcem[1] + vdtcem[2] << "\t"
    << vdtcem[3] + vdtcem[4] + vdtcem[5] << "\t"
    << vdtcmu[0] << "\t"
    << vdtcmu[1] + vdtcmu[2] << "\t"
    << vdtcmu[3] + vdtcmu[4] + vdtcmu[5] << "\t"
    << vdtcmx[0] << "\t"
    << vdtcmx[1] + vdtcmx[2] << "\t"
    << vdtcmx[3] + vdtcmx[4] + vdtcmx[5] << "\t"
    << vdacem[0] << "\t"
    << vdacem[1] + vdacem[2] << "\t"
    << vdacem[3] + vdacem[4] + vdacem[5] << "\t"
    << vdacmu[0] << "\t"
    << vdacmu[1] + vdacmu[2] << "\t"
    << vdacmu[3] + vdacmu[4] + vdacmu[5] << "\t"
    << vdacmx[0] << "\t"
    << vdacmx[1] + vdacmx[2] << "\t"
    << vdacmx[3] + vdacmx[4] + vdacmx[5] << endl;
    
cout << endl;
    
cout << "- Radial slices/m ["<<id1<<","<<id2<<","<<id3<<","<<id4<<","<<id5<<","<<id6<<","<<id7<<"]:"  << endl;
cout << "mclem rdistribution: " <<
    vdccem[0] << " " <<
    vdccem[1] << " " <<
    vdccem[2] << " " <<
    vdccem[3] << " " <<
    vdccem[4] << " " <<
    vdccem[5] << " " <<
    vdccem[6] << " " <<
    vdccem[7] << endl;
    cout << "mclmu rdistribution: " <<
    vdccmu[0] << " " <<
    vdccmu[1] << " " <<
    vdccmu[2] << " " <<
    vdccmu[3] << " " <<
    vdccmu[4] << " " <<
    vdccmu[5] << " " <<
    vdccmu[6] << " " <<
    vdccmu[7] << endl;
    cout << "mclmx rdistribution: " <<
    vdccmx[0] << " " <<
    vdccmx[1] << " " <<
    vdccmx[2] << " " <<
    vdccmx[3] << " " <<
    vdccmx[4] << " " <<
    vdccmx[5] << " " <<
    vdccmx[6] << " " <<
    vdccmx[7] << endl;
    cout << endl;
    
    cout << "- dTime intervals/ns ["<<it1<<","<<it2<<","<<it3<<","<<it4<<","<<it5 <<"]:" << endl;
   cout << "dtem distribution: " <<
   vdtcem[0] << " " <<
   vdtcem[1] << " " <<
   vdtcem[2] << " " <<
   vdtcem[3] << " " <<
   vdtcem[4] << " " <<
   vdtcem[5] << endl;
   cout << "dtmu distribution: " <<
   vdtcmu[0] << " " <<
   vdtcmu[1] << " " <<
   vdtcmu[2] << " " <<
   vdtcmu[3] << " " <<
   vdtcmu[4] << " " <<
   vdtcmu[5] << endl;
   cout << "dtmx distribution: " <<
   vdtcmx[0] << " " <<
   vdtcmx[1] << " " <<
   vdtcmx[2] << " " <<
   vdtcmx[3] << " " <<
   vdtcmx[4] << " " <<
   vdtcmx[5] << endl;
    cout << endl;
    
    cout << "- dAngle intervals/º ["<<ia1<<","<<ia2<<","<<ia3<<","<<ia4<<","<<ia5 <<"]:" << endl;
   cout << "daem distribution: " <<
   vdacem[0] << " " <<
   vdacem[1] << " " <<
   vdacem[2] << " " <<
   vdacem[3] << " " <<
   vdacem[4] << " " <<
   vdacem[5] << endl;
   cout << "damu distribution: " <<
   vdacmu[0] << " " <<
   vdacmu[1] << " " <<
   vdacmu[2] << " " <<
   vdacmu[3] << " " <<
   vdacmu[4] << " " <<
   vdacmu[5] << endl;
   cout << "damx distribution: " <<
   vdacmx[0] << " " <<
   vdacmx[1] << " " <<
   vdacmx[2] << " " <<
   vdacmx[3] << " " <<
   vdacmx[4] << " " <<
   vdacmx[5] << endl;
       
cout << endl;
printf("******  Proceso completado  ******");
file1.close();
file2.close();
file3.close();
file4.close();
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

