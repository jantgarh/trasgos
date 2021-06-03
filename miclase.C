#define miclase_cxx
#include "miclase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void miclase::Loop()
{
//   In a ROOT session, you can do:
//      root> .L miclase.C
//      root> miclase t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
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

   Long64_t nentries = fChain->GetEntriesFast();
   Float_t w = 1./nentries;
   cout << w << "   " << nentries << endl;

   Float_t sfh, sth, sph;  
   Int_t id;
   Float_t dx, dy, dr, dt;
   Float_t px, py, pz, pm, pt; 
   Float_t tx, ty, tr, th, ph;

   //Defining some histograms
   TH1F *h1d_mu_rad     = new TH1F("h1d_radial_mu","Mu radial distribution",100,1,1000000);
   TH1F *h1d_ele_rad    = new TH1F("h1d_radial_els","Ele radial distribution",100,1,1000000);
   TH1F *h1d_mu_tanr    = new TH1F("h1d_tanr_mu","Mu: tanr distribution",100,0,5);
   TH1F *h1d_ele_tanr   = new TH1F("h1d_tanr_ele","Ele: tanr distribution",100,0,5);
   TH2F *h2d_mu_rtanr   = new TH2F("h2d tanr_vs_r_mu ","Mu: tanr vs r",100,0,5,100,1,1000000);
   TH2F *h2d_ele_rtanr  = new TH2F("h2d tanr_vs_r_ele","Ele: tanr vs r",100,0,5,100,1,1000000);

   Long64_t nbytes = 0, nb = 0;
   //nentries = 10;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      sfh  = shower_FirstHeight/100;   // meters
      sth  = shower_Theta * (180/TMath::Pi()) ;   // meters
      sph  = shower_Phi   * (180/TMath::Pi()) ;   // meters

      //cout << particle__ << endl;

      for(Int_t ip=0; ip < particle__; ip++) {

	id  = particle__ParticleID[ip];
	dx  = particle__x[ip];
	dy  = particle__y[ip];
	dr  = sqrt(dx*dx + dy*dy);
        dt  = particle__Time[ip];
        px  = particle__Px[ip];
        py  = particle__Py[ip];
        pz  = particle__Pz[ip];
        pm  = sqrt(px*px + py*py + pz*pz);
        pt  = sqrt(px*px + py*py);  
	tx  = px/pz;
	ty  = py/pz; 
	tr  = pt/pz;                     
        th  = acos(pz/pm)*(180/TMath::Pi());
        ph  = atan(px/py)*(180/TMath::Pi());

	//cout << "*** " << id << endl;

        if( (id==2) || (id==3) ){ // electrons

	   //cout <<  "dr " << dr << " tr " << tr << endl;

	   h1d_ele_rad->Fill(dr,w);
      	   h1d_ele_tanr->Fill(tr,w); 
      	   h2d_ele_rtanr->Fill(tr,dr,w);
	
	   }

        if( (id==5) || (id==6) ){ // muons

	   //cout <<  "dr " << dr << " tr " << tr << endl;

	   h1d_mu_rad->Fill(dr,w);
      	   h1d_mu_tanr->Fill(tr,w); 
      	   h2d_mu_rtanr->Fill(tr,dr,w);

	   }
	       
      }

   }

   TCanvas *can1 = new TCanvas("can1","can1");
   can1->Divide(1,3);
   can1->cd(1);
   h1d_ele_rad->Draw();
   can1->cd(2);
   h1d_ele_tanr->Draw();
   can1->cd(3);
   h2d_ele_rtanr->Draw("colz");

   TCanvas *can2 = new TCanvas("can2","can2");
   can2->Divide(1,3);
   can2->cd(1);
   h1d_mu_rad->Draw();
   can2->cd(2);
   h1d_mu_tanr->Draw();
   can2->cd(3);
   h2d_mu_rtanr->Draw("colz");

}
