#define L1RateDisplacedPU_cxx
#include "L1RateDisplacedPU.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<TLorentzVector.h>
TH1D*hGenJetHT=new TH1D("hGenJetHT", "Gen H_{T}",40, 0, 800);
TH1D*hGenJetHT2LTkPass=new TH1D("hGenJetHT2LTkPass", "Gen H_{T}",40, 0, 800);
TH1D*hGenJetHTFastJetTkPass=new TH1D("hGenJetHTFastJetTkPass", "Gen H_{T}",40, 0, 800);
TH1D*hGenJetMHT=new TH1D("hGenJetMHT", "Gen Missing H_{T}",100, 0, 2000);
TH1D*hGenJetMHT2LTkPass=new TH1D("hGenJetMHT2LTkPass", "Gen H_{T}",100, 0, 2000);

TH1D*hTwoLJetHT=new TH1D("hTwoLJetHT", "TwoL Jet H_{T}",200, 0, 1000);
TH1D*hFastJetHT=new TH1D("hFastJetHT", "TwoL Jet H_{T}",200, 0, 1000);
TH1D*hTwoLJetHTDisp=new TH1D("hTwoLJetHTDisp", "TwoL Jet H_{T}",200, 0, 1000);
TH1D*hTwoLJetMHT=new TH1D("hTwoLJetMHT", "TwoL Jet H_{T}",200, 0, 500);

void L1RateDisplacedPU::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L L1RateDisplacedPU.C
//      Root > L1RateDisplacedPU t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       //if(MC_lep->at(0)>0)continue; 
       float GenHT=0;
              TLorentzVector tempGenMHT;
       for(unsigned int g=0; g<genjetak4_pt->size(); ++g){
           if(genjetak4_pt->at(g)<30)continue;
           if(fabs(genjetak4_eta->at(g))>2.4)continue;
           GenHT=GenHT+genjetak4_pt->at(g);
           TLorentzVector tempAdd ;
	   tempAdd.SetPtEtaPhiE(genjetak4_pt->at(g),genjetak4_eta->at(g),genjetak4_phi->at(g), genjetak4_p->at(g) );
	   tempGenMHT=tempGenMHT+tempAdd;
       }
       hGenJetHT->Fill(GenHT);
       hGenJetMHT->Fill(tempGenMHT.Pt());
       float TwoLJetHT=0;
       float TwoLJetHTDisp=0;
       int DispCounter=0;
       TLorentzVector tempMHT ;
       for(unsigned int j=0; j<Twoltrkjet_pt->size();++j){
          // if(fabs(pv_L1reco->at(0)-Twoltrkjet_vz->at(j))>3)continue;
           if(Twoltrkjet_ntracks->at(j)<2 && Twoltrkjet_pt->at(j)>50)continue;
           if(Twoltrkjet_ntracks->at(j)<3 && Twoltrkjet_pt->at(j)>100)continue;
           if(Twoltrkjet_pt->at(j)<5)continue;
           TLorentzVector tempAdd ;
	   tempAdd.SetPtEtaPhiE(Twoltrkjet_pt->at(j),Twoltrkjet_eta->at(j),Twoltrkjet_phi->at(j),Twoltrkjet_pt->at(j)*cosh(Twoltrkjet_eta->at(j)));
	   tempMHT=tempMHT+tempAdd;
           TwoLJetHT=TwoLJetHT+Twoltrkjet_pt->at(j);
           if(Twoltrkjet_nTightDisplaced->at(j)>2 )++DispCounter;//at least two displaced tracks
              }
	      //std::cout<<"MHT CrossCheck "<<
       	      hTwoLJetMHT->Fill(tempMHT.Pt());
       	      //hTwoLJetMHT->Fill(tkMHT);
              hTwoLJetHT->Fill(TwoLJetHT);
              if(TwoLJetHT>165)hGenJetHT2LTkPass->Fill(GenHT);
              if(tempMHT.Pt()>67.5)hGenJetMHT2LTkPass->Fill(tempGenMHT.Pt());
        for(unsigned int j=0; j<Twoltrkjet_pt->size();++j){
         //   if(fabs(pv_L1reco->at(0)-Twoltrkjet_vz->at(j))>3)continue;
            if(Twoltrkjet_ntracks->at(j)<2 && Twoltrkjet_pt->at(j)>50)continue;
            if(Twoltrkjet_ntracks->at(j)<3 && Twoltrkjet_pt->at(j)>100)continue;
            if(Twoltrkjet_pt->at(j)<5)continue;
               if(DispCounter<1)continue;//at least one displaced jet
            //if(DispCounter<2)
                TwoLJetHTDisp=TwoLJetHTDisp+Twoltrkjet_pt->at(j);

                     }
       hTwoLJetHTDisp->Fill(TwoLJetHTDisp);
              float FastJetHT=0;
              for(unsigned int j=0; j<trkjet_pt->size();++j){
                  if(trkjet_ntracks->at(j)<2 && trkjet_pt->at(j)>50)continue;
                  if(trkjet_ntracks->at(j)<3 && trkjet_pt->at(j)>100)continue;
                  if(trkjet_pt->at(j)<5)continue;
                  FastJetHT=FastJetHT+trkjet_pt->at(j);
              }
              hFastJetHT->Fill(FastJetHT);
              if(FastJetHT>165)hGenJetHTFastJetTkPass->Fill(GenHT);
              
              

   }
              TFile*fout=new TFile("TestRates6cmDz.root", "RECREATE");
	      hGenJetMHT2LTkPass->Write("hGenJetMHT2LTkPass");
	      hGenJetMHT->Write("GenMHT");
              hGenJetHT->Write("GenHT");
              hGenJetHT2LTkPass->Write("hGenJetHT2LTkPass");
              hGenJetHTFastJetTkPass->Write("hGenJetHTFastJetTkPass");
              hTwoLJetHT->Scale(1.0/hTwoLJetHT->Integral());
              hFastJetHT->Scale(1.0/hFastJetHT->Integral());
              hTwoLJetHTDisp->Scale(1.0/hTwoLJetHTDisp->Integral());
              hTwoLJetHT->Write("hTwoLJetHT");
              hTwoLJetHTDisp->Write("hTwoLJetHTDisp");
              hFastJetHT->Write("hFastJetHT");
    	      hTwoLJetMHT->Scale(1.0/hTwoLJetMHT->Integral());
	      hTwoLJetMHT->Write("TwoLMHT");
	    
}
