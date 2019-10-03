#define JetEfficiency_cxx
#include "JetEfficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
void JetEfficiency::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L JetEfficiency.C
//      Root > JetEfficiency t
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
    TH1D*GenPt=new TH1D("GenPt", "Gen #phi pT", 60, 0, 300);
    TH1D*GenEta=new TH1D("GenEta", "Gen #phi pT", 70, -3.5, 3.5);

    TH1D*GenPtEff=new TH1D("GenPtEff", "Gen #phi pT", 60, 0, 300);
    TH1D*GenEtaEff=new TH1D("GenEtaEff", "Gen #phi pT", 70, -3.5, 3.5);
    
    TH1D*GenPtFastEff=new TH1D("GenPtFastEff", "Gen #phi pT", 60, 0, 300);
    TH1D*GenEtaFastEff=new TH1D("GenEtaFastEff", "Gen #phi pT", 70, -3.5, 3.5);
    
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

       nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
       //if(MC_lep->at(0)>0)continue;
       if(fabs(pv_MC->at(0)-pv_L1reco->at(0))>1)continue;
       float GenHT=0;
       for(unsigned int g=0; g<genjetak4_pt->size(); ++g){
           if(genjetak4_pt->at(g)<30)continue;
           if(fabs(genjetak4_eta->at(g))>2.4)continue;
           GenHT=GenHT+genjetak4_pt->at(g);
           GenPt->Fill(genjetak4_pt->at(g));
           GenEta->Fill(genjetak4_eta->at(g));
           bool JetMatch=false;
           for(unsigned int j=0; j<Twoltrkjet_pt->size(); ++j){
               float deta=genjetak4_eta->at(g)-Twoltrkjet_eta->at(j);
               float dphi=genjetak4_phi->at(g)-Twoltrkjet_phi->at(j);
               float dR=sqrt((deta*deta)+(dphi*dphi));
               if(dR<0.4){
                   JetMatch=true;
               }
           }
               if(JetMatch){
                   GenPtEff->Fill(genjetak4_pt->at(g));
                   GenEtaEff->Fill(genjetak4_eta->at(g));

               }
           JetMatch=false;
           for(unsigned int j=0; j<trkjet_pt->size(); ++j){
               float deta=genjetak4_eta->at(g)-trkjet_eta->at(j);
               float dphi=genjetak4_phi->at(g)-trkjet_phi->at(j);
               float dR=sqrt((deta*deta)+(dphi*dphi));
               if(dR<0.4){
                   JetMatch=true;
               }
           }
           if(JetMatch){
               GenPtFastEff->Fill(genjetak4_pt->at(g));
               GenEtaFastEff->Fill(genjetak4_eta->at(g));
               
           }
       }
   }
       TFile*fout=new TFile("JetFindingEffH250mPhi30PU.root","RECREATE");
       GenPt->Write("GenPt");
    GenEta->Write("GenEta");

       GenEtaEff->Write("GenEtaEff");
    GenPtEff->Write("GenPtEff");
    
    GenEtaFastEff->Write("GenEtaFastEff");
    GenPtFastEff->Write("GenPtFastEff");

}
