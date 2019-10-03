//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 26 15:21:13 2019 by ROOT version 5.34/37
// from TTree eventTree/Event tree
// found on file: MinBiasD41Extended.root
//////////////////////////////////////////////////////////

#ifndef L1RateDisplacedPU_h
#define L1RateDisplacedPU_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class L1RateDisplacedPU {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString basename="H250mPhi60ctau0";
   // Declaration of leaf types
   vector<float>   *trk_p;
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_d0;
   vector<float>   *trk_z0;
   vector<float>   *trk_bconsist;
   vector<float>   *trk_sconsist;
   vector<float>   *trk_chi2;
   vector<int>     *trk_nstub;
   vector<int>     *trk_psnstub;
   vector<int>     *trk_genuine;
   vector<int>     *trk_loose;
   vector<int>     *trk_unknown;
   vector<int>     *trk_combinatoric;
   vector<int>     *trk_fake;
   vector<int>     *trk_matchtp_pdgid;
   vector<float>   *trk_matchtp_pt;
   vector<float>   *trk_matchtp_eta;
   vector<float>   *trk_matchtp_phi;
   vector<float>   *trk_matchtp_z0;
   vector<float>   *trk_matchtp_dxy;
   vector<float>   *tp_p;
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_nstublayers;
   vector<int>     *tp_eventid;
   vector<int>     *tp_charge;
   vector<float>   *matchtrk_p;
   vector<float>   *matchtrk_pt;
   vector<float>   *matchtrk_eta;
   vector<float>   *matchtrk_phi;
   vector<float>   *matchtrk_z0;
   vector<float>   *matchtrk_d0;
   vector<float>   *matchtrk_chi2;
   vector<int>     *matchtrk_nstub;
   Float_t         trueMET;
   Float_t         tkMET;
   Float_t         tkMHT;
   Float_t         tkHT;
   vector<float>   *pv_L1recotruesumpt;
   vector<float>   *pv_L1recosumpt;
   vector<float>   *pv_L1reco;
   vector<float>   *pv_L1TP;
   vector<float>   *pv_L1TPsumpt;
   vector<int>     *MC_lep;
   vector<float>   *pv_MCChgSumpT;
   vector<float>   *pv_MC;
   vector<float>   *tpjet_eta;
   vector<float>   *tpjet_vz;
   vector<float>   *tpjet_p;
   vector<float>   *tpjet_pt;
   vector<float>   *tpjet_phi;
   vector<int>     *tpjet_ntracks;
   vector<float>   *tpjet_tp_sumpt;
   vector<float>   *tpjet_truetp_sumpt;
   vector<float>   *Twol_eta;
   vector<float>   *Twoltrkjet_vz;
   vector<float>   *Twoltrkjet_p;
    vector<float>   *Twoltrkjet_eta;
   vector<float>   *Twoltrkjet_pt;
   vector<float>   *Twoltrkjet_phi;
   vector<int>     *Twoltrkjet_ntracks;
   vector<int>     *Twoltrkjet_nDisplaced;
    
   vector<int>     *Twoltrkjet_nTight;
    vector<int>     *Twoltrkjet_nTightDisplaced;
   vector<float>   *trkjet_eta;
   vector<float>   *trkjet_vz;
   vector<float>   *trkjet_p;
   vector<float>   *trkjet_pt;
   vector<float>   *trkjet_phi;
   vector<int>     *trkjet_ntracks;
   vector<float>   *trkjet_truetp_sumpt;
   vector<float>   *genjetak4_neufrac;
   vector<float>   *genjetak4_chgfrac;
   vector<float>   *genjetak4_metfrac;
   vector<float>   *genjetak4_eta;
   vector<float>   *genjetak4_phi;
   vector<float>   *genjetak4_p;
   vector<float>   *genjetak4_pt;
   vector<float>   *genjetchgak4_eta;
   vector<float>   *genjetchgak4_phi;
   vector<float>   *genjetchgak4_p;
   vector<float>   *genjetchgak4_z;
   vector<float>   *genjetchgak4_pt;

   // List of branches
   TBranch        *b_trk_p;   //!
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_bconsist;   //!
   TBranch        *b_trk_sconsist;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_nstub;   //!
   TBranch        *b_trk_psnstub;   //!
   TBranch        *b_trk_genuine;   //!
   TBranch        *b_trk_loose;   //!
   TBranch        *b_trk_unknown;   //!
   TBranch        *b_trk_combinatoric;   //!
   TBranch        *b_trk_fake;   //!
   TBranch        *b_trk_matchtp_pdgid;   //!
   TBranch        *b_trk_matchtp_pt;   //!
   TBranch        *b_trk_matchtp_eta;   //!
   TBranch        *b_trk_matchtp_phi;   //!
   TBranch        *b_trk_matchtp_z0;   //!
   TBranch        *b_trk_matchtp_dxy;   //!
   TBranch        *b_tp_p;   //!
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_nstublayers;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_charge;   //!
   TBranch        *b_matchtrk_p;   //!
   TBranch        *b_matchtrk_pt;   //!
   TBranch        *b_matchtrk_eta;   //!
   TBranch        *b_matchtrk_phi;   //!
   TBranch        *b_matchtrk_z0;   //!
   TBranch        *b_matchtrk_d0;   //!
   TBranch        *b_matchtrk_chi2;   //!
   TBranch        *b_matchtrk_nstub;   //!
   TBranch        *b_trueMET;   //!
   TBranch        *b_tkMET;   //!
   TBranch        *b_tkMHT;   //!
   TBranch        *b_tkHT;   //!
   TBranch        *b_pv_L1recotruesumpt;   //!
   TBranch        *b_pv_L1recosumpt;   //!
   TBranch        *b_pv_L1reco;   //!
   TBranch        *b_pv_L1TP;   //!
   TBranch        *b_pv_L1TPsumpt;   //!
   TBranch        *b_MC_lep;   //!
   TBranch        *b_pv_MCChgSumpT;   //!
   TBranch        *b_pv_MC;   //!
   TBranch        *b_tpjet_eta;   //!
   TBranch        *b_tpjet_vz;   //!
   TBranch        *b_tpjet_p;   //!
   TBranch        *b_tpjet_pt;   //!
   TBranch        *b_tpjet_phi;   //!
   TBranch        *b_tpjet_ntracks;   //!
   TBranch        *b_tpjet_tp_sumpt;   //!
   TBranch        *b_tpjet_truetp_sumpt;   //!
   TBranch        *b_2ltrkjet_eta;   //!
   TBranch        *b_2ltrkjet_vz;   //!
   TBranch        *b_2ltrkjet_p;   //!
   TBranch        *b_2ltrkjet_pt;   //!
   TBranch        *b_2ltrkjet_phi;   //!
   TBranch        *b_2ltrkjet_ntracks;   //!
   TBranch        *b_2ltrkjet_nDisplaced;   //!
    TBranch        *b_2ltrkjet_nTight;   //!

   TBranch        *b_2ltrkjet_nTightDisplaced;   //!
   TBranch        *b_trkjet_eta;   //!
   TBranch        *b_trkjet_vz;   //!
   TBranch        *b_trkjet_p;   //!
   TBranch        *b_trkjet_pt;   //!
   TBranch        *b_trkjet_phi;   //!
   TBranch        *b_trkjet_ntracks;   //!
   TBranch        *b_trkjet_truetp_sumpt;   //!
   TBranch        *b_genjetak4_neufrac;   //!
   TBranch        *b_genjetak4_chgfrac;   //!
   TBranch        *b_genjetak4_metfrac;   //!
   TBranch        *b_genjetak4_eta;   //!
   TBranch        *b_genjetak4_phi;   //!
   TBranch        *b_genjetak4_p;   //!
   TBranch        *b_genjetak4_pt;   //!
   TBranch        *b_genjetchgak4_eta;   //!
   TBranch        *b_genjetchgak4_phi;   //!
   TBranch        *b_genjetchgak4_p;   //!
   TBranch        *b_genjetchgak4_z;   //!
   TBranch        *b_genjetchgak4_pt;   //!

   L1RateDisplacedPU(TString name,TTree *tree=0);
   virtual ~L1RateDisplacedPU();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef L1RateDisplacedPU_cxx
L1RateDisplacedPU::L1RateDisplacedPU(TString name,TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   basename=name;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(basename+"HybridExtended.root");
      if (!f || !f->IsOpen()) {
         f = new TFile(basename+"HybridExtended.root");
      }
      TDirectory * dir = (TDirectory*)f->Get(basename+"HybridExtended.root:/L1TrackJetsFast");
      dir->GetObject("eventTree",tree);

   }
   Init(tree);
}

L1RateDisplacedPU::~L1RateDisplacedPU()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t L1RateDisplacedPU::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t L1RateDisplacedPU::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void L1RateDisplacedPU::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_p = 0;
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_d0 = 0;
   trk_z0 = 0;
   trk_bconsist = 0;
   trk_sconsist = 0;
   trk_chi2 = 0;
   trk_nstub = 0;
   trk_psnstub = 0;
   trk_genuine = 0;
   trk_loose = 0;
   trk_unknown = 0;
   trk_combinatoric = 0;
   trk_fake = 0;
   trk_matchtp_pdgid = 0;
   trk_matchtp_pt = 0;
   trk_matchtp_eta = 0;
   trk_matchtp_phi = 0;
   trk_matchtp_z0 = 0;
   trk_matchtp_dxy = 0;
   tp_p = 0;
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_pdgid = 0;
   tp_nmatch = 0;
   tp_nstub = 0;
   tp_nstublayers = 0;
   tp_eventid = 0;
   tp_charge = 0;
   matchtrk_p = 0;
   matchtrk_pt = 0;
   matchtrk_eta = 0;
   matchtrk_phi = 0;
   matchtrk_z0 = 0;
   matchtrk_d0 = 0;
   matchtrk_chi2 = 0;
   matchtrk_nstub = 0;
   pv_L1recotruesumpt = 0;
   pv_L1recosumpt = 0;
   pv_L1reco = 0;
   pv_L1TP = 0;
   pv_L1TPsumpt = 0;
   MC_lep = 0;
   pv_MCChgSumpT = 0;
   pv_MC = 0;
   tpjet_eta = 0;
   tpjet_vz = 0;
   tpjet_p = 0;
   tpjet_pt = 0;
   tpjet_phi = 0;
   tpjet_ntracks = 0;
   tpjet_tp_sumpt = 0;
   tpjet_truetp_sumpt = 0;
   Twoltrkjet_eta = 0;
   Twoltrkjet_vz = 0;
   Twoltrkjet_p = 0;
   Twoltrkjet_pt = 0;
   Twoltrkjet_phi = 0;
   Twoltrkjet_ntracks = 0;
   Twoltrkjet_nDisplaced = 0;
   Twoltrkjet_nTight = 0;
    Twoltrkjet_nTightDisplaced = 0;

   trkjet_eta = 0;
   trkjet_vz = 0;
   trkjet_p = 0;
   trkjet_pt = 0;
   trkjet_phi = 0;
   trkjet_ntracks = 0;
   trkjet_truetp_sumpt = 0;
   genjetak4_neufrac = 0;
   genjetak4_chgfrac = 0;
   genjetak4_metfrac = 0;
   genjetak4_eta = 0;
   genjetak4_phi = 0;
   genjetak4_p = 0;
   genjetak4_pt = 0;
   genjetchgak4_eta = 0;
   genjetchgak4_phi = 0;
   genjetchgak4_p = 0;
   genjetchgak4_z = 0;
   genjetchgak4_pt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trk_p", &trk_p, &b_trk_p);
   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_bconsist", &trk_bconsist, &b_trk_bconsist);
   fChain->SetBranchAddress("trk_sconsist", &trk_sconsist, &b_trk_sconsist);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
   fChain->SetBranchAddress("trk_psnstub", &trk_psnstub, &b_trk_psnstub);
   fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
   fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
   fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
   fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
   fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
   fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
   fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
   fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
   fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
   fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
   fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
   fChain->SetBranchAddress("tp_p", &tp_p, &b_tp_p);
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_nstublayers", &tp_nstublayers, &b_tp_nstublayers);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);
   fChain->SetBranchAddress("matchtrk_p", &matchtrk_p, &b_matchtrk_p);
   fChain->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
   fChain->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
   fChain->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
   fChain->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
   fChain->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
   fChain->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
   fChain->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
   fChain->SetBranchAddress("trueMET", &trueMET, &b_trueMET);
   fChain->SetBranchAddress("tkMET", &tkMET, &b_tkMET);
   fChain->SetBranchAddress("tkMHT", &tkMHT, &b_tkMHT);
   fChain->SetBranchAddress("tkHT", &tkHT, &b_tkHT);
   fChain->SetBranchAddress("pv_L1recotruesumpt", &pv_L1recotruesumpt, &b_pv_L1recotruesumpt);
   fChain->SetBranchAddress("pv_L1recosumpt", &pv_L1recosumpt, &b_pv_L1recosumpt);
   fChain->SetBranchAddress("pv_L1reco", &pv_L1reco, &b_pv_L1reco);
   fChain->SetBranchAddress("pv_L1TP", &pv_L1TP, &b_pv_L1TP);
   fChain->SetBranchAddress("pv_L1TPsumpt", &pv_L1TPsumpt, &b_pv_L1TPsumpt);
   fChain->SetBranchAddress("MC_lep", &MC_lep, &b_MC_lep);
   fChain->SetBranchAddress("pv_MCChgSumpT", &pv_MCChgSumpT, &b_pv_MCChgSumpT);
   fChain->SetBranchAddress("pv_MC", &pv_MC, &b_pv_MC);
   fChain->SetBranchAddress("tpjet_eta", &tpjet_eta, &b_tpjet_eta);
   fChain->SetBranchAddress("tpjet_vz", &tpjet_vz, &b_tpjet_vz);
   fChain->SetBranchAddress("tpjet_p", &tpjet_p, &b_tpjet_p);
   fChain->SetBranchAddress("tpjet_pt", &tpjet_pt, &b_tpjet_pt);
   fChain->SetBranchAddress("tpjet_phi", &tpjet_phi, &b_tpjet_phi);
   fChain->SetBranchAddress("tpjet_ntracks", &tpjet_ntracks, &b_tpjet_ntracks);
   fChain->SetBranchAddress("tpjet_tp_sumpt", &tpjet_tp_sumpt, &b_tpjet_tp_sumpt);
   fChain->SetBranchAddress("tpjet_truetp_sumpt", &tpjet_truetp_sumpt, &b_tpjet_truetp_sumpt);
   fChain->SetBranchAddress("2ltrkjet_eta", &Twoltrkjet_eta, &b_2ltrkjet_eta);
   fChain->SetBranchAddress("2ltrkjet_vz", &Twoltrkjet_vz, &b_2ltrkjet_vz);
   fChain->SetBranchAddress("2ltrkjet_p", &Twoltrkjet_p, &b_2ltrkjet_p);
   fChain->SetBranchAddress("2ltrkjet_pt", &Twoltrkjet_pt, &b_2ltrkjet_pt);
   fChain->SetBranchAddress("2ltrkjet_phi", &Twoltrkjet_phi, &b_2ltrkjet_phi);
   fChain->SetBranchAddress("2ltrkjet_ntracks", &Twoltrkjet_ntracks, &b_2ltrkjet_ntracks);
   fChain->SetBranchAddress("2ltrkjet_nDisplaced", &Twoltrkjet_nDisplaced, &b_2ltrkjet_nDisplaced);
   fChain->SetBranchAddress("2ltrkjet_nTight", &Twoltrkjet_nTight, &b_2ltrkjet_nTight);
    fChain->SetBranchAddress("2ltrkjet_nTightDisplaced", &Twoltrkjet_nTightDisplaced, &b_2ltrkjet_nTightDisplaced);

   fChain->SetBranchAddress("trkjet_eta", &trkjet_eta, &b_trkjet_eta);
   fChain->SetBranchAddress("trkjet_vz", &trkjet_vz, &b_trkjet_vz);
   fChain->SetBranchAddress("trkjet_p", &trkjet_p, &b_trkjet_p);
   fChain->SetBranchAddress("trkjet_pt", &trkjet_pt, &b_trkjet_pt);
   fChain->SetBranchAddress("trkjet_phi", &trkjet_phi, &b_trkjet_phi);
   fChain->SetBranchAddress("trkjet_ntracks", &trkjet_ntracks, &b_trkjet_ntracks);
   fChain->SetBranchAddress("trkjet_truetp_sumpt", &trkjet_truetp_sumpt, &b_trkjet_truetp_sumpt);
   fChain->SetBranchAddress("genjetak4_neufrac", &genjetak4_neufrac, &b_genjetak4_neufrac);
   fChain->SetBranchAddress("genjetak4_chgfrac", &genjetak4_chgfrac, &b_genjetak4_chgfrac);
   fChain->SetBranchAddress("genjetak4_metfrac", &genjetak4_metfrac, &b_genjetak4_metfrac);
   fChain->SetBranchAddress("genjetak4_eta", &genjetak4_eta, &b_genjetak4_eta);
   fChain->SetBranchAddress("genjetak4_phi", &genjetak4_phi, &b_genjetak4_phi);
   fChain->SetBranchAddress("genjetak4_p", &genjetak4_p, &b_genjetak4_p);
   fChain->SetBranchAddress("genjetak4_pt", &genjetak4_pt, &b_genjetak4_pt);
   fChain->SetBranchAddress("genjetchgak4_eta", &genjetchgak4_eta, &b_genjetchgak4_eta);
   fChain->SetBranchAddress("genjetchgak4_phi", &genjetchgak4_phi, &b_genjetchgak4_phi);
   fChain->SetBranchAddress("genjetchgak4_p", &genjetchgak4_p, &b_genjetchgak4_p);
   fChain->SetBranchAddress("genjetchgak4_z", &genjetchgak4_z, &b_genjetchgak4_z);
   fChain->SetBranchAddress("genjetchgak4_pt", &genjetchgak4_pt, &b_genjetchgak4_pt);
   Notify();
}

Bool_t L1RateDisplacedPU::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void L1RateDisplacedPU::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t L1RateDisplacedPU::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef L1RateDisplacedPU_cxx
