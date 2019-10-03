#include "../../PlottingTemplate/PlotTemplate.C"
void QuickEffPlot(){
gROOT->ProcessLine(".L ../../tdrstyle.C");
setTDRStyle();
    
    //TCanvas*c1=new TCanvas("c1","", 800,800);
    TCanvas* c1 = CreateCanvas("JetEff", false, false);

   /*
    gROOT->LoadMacro("CMS_lumi.C");
    writeExtraText = true;       // if extra text
    extraText  = "";  // default extra text is "Preliminary"
    lumi_13TeV  = "50 fb^{-1}"; // default is "19.7 fb^{-1}"
    lumi_sqrtS = "13 TeV";
    */

TFile*fin=new TFile("JetFindingEffTTBarNoPUNewVersion.root", "READ");
TFile*fin2=new TFile("JetFindingEffStopSamplePU.root", "READ");
FastJetEffTTBar=(TH1D*)fin->Get("GenEtaFastEff");
TwoLayerEffTTBar=(TH1D*)fin->Get("GenEtaEff");
    EtaDenomin=(TH1D*)fin->Get("GenEta");
    
    FastJetPtEffTTBar=(TH1D*)fin->Get("GenPtFastEff");
    TwoLayerPtEffTTBar=(TH1D*)fin->Get("GenPtEff");
    PtDenomin=(TH1D*)fin->Get("GenPt");
    
    
FastJetEffTTBar->SetLineWidth(2.0);
TwoLayerEffTTBar->SetLineWidth(2.0);

FastJetEffTTBar->SetLineColor(kBlue);
TwoLayerEffTTBar->SetLineColor(kCyan);

    c1->cd();
    TGraphAsymmErrors*grFastEta=new TGraphAsymmErrors(FastJetEffTTBar, EtaDenomin);
    TGraphAsymmErrors*grTwoLEta=new TGraphAsymmErrors(TwoLayerEffTTBar, EtaDenomin);
    TGraphAsymmErrors*grFastPt=new TGraphAsymmErrors(FastJetPtEffTTBar, PtDenomin);
    TGraphAsymmErrors*grTwoLPt=new TGraphAsymmErrors(TwoLayerPtEffTTBar, PtDenomin);

    
    grFastEta->GetXaxis()->SetRangeUser(-2.5,2.5);
    grFastEta->GetYaxis()->SetRangeUser(0.85,1.05);

    grFastEta->SetTitle("Efficiency vs. Gen Jet #eta; Gen Jet #eta [GeV]; Efficiency");
    
//FastJetEffTTBar->Draw();
    grFastEta->SetMarkerStyle(kOpenCircle);grFastEta->SetMarkerColor(kBlack);grFastEta->SetLineColor(kBlack);
    grFastEta->Draw("APE");
    grTwoLEta->SetMarkerStyle(kOpenCircle);grTwoLEta->SetMarkerColor(kRed);grTwoLEta->SetLineColor(kRed);

    grTwoLEta->Draw("PE");
    
    gPad->Update();
   // auto graph = FastJetEffTTBar->GetPaintedGraph();
   // graph->SetMinimum(0);
   // graph->SetMaximum(1);
    gPad->Update();
//TwoLayerEffTTBar->Draw("same");
    //return;
  TLegend* leg1=new TLegend(0.4423559,0.1909677,0.7919799,0.4206452);
 leg1->SetFillStyle(0);
 leg1->SetBorderSize(0);
 leg1->SetTextSize(0.04);
 leg1->SetTextFont(42);
 // leg1->SetHeader("Simulated Hadronic t#bar{t} Events PU=0");
    leg1->SetHeader("Simulated SUSY Stop m_{#tilde{t}}=1000 m_{LSP}=775");

  leg1->AddEntry(grFastEta, "FastJet anti-k_{T}, R=0.3", "PE");
  leg1->AddEntry(grTwoLEta, "L1 Track Jets ", "PE");
  leg1->Draw();
   // writeExtraText = true;       // if extra text
  //  extraText  = " Phase-2 Simulation ";
   // lumi_sqrtS = "PU 0 (14 TeV)";
//lumi_sqrtS = "PU 200 (14 TeV)";
    //lumi_sqrtS = "PU 0 (14 TeV)";
   // CMS_lumi( c1, 0, 1 );
    c1->Update();
    c1->Print("jet_eta_eff.pdf");
    //return;

    grFastPt->SetTitle("Efficiency vs. Gen Jet #eta; Gen Jet p_{T} [GeV]; Efficiency");

    grFastPt->SetMarkerStyle(kOpenCircle);grFastPt->SetMarkerColor(kBlack);grFastPt->SetLineColor(kBlack);
    grFastPt->GetXaxis()->SetRangeUser(30,300);

    grFastPt->GetYaxis()->SetRangeUser(0.85,1.05);

    grFastPt->Draw("APE");
    grTwoLPt->SetMarkerStyle(kOpenCircle);grTwoLPt->SetMarkerColor(kRed);grTwoLPt->SetLineColor(kRed);

    grTwoLPt->Draw("PE");
    leg1->Draw();
      DrawPrelimLabel(c1);
    DrawLumiLabel(c1,"");
    //CMS_lumi( c1, 0, 1 );
    c1->Update();
    return;

    TFile*fin3=new TFile("TTBarTurnOn.root", "READ");

TwoLayerTurnOnHT=(TEfficiency*)fin2->Get("RecoThresh140HT");
TwoLayerTurnOnQuad=(TEfficiency*)fin2->Get("RecoThreshQuadPt15");
    FastJetHT=(TH1D*)fin3->Get("HTTkJetPass");
    FastJetQuad=(TH1D*)fin3->Get("QuadTkJetPass");

    FastJetHTDen=(TH1D*)fin3->Get("GenHT");
    FastJetQuadDen=(TH1D*)fin3->Get("GenQuadJetpT");
    TGraphAsymmErrors*fjTurnOnHT= new TGraphAsymmErrors();
    TGraphAsymmErrors*fjTurnOnQuad= new TGraphAsymmErrors();

    fjTurnOnHT->BayesDivide(FastJetHT,FastJetHTDen);
    fjTurnOnQuad->BayesDivide(FastJetQuad,FastJetQuadDen);
    fjTurnOnHT->SetMarkerColor(kBlue);
    fjTurnOnHT->SetLineColor(kBlue); fjTurnOnHT->SetLineWidth(3.0);
    TwoLayerTurnOnHT->SetLineColor(kCyan);TwoLayerTurnOnHT->SetLineWidth(3.0);TwoLayerTurnOnHT->SetMarkerColor(kCyan);
    
    TCanvas*c1=new TCanvas("c1", "", 800, 800);
    fjTurnOnHT->SetTitle("HT Turn On TTBar (Hadronic) PU 200; Gen Level H_{T}; Efficiency;");
    fjTurnOnHT->Draw();
    TwoLayerTurnOnHT->Draw("same");
    leg1->Draw();
    
    
    fjTurnOnQuad->SetTitle("Quad Turn On TTBar (Hadronic) PU 200; Gen Level Quad Jet p_{T}; Efficiency;");
    fjTurnOnQuad->SetMarkerColor(kBlue);
    fjTurnOnQuad->SetLineColor(kBlue); fjTurnOnQuad->SetLineWidth(3.0);
    TwoLayerTurnOnQuad->SetLineColor(kCyan);TwoLayerTurnOnQuad->SetLineWidth(3.0);TwoLayerTurnOnQuad->SetMarkerColor(kCyan);
    fjTurnOnQuad->Draw();
    TwoLayerTurnOnQuad->Draw("same");
    leg1->Draw();

}

