#include "PlotTemplate.C"
#include<iostream>
void QuickRateEff( int mH=125, int mPhi=15){
    //std::cout<<signame<<std::endl;
   gROOT->ProcessLine(".L ../tdrStyle.C");
   setTDRStyle();
    TFile*fin=new TFile("HybridPromptTestFinalCuts.root", "READ");
    TFile*finExt=new TFile("HybridExtendedFinalCuts.root", "READ");


   TwoLMinBiasHT=(TH1F*)fin->Get("hTwoLJetHT");
   L1RateTwoLMinBiasHT=(TH1F*)TwoLMinBiasHT->Clone("L1RateTwoLMinBiasHT");
    
    
    TwoLMinBiasExtHT=(TH1F*)finExt->Get("hTwoLJetHTDisp");
    L1RateTwoLMinBiasExtHT=(TH1F*)TwoLMinBiasExtHT->Clone("L1RateTwoLMinBiasExtHT");
    
    for(unsigned int i=1; i<=TwoLMinBiasHT->GetNbinsX(); ++i){
        //need to check for zero bins
     
        L1RateTwoLMinBiasHT->SetBinContent(i,2760*11.246*TwoLMinBiasHT->Integral(i,200));
        L1RateTwoLMinBiasExtHT->SetBinContent(i,2760*11.246*TwoLMinBiasExtHT->Integral(i,200));


    }
    TCanvas*c1=new TCanvas("c1", "", 800,800);
    leg=new TLegend(0.5263158,0.7260813,0.9423559,0.9174312,NULL,"brNDC");
    writeExtraText = true;       // if extra text
    extraText  = " Phase-2 Simulation";
    lumi_sqrtS = "PU 200 (14 TeV) ";

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->AddEntry(L1RateTwoLMinBiasHT, "L1 Track Jets ", "L");

    c1->SetLogy();
 


    L1RateTwoLMinBiasHT->GetXaxis()->SetRangeUser(0,500);
    L1RateTwoLMinBiasHT->SetTitle("MinBias Events PU 200; Trk H_{T} [GeV]; L1 Rate (kHz)");

    L1RateTwoLMinBiasHT->GetYaxis()->SetTitleOffset(1.0);
    L1RateTwoLMinBiasHT->SetLineColor(kBlack);L1RateTwoLMinBiasHT->SetLineWidth(3.0);
   
    L1RateTwoLMinBiasHT->Draw("hist");


   // return;
    std::cout<<"HT Threshold at 25kHz Fast Jets "<<L1RateTwoLMinBiasHT->GetBinLowEdge(L1RateTwoLMinBiasHT->FindLastBinAbove(25))<<std::endl;
    c1->SetGridx();
    c1->SetGridy();
    L1RateTwoLMinBiasHT->GetYaxis()->CenterTitle(true);
    L1RateTwoLMinBiasHT->GetYaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->GetXaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->SetLineColor(kBlack);

    leg->Draw("hist");
  //  CMS_lumi( c1, 0, 1 );
    c1->Update();
    c1->Print("rates_HT.pdf");
    
    
   // return;

    leg->Draw();
    
    // return;
    c1->SetGridx();
    c1->SetGridy();
    L1RateTwoLMinBiasHT->GetYaxis()->CenterTitle(true);
    L1RateTwoLMinBiasHT->GetYaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->GetXaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->SetLineColor(kBlack);
    
    leg->Draw("hist");
  //  CMS_lumi( c1, 0, 1 );
    c1->Update();
    c1->Print("rates_HT.pdf");
    TString fname=TString::Format("H%dmPhi%d", mH, mPhi);
    TFile*fsig=new TFile(fname+"ctau0TurnOnRate.root", "READ");
    //TFile*fsig=new TFile(TString::Format("PUJetsOutputHiddenHiggs%sctau0TDRTree.root",signame), "READ");
    TFile*fsig2=new TFile(fname+"ctau10TurnOnRate.root", "READ");
    TFile*fsig3=new TFile(fname+"ctau50TurnOnRate.root", "READ");
    


    HTSigEff=(TH1D*)fsig->Get("hTwoLJetHT");
    HTSigEffCtau10=(TH1D*)fsig2->Get("hTwoLJetHT");
    HTSigEffCtau50=(TH1D*)fsig3->Get("hTwoLJetHT");

    HTSigEff->Scale(1.0/HTSigEff->Integral());
    HTSigEffCtau10->Scale(1.0/HTSigEffCtau10->Integral());
    HTSigEffCtau50->Scale(1.0/HTSigEffCtau50->Integral());

    TGraph*gr=new TGraph();
    TGraph*gr2=new TGraph();
    TGraph*gr3=new TGraph();

    for(unsigned int i=1; i<L1RateTwoLMinBiasHT->GetNbinsX(); ++i){
        if(L1RateTwoLMinBiasHT->GetBinContent(i)>50)continue;
       // std::cout<<"Signal Eff "<<HTSigEff->Integral(i,200)<<std::endl;
                gr->SetPoint(gr->GetN(), HTSigEff->Integral(i,200), L1RateTwoLMinBiasHT->GetBinContent(i) );
                        gr2->SetPoint(gr2->GetN(), HTSigEffCtau10->Integral(i,200), L1RateTwoLMinBiasHT->GetBinContent(i) );
        gr3->SetPoint(gr3->GetN(), HTSigEffCtau50->Integral(i,200), L1RateTwoLMinBiasHT->GetBinContent(i) );

    }
   // TCanvas*c2=new TCanvas("c2", "", 800,800);
   // c2->SetLogy();
    TCanvas* c2 = CreateCanvas("RatevEff", true, false);
    writeExtraText = true;       // if extra text
    extraText  = " Phase-2 Simulation";
    lumi_sqrtS = "PU 200 (14 TeV) ";
    c1->SetLogy();
    gr->SetLineStyle(kDashed);gr->SetLineColor(kBlue);gr->SetLineWidth(4.0);
    gr2->SetLineStyle(kDashed);gr2->SetLineColor(kRed);gr2->SetLineWidth(4.0);
    gr3->SetLineStyle(kDashed);gr3->SetLineColor(kGreen);gr3->SetLineWidth(4.0);

    gr->GetXaxis()->SetRangeUser(0.0, 1.0);
    gr->GetYaxis()->SetRangeUser(1.0, 50);

 
    //gr3->Draw("LSame");
    //HT plot

    
    TString fname=TString::Format("H%dmPhi%d", mH, mPhi);

    TFile*fExtsig=new TFile(fname+"ctau0ExtendedTurnOnRate.root", "READ");
    TFile*fExtsig2=new TFile(fname+"ctau10ExtendedTurnOnRate.root", "READ");
    TFile*fExtsig3=new TFile(fname+"ctau50ExtendedTurnOnRate.root", "READ");
    HTExtSigEff=(TH1D*)fExtsig->Get("hTwoLJetHTDisp");
    HTExtSigEffCtau10=(TH1D*)fExtsig2->Get("hTwoLJetHTDisp");
    HTExtSigEffCtau50=(TH1D*)fExtsig3->Get("hTwoLJetHTDisp");
    TGraph*grExt=new TGraph();
    TGraph*gr2Ext=new TGraph();
    TGraph*gr3Ext=new TGraph();
    for(unsigned int i=1; i<L1RateTwoLMinBiasExtHT->GetNbinsX(); ++i){
        if(L1RateTwoLMinBiasExtHT->GetBinContent(i)>50)continue;
        // std::cout<<"Signal Eff "<<HTSigEff->Integral(i,200)<<std::endl;
        grExt->SetPoint(grExt->GetN(), HTExtSigEff->Integral(i,200), L1RateTwoLMinBiasExtHT->GetBinContent(i) );
        gr2Ext->SetPoint(gr2Ext->GetN(), HTExtSigEffCtau10->Integral(i,200), L1RateTwoLMinBiasExtHT->GetBinContent(i) );
        gr3Ext->SetPoint(gr3Ext->GetN(), HTExtSigEffCtau50->Integral(i,200), L1RateTwoLMinBiasExtHT->GetBinContent(i) );
        
    }
    
    

   grExt->SetLineColor(kBlue);grExt->SetLineWidth(4.0);
    gr2Ext->SetLineColor(kRed);gr2Ext->SetLineWidth(4.0);
    gr3Ext->SetLineColor(kGreen);gr3Ext->SetLineWidth(4.0);
    gr2Ext->SetTitle("; Signal Efficiency; L1 Track H_{T} Rate ( kHz )");

    gr2Ext->GetXaxis()->SetLimits(0,0.2);
    gr2Ext->GetYaxis()->SetRangeUser(1.0,40);
    gr2Ext->GetXaxis()->SetNdivisions(509);
    gr2Ext->Draw("AL");
  
    c2->Update();
    grExt->Draw("LSame");
    gr3Ext->Draw("LSame");
    gr->Draw("LSame");
    gr2->Draw("LSame");
    gr3->Draw("LSame");
    
    TLegend* leg1=new TLegend(0.6328321,0.323871,0.9974937,0.4941935);
    //  TLegend* leg1=new TLegend(0.6240602,0.3264516,0.9736842,0.4967742);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->SetTextFont(42);
    leg1->SetHeader("Prompt Tracks");
    leg1->AddEntry(gr,"c#tau = 0 cm" ,"L");
    leg1->AddEntry(gr2,"c#tau = 1 cm" ,"L");
    leg1->AddEntry(gr3,"c#tau = 5 cm" ,"L");
    leg1->Draw();
    TLegend* leg2=new TLegend(0.6328321,0.1677419,0.9824561,0.3251613);
    //TLegend* leg2=new TLegend(0.6253133,0.1677419,0.9749373,0.3251613);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.03);
    leg2->SetTextFont(42);
    leg2->SetHeader("Extended Displaced Tracks");
    leg2->AddEntry(grExt,"c#tau = 0 cm" ,"L");
    leg2->AddEntry(gr2Ext,"c#tau= 1 cm " ,"L");
    leg2->AddEntry(gr3Ext,"c#tau = 5 cm" ,"L");
    leg2->Draw();
    
    TLegend* leg4=new TLegend(0.75,0.4671916,0.9974937,0.6377953);
    //  TLegend* leg3=new TLegend(0.6240602,0.4348387,0.9736842,0.6051613);
    
    leg4->SetFillStyle(0);
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.033);
    leg4->SetTextFont(42);
    leg4->SetHeader("h#rightarrow #phi#phi");
    leg4->Draw();

    TLegend* leg3=new TLegend(0.6027569,0.4225722,0.952381,0.5931759);
    //  TLegend* leg3=new TLegend(0.6240602,0.4348387,0.9736842,0.6051613);
    
    leg3->SetFillStyle(0);
    leg3->SetBorderSize(0);
    leg3->SetTextSize(0.033);
    leg3->SetTextFont(42);
    leg3->SetHeader(TString::Format("m_{H} = %d GeV, m_{#phi} = %d GeV",mH,mPhi) );
    // leg3->SetHeader("h#rightarrow #phi#phi : m_{H} = 250 GeV, m_{#phi} = 30 GeV");
    leg3->Draw();
   // CMS_lumi( c2, 0, 1 );
    c2->Update();
    DrawPrelimLabel(c2);
    DrawLumiLabel(c2, "");
               
    c2->Print("Results_effRate"+fname+".pdf");
    
    TFile*output=new TFile(TString::Format("output_hist_H%d_mPhi%d.root",mH,mPhi), "RECREATE");
    output->cd();
    gr->Write("PromptCtau0");
    gr2->Write("PromptCtau10");
    gr3->Write("PromptCtau50");
    grExt->Write("DispCtau0");
    gr2Ext->Write("DispCtau10");
    gr3Ext->Write("DispCtau50");
}

