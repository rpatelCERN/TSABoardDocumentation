#include "PlotTemplate.C"

void QuickRate(){
   gROOT->ProcessLine(".L ../../tdrStyle.C");
   setTDRStyle();
    //TFile*fin=new TFile("PUJetsOutputMinBias10XTestNewCutsPrompt.root", "READ");
    TFile*fin=new TFile("HybridPromptTestTighterCuts.root", "READ");
    TFile*fin2=new TFile("PUJetsOutputMinBiasNoQualityCutsTDRTree.root", "READ");
    TFile*fin3=new TFile("TestRates1cmDz1cmPromptTrackletEmulation.root", "READ");

    /*
    TFile*fin2=new TFile("PUJetsOutputMinBias2cmTDRTree.root", "READ");
    TFile*fin3=new TFile("PUJetsOutputMinBias3cmTDRTree.root", "READ");
    TFile*fin4=new TFile("PUJetsOutputMinBias6cmTDRTree.root", "READ");
*/
    

   TwoLMinBiasHT=(TH1F*)fin->Get("hTwoLJetHT");
   L1RateTwoLMinBiasHT=(TH1F*)TwoLMinBiasHT->Clone("L1RateTwoLMinBiasHT");
    
    TwoLMinBiasMHT=(TH1F*)fin->Get("TwoLMHT");
    L1RateTwoLMinBiasMHT=(TH1F*)TwoLMinBiasMHT->Clone("L1RateTwoLMinBiasMHT");
    
    
    //TwoLMinBiasMET=(TH1F*)fin->Get("TwoLMET");
    //L1RateTwoLMinBiasMET=(TH1F*)TwoLMinBiasMET->Clone("L1RateTwoLMinBiasMET");
    
    TwoLMinBiasHTNoCuts=(TH1F*)fin2->Get("TwoLHT");
    L1RateTwoLMinBiasHTNoCuts=(TH1F*)TwoLMinBiasHT->Clone("L1RateTwoLMinBiasHTNoCuts");
    
    TwoLMinBiasMHTNoCuts=(TH1F*)fin2->Get("TwoLMHT");
    L1RateTwoLMinBiasMHTNoCuts=(TH1F*)TwoLMinBiasMHT->Clone("L1RateTwoLMinBiasMHTNoCuts");
    
    for(unsigned int i=1; i<=TwoLMinBiasHT->GetNbinsX(); ++i){
        //need to check for zero bins
     
        L1RateTwoLMinBiasHT->SetBinContent(i,2760*11.246*TwoLMinBiasHT->Integral(i,200));
        L1RateTwoLMinBiasHTNoCuts->SetBinContent(i,2760*11.246*TwoLMinBiasHTNoCuts->Integral(i,200));

        L1RateTwoLMinBiasMHT->SetBinContent(i,2760*11.246*TwoLMinBiasMHT->Integral(i,200));
        //L1RateTwoLMinBiasMET->SetBinContent(i,2760*11.246*TwoLMinBiasMET->Integral(i,200));
 L1RateTwoLMinBiasMHTNoCuts->SetBinContent(i,2760*11.246*TwoLMinBiasMHTNoCuts->Integral(i,200));

 
    }
    TCanvas*c1=new TCanvas("c1", "", 800,800);
    leg=new TLegend(0.5263158,0.7260813,0.9423559,0.9174312,NULL,"brNDC");
    writeExtraText = true;       // if extra text
    extraText  = " Phase-II Simulation";
    lumi_sqrtS = "PU 200 (14 TeV) ";

  //  leg->AddEntry(L1RateCaloMinBiasSingleJet, "L1 Calo Jets ", "L");
    //leg->AddEntry(L1RateCaloTkMinBiasSingleJet, "L1 Calo Jets +Tk ", "L");
   // L1RateFastJetMinBiasSingleJet->SetLineStyle(kDashed);
  // leg->AddEntry(L1RateFastJetMinBiasSingleJet, " L1 TwoLayer Jets (Old cuts)", "L");
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
    
    
//    return;

    L1RateTwoLMinBiasMHT->GetXaxis()->SetRangeUser(0,200);
    L1RateTwoLMinBiasMHT->SetTitle("MinBias Events PU 200; Trk H^{miss}_{T} [GeV]; L1 Rate (kHz)");
    
    L1RateTwoLMinBiasMHT->GetYaxis()->SetTitleOffset(1.0);
    L1RateTwoLMinBiasMHT->SetLineColor(kBlack);L1RateTwoLMinBiasMHT->SetLineWidth(3.0);
    
    L1RateTwoLMinBiasMHT->Draw("hist");
    std::cout<<"MHT Threshold at 35kHz Fast Jets "<<L1RateTwoLMinBiasMHT->GetBinLowEdge(L1RateTwoLMinBiasMHT->FindLastBinAbove(35))<<std::endl;
   // return;
    leg->Draw();
    
    // return;
    std::cout<<"HT Threshold at 25kHz Fast Jets "<<L1RateTwoLMinBiasHTNoCuts->GetBinLowEdge(L1RateTwoLMinBiasHTNoCuts->FindLastBinAbove(25))<<std::endl;
    c1->SetGridx();
    c1->SetGridy();
    L1RateTwoLMinBiasHT->GetYaxis()->CenterTitle(true);
    L1RateTwoLMinBiasHT->GetYaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->GetXaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasHT->SetLineColor(kBlack);
    
    leg->Draw("hist");
    //CMS_lumi( c1, 0, 1 );
    c1->Update();
    c1->Print("rates_HT.pdf");
    
    TFile*fsig=new TFile("StopTurnOn.root", "READ");
    TFile*fsig2=new TFile("PUJetsOutputStopSampleNoQualityTDRTree.root", "READ");

    HTSigEff=(TH1D*)fsig->Get("hTwoLJetHT");
    HTSigEffNoCuts=(TH1D*)fsig2->Get("TwoLHT");
    MHTSigEff=(TH1D*)fsig->Get("TwoLMHT");
    MHTSigEffNoCuts=(TH1D*)fsig2->Get("TwoLMHT");

    HTSigEff->Scale(1.0/HTSigEff->Integral());
    HTSigEffNoCuts->Scale(1.0/HTSigEffNoCuts->Integral());
   // MHTSigEff->Scale(1.0/MHTSigEff->Integral());
    TGraph*gr=new TGraph();
    TGraph*gr2=new TGraph();

    TGraph*gr3=new TGraph();
    TGraph*gr4=new TGraph();

    for(unsigned int i=1; i<L1RateTwoLMinBiasHT->GetNbinsX(); ++i){
        if(L1RateTwoLMinBiasHT->GetBinContent(i)>50)continue;
       // std::cout<<"Signal Eff "<<HTSigEff->Integral(i,200)<<std::endl;
                gr->SetPoint(gr->GetN(), HTSigEff->Integral(i,200), L1RateTwoLMinBiasHT->GetBinContent(i) );
    }
  
    for(unsigned int i=1; i<L1RateTwoLMinBiasMHT->GetNbinsX(); ++i){
        if(L1RateTwoLMinBiasMHT->GetBinContent(i)>=50)continue;
        // std::cout<<"Signal Eff "<<HTSigEff->Integral(i,200)<<std::endl;
        gr2->SetPoint(gr2->GetN(), MHTSigEff->Integral(i,200), L1RateTwoLMinBiasMHT->GetBinContent(i) );
    }
    
    for(unsigned int i=1; i<L1RateTwoLMinBiasMHTNoCuts->GetNbinsX(); ++i){
        if(L1RateTwoLMinBiasMHTNoCuts->GetBinContent(i)>=50)continue;
        // std::cout<<"Signal Eff "<<HTSigEff->Integral(i,200)<<std::endl;
        gr4->SetPoint(gr4->GetN(), MHTSigEffNoCuts->Integral(i,200), L1RateTwoLMinBiasMHTNoCuts->GetBinContent(i) );
    }
   
    for(unsigned int i=1; i<L1RateTwoLMinBiasHTNoCuts->GetNbinsX(); i++){
        if(L1RateTwoLMinBiasHTNoCuts->GetBinContent(i)>50)continue;
        gr3->SetPoint(gr3->GetN(), HTSigEffNoCuts->Integral(i,200), L1RateTwoLMinBiasHTNoCuts->GetBinContent(i) );

    }
   /// TCanvas*c2=new TCanvas("c2", "", 800,800);
    TCanvas* c2 = CreateCanvas("c2", false, false);

    c2->SetLogy();
    writeExtraText = true;       // if extra text
    extraText  = " Phase-II Simulation";
    lumi_sqrtS = "PU 200 (14 TeV) ";
    c1->SetLogy();
    gr->SetTitle("; Signal Efficiency; L1 Track H_{T} Rate ( kHz )");
    gr->SetLineWidth(4.0);;gr->SetLineColor(kBlack);gr->SetLineStyle(kBlack);
    gr3->SetLineWidth(2.0);;gr3->SetLineColor(kBlack);gr3->SetLineStyle(kDashed);

    gr->GetXaxis()->SetRangeUser(0.0, 1.0);
    gr->GetYaxis()->SetRangeUser(1.0, 50);

    gr->Draw("AL");
    gr3->Draw("LSame");
    //HT plot
    DrawPrelimLabel(c2);
    DrawLumiLabel(c2,"");
//    CMS_lumi( c2, 0, 1 );
    c2->Update();
    //gr4->Draw("LSame");
    leg2=new TLegend(0.575188,0.1845161,0.820802,0.3767742,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.025);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(gr3, "All Tracks p_{T}>2 GeV", "L");
    leg2->AddEntry(gr, "Tracks passing quality criteria", "L");
    
    leg2->Draw();
   // return;
    gr2->SetTitle("; Signal Efficiency; L1 Track H^{miss}_{T} Rate ( kHz )");

    gr2->SetLineWidth(4.0);;gr2->SetLineColor(kBlack);gr2->SetLineStyle(kBlack);
    gr4->SetLineWidth(2.0);;gr4->SetLineColor(kBlack);gr4->SetLineStyle(kDashed);
    gr2->GetXaxis()->SetRangeUser(0.0, 1.0);
    gr2->GetYaxis()->SetRangeUser(1.0, 50);
    gr2->Draw("AL");
    gr4->Draw("LSame");
   // CMS_lumi( c2, 0, 1 );
    DrawPrelimLabel(c2);
    DrawLumiLabel(c2,"");
    c2->Update();
    leg2->Draw();

    
    
    
    
    return;
    
    leg2=new TLegend(0.5263158,0.7260813,0.9423559,0.9174312,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    TLatex *t1=new TLatex();
    t1->SetNDC();
    t1->SetTextAlign(22);
    t1->SetTextFont(63);
    t1->SetTextSizePixels(22);
    t1->DrawLatex(0.4,0.9,"SUSY Stop m_{#tilde{t}}=1000 m_{LSP}=775");

   // leg->SetTItle("Stop: m_{#tilde{t}}=1000 m_{LSP}=775");
    //eg2->AddEntry(gr, "Stop Sample: m_{#tilde{t}}=1000 m_{LSP}=775","L");
    //leg2->AddEntry(gr3, "H(250) m#phi(30) c#tau=5cm (Old 93X sample)","L");
    //leg2->AddEntry(gr4, "H(250) m#phi(30) c#tau=5cm (New 10X sample)","L");
    leg2->Draw();
    L1RateTwoLMinBiasMHT->SetTitle("MinBias Events PU 200; Missing H_{T} [GeV]; L1 Rate (kHz)");
    L1RateTwoLMinBiasMHT->GetXaxis()->SetRangeUser(0,200);
    L1RateTwoLMinBiasMHT->SetLineColor(kBlack);L1RateTwoLMinBiasMHT->SetLineWidth(3.0);
    

    //L1RateTwoLMinBiasMET->SetLineColor(kBlack);L1RateTwoLMinBiasMET->SetLineWidth(3.0);

    L1RateTwoLMinBiasMHT->GetYaxis()->CenterTitle(true);
    L1RateTwoLMinBiasMHT->GetYaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasMHT->GetXaxis()->SetTitleSize(0.05);
    L1RateTwoLMinBiasMHT->SetLineColor(kBlack);
    L1RateTwoLMinBiasMHT->Draw("hist");
   
   // L1RateTwoLMinBiasMET->Draw("same");

    
//   leg2->AddEntry(L1RateTwoLMinBiasMET, "L1 Track MET ", "L");
    leg->Draw();
    CMS_lumi( c1, 0, 1 );
    c1->Update();

    c1->Print("rates_MHT.pdf");
    std::cout<<"MHT Threshold at 35kHz Fast Jets "<<L1RateTwoLMinBiasMHT->GetBinLowEdge(L1RateTwoLMinBiasMHT->FindLastBinAbove(35));
   // <<" MET Threshold at 35kHz "<< L1RateTwoLMinBiasMET->GetBinLowEdge(L1RateTwoLMinBiasMET->FindLastBinAbove(35))<<std::endl;

  //  L1RateTwoLMinBiasMET->Draw();


}

