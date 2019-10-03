#include "../tdrstyle.C"
#include "../CMS_lumi.C"
#include "PlotTemplate.C"
Double_t errorFun(Double_t *x, Double_t *par) {

  return 0.5*par[0]*(1. + TMath::Erf( (x - par[1]) / (sqrt(2.)*par[2]) )) ;
  // return par[3] + par[0] * TMath::Freq( (x[0] - par[1]) / par[2] ) * TMath::Freq( (x[0] - par[4] ) / par[5] );

}

Double_t CrystalBall(Double_t x,Double_t *par) {

  Double_t t = (x-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    //return exp(-0.5*t*t);
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha;

    // return (a/TMath::Power(b - t, par[1]));
    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}

Double_t cbgaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      //   mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[7];
      xupp = x[0] + sc * par[7];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
	 fland = CrystalBall(xx,par);
	 //  fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
	 fland = CrystalBall(xx,par);
	 //     fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[8]*step * sum ); //* invsq2pi / par[3]);
}

Double_t EMG(Double_t *x, Double_t *par) {

  Double_t fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*x[0]) ); 
  return par[0] * fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-x[0])/(sqrt(2.)*par[2]) );

}

// Exponentially modified Gaussian
Double_t expgaus(Double_t *x, Double_t *par) {

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 50.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      bool x_initial = false;

      // MP shift correction
      //   mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      //      xlow = x[0] - sc * par[2];
      //xupp = x[0] + sc * par[2];
      xlow = par[4];
      xupp=x[0];
      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;

         //cout << xx << endl;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );
	   //errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;


         //cout << xx << endl;
         //cout << "=========" << endl;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );
         //cout << "Par 3  " << par[3] << " x value " << *x << endl;
      }



      return (par[0] *step * sum ); //* invsq2pi / par[3]);
}

/*
// Exponentially modified Gaussian with derived fit parameters
Double_t derivedexpgaus(Double_t *x, Double_t *par) {

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;
      Double_t lamda;

      // MP shift correction
      //   mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      //      xlow = x[0] - sc * par[2];
      //xupp = x[0] + sc * par[2];
      xlow=-50.;
      xupp=x[0];

      step = (xupp-xlow) / np;
      //lambda = fabs((par[5]*(*x)) + par[4]);
      //cout <<"x is " << *x <<  " lambda is " << lambda << endl;
      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         //par[3] = fabs((par[5]*(*x)) + par[4]);
         par[3] = 1.02;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );
	   //errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
	 fland = (par[3]/2.)*TMath::Exp( (par[3]/2.)*(2.*par[1]+par[3]*par[2]*par[2]-2.*xx) ); 
         sum += fland * TMath::Erfc( (par[1]+par[3]*par[2]*par[2]-xx)/(sqrt(2.)*par[2]) );

      }


      

      return (par[0] *step * sum ); //* invsq2pi / par[3]);
}
*/

Double_t expgaus2(Double_t *x, Double_t *par) {

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      //   mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[1];
      xupp = x[0] + sc * par[1];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
	 fland = (par[2]/2.)*TMath::Exp( (par[2]/2.)*(2.*x[0]+par[2]*par[1]*par[1]-2.*xx) ); 
         sum += fland * TMath::Erfc( (x[0]+par[2]*par[1]*par[1]-xx)/(sqrt(2.)*par[1]) );
	   //errorFun(xx,&par[5]); //TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
	 fland = (par[2]/2.)*TMath::Exp( (par[2]/2.)*(2.*x[0]+par[2]*par[1]*par[1]-2.*xx) ); 
         sum += fland * TMath::Erfc( (x[0]+par[2]*par[1]*par[1]-xx)/(sqrt(2.)*par[1]) );
      }

      return (par[0] *step * sum ); //* invsq2pi / par[3]);
}

TLegend *legend() {
  
 TLegend *leg2 = new TLegend(0.17,0.67,0.52,0.90);
 leg2->SetFillStyle(0);
 leg2->SetBorderSize(0);
 leg2->SetTextSize(0.04);
 leg2->SetTextFont(42); 
 
 return leg2;

}
TH2* readHist2D(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);
 // file->ls();
 // TDirectory* dir = (TDirectory*)file->Get(Dirname);
 // dir->ls();
 TH2* hist = (TH2*)file->Get(nameHist);
 // hist->SetLineWidth(2);
 if(rebin>0) hist->RebinX(rebin);
 // hist->RebinY(2);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);

 return hist;

}

TH1* readHistRef(TString nameHist,TString nameFile, int rebin)
{
 TFile* file = new TFile(nameFile);
 // file->ls();
 // TDirectory* dir = (TDirectory*)file->Get(Dirname);
 // dir->ls();
 TH1* hist = (TH1*)file->Get(nameHist);
 hist->GetSumw2();
 // hist->SetLineWidth(2);
 if(rebin>0) hist->Rebin(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 hist->SetStats(kFALSE);
 return hist;
}

TH1* readHistSel(TString nameHist,TString nameFile, int rebin)
{
 TFile* file =  new TFile(nameFile);
 // file->ls();
 // TDirectory* dir = (TDirectory*)file->Get(Dirname);
 // dir->ls();
 TH1* hist = (TH1*)file->Get(nameHist);
 // hist->SetLineWidth(2);
 hist->GetSumw2();
 hist->SetLineColor(9); hist->SetFillColor(9);
 hist->SetFillStyle(3035);

 if(rebin>0) hist->Rebin(rebin);
 hist->GetXaxis()->SetTitleSize(.055);
 hist->GetYaxis()->SetTitleSize(.055);
 hist->GetXaxis()->SetLabelSize(.05);
 hist->GetYaxis()->SetLabelSize(.05);
 // hist->SetStats(kFALSE);
 

 return hist;
}


TCanvas* getaCanvas(TString name)
{

  TCanvas* aCanvas = new TCanvas(name);

  aCanvas->SetFillColor(0);
  aCanvas->SetBottomMargin(0.125);
  aCanvas->SetLeftMargin(0.125);
  aCanvas->SetFrameFillColor(0);
  aCanvas->SetFrameBorderMode(0);
  aCanvas->SetFrameLineWidth(2);
  return aCanvas;
}

void drawPair(TString hname1, TString hname2, TString filename, TString entry1, TString entry2, TString xtitle, int rebin, float xhigh) {

  TH1D *ref = readHistRef(hname1, filename, rebin);
  TH1D *sel = readHistSel(hname2, filename, rebin);

  TCanvas *c=getaCanvas(hname1+hname2);

  TLegend *leg = legend();
  leg->AddEntry(ref,entry1,"LF");
  leg->AddEntry(sel,entry2,"LF");

  //gPad->SetLogy(); //RefEt->SetStats(kFALSE);
  ref->Draw("E1HIST"); sel->Draw("HISTE1SAMES");

  ref->SetFillStyle(3035); ref->SetFillColor(1);
  sel->SetFillColor(9); sel->SetLineWidth(2);
  ref->GetXaxis()->SetTitle(xtitle);
  ref->GetXaxis()->SetRangeUser(0.,xhigh);

  leg->Draw("SAME");

}

TGraphAsymmErrors* drawEff(TString hname1, TString hname2, TString filename, TString header, TString xtitle, int rebin, bool asymBin, float xlow, float xhigh, float xmin, float xmax, int icol, int imark, TString draw, double mup,int lstyle){//,double y, double c) {

  // gStyle->SetOptStat(1);
  // gStyle->SetOptFit(1);
  cout << asymBin << " Asym bin " << endl;
  TH1D *ref = readHistRef(hname1, filename, rebin);
  TH1D *sel = readHistSel(hname2, filename, rebin);

  ref->SetName(hname1+header);
  sel->SetName(hname2+header);

  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  if (asymBin) {
    Double_t xbins[14]={0.,4.,6.,8.,10.,12.,14.,18.
		      22.,28.,34.,42.,52.,100.};

    sel->Rebin(14,"nsel",xbins);
    ref->Rebin(14,"nref",xbins);
    Eff->BayesDivide(nsel, nref);
  } else {
    Eff->BayesDivide(sel, ref);
    //Eff->Divide(sel,ref,"cl=0.683 b(1,1) mode");
    //Eff->Divide(sel,ref,"cl=0.683 b(1,1) mode v n");

  }

  Double_t Nbins=(double)ref->GetNbinsX();
  Double_t Nconv=(xhigh-xlow)/100.;

  if (Nbins>=Nconv) {
    cout << "Histo bins >= Nconvolution steps; OK" << endl;
  } else {
    cout << "Histo bins < Nconvolution steps!!! Not OK" << endl;
  }


  // gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  TF1 *fermiFunction = new TF1("fermiFunction",expgaus,xlow,xhigh,5);
  //TF1 *derivedfermiFunction = new TF1("derivedfermiFunction",derivedexpgaus,xlow,xhigh,6);

  //  Double_t params[9] = {-2.,10.,-20.,10.,10.,-1.,1.,10.,10.};  // Parameters for cbgaufun          
  Double_t params[5] = {1.,mup,20.,0.05,xlow};
  //TF1 *fermiFunction = new TF1("fermiFunction",errorFun,xlow,xhigh,3);
  //Double_t params[3] = {1.,mup,20.}; // for errorFun

  fermiFunction->SetParameters(params);
  fermiFunction->SetParNames("#epsilon","#mu","#sigma","#lambda","xlow");
  if (hname2 != "L1JetAnalysis/SumEt60"  ){ 
    fermiFunction->FixParameter(4,xlow);
    // fermiFunction->SetParLimits(1,100.,150.);
  }
  //  fermiFunction->SetParLimits(1,100.,150.);
  fermiFunction->SetParLimits(0,0.9,1.0);
  fermiFunction->SetLineColor(icol);
  fermiFunction->SetLineWidth(3);
    fermiFunction->SetLineStyle(lstyle);

  // fermiFunction->Draw();
  //derivedfermiFunction->SetParameters(derivedparams);
  //derivedfermiFunction->SetParNames("#epsilon","#mu","#sigma","#lambda","y","c");
 
  
  for(int i=0; i != 1 ; i++){  
  
   Eff->Fit("fermiFunction","RS+");
   //Eff->Fit("fermiFunction","R");
   //TFitResultPtr r = Eff->Fit("fermiFunction","RS+");

  }
  //Eff->Fit("derivedfermiFunction","R");

  Eff->SetMarkerColor(icol);// EtEff0->SetMarkerColor(1);
  Eff->SetLineColor(icol); Eff->SetLineWidth(3);
  Eff->SetMarkerStyle(imark); Eff->SetMarkerSize(1.);
  Eff->GetXaxis()->SetTitle(xtitle);
  Eff->GetYaxis()->SetTitle("Efficiency");
  Eff->GetXaxis()->SetLabelSize(0.03);
  Eff->GetYaxis()->SetLabelSize(0.03);
  Eff->GetXaxis()->SetTitleSize(0.05);
  Eff->GetYaxis()->SetTitleSize(0.05);
  Eff->GetXaxis()->SetTitleOffset(1.);

  Eff->GetYaxis()->SetRangeUser(0.,1.1);
  Eff->GetXaxis()->SetRangeUser((xmin+0.01),xmax);
  Eff->GetXaxis()->SetMoreLogLabels();

  Double_t xpoint;
  Double_t ypoint;

  TF1 *fit = Eff->GetFunction("fermiFunction");


  /*for(int i =0; i != Eff->GetN() ; i++){  

    Eff->GetPoint(i,xpoint,ypoint);
    cout << xpoint << "   "  << fit->Eval(xpoint) << "    " <<   ypoint << endl;
  }*/

  if (draw=="") {
    Eff->Draw("AP");
  } else {
    Eff->Draw("P");
  }
  
  return Eff;

}

void plot_EtaPhiEff(TString s,TString hname1, TString hname2, TString nameFile, TString hname, TString header, double ETthr, TString same, int icol) 
{
  TFile* file = new TFile(nameFile);

  TH3D *href=(TH3D*)file->Get(hname2);
  TH3D *hsel=(TH3D*)file->Get(hname1);

  int Zbin=hsel->GetZaxis()->FindBin(ETthr);
  //  std::cout << "Threshold bin is : " << Zbin << " out of " << hsel->GetNbinsZ() << " bins" << std::endl;


  TH1D *ref, *sel;
  //  double alpha, beta;

//   alpha=hsel->ProjectionZ(hname+"_sel",0,-1,0,-1)->Integral(Zbin,-1);
//   beta=href->ProjectionZ(hname+"_ref",0,-1,0,-1)->Integral(Zbin,-1);
//      std::cout << " eff = " << alpha/beta << std::endl;


  //  alpha=hsel->ProjectionZ(hname+"_sel",0,-1,0,-1)->Integral(Zbin,-1);
  // beta=href->ProjectionZ(hname+"_ref",0,-1,0,-1)->Integral(Zbin,-1);
  // std::cout << " eff = " << heff->ProjectionZ(hname+"_eff",0,-1,0,-1)->Integral(Zbin,0) << std::endl;

  if (s=="eta") {
    ref=href->ProjectionX(hname+"_ref",0,-1,-1,-1);
    sel=hsel->ProjectionX(hname+"_sel",0,-1,-1,-1);
  } else if (s=="phi") {
    ref=href->ProjectionY(hname+"_ref",0,-1,-1,-1);
    sel=hsel->ProjectionY(hname+"_sel",0,-1,-1,-1);
    //    ref->Rebin(2); sel->Rebin(2);
  }

 //  std::cout << "Sel bins = " << sel->GetEntries() << std::endl;
//   std::cout << "Ref bins = " << ref->GetEntries() << std::endl;


 //  TCanvas *c;
//   if (same=="") {
//     c=getaCanvas(hname1+hname2+s+"_etaeff");
//     gPad->SetGridy();  gPad->SetGridx();
//   }
  TLegend *leg = legend();
  leg->SetHeader(header);

  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  Eff->BayesDivide(sel, ref);

  Eff->GetYaxis()->SetTitle("Efficiency");
  Eff->GetYaxis()->SetRangeUser(0.,1.04);
  Eff->SetMarkerSize(1.); Eff->SetLineWidth(2.);
  if (same=="") {
    Eff->Draw("AP"); } else {
    Eff->Draw("P");
    Eff->SetMarkerColor(icol);
    Eff->SetLineColor(icol);
  }

  if (s=="eta") {
    Eff->GetXaxis()->SetTitle("RecoJet #eta");
  } else if (s=="phi") {
    Eff->GetXaxis()->SetTitle("RecoJet #phi");
  }

  leg->Draw("SAME");


}

void drawEff2(TString hname1, TString hname2, TString filename, TString header, TString xtitle, int rebin, float xlow, float xhigh) {

  // gStyle->SetOptStat(1);
  // gStyle->SetOptFit(1);

  TH1D *ref = readHistRef(hname1, filename, rebin);
  TH1D *sel = readHistSel(hname2, filename, rebin);

  ref->SetName(hname1+header);
  sel->SetName(hname2+header);

  TCanvas *c=getaCanvas(hname1+hname2+"_eff");
  gPad->SetGridy();  gPad->SetGridx();
 
  TLegend *leg = legend();
  leg->SetHeader(header);

  TGraphAsymmErrors *Eff = new TGraphAsymmErrors();
  Eff->BayesDivide(sel, ref);
  
  Eff->GetYaxis()->SetRangeUser(0.5,1.04);
  Eff->Draw("AP");

  leg->Draw("SAME");

}
void draw2D(TString hname, TString filename, TString header, TString xtitle,TString ytitle, int rebin, 
	    bool col, float xlow,float xhigh, float ylow, float yhigh) {

  TH2D *ref = readHist2D(hname, filename, rebin);

  TCanvas *c=getaCanvas(hname);
  gPad->SetGridx(); gPad->SetGridy();

  TLegend *leg = legend();
  leg->SetHeader(header);

  gStyle->SetPalette(1);

  if (col) {
    ref->Draw("COLZ"); 
  } else {
    ref->Draw();
  }
  ref->GetXaxis()->SetTitle(xtitle);
  ref->GetYaxis()->SetTitle(ytitle);
  ref->GetXaxis()->SetRangeUser(xlow,xhigh);
  ref->GetYaxis()->SetRangeUser(ylow,yhigh);

  //  pt->SetLabel("#sqrt{s}=900 GeV");
  //  pt->Draw("SAME");
  leg->Draw("SAME");

}




void effj() {

    TString ifile="StopTurnOn.root";

  TString ifile2="PUJetsOutputStopSampleNoQualityTDRTree.root";

//  TString ifile2="../../Efficiency/PUJetsOutputStopSampleTDRTree.root";

  TString dir="";
  
  int rbin=2;

  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

 // TCanvas *cht=getaCanvas("_hteff");
    TCanvas* cht = CreateCanvas("TurnOn", false, false);


    
  cht->cd();
  //gPad->SetGridy();  gPad->SetGridx();

    TLatex *t1=new TLatex();
    t1->SetNDC();
    t1->SetTextAlign(22);
    t1->SetTextFont(63);
    t1->SetTextSizePixels(22);
    t1->DrawLatex(0.4,0.9,"SUSY Stop m_{#tilde{t}}=1000 m_{LSP}=775");
    Int_t n=5;
    Double_t xx[5], genxx[5];
    Double_t yy[5]={50.,75.,100.,150.,200.};

  //TGraphAsymmErrors *effJETCalo=drawEff("GenLeadJetpT","SingleCaloJetPass",ifile,"E_{T}>130GeV","GEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kRed+2,20,"",75.);
  //TGraphAsymmErrors *effJETTkCalo=drawEff("GenLeadJetpT","SingleCaloTkJetPass",ifile,"E_{T}>130GeV","GEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kRed,20,"same",75.);
/*
    TGraphAsymmErrors *effJETTk=drawEff("GenLeadJetpT","Single2lTkJetPass",ifile,"E_{T}>65GeV","GEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kBlue,20,"",75.);
    TGraphAsymmErrors *effJETTkStubPt=drawEff("GenLeadJetpT","SinglePuppiJetPass",ifile,"E_{T}>55GeV","GEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kCyan,20,"psame",75.);

  


  TLegend* leg1=new TLegend(0.6137124,0.1982609,0.9632107,0.4278261,NULL,"brNDC");
  leg1->SetHeader("SingleJet triggers");
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.04);
    leg1->SetLineColor(1);
    leg1->SetLineStyle(1);
    leg1->SetLineWidth(1);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(0);
 // leg1->AddEntry(effJETCalo,"Calo Jet E_{T}^{L1}>140 GeV","LP");
   // leg1->AddEntry(effJETTkCalo,"Calo Tk Jet E_{T}^{L1}>140 GeV","LP");
    leg1->AddEntry(effJETTk,"Track Jets ","LP");
    leg1->AddEntry(effJETTkStubPt,"PF Puppi Jets","LP");

  leg1->Draw("");
  gPad->SetGridy();  gPad->SetGridx();
  gStyle->SetOptFit(1);
    

//  TF1 *f31=effJETCalo->GetFunction("fermiFunction"); xx[0]=f31->GetX(0.95);
//   TF1*f32=effJETTkCalo->GetFunction("fermiFunction"); xx[1]=f32->GetX(0.95);
    TF1*f33=effJETTk->GetFunction("fermiFunction"); xx[2]=f33->GetX(0.95);
    TF1*f34=effJETTkStubPt->GetFunction("fermiFunction"); xx[3]=f34->GetX(0.95);

  cout << "" << endl; cout << "" << endl;
 //   cout << "Online Calo Jet at 140 (GeV) = " << xx[0] << " GeV (offline)" <<  endl;
 //   cout << "Online Calo TkJet at 140 (GeV) = " << xx[1] << " GeV (offline)" <<  endl;
    cout << "Online TkJet at 57.5 (GeV) = " << xx[2] << " GeV (offline)" <<  endl;
    cout << "Online Puppi at 180 (GeV) = " << xx[3] << " GeV (offline)" <<  endl;
  //  return;
    TGraphAsymmErrors *effSubJETTk=drawEff("GenSubLeadJetpT","Di2lTkJetPass",ifile,"E_{T}>65GeV","SubleadGEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kBlue,20,"",75.);
    TGraphAsymmErrors *effSubJETTkStubPt=drawEff("GenSubLeadJetpT","DiPuppiJetPass",ifile,"E_{T}>30GeV","SubleadGEN-JET E_{T} (GeV)",2,0,40.,400.,0.,700.,kCyan,20,"same",75.);

    TF1*f33=effSubJETTk->GetFunction("fermiFunction"); xx[2]=f33->GetX(0.95);
    TF1*f34=effSubJETTkStubPt->GetFunction("fermiFunction"); xx[3]=f34->GetX(0.95);

    cout << "Online  Sub Lead TkJet at 35 (GeV) = " << xx[2] << " GeV (offline)" <<  endl;
    cout << "Online  Sub Lead Puppi at 127.5 (GeV) = " << xx[3] << " GeV (offline)" <<  endl;

    TLegend* leg2=new TLegend(0.6137124,0.1982609,0.9632107,0.4278261,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetHeader("DiJet triggers");
    leg2->AddEntry(effJETTk,"Track Jets ","LP");
    leg2->AddEntry(effJETTkStubPt,"PF Puppi Jets","LP");

    leg2->Draw("");
   // return;

    TGraphAsymmErrors *effQuadJETTk=drawEff("GenQuadJetpT","Quad2lTkJetPass",ifile,"E_{T}>65GeV","QuadGEN-JET E_{T} (GeV)",2,0,0.,400.,0.,700.,kBlue,20,"",75.);
    TGraphAsymmErrors *effQuadJETTkStubPt=drawEff("GenQuadJetpT","QuadPuppiJetPass",ifile,"E_{T}>5GeV","QuadGEN-JET E_{T} (GeV)",2,0,0.,400.,0.,700.,kCyan,20,"same",75.);

    TF1*f43=effQuadJETTk->GetFunction("fermiFunction"); xx[2]=f43->GetX(0.95);
    TF1*f44=effQuadJETTkStubPt->GetFunction("fermiFunction"); xx[3]=f44->GetX(0.95);

    cout << "Online  Quad Lead TkJet at 14 (GeV) = " << xx[2] << " GeV (offline)" <<  endl;
    cout << "Online  Quad Lead PuppiJet at 73 (GeV) = " << xx[3] << " GeV (offline)" <<  endl;

    TLegend* leg4=new TLegend(0.6137124,0.1982609,0.9632107,0.4278261,NULL,"brNDC");
    leg4->SetBorderSize(0);
    leg4->SetTextSize(0.04);
    leg4->SetLineColor(1);
    leg4->SetLineStyle(1);
    leg4->SetLineWidth(1);
    leg4->SetFillColor(0);
    leg4->SetFillStyle(0);
    leg4->SetHeader("QuadJet triggers");
    leg4->AddEntry(effJETTk,"Track Jets ","LP");
    leg4->AddEntry(effJETTkStubPt,"PF Puppi Jets","LP");

    leg4->Draw("");
    //return;
*/
    TGraphAsymmErrors *effHTJETTk=drawEff("GenHT","hGenJetHT2LTkPass",ifile,"E_{T}>65GeV","Gen Jet H_{T} (GeV)",2,0,0.,800.,0.,800.,kBlack,20,"",75., 1);
    TGraphAsymmErrors *effNoQualHTJETTk=drawEff("GenHT","HT2lTkJetPass",ifile2,"E_{T}>65GeV","GEN-JET H_{T} (GeV)",2,0,300.,800.,300.,800.,kBlack,24,"same",75.,7);

//    TGraphAsymmErrors *effHTJETTkStubPt=drawEff("GenHT","HTPuppiJetPass",ifile,"E_{T}>65GeV","GEN-JET H_{T} (GeV)",2,0,0.,800.,0.,800.,kBlack,20,"same",75.);

    TF1*f43=effHTJETTk->GetFunction("fermiFunction"); xx[2]=f43->GetX(0.95);
    TF1*f44=effNoQualHTJETTk->GetFunction("fermiFunction"); xx[3]=f44->GetX(0.95);

    cout << "Online  HT TkJet at 150 (GeV) = " << xx[2] << " GeV (offline)" <<  endl;
    cout << "Online  HT TkJet at 670 (GeV) = " << xx[3] << " GeV (offline)" <<  endl;
   // return;
    TLegend* leg5=new TLegend(0.6137124,0.1982609,0.9632107,0.4278261,NULL,"brNDC");
    leg5->SetBorderSize(0);
    leg5->SetTextSize(0.025);
    leg5->SetLineColor(1);
    leg5->SetLineStyle(1);
    leg5->SetLineWidth(1);
    leg5->SetFillColor(0);
    leg5->SetFillStyle(0);
    leg5->SetHeader("Track-based H_{T} 25 kHz");
    leg5->AddEntry(f43,"All Track Quality cuts","L");
    leg5->AddEntry(f44,"Track p_{T}>2 GeV","L");

    t1->DrawLatex(0.4,0.9,"SUSY Stop m_{#tilde{t}}=1000 m_{LSP}=775");

   // leg5->AddEntry(effHTJETTkStubPt,"PF Puppi Jets","LP");
    
    leg5->Draw("");
  //  writeExtraText = true;       // if extra text
  //  extraText  = " Phase-II Simulation";
  //  lumi_sqrtS = "PU 200 (14 TeV) ";
    //CMS_lumi( cht, 0, 1 );
    DrawPrelimLabel(cht);
    DrawLumiLabel(cht,"");
    cht->Update();
   return;
    
    TGraphAsymmErrors *effMHTJETTk=drawEff("GenMHT","hGenJetMHT2LTkPass",ifile,"E_{T}>65GeV","Gen Jet Missing H_{T} (GeV)",2,0,0.,800.,0.,800.,kBlack,20,"",75.,1);
    TGraphAsymmErrors *effNoQualMHTJETTk=drawEff("GenMHT","MHT2lTkJetPass",ifile2,"E_{T}>65GeV","GEN-JET Missing H_{T} (GeV)",2,0,0.,800.,0.,800.,kBlack,24,"same",75., 7);

    //TGraphAsymmErrors *effMHTJETTkStubPt=drawEff("GenMHT","MHTPuppiJetPass",ifile,"E_{T}>65GeV","GEN-JET H_{T} (GeV)",2,0,30.,800.,0.,800.,kBlack,20,"same",75.);
    
    TF1*f43=effMHTJETTk->GetFunction("fermiFunction"); xx[2]=f43->GetX(0.95);
    TF1*f44=effNoQualMHTJETTk->GetFunction("fermiFunction"); xx[3]=f44->GetX(0.95);
    
    cout << "Online  MHT TkJet at 160 (GeV) = " << xx[2] << " GeV (offline)" <<  endl;
    cout << "Online  MHT Puppi at 263 (GeV) = " << xx[3] << " GeV (offline)" <<  endl;
    
    TLegend* leg6=new TLegend(0.6137124,0.1982609,0.9632107,0.4278261,NULL,"brNDC");
    leg6->SetBorderSize(0);
    leg6->SetTextSize(0.025);
    leg6->SetLineColor(1);
    leg6->SetLineStyle(1);
    leg6->SetLineWidth(1);
    leg6->SetFillColor(0);
    leg6->SetFillStyle(0);

    leg6->SetHeader("Track-based Missing H_{T} 35 kHz");
    leg6->AddEntry(f43,"All Track Quality cuts ","L");
    leg6->AddEntry(f44,"Track p_{T}>2 GeV","L");

   // leg6->AddEntry(effMHTJETTkStubPt,"PF Puppi Jets","LP");
    t1->DrawLatex(0.4,0.9,"SUSY Stop m_{#tilde{t}}=1000 m_{LSP}=775");

    leg6->Draw("");
    CMS_lumi( cht, 0, 1 );
    cht->Update();
 
}




