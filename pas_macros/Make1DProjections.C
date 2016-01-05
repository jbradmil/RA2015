
#include <iostream>
#include <vector>
#include <math.h> 
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "THStack.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TFileCollection.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TText.h"
#include "TGraphAsymmErrors.h"
//#include "<ROOT/Math.h>"
#include "inc/jack_style.h"
#include "inc/CMS_lumi.C"
#include "Math/QuantFuncMathCore.h"

using namespace std;

const double alpha = 1 - 0.6827;


TString plotdir = "plots/";
//const double predSF = 1.63947463;

TGraphAsymmErrors* GetGDataObs(TH1D* hdata_obs, bool suppress_zeroes=false, bool skip1=false) {

  int nbinsx=0;
  if (!suppress_zeroes) nbinsx=hdata_obs->GetNbinsX();
  else {
    for (int bin(0); bin<hdata_obs->GetNbinsX(); bin++) {
      if (hdata_obs->GetBinContent(bin+1)>0) nbinsx++;
    }
  }
    
  const int nbins=nbinsx;
  
  Double_t x[nbins];
  Double_t xl[nbins];
  Double_t xh[nbins];

  
  Double_t data_cv[nbins];
  Double_t data_pois_up[nbins];
  Double_t data_pois_down[nbins];

  if (skip1) {
    for (int bin(1); bin<hdata_obs->GetNbinsX(); bin++) {
      x[bin] = hdata_obs->GetBinCenter(bin+1);
      xl[bin]=hdata_obs->GetBinWidth(bin+1)/1000.;
      xh[bin]=hdata_obs->GetBinWidth(bin+1)/1000.;
      data_cv[bin]=hdata_obs->GetBinContent(bin+1);
      double N=data_cv[bin];
      double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U =  (N==0) ? 0  : ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      data_pois_up[bin]=(U-N);
      data_pois_down[bin]=(N-L);
    }
  } else {
    for (int bin(0); bin<nbins; bin++) {
      x[bin] = hdata_obs->GetBinCenter(bin+1);
      xl[bin]=hdata_obs->GetBinWidth(bin+1)/1000.;
      xh[bin]=hdata_obs->GetBinWidth(bin+1)/1000.;
      data_cv[bin]=hdata_obs->GetBinContent(bin+1);
      double N=data_cv[bin];
      double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
      data_pois_up[bin]=(U-N);
      data_pois_down[bin]=(N-L);
    }
  }

  TGraphAsymmErrors* gdata_obs = new TGraphAsymmErrors(nbins, x, data_cv, xl, xh, data_pois_down, data_pois_up);

  return gdata_obs;

}


TGraphAsymmErrors* GetBGErr(TString graph_name, TH1D* hdata_obs, TFile* f_lostlep, TFile* f_hadtau, TFile* f_qcd, TFile* f_znn) {

  TString full_name = "g"+graph_name+"Full";
  TString stat_name = "g"+graph_name+"Stat";
  TString syst_name = "g"+graph_name+"Syst";

  TH1D* hqcd = (TH1D*) f_qcd->Get("hPred"+graph_name);
  TH1D* hlostlep = (TH1D*) f_lostlep->Get("hPred"+graph_name);
  TH1D* hhadtau = (TH1D*) f_hadtau->Get("hPred"+graph_name);
  TH1D* hznn = (TH1D*) f_znn->Get("hPred"+graph_name);

  TGraphAsymmErrors* glostlepstat = (TGraphAsymmErrors*) f_lostlep->Get(stat_name);
  TGraphAsymmErrors* ghadtaustat = (TGraphAsymmErrors*) f_hadtau->Get(stat_name);
  TGraphAsymmErrors* glostlepsyst = (TGraphAsymmErrors*) f_lostlep->Get(syst_name);
  TGraphAsymmErrors* ghadtausyst = (TGraphAsymmErrors*) f_hadtau->Get(syst_name);
  TGraphAsymmErrors* gqcdstat = (TGraphAsymmErrors*) f_qcd->Get(stat_name);
  TGraphAsymmErrors* gznnstat = (TGraphAsymmErrors*) f_znn->Get(stat_name);
  TGraphAsymmErrors* gqcdsyst = (TGraphAsymmErrors*) f_qcd->Get(syst_name);
  TGraphAsymmErrors* gznnsyst = (TGraphAsymmErrors*) f_znn->Get(syst_name);


  const int nbins = hdata_obs->GetNbinsX();

  Double_t x[nbins];
  Double_t xl[nbins];
  Double_t xh[nbins];

  
  Double_t pred_cv[nbins];
  Double_t full_stat_up[nbins];
  Double_t full_stat_down[nbins];
  Double_t full_syst_up[nbins];
  Double_t full_syst_down[nbins];
  Double_t full_err_up[nbins];
  Double_t full_err_down[nbins];

  for (int bin(0); bin<nbins; bin++) {
    x[bin] = hdata_obs->GetBinCenter(bin+1);
    xl[bin]=hdata_obs->GetBinWidth(bin+1)/2.;
    xh[bin]=hdata_obs->GetBinWidth(bin+1)/2.;
    pred_cv[bin]=hlostlep->GetBinContent(bin+1) + hhadtau->GetBinContent(bin+1) + hqcd->GetBinContent(bin+1)+ hznn->GetBinContent(bin+1);

    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    // double wtop_syst_up = sqrt(pow(glostlepsyst->GetErrorYhigh(bin)+ghadtausyst->GetErrorYhigh(bin),2.));
    // double wtop_syst_down = sqrt(pow(glostlepsyst->GetErrorYlow(bin)+ghadtausyst->GetErrorYlow(bin),2.));
    full_stat_up[bin] = sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.));
    full_stat_down[bin] = sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.));
    // printf("Bin %d, LL stat err: + %3.3f - %3.3f\n",bin+1, glostlepstat->GetErrorYhigh(bin), glostlepstat->GetErrorYlow(bin));
    // printf("Bin %d, LL syst err: + %3.3f - %3.3f\n",bin+1, glostlepsyst->GetErrorYhigh(bin), glostlepsyst->GetErrorYlow(bin));
    // printf("Bin %d, wtop err: + %3.3f - %3.3f\n",bin+1,wtop_err_up, wtop_err_down);
    full_syst_up[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.));
    full_syst_down[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.));
    // full_err_up[bin] = sqrt(pow(wtop_err_up,2.)+pow(gqcd->GetErrorYhigh(bin),2.)+pow(gznn->GetErrorYhigh(bin),2.));
    // full_err_down[bin] = sqrt(pow(wtop_err_down,2.)+pow(gqcd->GetErrorYlow(bin),2.)+pow(gznn->GetErrorYlow(bin),2.));
    full_err_up[bin] = sqrt(pow(full_stat_up[bin], 2.)+pow(full_syst_up[bin], 2.));
    full_err_down[bin] = sqrt(pow(full_stat_down[bin], 2.)+pow(full_syst_down[bin], 2.));
  }

  TGraphAsymmErrors* gBGErr = new TGraphAsymmErrors(nbins, x, pred_cv, xl, xh, full_err_down, full_err_up);

  return gBGErr;

}

void MakePlot(TString plot_title, TGraphAsymmErrors* gdata_obs, TGraphAsymmErrors* gerr, TH1D* hdata_obs, TH1D* hlostlep, TH1D* hhadtau, TH1D* hqcd, TH1D* hznn, TH1D* hsig1, TH1D* hsig2, TString siglabel, TString siglabel2, /*TH1D* hsig2, TH1D* hsig3,*/ bool logy=false, double manual_max=-1., TString cut1="", TString cut2="", TString cut3="")
{

  hdata_obs->Sumw2();
  gStyle->SetEndErrorSize(0);

  set_style_lite(hlostlep, "lost_lep");
  set_style_lite(hhadtau, "had_tau");
  set_style_lite(hqcd, "qcd");
  set_style_lite(hznn, "znn");
  set_style(hdata_obs, "data_obs");
  set_style_lite(hsig1, "single_top");
  set_style_lite(hsig2, "single_top");


  hsig1->SetLineWidth(3);
  hsig1->SetLineColor(4);
  hsig1->SetFillColor(0);

  hsig2->SetLineWidth(4);
  hsig2->SetLineStyle(7);
  hsig2->SetLineColor(4);
  hsig2->SetFillColor(0);

  //  cout << "Sum up the BGs" << endl;
  TH1D * hbg_pred = (TH1D*)hlostlep->Clone("bg_pred");
  hbg_pred->Reset();
  hbg_pred->SetTitle("");
  hbg_pred->GetYaxis()->SetTitle("Events");




  hbg_pred->Add(hlostlep);
  hbg_pred->Add(hhadtau);
  hbg_pred->Add(hqcd);
  hbg_pred->Add(hznn);


  //hbg_pred->Print("all");
  gdata_obs->Print("all");
  gerr->Print("all");


  hbg_pred->SetMarkerSize(0);
  hbg_pred->SetMarkerColor(0);
  hbg_pred->SetLineWidth(0);
  hbg_pred->SetLineColor(0);


 
  THStack * hs = new THStack("hs", "");
  hs->Add(hqcd);
  hs->Add(hznn);
  hs->Add(hhadtau);
  hs->Add(hlostlep);

  gdata_obs->SetLineWidth(1);
  gdata_obs->SetMarkerStyle(20);
  gdata_obs->SetMarkerSize(1.75);
  gdata_obs->SetLineColor(1);
  
  
  gerr->SetFillColor(14);
  gerr->SetMarkerSize(0);
  gerr->SetLineWidth(0);
  gerr->SetLineColor(0);
  gerr->SetFillStyle(3445);

  hdata_obs->SetBinErrorOption(TH1::kPoisson);
  
  //cout << "Compute ratio hist..." << endl;
  TH1D * ratio = (TH1D *) hdata_obs->Clone("ratio");
  TH1D * hratiogerr = (TH1D *) hdata_obs->Clone("hratiogerr");
  TLine* ratiounity = new TLine(hbg_pred->GetBinLowEdge(1),0,hbg_pred->GetBinLowEdge(hbg_pred->GetNbinsX()+1),0);

  TGraphAsymmErrors* ratioderr = new TGraphAsymmErrors(gerr->GetN(), gerr->GetX(), gerr->GetY(), gerr->GetEXlow(), gerr->GetEXhigh(), gerr->GetEYlow(), gerr->GetEYhigh());
  TGraphAsymmErrors* ratiogerr = new TGraphAsymmErrors(gerr->GetN(), gerr->GetX(), gerr->GetY(), gerr->GetEXlow(), gerr->GetEXhigh(), gerr->GetEYlow(), gerr->GetEYhigh());
  for (Int_t i = 0; i < gerr->GetN(); i++) {
    ratiogerr->SetPoint(i, ratio->GetBinCenter(i+1), 0.);
    if (hbg_pred->GetBinContent(i+1)>0) {
      ratiogerr->SetPointError(i, ratio->GetBinWidth(i+1)/2., ratio->GetBinWidth(i+1)/2., ratiogerr->GetErrorYlow(i)/hbg_pred->GetBinContent(i+1), ratiogerr->GetErrorYhigh(i)/hbg_pred->GetBinContent(i+1));
      if (hdata_obs->GetBinContent(i+1)>0) {
	ratio->SetBinContent(i+1, (hdata_obs->GetBinContent(i+1)-hbg_pred->GetBinContent(i+1))/hbg_pred->GetBinContent(i+1));
	ratioderr->SetPoint(i, ratio->GetBinCenter(i+1), (hdata_obs->GetBinContent(i+1)-hbg_pred->GetBinContent(i+1))/hbg_pred->GetBinContent(i+1));
	ratioderr->SetPointError(i, ratio->GetBinWidth(i+1)/1000., ratio->GetBinWidth(i+1)/1000., gdata_obs->GetErrorYlow(i)/hbg_pred->GetBinContent(i+1), gdata_obs->GetErrorYhigh(i)/hbg_pred->GetBinContent(i+1));
      }
      else {
	ratioderr->SetPoint(i, ratio->GetBinCenter(i+1), -999.);
      }
    }
    else {
      ratiogerr->SetPointError(i, ratio->GetBinWidth(i+1)/2., ratio->GetBinWidth(i+1)/2., 0., 1.);
      ratioderr->SetPoint(i, ratio->GetBinCenter(i+1), -999.);
    }
    //    printf("Bin %d: pred = %3.2f + %3.2f - %3.2f; ratiogerr: + %3.2f - %3.2f\n", i+1, hbg_pred->GetBinContent(i+1), gerr->GetErrorYhigh(i), gerr->GetErrorYlow(i), ratiogerr->GetErrorYhigh(i), ratiogerr->GetErrorYlow(i));
  }
  set_style(ratio, "data_obs");
  ratioderr->SetLineWidth(1);
  ratioderr->SetMarkerStyle(20);
  ratioderr->SetMarkerSize(1.75);
  ratioderr->SetLineColor(1);
  hratiogerr->SetStats(0);
  ratio->SetTitle(hbg_pred->GetTitle());
  ratio->GetYaxis()->SetTitle("#frac{(Obs.-Exp.)}{Exp.}");
  ratio->SetMaximum(1.15);
  ratio->SetMinimum(-1.15);
  ratio->SetMarkerSize(0);
  hratiogerr->SetFillColor(0);
  hratiogerr->SetFillStyle(0);

  ratiogerr->SetFillColor(14);
  ratiogerr->SetMarkerSize(0);
  ratiogerr->SetLineWidth(1504);
  ratiogerr->SetLineColor(0);
  //ratiogerr->SetFillStyle(3004);
  ratiogerr->SetFillStyle(3445);
    
  ratio->GetXaxis()->SetLabelSize(0.15);
  ratio->GetXaxis()->SetLabelOffset(0.03);
  ratio->GetXaxis()->SetTitleSize(0.22);
  ratio->GetXaxis()->SetTitleOffset(0.95);
  ratio->GetYaxis()->SetLabelSize(0.10);
  ratio->GetYaxis()->SetNdivisions(505);
  ratiounity->SetLineStyle(2);
 


    
  // // Setup legends
  TLegend * legdata;
  TLegend * leg1;
  TLegend * leg2;
  TLegend * leg3;
  TLegend * leg4;
  TLegend * legsig;

  legdata = new TLegend(0.25-0.08, 0.75, 0.45-0.12, 0.85);
  leg1 = new TLegend(0.29, 0.756, 0.58, 0.85);
  leg2 = new TLegend(0.46, 0.756, 0.75, 0.85);
  leg3 = new TLegend(0.67, 0.756, 0.96, 0.85);
  leg4 = new TLegend(0.83, 0.756, 1.12, 0.85);

  
  set_style(legdata,0.04);
  set_style(leg1,0.04);
  set_style(leg2,0.04);
  set_style(leg3,0.04);
  set_style(leg4,0.04);
  legdata->AddEntry(gdata_obs, "Data", "pe");
  leg1->AddEntry(hlostlep, "#splitline{Lost}{lepton}", "f");
  leg2->AddEntry(hhadtau, "#splitline{Hadronic}{#tau lepton}", "f");
  leg3->AddEntry(hznn, "Z#rightarrow#nu#bar{#nu}", "f");
  leg4->AddEntry(hqcd, "QCD", "f");
  legsig = new TLegend(0.17, 0.55, 0.65, 0.75);
  set_style(legsig,0.04);
  legsig->SetMargin(0.2);
  legsig->AddEntry(hsig1, siglabel, "l");
  legsig->AddEntry(hsig2, siglabel2, "l");

 
  double ymax = hbg_pred->GetMaximum();
  //  double smax = hsig1->GetMaximum();
  if (hdata_obs->GetMaximum()>ymax) ymax=hdata_obs->GetMaximum();

  if (manual_max<0) {
    if(logy) {
      hbg_pred->SetMaximum(1000*ymax);
      hbg_pred->SetMinimum(0.07);
    }
    else {
      hbg_pred->SetMinimum(0.00000001);
      hbg_pred->SetMinimum(0);
      hbg_pred->SetMaximum(1.5*ymax);
    }
  }
  else {
    //    hbg_pred->SetMinimum(0.00000001);
    hbg_pred->SetMinimum(0);
    hbg_pred->SetMaximum(manual_max);
  }


  // Setup canvas and pads

  int W = 800;
  int H = 800;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 800; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.1*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TCanvas* canv = new TCanvas("canvName","canvName",50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  double up_height     = 0.8;  // please tune so that the upper figures size will meet your requirement
  double dw_correction = 1.30; // please tune so that the smaller canvas size will work in your environment
  double font_size_dw  = 0.2;  // please tune the font size parameter for bottom figure
  double dw_height     = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02; // KH, added to put the bottom one closer to the top panel

  
  TPad * pad1 = new TPad("pad1", "top pad" , 0.0, 0.4, 1.0, 1.0);
  TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.2);

  pad1->SetTickx(0);
  pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height,    1., 1.00);
  //
  pad1->SetFrameFillColor(0);
  pad1->SetFillColor(0);
  pad1->SetTopMargin(0.12);
  pad1->SetLeftMargin(0.15);
  // pad1->SetBottomMargin(0.0);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0);
  pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.38);
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy(logy);
 
  // // Draw hists

  hbg_pred->Draw();

  hbg_pred->GetYaxis()->SetLabelSize(0.035*1.4);
  hbg_pred->GetYaxis()->SetTitleSize(0.04*1.34);
  hbg_pred->GetYaxis()->SetTitleOffset(0.9);
  hbg_pred->GetYaxis()->SetTitleFont(42);
  
  hs->Draw("hist,same");

  hbg_pred->GetXaxis()->SetLabelSize(0);
  
  hdata_obs->SetMarkerStyle(20);
  set_style(hdata_obs, "data_obs");
  gdata_obs->SetMarkerStyle(20);
  gerr->Draw("2 same");

  hsig1->Draw("hist, same");
  hsig2->Draw("hist, same");

  gdata_obs->Draw("p same");

  // // Draw legends
  legdata->Draw();
  leg1->Draw();
  leg2->Draw();
  leg3->Draw();
  leg4->Draw();
  legsig->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.06);

  TString cut = cut1;
  if (cut2!="") cut+=(", "+cut2);
  if (cut3!="") cut+=(", "+cut3);

  cout << "Cut length: " << cut.Length() << endl;
  if (cut.Length()<20) latex->DrawLatex(0.77, 0.49, cut);
  else if (cut.Length()<35) latex->DrawLatex(0.6, 0.49, cut);
  else if (cut.Length()<40) latex->DrawLatex(0.44, 0.49, cut);
  else latex->DrawLatex(0.27, 0.49, cut);
  latex->SetTextSize(0.052);


  // Luminosity information for scaling
  double lumi     = 2.153738; // normaliza to this lumi (fb-1)
  //  double lumi_ref = 2.109271; // normaliza to 3 (fb-1)

  char tempname[200];
  TString line = "";
  sprintf(tempname,"%8.1f",lumi);
  line+=tempname;
  line+=" fb^{-1} (13 TeV)";
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=0;
    
  writeExtraText = false;
  extraText   = "        Preliminary";
  //TString lumi_sqrtS = line;
  
  TPaveText * pave = new TPaveText(0.18, 0.86, 0.4, 0.96, "brNDC");
  //  TText * text = NULL; 
  TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
  
  pad2->cd();
  pad2->SetGridy(0);
   
  //  hratiogerr->Draw("e2");
  set_style(ratio, "data_obs");
  ratio->Draw("axis");

  ratio->GetXaxis()->SetLabelSize(font_size_dw);
  ratio->GetXaxis()->SetTitleSize(font_size_dw);
  ratio->GetYaxis()->SetLabelSize(font_size_dw);
  ratio->GetYaxis()->SetTitleSize(font_size_dw);

  ratio->GetXaxis()->SetLabelSize(0.14*1.15);
  ratio->GetXaxis()->SetTitleSize(0.14*1.2);
  ratio->GetXaxis()->SetTitleOffset(1.);
  ratio->GetXaxis()->SetTitleFont(42);
  ratio->GetYaxis()->SetTitleSize(0.12*1.18);
  ratio->GetYaxis()->SetLabelSize(0.12*1.18);
  ratio->GetYaxis()->SetTitleOffset(0.46);
  ratio->GetYaxis()->SetTitleFont(42);
  //ratio->GetXaxis()->SetTitle("Search region bin number");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTickLength(0.015);
  ratio->GetXaxis()->SetTickLength(0.08);
  // ratio->GetXaxis()->SetTitleSize(0.12);
  // ratio->GetXaxis()->SetLabelSize(0.12);
  hbg_pred->GetXaxis()->SetTitleSize(0.2);
  //    ratio->GetXaxis()->SetLabelSize(0.12);
  // //ratiosysterr->Draw("e2 same");
  //ratio->Draw("e2 same");
  ratiounity->Draw();
  ratiogerr->Draw("2 same");
  ratioderr->Draw("p same");
  

  pave->SetLineColor(0);
  pave->SetLineWidth(0);
  pave->SetFillStyle(4000);
  pave->SetShadowColor(0);
  pave->SetBorderSize(1);
  // double nchisq = hdata_obs->Chi2Test(hbg_pred, "UWCHI2/NDF"); // MC uncert. (stat)
  //  double p_value = hdata_obs->Chi2Test(hbg_pred, "UW"); // MC uncert. (stat)
  // //double kolprob = hdata_obs->KolmogorovTest(hbg_pred); // MC uncert. (stat)
  // //TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
  //text = pave->AddText(Form("#chi_{#nu}^{2}/ndf = %.3f, p = %.3f", nchisq, p_value));
  // text = pave->AddText(Form(""));
  // text->SetTextFont(62);
  // text->SetTextSize(0.07);
  // text->SetTextSize(0.06);
  // pave->Draw();
  
  set_style(ratioleg);
  ratioleg->SetTextSize(0.07);
  ratioleg->AddEntry(gerr, "Pred. uncert. (stat#oplussyst)", "f");  
  //ratioleg->Draw();
  pad1->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  pad2->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();

  
  canv->cd();
  CMS_lumi(canv, iPeriod, iPos, line);

 
  gPad->Print(plotdir+plot_title+".pdf");
  gPad->Print(plotdir+plot_title+".png");


  // Clean up
  // delete hqcd;
  // delete hznn;
  // delete hdata_obs;
  delete hbg_pred;

  // delete gerr;
  delete ratio;
  delete hratiogerr;
  delete ratiogerr;
  delete ratioderr;
  delete ratioleg;
  delete pad1;
  delete pad2;
  // delete ratiosysterr;
  delete hs;
  delete legdata;
  delete leg1;
  delete leg2;
  delete legsig;
  delete latex;
  delete pave;
  delete canv;

  cout << "SaveHist(): DONE!" << endl;


  return;
}



void Make1DProjections() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);


  if (gSystem->AccessPathName(plotdir))
    gSystem->mkdir(plotdir);
  // gInterpreter->GenerateDictionary("vector<TLorentzVector>","TLorentzVector.h;vector");

  // Setup style
  cout << "Setting tdr style."  << endl;
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  setTDRStyle(tdrStyle);
  tdrStyle->cd();

  gStyle->SetHatchesLineWidth(1);

  gStyle->SetHatchesSpacing(1);
  //  gStyle->SetLineScalePS(4);

  
  TFile* f_lostlep = new TFile("bg_hists/lostlep_hists_fine.root", "read");
  TFile* f_hadtau = new TFile("bg_hists/hadtau_hists_fine.root", "read");
  TFile* f_qcd = new TFile("bg_hists/qcd_hists_fine.root", "read");
  TFile* f_znn = new TFile("bg_hists/znn_hists_fine.root", "read");
  TFile* f_data_obs = new TFile("data_hists/data_hists.root", "read");
  TFile* f_sig = new TFile("signal_hists/fastsim_signal.root", "read");

 
  TH1D* hdata_obs_mht_3b = (TH1D*) f_data_obs->Get("hObsMHT_3b");
  TH1D* hlostlep_mht_3b = (TH1D*) f_lostlep->Get("hPredMHT_3b");
  TH1D* hhadtau_mht_3b = (TH1D*) f_hadtau->Get("hPredMHT_3b");
  TH1D* hqcd_mht_3b = (TH1D*) f_qcd->Get("hPredMHT_3b");
  TH1D* hznn_mht_3b = (TH1D*) f_znn->Get("hPredMHT_3b");
  TH1D* ht1bbbb_1500_100_mht_3b = (TH1D*) f_sig->Get("hPredMHT_nj-1-5_nb-4-4_RA2bin_T1bbbb_1400_50_fast");
  TH1D* ht1bbbb_1000_900_mht_3b = (TH1D*) f_sig->Get("hPredMHT_nj-1-5_nb-4-4_RA2bin_T1bbbb_900_700_fast");
  TGraphAsymmErrors* gbg_mht_3b = GetBGErr("MHT_3b", hdata_obs_mht_3b, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_mht_3b = GetGDataObs(hdata_obs_mht_3b);
  MakePlot("T1bbbb-projection", gdata_obs_mht_3b, gbg_mht_3b, hdata_obs_mht_3b, hlostlep_mht_3b, hhadtau_mht_3b, hqcd_mht_3b, hznn_mht_3b,
	   ht1bbbb_1500_100_mht_3b, ht1bbbb_1000_900_mht_3b,
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)",
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 900 GeV, m_{#tilde{#chi}_{1}^{0}} = 700 GeV)",
	   false, 70, "N_{b-jet} #geq 3");
 
  TH1D* hdata_obs_mht_9j_2b = (TH1D*) f_data_obs->Get("hObsMHT_9j_2b");
  TH1D* hlostlep_mht_9j_2b = (TH1D*) f_lostlep->Get("hPredMHT_9j_2b");
  TH1D* hhadtau_mht_9j_2b = (TH1D*) f_hadtau->Get("hPredMHT_9j_2b");
  TH1D* hqcd_mht_9j_2b = (TH1D*) f_qcd->Get("hPredMHT_9j_2b");
  TH1D* hznn_mht_9j_2b = (TH1D*) f_znn->Get("hPredMHT_9j_2b");
  TH1D* ht1tttt_1500_100_mht_9j_2b = (TH1D*) f_sig->Get("hPredMHT_nj-5-5_nb-3-4_RA2bin_T1tttt_1350_50_fast");
  TH1D* ht1tttt_1200_800_mht_9j_2b = (TH1D*) f_sig->Get("hPredMHT_nj-5-5_nb-3-4_RA2bin_T1tttt_1025_600_fast");
  TGraphAsymmErrors* gbg_mht_9j_2b = GetBGErr("MHT_9j_2b", hdata_obs_mht_9j_2b, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_mht_9j_2b = GetGDataObs(hdata_obs_mht_9j_2b);
  MakePlot("T1tttt-projection", gdata_obs_mht_9j_2b, gbg_mht_9j_2b, hdata_obs_mht_9j_2b, hlostlep_mht_9j_2b, hhadtau_mht_9j_2b, hqcd_mht_9j_2b, hznn_mht_9j_2b,
  	   ht1tttt_1500_100_mht_9j_2b,  ht1tttt_1200_800_mht_9j_2b,
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1350 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)",
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1025 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)",
  	   false, 30, "N_{jet} #geq 9", "N_{b-jet} #geq 2");

  TH1D* hdata_obs_ht_MHT500 = (TH1D*) f_data_obs->Get("hObsHT_MHT500");
  TH1D* hlostlep_ht_MHT500 = (TH1D*) f_lostlep->Get("hPredHT_MHT500");
  TH1D* hhadtau_ht_MHT500 = (TH1D*) f_hadtau->Get("hPredHT_MHT500");
  TH1D* hqcd_ht_MHT500 = (TH1D*) f_qcd->Get("hPredHT_MHT500");
  TH1D* hznn_ht_MHT500 = (TH1D*) f_znn->Get("hPredHT_MHT500");
  TH1D* ht1qqqq_1400_100_ht_MHT500 = (TH1D*) f_sig->Get("hPredHT_nj-3-5_nb-1-1_mht-3-4_RA2bin_T1qqqq_1300_250_fast");
  TH1D* ht1qqqq_1000_800_ht_MHT500 = (TH1D*) f_sig->Get("hPredHT_nj-3-5_nb-1-1_mht-3-4_RA2bin_T1qqqq_750_600_fast");
  TGraphAsymmErrors* gbg_ht_MHT500 = GetBGErr("HT_MHT500", hdata_obs_ht_MHT500, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_ht_MHT500 = GetGDataObs(hdata_obs_ht_MHT500);
  MakePlot("T1qqqq-projection", gdata_obs_ht_MHT500, gbg_ht_MHT500, hdata_obs_ht_MHT500, hlostlep_ht_MHT500, hhadtau_ht_MHT500, hqcd_ht_MHT500, hznn_ht_MHT500,
  	   ht1qqqq_1400_100_ht_MHT500, ht1qqqq_1000_800_ht_MHT500,
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 250 GeV)",
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q} #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 750 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)",
  	   false, 60, "N_{jet} #geq 6", "N_{b-jet} = 0", "H_{T}^{miss} > 500 GeV");
  
  TH1D* hdata_obs_ht_MHT500_0b_7j = (TH1D*) f_data_obs->Get("hObsHT_MHT500_0b_7j");
  TH1D* hlostlep_ht_MHT500_0b_7j = (TH1D*) f_lostlep->Get("hPredHT_MHT500_0b_7j");
  TH1D* hhadtau_ht_MHT500_0b_7j = (TH1D*) f_hadtau->Get("hPredHT_MHT500_0b_7j");
  TH1D* hqcd_ht_MHT500_0b_7j = (TH1D*) f_qcd->Get("hPredHT_MHT500_0b_7j");
  TH1D* hznn_ht_MHT500_0b_7j = (TH1D*) f_znn->Get("hPredHT_MHT500_0b_7j");
  TH1D* ht5qqqqVV_1300_50_ht_MHT500_0b_7j = (TH1D*) f_sig->Get("hPredHT_nj-4-5_nb-1-1_mht-3-4_RA2bin_T5qqqqVV_1300_50_fast");
  TH1D* ht5qqqqVV_1000_800_ht_MHT500_0b_7j = (TH1D*) f_sig->Get("hPredHT_nj-4-5_nb-1-1_mht-3-4_RA2bin_T5qqqqVV_750_600_fast");
  TGraphAsymmErrors* gbg_ht_MHT500_0b_7j = GetBGErr("HT_MHT500_0b_7j", hdata_obs_ht_MHT500_0b_7j, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_ht_MHT500_0b_7j = GetGDataObs(hdata_obs_ht_MHT500_0b_7j);
  MakePlot("T5qqqqVV-projection", gdata_obs_ht_MHT500_0b_7j, gbg_ht_MHT500_0b_7j, hdata_obs_ht_MHT500_0b_7j,
	   hlostlep_ht_MHT500_0b_7j, hhadtau_ht_MHT500_0b_7j, hqcd_ht_MHT500_0b_7j, hznn_ht_MHT500_0b_7j,
	   ht5qqqqVV_1300_50_ht_MHT500_0b_7j, ht5qqqqVV_1000_800_ht_MHT500_0b_7j,
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)",
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 750 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)",
	   false, 35, "N_{jet} #geq 7", "N_{b-jet} = 0", "H_{T}^{miss} > 500 GeV"); 

  TH1D* hdata_obs_ht_7j_0b = (TH1D*) f_data_obs->Get("hObsHT_7j_0b");
  TH1D* hlostlep_ht_7j_0b = (TH1D*) f_lostlep->Get("hPredHT_7j_0b");
  TH1D* hhadtau_ht_7j_0b = (TH1D*) f_hadtau->Get("hPredHT_7j_0b");
  TH1D* hqcd_ht_7j_0b = (TH1D*) f_qcd->Get("hPredHT_7j_0b");
  TH1D* hznn_ht_7j_0b = (TH1D*) f_znn->Get("hPredHT_7j_0b");
  TH1D* ht5qqqqVV_1300_50_ht_7j_0b = (TH1D*) f_sig->Get("hPredHT_nj-4-5_nb-1-1_mht-2-4_RA2bin_T5qqqqVV_1250_300_fast");
  TH1D* ht5qqqqVV_1000_800_ht_7j_0b = (TH1D*) f_sig->Get("hPredHT_nj-4-5_nb-1-1_mht-2-4_RA2bin_T5qqqqVV_750_600_fast");
  TGraphAsymmErrors* gbg_ht_7j_0b = GetBGErr("HT_7j_0b", hdata_obs_ht_7j_0b, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_ht_7j_0b = GetGDataObs(hdata_obs_ht_7j_0b);
  MakePlot("T5qqqqVV-projection-v2", gdata_obs_ht_7j_0b, gbg_ht_7j_0b, hdata_obs_ht_7j_0b, hlostlep_ht_7j_0b, hhadtau_ht_7j_0b, hqcd_ht_7j_0b, hznn_ht_7j_0b,
  	   ht5qqqqVV_1300_50_ht_7j_0b, ht5qqqqVV_1000_800_ht_7j_0b,
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)",
  	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 750 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)",
  	   false, 80, "N_{jet} #geq 7", "N_{b-jet} = 0", "H_{T}^{miss} > 300 GeV");

  TH1D* hdata_obs_NJets_0b = (TH1D*) f_data_obs->Get("hObsNJets_0b");
  TH1D* hlostlep_NJets_0b = (TH1D*) f_lostlep->Get("hPredNJets_0b");
  TH1D* hhadtau_NJets_0b = (TH1D*) f_hadtau->Get("hPredNJets_0b");
  TH1D* hqcd_NJets_0b = (TH1D*) f_qcd->Get("hPredNJets_0b");
  TH1D* hznn_NJets_0b = (TH1D*) f_znn->Get("hPredNJets_0b");
  TH1D* ht5qqqqVV_1300_50_NJets_0b = (TH1D*) f_sig->Get("hPredNJets_nb-1-1_box-7-11_RA2bin_T5qqqqVV_1300_50_fast");
  TH1D* ht5qqqqVV_1000_800_NJets_0b = (TH1D*) f_sig->Get("hPredNJets_nb-1-1_box-7-11_RA2bin_T5qqqqVV_750_600_fast");
  TGraphAsymmErrors* gbg_NJets_0b = GetBGErr("NJets_0b", hdata_obs_NJets_0b, f_lostlep, f_hadtau, f_qcd, f_znn);
  TGraphAsymmErrors* gdata_obs_NJets_0b = GetGDataObs(hdata_obs_NJets_0b);
  MakePlot("T5qqqqVV-projection-v3", gdata_obs_NJets_0b, gbg_NJets_0b, hdata_obs_NJets_0b,
	   hlostlep_NJets_0b, hhadtau_NJets_0b, hqcd_NJets_0b, hznn_NJets_0b,
	   ht5qqqqVV_1300_50_NJets_0b, ht5qqqqVV_1000_800_NJets_0b,
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 1300 GeV, m_{#tilde{#chi}_{1}^{0}} = 50 GeV)",
	   "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q}V #tilde{#chi}_{1}^{0} (m_{#tilde{g}} = 750 GeV, m_{#tilde{#chi}_{1}^{0}} = 600 GeV)",
	   false, 200, "N_{b-jet} = 0", "H_{T}^{miss} > 500 GeV"); 
  
  return;
  
}

