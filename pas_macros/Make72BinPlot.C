
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
#include "Math/QuantFuncMathCore.h"
#include "inc/jack_style.h"
#include "inc/CMS_lumi.C"

using namespace std;

TString plotdir = "plots/";
//const double predSF = 1.63947463;

const double logmin=0.09;

void MakePlot(TString plot_title, TGraphAsymmErrors* gdata_obs, TGraphAsymmErrors* gerr, TH1D* hdata_obs, TH1D* hlostlep, TH1D* hhadtau, TH1D* hqcd, TH1D* hznn, TFile* outfile, bool logy=false)
{

  //  cout << "Set Sumw2()" << endl;

    
  hdata_obs->Sumw2();
  gStyle->SetEndErrorSize(0);
  //gStyle->SetErrorX(0);

  set_style_lite(hlostlep, "lost_lep");
  set_style_lite(hhadtau, "had_tau");
  set_style_lite(hqcd, "qcd");
  set_style_lite(hznn, "znn");
  set_style(hdata_obs, "data_obs");
    

  //  cout << "Sum up the BGs" << endl;
  TH1D * hbg_pred = (TH1D*)hdata_obs->Clone("bg_pred");
  hbg_pred->Reset();
  hbg_pred->SetTitle("");
  hbg_pred->GetYaxis()->SetTitle("Events");

  
  for (int bin (0); bin<hbg_pred->GetNbinsX(); bin++) {
    hbg_pred->SetBinContent(bin+1, hlostlep->GetBinContent(bin+1)+hhadtau->GetBinContent(bin+1));
  }
  hbg_pred->Add(hqcd);
  hbg_pred->Add(hznn);


  
  hbg_pred->SetMarkerSize(0);
  hbg_pred->SetMarkerColor(0);
  hbg_pred->SetLineWidth(0);
  hbg_pred->SetLineColor(0);


 
  THStack * hs = new THStack("hs", "");
  hs->Add(hqcd);
  hs->Add(hznn);
  hs->Add(hhadtau);
  hs->Add(hlostlep);


  gdata_obs->SetMarkerSize(1);
  gdata_obs->SetLineWidth(1);
  gdata_obs->SetMarkerStyle(20);
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
    ratiogerr->SetPoint(i, i+1, 0.);
    if (hbg_pred->GetBinContent(i+1)>0) {
      ratiogerr->SetPointError(i, 0.5, 0.5, ratiogerr->GetErrorYlow(i)/hbg_pred->GetBinContent(i+1), ratiogerr->GetErrorYhigh(i)/hbg_pred->GetBinContent(i+1));
      if (hdata_obs->GetBinContent(i+1)>0) {
	ratioderr->SetPoint(i, i+1, (hdata_obs->GetBinContent(i+1)-hbg_pred->GetBinContent(i+1))/hbg_pred->GetBinContent(i+1));
	ratioderr->SetPointError(i, 0.5, 0.5, gdata_obs->GetErrorYlow(i)/hbg_pred->GetBinContent(i+1), gdata_obs->GetErrorYhigh(i)/hbg_pred->GetBinContent(i+1));
      }
      else {
	ratioderr->SetPoint(i, i+1, -999.);
      }
    }
    else {
      ratiogerr->SetPointError(i, 0.5, 0.5, 0., 5.);
      ratioderr->SetPoint(i, i+1, -999.);
    }
  }
  set_style(ratio, "data_obs");
  ratioderr->SetMarkerSize(1);
  ratioderr->SetLineWidth(1);
  ratioderr->SetMarkerStyle(20);
  ratioderr->SetLineColor(1);
  hratiogerr->SetStats(0);
  ratio->SetTitle(hbg_pred->GetTitle());
  ratio->GetYaxis()->SetTitle("#frac{(Obs.-Exp.)}{Exp.}");
  ratio->SetMaximum(2.3);
  ratio->SetMinimum(-2.3);
  ratio->SetMarkerSize(0);
  hratiogerr->SetFillColor(0);
  hratiogerr->SetFillStyle(0);

  ratiogerr->SetFillColor(14);
  ratiogerr->SetMarkerSize(0);
  ratiogerr->SetLineWidth(0);
  ratiogerr->SetLineColor(0);
  ratiogerr->SetFillStyle(3445);
      
  ratio->GetXaxis()->SetLabelSize(0.15);
  ratio->GetXaxis()->SetLabelOffset(0.03);
  ratio->GetXaxis()->SetTitleSize(0.18);
  ratio->GetXaxis()->SetTitleOffset(0.95);
  ratio->GetYaxis()->SetLabelSize(0.10);
  ratio->GetYaxis()->SetTitleSize(0.16);
  ratio->GetYaxis()->SetTitleOffset(0.4);
  ratio->GetYaxis()->SetNdivisions(505);
  ratiounity->SetLineStyle(2);
 
  
  // // Setup legends
  TLegend * leg1 = new TLegend(0.7, 0.45, 0.945, 0.77);
  //if (wide) set_style(leg1,0.045);
  set_style(leg1,0.035);
  leg1->AddEntry(hdata_obs, "Data", "pes");
  leg1->AddEntry(hlostlep, "#splitline{Lost}{lepton}", "f");
  leg1->AddEntry(hhadtau, "#splitline{Hadronic}{#tau lepton}", "f");
  TLegend * leg2 = new TLegend(0.855, 0.45, 1.1, 0.77);
  set_style(leg2,0.035);
  leg2->AddEntry(hbg_pred, "", "f");
  leg2->AddEntry(hznn, "Z#rightarrow#nu#bar{#nu}", "f");
  leg2->AddEntry(hqcd, "QCD", "f");

 
  double ymax = hbg_pred->GetMaximum();
  if (hdata_obs->GetMaximum()>ymax) ymax=hdata_obs->GetMaximum();

  if(logy) {
    hbg_pred->SetMaximum(500*ymax);
    hbg_pred->SetMinimum(logmin);
  }
  else {
    hbg_pred->SetMinimum(0);
    hbg_pred->SetMaximum(10);
  }
 
  // Setup canvas and pads

  int W = 800;
  int H = 600;

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
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
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
  double font_size_dw  = 0.1;  // please tune the font size parameter for bottom figure
  double dw_height     = (1. - up_height) * dw_correction;
  double dw_height_offset = 0.02; // KH, added to put the bottom one closer to the top panel

  
  TPad * pad1 = new TPad("pad1", "top pad" , 0.0, 0.3, 1.0, 1.0);
  TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
  
  pad1->SetTickx(0);
  pad1->SetTicky(0);
  pad1->SetPad(0., 1 - up_height,    1., 1.00);
  //
  pad1->SetFrameFillColor(0);
  pad1->SetFillColor(0);
  pad1->SetTopMargin(0.12);
  pad1->SetLeftMargin(0.1);
  pad1->Draw();

  pad2->SetPad(0., 0., 1., dw_height+dw_height_offset);
  pad2->SetFillColor(0);
  pad2->SetFrameFillColor(0);
  pad2->SetBottomMargin(0.35);
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.1);
  pad2->Draw();
  pad1->cd();
  pad1->SetLogy(logy);
 
  // // Draw hists

  hbg_pred->Draw();

  hbg_pred->GetYaxis()->SetLabelSize(0.035*1.18);
  hbg_pred->GetYaxis()->SetTitleSize(0.045*1.18);
  hbg_pred->GetYaxis()->SetTitleOffset(0.75);
  hbg_pred->GetYaxis()->SetTitleFont(42);
  
  hs->Draw("hist,same");

  float ymax_top = hbg_pred->GetMaximum();
  float ymin_top = logmin;

  float ymax2_top = 1000.;
  float ymax3_top = 200.;
  float ymax4_top = 30.;

  float ymax_bottom = 1.99;
  float ymin_bottom = 0.01;

  float ymax2_bottom = 2.15;
  float ymax3_bottom = 2.15;
  float ymax4_bottom = 2.15;

  
  // Njet separation lines
  TLine *tl_njet = new TLine();
  tl_njet->SetLineStyle(2);
  tl_njet->DrawLine(25.-0.5,ymin_top,25.-0.5,hbg_pred->GetMaximum()); 
  tl_njet->DrawLine(49.-0.5,ymin_top,49.-0.5,hbg_pred->GetMaximum()); 

  // Njet labels
  TLatex * ttext_njet = new TLatex();
  ttext_njet->SetTextFont(42);
  ttext_njet->SetTextSize(0.060);
  ttext_njet->SetTextAlign(22);
  ttext_njet->DrawLatex(13.-0.5 , ymax_top/4. , "4 #leq N_{#scale[0.2]{ }jet} #leq 6");
  ttext_njet->DrawLatex(37.-0.5 , ymax_top/4. , "7 #leq N_{#scale[0.2]{ }jet} #leq 8");
  ttext_njet->DrawLatex(61.-0.5 , ymax_top/4. , "N_{#scale[0.2]{ }jet} #geq 9");

  // Nb separation lines
  TLine *tl_nb = new TLine();
  tl_nb->SetLineStyle(3);
  tl_nb->DrawLine( 7.-0.5,ymin_top, 7.-0.5,ymax2_top); 
  tl_nb->DrawLine(13.-0.5,ymin_top,13.-0.5,ymax2_top); 
  tl_nb->DrawLine(19.-0.5,ymin_top,19.-0.5,ymax2_top); 
  tl_nb->DrawLine(31.-0.5,ymin_top,31.-0.5,ymax3_top); 
  tl_nb->DrawLine(37.-0.5,ymin_top,37.-0.5,ymax3_top); 
  tl_nb->DrawLine(43.-0.5,ymin_top,43.-0.5,ymax3_top); 
  tl_nb->DrawLine(55.-0.5,ymin_top,55.-0.5,ymax4_top); 
  tl_nb->DrawLine(61.-0.5,ymin_top,61.-0.5,ymax4_top); 
  tl_nb->DrawLine(67.-0.5,ymin_top,67.-0.5,ymax4_top); 
    
  // Nb labels
  TLatex * ttext_nb = new TLatex();
  ttext_nb->SetTextFont(42);
  ttext_nb->SetTextSize(0.050);
  ttext_nb->SetTextAlign(22);
    
  ttext_nb->DrawLatex( 4.75 , ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
  ttext_nb->DrawLatex( 4.-0.5 , ymax_top/40. , "0");
  ttext_nb->DrawLatex(10.-0.5 , ymax_top/40. , "1");
  ttext_nb->DrawLatex(16.-0.5 , ymax_top/40. , "2");
  ttext_nb->DrawLatex(22.-0.5 , ymax_top/40. , "#geq 3");

  

  hdata_obs->SetMarkerStyle(20);
  set_style(hdata_obs, "data_obs");
  gdata_obs->SetMarkerStyle(20);
  hbg_pred->GetXaxis()->SetLabelSize(0);
  
  gerr->Draw("2 same");
  gdata_obs->Draw("p same");
  
  // // Draw legends
  leg1->Draw();
  leg2->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.052);
  latex->SetTextSize(0.04);
  latex->SetTextSize(0.052);

  // Luminosity information for scaling
  double lumi     = 2.262; // normaliza to this lumi (fb-1)
  double lumi_ref = 2.262; // normaliza to 3 (fb-1)

  char tempname[200];
  TString line = "";
  sprintf(tempname,"%8.1f",lumi);
  line+=tempname;
  line+=" fb^{-1} (13 TeV)";
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=0;
    
  writeExtraText = false;
  extraText   = "       Preliminary";
  TString lumi_sqrtS = line;
  
  TPaveText * pave = new TPaveText(0.18, 0.86, 0.4, 0.96, "brNDC");
  //  TText * text = NULL; 
  TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
  
  pad2->cd();
  pad2->SetGridy(0);




  set_style(ratio, "data_obs");
  ratio->Draw("axis");

  ratio->GetXaxis()->SetLabelSize(font_size_dw);
  ratio->GetXaxis()->SetTitleSize(font_size_dw);
  ratio->GetYaxis()->SetLabelSize(font_size_dw);
  ratio->GetYaxis()->SetTitleSize(font_size_dw);

  ratio->GetXaxis()->SetLabelSize(0.14);
  ratio->GetXaxis()->SetTitleSize(0.14);
  ratio->GetXaxis()->SetTitleOffset(1.1);
  ratio->GetXaxis()->SetTitleFont(42);
  ratio->GetYaxis()->SetLabelSize(0.13);
  ratio->GetYaxis()->SetTitleSize(0.13);
  ratio->GetYaxis()->SetTitleOffset(0.32);
  ratio->GetYaxis()->SetTitleFont(42);
  ratio->GetXaxis()->SetTitle("Search region bin number");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTickLength(0.015);
  ratio->GetXaxis()->SetTickLength(0.08);
  hbg_pred->GetXaxis()->SetTitleSize(0.14);

  ratiounity->Draw();
  ratiogerr->Draw("e2 same");
  ratioderr->Draw("p same");

    // tl_njet->DrawLine(25.-0.5,ymin_bottom,25.-0.5,hbg_pred->GetMaximum()); 
  // tl_njet->DrawLine(49.-0.5,ymin_bottom,49.-0.5,hbg_pred->GetMaximum());
  tl_nb->DrawLine( 7.-0.5,-2.3, 7.-0.5,2.3); 
  tl_nb->DrawLine(13.-0.5,-2.3,13.-0.5,2.3); 
  tl_nb->DrawLine(19.-0.5,-2.3,19.-0.5,2.3); 
  tl_njet->DrawLine(25.-0.5,-2.3,25.-0.5,2.3); 
  tl_nb->DrawLine(31.-0.5,-2.3,31.-0.5,2.3); 
  tl_nb->DrawLine(37.-0.5,-2.3,37.-0.5,2.3); 
  tl_nb->DrawLine(43.-0.5,-2.3,43.-0.5,2.3); 
  tl_njet->DrawLine(49.-0.5,-2.3,49.-0.5,2.3); 
  tl_nb->DrawLine(55.-0.5,-2.3,55.-0.5,2.3); 
  tl_nb->DrawLine(61.-0.5,-2.3,61.-0.5,2.3); 
  tl_nb->DrawLine(67.-0.5,-2.3,67.-0.5,2.3);

  ratioderr->Print("all");
  

  pave->SetLineColor(0);
  pave->SetLineWidth(0);
  pave->SetFillStyle(4000);
  pave->SetShadowColor(0);
  pave->SetBorderSize(1);
  
  set_style(ratioleg);
  ratioleg->SetTextSize(0.07);
  ratioleg->AddEntry(gerr, "Pred. uncert. (stat#oplussyst)", "f");  
  pad1->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();
  pad2->cd();
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();

  
  canv->cd();
  CMS_lumi(canv, iPeriod, iPos, lumi_sqrtS);

 
  gPad->Print(plotdir+plot_title+".pdf");
  gPad->Print(plotdir+plot_title+".png");

  outfile->cd();
  canv->Write();

  // Clean up
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
  delete leg1;
  delete leg2;
  delete latex;
  delete pave;
  delete canv;

  cout << "SaveHist(): DONE!" << endl;


  return;
}

void Make72BinPlot() {

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
  gStyle->SetHatchesSpacing(1.);


  TFile* f_lostlep = new TFile("bg_hists/lostlep_hists.root", "read");
  TFile* f_hadtau = new TFile("bg_hists/hadtau_hists.root", "read");
  TFile* f_qcd = new TFile("bg_hists/qcd_hists.root", "read");
  TFile* f_znn = new TFile("bg_hists/znn_hists.root", "read");
  TFile* f_data_obs = new TFile("data_hists/data_hists.root", "read");

  TH1D* hdata_obs = (TH1D*) f_data_obs->Get("hObsAllBins");
  TH1D* hqcd = (TH1D*) f_qcd->Get("hPredAllBins");
  TH1D* hlostlep = (TH1D*) f_lostlep->Get("hPredAllBins");
  TH1D* hhadtau = (TH1D*) f_hadtau->Get("hPredAllBins");
  TH1D* hznn = (TH1D*) f_znn->Get("hPredAllBins");
  
  hlostlep->Sumw2();
  hhadtau->Sumw2();
  hqcd->Sumw2();
  hznn->Sumw2();



  TGraphAsymmErrors* glostlepstat = (TGraphAsymmErrors*) f_lostlep->Get("Graph;2");
  TGraphAsymmErrors* ghadtaustat = (TGraphAsymmErrors*) f_hadtau->Get("Graph;2");
  TGraphAsymmErrors* glostlepsyst = (TGraphAsymmErrors*) f_lostlep->Get("Graph;3");
  TGraphAsymmErrors* ghadtausyst = (TGraphAsymmErrors*) f_hadtau->Get("Graph;3");
  TGraphAsymmErrors* gqcdstat = (TGraphAsymmErrors*) f_qcd->Get("Graph;2");
  TGraphAsymmErrors* gznnstat = (TGraphAsymmErrors*) f_znn->Get("Graph;2");
  TGraphAsymmErrors* gqcdsyst = (TGraphAsymmErrors*) f_qcd->Get("Graph;3");
  TGraphAsymmErrors* gznnsyst = (TGraphAsymmErrors*) f_znn->Get("Graph;3");


  const double alpha = 1 - 0.6827;
  
  Double_t x[72];
  Double_t xl[72];
  Double_t xh[72];
  Double_t xld[72];
  Double_t xhd[72];
  
  Double_t data_cv[72];
  Double_t data_pois_up[72];
  Double_t data_pois_down[72];
  
  Double_t pred_cv[72];
  Double_t full_stat_up[72];
  Double_t full_stat_down[72];
  Double_t full_syst_up[72];
  Double_t full_syst_down[72];
  Double_t full_err_up[72];
  Double_t full_err_down[72];

  for (unsigned int bin(0); bin<72; bin++) {
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    xld[bin]=0.1;
    xhd[bin]=0.1;

    data_cv[bin]=hdata_obs->GetBinContent(bin+1);
    double N=data_cv[bin];
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? 0  : ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    data_pois_up[bin]=(U-N);
    data_pois_down[bin]=(N-L);
    printf("Bin %d: Nobs=%d, estat=+%3.3f-%3.3f\n", bin+1, (int)hdata_obs->GetBinContent(bin+1), data_pois_up[bin], data_pois_down[bin]);
  
    pred_cv[bin]=glostlepstat->Eval(bin+1)+ghadtaustat->Eval(bin+1)+gqcdstat->Eval(bin+1)+gznnstat->Eval(bin+1);
    if (pred_cv[bin]<logmin&&pred_cv[bin]>0) pred_cv[bin]=logmin+0.005; 
    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    full_stat_up[bin] = sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.));
    full_stat_down[bin] = sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.));
    full_syst_up[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.));
    full_syst_down[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.));
    full_err_up[bin] = sqrt(pow(full_stat_up[bin], 2.)+pow(full_syst_up[bin], 2.));
    full_err_down[bin] = sqrt(pow(full_stat_down[bin], 2.)+pow(full_syst_down[bin], 2.));
  }
  TGraphAsymmErrors* gdata_obs = new TGraphAsymmErrors(72, x, data_cv, xld, xhd, data_pois_down, data_pois_up);
  TGraphAsymmErrors* gbg = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, full_err_down, full_err_up);
  
  TFile* outfile = new TFile("plots/summary_plots.root", "recreate");
  cout << "Saving output to " << outfile->GetName() << endl;
  outfile->cd();


  cout << "Make plots..." << endl;
  MakePlot("results-plot-prefit", gdata_obs, gbg, hdata_obs, hlostlep, hhadtau, hqcd, hznn, outfile, true);


  outfile->Close();

  
  return;
  
}

