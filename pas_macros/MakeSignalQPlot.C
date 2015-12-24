
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

double GetQ(double S, double B) {
  return 2*(sqrt(S+B)-sqrt(B));
}

TFile* outfile;


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
    full_stat_up[bin] = sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.));
    full_stat_down[bin] = sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.));
    full_syst_up[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.));
    full_syst_down[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.));
    full_err_up[bin] = sqrt(pow(full_stat_up[bin], 2.)+pow(full_syst_up[bin], 2.));
    full_err_down[bin] = sqrt(pow(full_stat_down[bin], 2.)+pow(full_syst_down[bin], 2.));
  }

  TGraphAsymmErrors* gBGErr = new TGraphAsymmErrors(nbins, x, pred_cv, xl, xh, full_err_down, full_err_up);

  return gBGErr;

}

void MakePlot(TString plot_title, TGraphAsymmErrors* gerr, TH1D* hlostlep, TH1D* hhadtau, TH1D* hqcd, TH1D* hznn, 
	      TH1D* ht1tttt_1500_100, TH1D* ht1tttt_1200_800, TH1D* ht1bbbb_1500_100, TH1D* ht1bbbb_1000_900, TH1D* ht1qqqq_1400_100, TH1D* ht1qqqq_1000_800,
	      bool logy=false)
{


    
  gStyle->SetEndErrorSize(0);
  gerr->SetFillColor(14);
  gerr->SetMarkerSize(0);
  gerr->SetLineWidth(0);
  gerr->SetLineColor(0);
  gerr->SetFillStyle(3445);
  


    

  //  cout << "Sum up the BGs" << endl;
  TH1D * hbg_pred = (TH1D*)hlostlep->Clone("bg_pred");
  hbg_pred->Reset();
  hbg_pred->SetTitle("");
  hbg_pred->GetYaxis()->SetTitle("Events");

  TH1D* htemp = (TH1D*)hlostlep->Clone("temp");
  htemp->Reset();
  htemp->SetTitle("");
  htemp->GetXaxis()->SetTitle("Search region bin number");
  
  for (int bin (0); bin<hbg_pred->GetNbinsX(); bin++) {
    hbg_pred->SetBinContent(bin+1, hlostlep->GetBinContent(bin+1)+hhadtau->GetBinContent(bin+1));
  }


  hbg_pred->Add(hqcd);
  hbg_pred->Add(hznn);

  set_style(hbg_pred, "lost_lep");
  hbg_pred->SetFillColor(3002);


  set_style(ht1tttt_1500_100,"data_obs");
  ht1tttt_1500_100->SetLineColor(5000);
  ht1tttt_1500_100->SetMarkerColor(5000);
  ht1tttt_1500_100->SetMarkerStyle(20);
  set_style(ht1tttt_1200_800,"data_obs");
  ht1tttt_1200_800->SetLineColor(5000);
  ht1tttt_1200_800->SetMarkerColor(5000);
  ht1tttt_1200_800->SetMarkerStyle(24);

  set_style(ht1bbbb_1500_100,"data_obs");
  ht1bbbb_1500_100->SetLineColor(5002);
  ht1bbbb_1500_100->SetMarkerColor(5002);
  ht1bbbb_1500_100->SetMarkerStyle(21);
  set_style(ht1bbbb_1000_900,"data_obs");
  ht1bbbb_1000_900->SetLineColor(5002);
  ht1bbbb_1000_900->SetMarkerColor(5002);
  ht1bbbb_1000_900->SetMarkerStyle(25);

  set_style(ht1qqqq_1400_100,"data_obs");
  ht1qqqq_1400_100->SetLineColor(5004);
  ht1qqqq_1400_100->SetMarkerColor(5004);
  ht1qqqq_1400_100->SetMarkerStyle(22);
  set_style(ht1qqqq_1000_800,"data_obs");
  ht1qqqq_1000_800->SetLineColor(5004);
  ht1qqqq_1000_800->SetMarkerColor(5004);
  ht1qqqq_1000_800->SetMarkerStyle(26);

  ht1tttt_1500_100->Scale(2.153738/3);;
  ht1tttt_1200_800->Scale(2.153738/3);;
  ht1bbbb_1500_100->Scale(2.153738/3);;
  ht1bbbb_1000_900->Scale(2.153738/3);;
  ht1qqqq_1400_100->Scale(2.153738/3);;
  ht1qqqq_1000_800->Scale(2.153738/3);;

  TH1D * q_t1tttt_1500_100 = (TH1D *) ht1tttt_1500_100->Clone("q_t1tttt_1500_100");
  TH1D * q_t1tttt_1200_800 = (TH1D *) ht1tttt_1200_800->Clone("q_t1tttt_1200_800");
  TH1D * q_t1bbbb_1500_100 = (TH1D *) ht1bbbb_1500_100->Clone("q_t1bbbb_1500_100");
  TH1D * q_t1bbbb_1000_900 = (TH1D *) ht1bbbb_1000_900->Clone("q_t1bbbb_1000_900");
  TH1D * q_t1qqqq_1400_100 = (TH1D *) ht1qqqq_1400_100->Clone("q_t1qqqq_1400_100");
  TH1D * q_t1qqqq_1000_800 = (TH1D *) ht1qqqq_1000_800->Clone("q_t1qqqq_1000_800");

  set_style(q_t1tttt_1500_100,"data_obs");
  q_t1tttt_1500_100->SetLineColor(5000);
  q_t1tttt_1500_100->SetMarkerColor(5000);
  q_t1tttt_1500_100->SetMarkerStyle(20);
  set_style(q_t1tttt_1200_800,"data_obs");
  q_t1tttt_1200_800->SetLineColor(5000);
  q_t1tttt_1200_800->SetMarkerColor(5000);
  q_t1tttt_1200_800->SetMarkerStyle(24);

  set_style(q_t1bbbb_1500_100,"data_obs");
  q_t1bbbb_1500_100->SetLineColor(5002);
  q_t1bbbb_1500_100->SetMarkerColor(5002);
  q_t1bbbb_1500_100->SetMarkerStyle(21);
  set_style(q_t1bbbb_1000_900,"data_obs");
  q_t1bbbb_1000_900->SetLineColor(5002);
  q_t1bbbb_1000_900->SetMarkerColor(5002);
  q_t1bbbb_1000_900->SetMarkerStyle(25);

  set_style(q_t1qqqq_1400_100,"data_obs");
  q_t1qqqq_1400_100->SetLineColor(5004);
  q_t1qqqq_1400_100->SetMarkerColor(5004);
  q_t1qqqq_1400_100->SetMarkerStyle(22);
  set_style(q_t1qqqq_1000_800,"data_obs");
  q_t1qqqq_1000_800->SetLineColor(5004);
  q_t1qqqq_1000_800->SetMarkerColor(5004);
  q_t1qqqq_1000_800->SetMarkerStyle(26);

  htemp->SetStats(0);
  htemp->GetYaxis()->SetTitle("Q = 2[#sqrt{S+B}-#sqrt{B}]");

  htemp->GetXaxis()->SetLabelSize(0.15);
  htemp->GetXaxis()->SetLabelOffset(0.03);
  htemp->GetXaxis()->SetTitleSize(0.14);
  htemp->GetXaxis()->SetTitleOffset(1.2);
  htemp->GetYaxis()->SetLabelSize(0.10);
  htemp->GetYaxis()->SetTitleSize(0.12);
  htemp->GetYaxis()->SetTitleOffset(0.35);
  htemp->GetYaxis()->SetNdivisions(505);
  TLine* qp1 = new TLine(hbg_pred->GetBinLowEdge(1),1,hbg_pred->GetBinLowEdge(hbg_pred->GetNbinsX()+1),1);
  TLine* qp2 = new TLine(hbg_pred->GetBinLowEdge(1),2,hbg_pred->GetBinLowEdge(hbg_pred->GetNbinsX()+1),2);
  TLine* qp3 = new TLine(hbg_pred->GetBinLowEdge(1),3,hbg_pred->GetBinLowEdge(hbg_pred->GetNbinsX()+1),3);
  qp1->SetLineStyle(2);
  qp2->SetLineStyle(2);
  qp3->SetLineStyle(2);

  for (Int_t bin = 0; bin < hbg_pred->GetNbinsX(); bin++) {
    q_t1tttt_1500_100->SetBinContent(bin+1, GetQ(ht1tttt_1500_100->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1tttt_1500_100->SetBinError(bin+1, 0);
    ht1tttt_1500_100->SetBinError(bin+1, 0);
    q_t1tttt_1200_800->SetBinContent(bin+1, GetQ(ht1tttt_1200_800->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1tttt_1200_800->SetBinError(bin+1, 0);
    ht1tttt_1200_800->SetBinError(bin+1, 0);
    q_t1bbbb_1500_100->SetBinContent(bin+1, GetQ(ht1bbbb_1500_100->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1bbbb_1500_100->SetBinError(bin+1, 0);
    ht1bbbb_1500_100->SetBinError(bin+1, 0);
    q_t1bbbb_1000_900->SetBinContent(bin+1, GetQ(ht1bbbb_1000_900->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1bbbb_1000_900->SetBinError(bin+1, 0);
    ht1bbbb_1000_900->SetBinError(bin+1, 0);
    q_t1qqqq_1400_100->SetBinContent(bin+1, GetQ(ht1qqqq_1400_100->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1qqqq_1400_100->SetBinError(bin+1, 0);
    ht1qqqq_1400_100->SetBinError(bin+1, 0);
    q_t1qqqq_1000_800->SetBinContent(bin+1, GetQ(ht1qqqq_1000_800->GetBinContent(bin+1), hbg_pred->GetBinContent(bin+1)));
    q_t1qqqq_1000_800->SetBinError(bin+1, 0);
    ht1qqqq_1000_800->SetBinError(bin+1, 0);
  }
 
  // Setup legends                                                                                                                                                                                                                                                         
  TLegend * leg1 = new TLegend(0.37, 0.5, 0.77, 0.77);
  set_style(leg1,0.025);
  leg1->AddEntry(ht1tttt_1500_100, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1500 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)}", "p");
  leg1->AddEntry(ht1bbbb_1500_100, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1500 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)}", "p");
  leg1->AddEntry(ht1qqqq_1400_100, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1400 GeV, m_{#tilde{#chi}_{1}^{0}} = 100 GeV)}", "p");

  TLegend * leg2 = new TLegend(0.7, 0.5, 0.94, 0.77);
  set_style(leg2,0.025);
  leg2->AddEntry(ht1tttt_1200_800, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1200 GeV, m_{#tilde{#chi}_{1}^{0}} = 800 GeV)}", "p");
  leg2->AddEntry(ht1bbbb_1000_900, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow b#bar{b} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 900 GeV)}", "p");
  leg2->AddEntry(ht1qqqq_1000_800, "#splitline{pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow q#bar{q} #tilde{#chi}_{1}^{0}}{(m_{#tilde{g}} = 1000 GeV, m_{#tilde{#chi}_{1}^{0}} = 800 GeV)}", "p");

  TLegend * leg3 = new TLegend(0.7, 0.4, 0.94, 0.5);
  set_style(leg3,0.035);
  leg3->AddEntry(hbg_pred, "Total BG", "f");

  
  double ymax = hbg_pred->GetMaximum();
  hbg_pred->SetMaximum(500*ymax);
  hbg_pred->SetMinimum(0.07);



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
  double dw_height_offset = 0.04; // KH, added to put the bottom one closer to the top panel

  
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

  hbg_pred->Draw("hist");
  gerr->Draw("2 same");
  hbg_pred->GetYaxis()->SetLabelSize(0.035*1.15);
  hbg_pred->GetYaxis()->SetTitleSize(0.045*1.15);
  hbg_pred->GetYaxis()->SetTitleOffset(1);
  hbg_pred->GetYaxis()->SetTitleFont(42);
  hbg_pred->GetXaxis()->SetLabelSize(0);
  cout << "Draw hists..." << endl;
  ht1tttt_1500_100->Draw("p,same");
  ht1tttt_1200_800->Draw("p,same");
  ht1bbbb_1500_100->Draw("p,same");
  ht1bbbb_1000_900->Draw("p,same");
  ht1qqqq_1400_100->Draw("p,same");
  ht1qqqq_1000_800->Draw("p,same");


  float ymax_top = hbg_pred->GetMaximum();
  float ymin_top = 0.015;

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
    
  ttext_nb->DrawLatex( 5. , ymax_top/12. , "N_{#scale[0.2]{ }b-jet}");
  ttext_nb->DrawLatex( 4.-0.5 , ymax_top/40. , "0");
  ttext_nb->DrawLatex(10.-0.5 , ymax_top/40. , "1");
  ttext_nb->DrawLatex(16.-0.5 , ymax_top/40. , "2");
  ttext_nb->DrawLatex(22.-0.5 , ymax_top/40. , "#geq 3");

  

  hbg_pred->GetXaxis()->SetLabelSize(0);
  
  // // Draw legends
  leg1->Draw();
  leg2->Draw();
  leg3->Draw();
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextFont(62);
  latex->SetTextSize(0.052);
  latex->SetTextSize(0.04);
  latex->SetTextSize(0.052);

  // Luminosity information for scaling
  double lumi     = 2.153738; // normaliza to this lumi (fb-1)
  double lumi_ref = 2.153738; // normaliza to 3 (fb-1)

  char tempname[200];
  TString line = "";
  sprintf(tempname,"%8.1f",lumi);
  line+=tempname;
  line+=" fb^{-1} (13 TeV)";
  
  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos=0;
    
  writeExtraText = true;
  extraText   = "       Preliminary";
  TString lumi_sqrtS = line;
  
  TPaveText * pave = new TPaveText(0.18, 0.86, 0.4, 0.96, "brNDC");
  //  TText * text = NULL; 
  TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
  
  pad2->cd();
  pad2->SetGridy(0);

  htemp->SetMaximum(1.85);
  htemp->SetMinimum(0);
    
  htemp->Draw("axis");
  htemp->GetXaxis()->SetTitleSize(0.12);
  htemp->GetXaxis()->SetLabelSize(0.12);
  hbg_pred->GetXaxis()->SetTitleSize(0.12);
  htemp->GetXaxis()->SetLabelSize(0.12);
  q_t1tttt_1500_100->Draw("p same");
  q_t1tttt_1200_800->Draw("p, same");
  q_t1bbbb_1500_100->Draw("p, same");
  q_t1bbbb_1000_900->Draw("p, same");
  q_t1qqqq_1400_100->Draw("p, same");
  q_t1qqqq_1000_800->Draw("p, same");
  qp1->Draw();
  qp2->Draw();
  qp3->Draw();

  // // hratiogerr->GetXaxis()->SetRangeUser(0,6);
  //  hratiogerr->Draw("e2");


  q_t1tttt_1500_100->GetXaxis()->SetLabelSize(font_size_dw);
  q_t1tttt_1500_100->GetXaxis()->SetTitleSize(font_size_dw);
  q_t1tttt_1500_100->GetYaxis()->SetLabelSize(font_size_dw);
  q_t1tttt_1500_100->GetYaxis()->SetTitleSize(font_size_dw);

  q_t1tttt_1500_100->GetXaxis()->SetTitleSize(0.12);
  q_t1tttt_1500_100->GetXaxis()->SetTitleOffset(1.1);
  q_t1tttt_1500_100->GetXaxis()->SetTitleFont(42);
  q_t1tttt_1500_100->GetYaxis()->SetTitleSize(0.13);
  q_t1tttt_1500_100->GetYaxis()->SetTitleOffset(0.32);
  q_t1tttt_1500_100->GetYaxis()->SetTitleFont(42);
  q_t1tttt_1500_100->GetXaxis()->SetTitle("Search region bin number");
  q_t1tttt_1500_100->GetYaxis()->SetNdivisions(505);
  q_t1tttt_1500_100->GetYaxis()->SetTickLength(0.015);
  q_t1tttt_1500_100->GetXaxis()->SetTickLength(0.08);
  hbg_pred->GetXaxis()->SetTitleSize(0.12);


  pave->SetLineColor(0);
  pave->SetLineWidth(0);
  pave->SetFillStyle(4000);
  pave->SetShadowColor(0);
  pave->SetBorderSize(1);
 
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
  gPad->Write();
  gerr->Write();

 
  delete hbg_pred;

  delete pad1;
  delete pad2;
  delete leg1;
  delete leg2;
  delete latex;
  delete pave;
  delete canv;

  cout << "SaveHist(): DONE!" << endl;


  return;
}

void MakeSignalQPlot() {

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

 
  TFile* f_lostlep = new TFile("bg_hists/lostlep_hists.root", "read");
  TFile* f_hadtau = new TFile("bg_hists/hadtau_hists.root", "read");
  TFile* f_qcd = new TFile("bg_hists/qcd_hists.root", "read");
  TFile* f_znn = new TFile("bg_hists/znn_hists.root", "read");
  TFile* f_sig = new TFile("signal_hists/RA2bin_signal.root", "read");


  TH1D* hqcd = (TH1D*) f_qcd->Get("hPredAllBins");
  TH1D* hlostlep = (TH1D*) f_lostlep->Get("hPredAllBins");
  TH1D* hhadtau = (TH1D*) f_hadtau->Get("hPredAllBins");
  TH1D* hznn = (TH1D*) f_znn->Get("hPredAllBins");
  
  hlostlep->Sumw2();
  hhadtau->Sumw2();
  hqcd->Sumw2();
  hznn->Sumw2();

  TH1D* ht1tttt_1500_100 = (TH1D*) f_sig->Get("RA2bin_T1tttt_1500_100_fast");
  TH1D* ht1tttt_1200_800 = (TH1D*) f_sig->Get("RA2bin_T1tttt_1200_800_fast");
  TH1D* ht1bbbb_1500_100 = (TH1D*) f_sig->Get("RA2bin_T1bbbb_1500_100_fast");
  TH1D* ht1bbbb_1000_900 = (TH1D*) f_sig->Get("RA2bin_T1bbbb_1000_900_fast");
  TH1D* ht1qqqq_1400_100 = (TH1D*) f_sig->Get("RA2bin_T1qqqq_1400_800_fast");
  TH1D* ht1qqqq_1000_800 = (TH1D*) f_sig->Get("RA2bin_T1qqqq_1000_800_fast");

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

    pred_cv[bin]=glostlepstat->Eval(bin+1)+ghadtaustat->Eval(bin+1)+gqcdstat->Eval(bin+1)+gznnstat->Eval(bin+1);
    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    full_stat_up[bin] = sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.));
    full_stat_down[bin] = sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.));
    full_syst_up[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.));
    full_syst_down[bin] = sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.));
    full_err_up[bin] = sqrt(pow(full_stat_up[bin], 2.)+pow(full_syst_up[bin], 2.));
    full_err_down[bin] = sqrt(pow(full_stat_down[bin], 2.)+pow(full_syst_down[bin], 2.));
  }
  TGraphAsymmErrors* gbg = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, full_err_down, full_err_up);
  
  outfile = new TFile("test.root","recreate");

  cout << "Make plots..." << endl;
  MakePlot("bins_all_signal_BG_data-v3", gbg, hlostlep, hhadtau, hqcd, hznn, ht1tttt_1500_100, ht1tttt_1200_800, ht1bbbb_1500_100, ht1bbbb_1000_900, ht1qqqq_1400_100, ht1qqqq_1000_800, true);

  cout << gStyle->GetHatchesSpacing() << endl;
  cout << gStyle->GetHatchesLineWidth() << endl;


  //  outfile->Close();

  
  return;
  
}

