
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TChain.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "TCut.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TExec.h"
#include "TPie.h"
#include "TFile.h"
#include <vector>


#include <iostream>
#include "inc/jack_style.h"
#include "inc/CMS_lumi.C"

using namespace std;

void DrawPieCharts()
{ 


  TH1::SetDefaultSumw2();

  TFile* f_lostlep = new TFile("bg_hists/lostlep_hists.root", "read");
  TFile* f_hadtau = new TFile("bg_hists/hadtau_hists.root", "read");
  TFile* f_qcd = new TFile("bg_hists/qcd_hists.root", "read");
  TFile* f_znn = new TFile("bg_hists/znn_hists.root", "read");
  TH1D* hlostlep = (TH1D*) f_lostlep->Get("hPredAllBins");
  TH1D* hhadtau = (TH1D*) f_hadtau->Get("hPredAllBins");
  TH1D* hqcd = (TH1D*) f_qcd->Get("hPredAllBins");
  TH1D* hznn = (TH1D*) f_znn->Get("hPredAllBins");
  
  
  vector<TPie*> bgPies;

  for (unsigned int box(0); box<12; box++) {
    cout << "Box " << box+1 << "...";
    vector<TString> labels;
    labels.push_back("Lost-e/#mu");
    labels.push_back("#tau#rightarrowhad.");
    labels.push_back("QCD");
    labels.push_back("Z+jets");
    vector<double> yields;
    yields.push_back(hlostlep->Integral(box*6+1,box*6+6));
    yields.push_back(hhadtau->Integral(box*6+1,box*6+6));
    yields.push_back(hqcd->Integral(box*6+1,box*6+6));
    yields.push_back(hznn->Integral(box*6+1,box*6+6));
    vector<int> colors;
    colors.push_back(2006);
    colors.push_back(2007);
    colors.push_back(2001);
    colors.push_back(2002);

    TPie* pie = new TPie(Form("pie%d",box+1), "", yields.size(), &yields[0]);
    pie->SetFillColors(&colors[0]);
    pie->SetLabelFormat("");
    pie->SetRadius(0.3);
    bgPies.push_back(pie);
  }
  

  TCanvas * thecanvas= new TCanvas("thecanvas","the canvas",800, 700);
  thecanvas->SetFillStyle(4000);

  set_style_lite(hlostlep, "lost_lep");
  set_style_lite(hhadtau, "had_tau");
  set_style_lite(hqcd, "qcd");
  set_style_lite(hznn, "znn");
  //  TLegend * leg1 = new TLegend(0.12, 0.815, 0.9, 0.865);
  TLegend * leg1 = new TLegend(0.12, 0.79, 0.9, 0.84);
  set_style(leg1,0.025);
  leg1->SetNColumns(4);
  leg1->AddEntry(hlostlep, "#splitline{Lost}{lepton}", "f");
  leg1->AddEntry(hhadtau, "#splitline{Hadronic #tau}{lepton}", "f");
  leg1->AddEntry(hqcd, "QCD   ", "f");
  leg1->AddEntry(hznn, "Z+jets", "f");


  //  TH2Poly * Hcomp = new TH2Poly("Hcomp","BG Composition;N_{b-jet} (p_{T} > 30 GeV);N_{jet} (p_{T} > 30 GeV)", -0.5, 3.5, 3.5, 12.5);
  TH2Poly * Hcomp = new TH2Poly("Hcomp","", 0, 4, 0, 3.5);
  Hcomp->SetStats(0);
  Hcomp->GetXaxis()->SetLabelSize(0.05);
  Hcomp->GetYaxis()->SetLabelSize(0.05);
  Hcomp->GetXaxis()->SetBinLabel(12, "0 b-jets");
  Hcomp->GetXaxis()->SetBinLabel(37, "1 b-jet");
  Hcomp->GetXaxis()->SetBinLabel(62, "2 b-jets");
  Hcomp->GetXaxis()->SetBinLabel(87, "3+ b-jets");
  Hcomp->GetYaxis()->SetBinLabel(14, "#splitline{4-6}{jets}");
  Hcomp->GetYaxis()->SetBinLabel(42, "#splitline{7-8}{jets}");
  Hcomp->GetYaxis()->SetBinLabel(71, "#splitline{ 9+}{jets}");
  Hcomp->GetXaxis()->LabelsOption("h");
  Hcomp->AddBin(0, 0, 1, 1);
  Hcomp->AddBin(1, 0, 2, 1);
  Hcomp->AddBin(2, 0, 3, 1);
  Hcomp->AddBin(3, 0, 4, 1);
  Hcomp->AddBin(0, 1, 1, 2);
  Hcomp->AddBin(1, 1, 2, 2);
  Hcomp->AddBin(2, 1, 3, 2);
  Hcomp->AddBin(3, 1, 4, 2);
  Hcomp->AddBin(0, 2, 1, 3);
  Hcomp->AddBin(1, 2, 2, 3);
  Hcomp->AddBin(2, 2, 3, 3);
  Hcomp->AddBin(3, 2, 4, 3);
  Hcomp->AddBin(0, 3, 4, 3.5);

  Hcomp->Draw();

  leg1->Draw();
  // leg2->Draw();
  // leg3->Draw();
 
  vector<TPad*> pads;
  TPad* pad1 = new TPad("pad1", "", 0.06, 0.05, 0.343, 0.370);
  TPad* pad2 = new TPad("pad2", "", 0.255, 0.05, 0.545, 0.370);
  TPad* pad3 = new TPad("pad3", "", 0.457, 0.05, 0.747, 0.370);
  TPad* pad4 = new TPad("pad4", "", 0.659, 0.05, 0.949, 0.370);
  TPad* pad5 = new TPad("pad5", "", 0.06, 0.282, 0.343, 0.603);
  TPad* pad6 = new TPad("pad6", "", 0.255, 0.282, 0.545, 0.603);
  TPad* pad7 = new TPad("pad7", "", 0.457, 0.282, 0.747, 0.603);
  TPad* pad8 = new TPad("pad8", "", 0.659, 0.282, 0.949, 0.603);
  TPad* pad9 = new TPad("pad9", "", 0.06, 0.514, 0.343, 0.835);
  TPad* pad10 = new TPad("pad10", "", 0.255, 0.514, 0.545, 0.835);
  TPad* pad11 = new TPad("pad11", "", 0.457, 0.514, 0.747, 0.835);
  TPad* pad12 = new TPad("pad12", "", 0.659, 0.514, 0.949, 0.835);

  pads.push_back(pad1);
  pads.push_back(pad2);
  pads.push_back(pad3);
  pads.push_back(pad4);
  pads.push_back(pad5);
  pads.push_back(pad6);
  pads.push_back(pad7);
  pads.push_back(pad8);
  pads.push_back(pad9);
  pads.push_back(pad10);
  pads.push_back(pad11);
  pads.push_back(pad12);
  for (unsigned int box(0); box<12; box++) {
    pads[box]->SetFillStyle(4000);
    thecanvas->cd();
    pads[box]->Draw();  
    pads[box]->cd();
    bgPies[box]->Draw();
  }



  // CMS_lumi(thecanvas, iPeriod, iPos, lumi_sqrtS);

  thecanvas->cd();
  TLatex * latex = new TLatex();
  latex->SetNDC();

  float t=thecanvas->GetTopMargin();
  float r=thecanvas->GetRightMargin();
  TString cmsText     = "CMS";
  float cmsTextFont   = 61;
  float cmsTextSize      = 0.55;
  latex->SetTextFont(cmsTextFont);
  latex->SetTextSize(cmsTextSize*t);
  latex->DrawLatex(0.12, 0.85, cmsText);

  TString extraText   = "            Supplementary";
  latex->SetTextFont(52);
  float extraOverCmsTextSize  = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;
  latex->SetTextSize(extraTextSize*t);
  float relExtraDY = 1.2;
  latex->DrawLatex(0.12, 0.85, extraText);
  
  latex->SetTextFont(42);
  latex->SetTextAlign(31);
  float lumiTextSize=0.6;
  float lumiTextOffset=0.2;

  latex->SetTextSize(lumiTextSize*t);    
  latex->DrawLatex(1-r,1-t+lumiTextOffset*thecanvas->GetTopMargin(), "(13 TeV)");

    latex->SetTextAlign(12);
  latex->SetTextFont(42);
  latex->SetTextColor(2000);
  latex->SetTextSize(0.035);
  latex->DrawLatex(0.5, 0.935, "arXiv:1602.06581");
  
  thecanvas->Print("plots/Data_BG_Pie_vs_NJets_NBJets.pdf");
  thecanvas->Print("plots/Data_BG_Pie_vs_NJets_NBJets.png");

  // TFile* ofile = new TFile("test_pie.root","recreate");
  // ofile->cd();
  // thecanvas->Write();  
  

  delete Hcomp;
  delete thecanvas;
  delete hlostlep;
  delete hqcd;
  delete hznn;
  delete hhadtau;

  //ofile->Close();

  
}
