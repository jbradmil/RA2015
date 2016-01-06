/* 
   Script for drawing a branch from reduec trees with a set of cuts.
*/

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
#include "TGraphAsymmErrors.h"
#include "inc/jack_style.h"

using namespace std;
//const double predSF = 1.63947463;



void fill_znn_hists() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile = new TFile("bg_hists/ZinvHistos.root", "read");
  TH1D* hin = (TH1D*)infile->Get("ZinvBGpred");
  TH1D* hin_systup = (TH1D*)infile->Get("ZinvBGsysUp");
  TH1D* hin_0evt_statup = (TH1D*)infile->Get("ZinvBG0EVsysUp");
  TH1D* hin_systdown = (TH1D*)infile->Get("ZinvBGsysLow");
  // TH1D* hin_0evt_statdown = (TH1D*)infile->Get("ZinvBG0EVsysLow");
  
  TFile* outfile = new TFile("bg_hists/znn_hists.root", "recreate");
  outfile->cd();

  Double_t pred_cv[72];
  Double_t stat_up[72];
  Double_t stat_down[72];
  Double_t syst_up[72];
  Double_t syst_down[72];
  Double_t full_err_up[72];
  Double_t full_err_down[72];
  Double_t x[72];
  Double_t xl[72];
  Double_t xh[72];

  TH1D* hPredAllBins = new TH1D("hPredAllBins", ";Search Bin;Events / Bin", 72, 0.5, 72.5);
  TH1D* hRZ0 = new TH1D("hRZ0", ";Search Bin;RZ0", 72, 0.5, 72.5);
  for (unsigned int bin(0); bin<72; bin++) {

    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;

    stat_up[bin]=hin->GetBinError(bin+1);
    stat_down[bin]=hin->GetBinError(bin+1);
    if (stat_up[bin]==0) stat_up[bin]= hin_0evt_statup->GetBinError(bin+1);
    syst_up[bin]=hin_systup->GetBinContent(bin+1);
    syst_down[bin]=hin_systdown->GetBinContent(bin+1);
   
    hRZ0->SetBinContent(bin+1, hin_0evt_statup->GetBinContent(bin+1));

    hPredAllBins->SetBinContent(bin+1, hin->GetBinContent(bin+1));
    pred_cv[bin] = hin->GetBinContent(bin+1);

   //  if (bin<24) {
   //    hPredAllBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
   //    pred_cv[bin] = 0.92*hin->GetBinContent(bin+1);
   //    stat_up[bin]*=0.92;
   //    stat_down[bin]*=0.92;
   //    syst_up[bin]*=0.92;
   //    syst_down[bin]*=0.92;
   // }
   //  else if (bin<48) {
   //    hPredAllBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
   //    pred_cv[bin] = 0.92*hin->GetBinContent(bin+1);
   //    stat_up[bin]*=0.92;
   //    stat_down[bin]*=0.92;
   //    syst_up[bin]*=0.92;
   //    syst_down[bin]*=0.92;
   //  }
   //  else {
   //    hPredAllBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
   //    pred_cv[bin] = 0.92*hin->GetBinContent(bin+1);
   //    stat_up[bin]*=0.92;
   //    stat_down[bin]*=0.92;
   //    syst_up[bin]*=0.92;
   //    syst_down[bin]*=0.92;
   //  }


    if (stat_down[bin]>pred_cv[bin]) stat_down[bin]=pred_cv[bin];
    if (stat_down[bin]+syst_down[bin]>=pred_cv[bin]) syst_down[bin]=pred_cv[bin]-stat_down[bin];
    
    full_err_up[bin]=sqrt(stat_up[bin]*stat_up[bin]+syst_up[bin]*syst_up[bin]);
    full_err_down[bin]=sqrt(stat_down[bin]*stat_down[bin]+syst_down[bin]*syst_down[bin]);
    printf("Bin %d: %3.5f + %3.5f + %3.5f - %3.5f - %3.5f\n", bin+1, pred_cv[bin], stat_up[bin], syst_up[bin], stat_down[bin], syst_down[bin]);
  }

  TGraphAsymmErrors* gFull = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, full_err_down, full_err_up);
  TGraphAsymmErrors* gStat = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, stat_down, stat_up);
  TGraphAsymmErrors* gSyst = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, syst_down, syst_up);
  gFull->Write();
  gStat->Write();
  gSyst->Write();
  hPredAllBins->Write();
  hRZ0->Write();
  outfile->Close();

  
  return;
  
}

