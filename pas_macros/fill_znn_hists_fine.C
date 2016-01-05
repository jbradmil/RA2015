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
#include "inc/make_1D_projections.h"

using namespace std;
//const double predSF = 1.63947463;



void fill_znn_hists_fine() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile_fine = new TFile("bg_hists/ZinvHistos_hdp.root", "read");
  TH1D* hin = (TH1D*)infile_fine->Get("ZinvBGpred");
  TH1D* hin_systup = (TH1D*)infile_fine->Get("ZinvBGsysUp");
  TH1D* hin_0evt_statup = (TH1D*)infile_fine->Get("ZinvBG0EVsysUp");
  TH1D* hin_systdown = (TH1D*)infile_fine->Get("ZinvBGsysLow");
  // TH1D* hin_0evt_statdown = (TH1D*)infile_fine->Get("ZinvBG0EVsysLow");

  TFile* infile = new TFile("bg_hists/znn_hists.root","read");
  TH1D* hPredAllBins = (TH1D*) infile->Get("hPredAllBins");
  TGraphAsymmErrors* gAllFull = (TGraphAsymmErrors*) infile->Get("Graph;1");
  TGraphAsymmErrors* gAllStat = (TGraphAsymmErrors*) infile->Get("Graph;2");
  TGraphAsymmErrors* gAllSyst = (TGraphAsymmErrors*) infile->Get("Graph;3");
  
  TFile* outfile = new TFile("bg_hists/znn_hists_fine.root", "recreate");
  outfile->cd();

  Double_t pred_cv[220];
  Double_t stat_up[220];
  Double_t stat_down[220];
  Double_t syst_up[220];
  Double_t syst_down[220];
  Double_t full_err_up[220];
  Double_t full_err_down[220];
  Double_t x[220];
  Double_t xl[220];
  Double_t xh[220];

  TH1D* hPredFineBins = new TH1D("hPredFineBins", ";Search Bin;Events / Bin", 220, 0.5, 220.5);
  TH1D* hRZ0 = new TH1D("hRZ0", ";Search Bin;RZ0", 220, 0.5, 220.5);
  for (unsigned int bin(0); bin<220; bin++) {

    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;

    stat_up[bin]=hin->GetBinError(bin+1);
    stat_down[bin]=hin->GetBinError(bin+1);
    if (stat_up[bin]==0) stat_up[bin]= hin_0evt_statup->GetBinError(bin+1);
    syst_up[bin]=hin_systup->GetBinContent(bin+1);
    syst_down[bin]=hin_systdown->GetBinContent(bin+1);
   
    hRZ0->SetBinContent(bin+1, hin_0evt_statup->GetBinContent(bin+1));

    hPredFineBins->SetBinContent(bin+1, hin->GetBinContent(bin+1));
    pred_cv[bin] = hin->GetBinContent(bin+1);

   //  if (bin<24) {
   //    hPredFineBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
   //    pred_cv[bin] = 0.92*hin->GetBinContent(bin+1);
   //    stat_up[bin]*=0.92;
   //    stat_down[bin]*=0.92;
   //    syst_up[bin]*=0.92;
   //    syst_down[bin]*=0.92;
   // }
   //  else if (bin<48) {
   //    hPredFineBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
   //    pred_cv[bin] = 0.92*hin->GetBinContent(bin+1);
   //    stat_up[bin]*=0.92;
   //    stat_down[bin]*=0.92;
   //    syst_up[bin]*=0.92;
   //    syst_down[bin]*=0.92;
   //  }
   //  else {
   //    hPredFineBins->SetBinContent(bin+1, 0.92*hin->GetBinContent(bin+1));
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

  TGraphAsymmErrors* gFull = new TGraphAsymmErrors(220, x, pred_cv, xl, xh, full_err_down, full_err_up);
  gFull->Write("gFull");
  TGraphAsymmErrors* gStat = new TGraphAsymmErrors(220, x, pred_cv, xl, xh, stat_down, stat_up);
  gStat->Write("gStat");
  TGraphAsymmErrors* gSyst = new TGraphAsymmErrors(220, x, pred_cv, xl, xh, syst_down, syst_up);
  gSyst->Write("gSyst");
  hPredFineBins->Write();
  hRZ0->Write();

TH1D* hPredNJets = MakeNJetsProjection(hPredFineBins);
  hPredNJets->Write("hPredNJets");

  TGraphAsymmErrors* gNJetsFull = MakeNJetsProjection(hPredFineBins, gFull);
  gNJetsFull->Write("gNJetsFull");
  TGraphAsymmErrors* gNJetsStat = MakeNJetsProjection(hPredFineBins, gStat);
  gNJetsStat->Write("gNJetsStat");
  TGraphAsymmErrors* gNJetsSyst = MakeNJetsProjection(hPredFineBins, gSyst);
  gNJetsSyst->Write("gNJetsSyst");

  TH1D* hPredNJets_2b = MakeNJetsProjection(hPredFineBins, 22, 43);
  hPredNJets_2b->Write("hPredNJets_2b");

  TGraphAsymmErrors* gNJets_2bFull = MakeNJetsProjection(hPredFineBins, gFull, 22, 43);
  gNJets_2bFull->Write("gNJets_2bFull");
  TGraphAsymmErrors* gNJets_2bStat = MakeNJetsProjection(hPredFineBins, gStat, 22, 43);
  gNJets_2bStat->Write("gNJets_2bStat");
  TGraphAsymmErrors* gNJets_2bSyst = MakeNJetsProjection(hPredFineBins, gSyst, 22, 43);
  gNJets_2bSyst->Write("gNJets_2bSyst");

  TH1D* hPredNJets_0b = MakeNJetsProjection(hPredFineBins, 6, 10);
  hPredNJets_0b->Write("hPredNJets_0b");

  TGraphAsymmErrors* gNJets_0bFull = MakeNJetsProjection(hPredFineBins, gFull, 6, 10);
  gNJets_0bFull->Write("gNJets_0bFull");
  TGraphAsymmErrors* gNJets_0bStat = MakeNJetsProjection(hPredFineBins, gStat, 6, 10);
  gNJets_0bStat->Write("gNJets_0bStat");
  TGraphAsymmErrors* gNJets_0bSyst = MakeNJetsProjection(hPredFineBins, gSyst, 6, 10);
  gNJets_0bSyst->Write("gNJets_0bSyst");

  TH1D* hPredNBJets = MakeNBJetsProjection(hPredFineBins);
  hPredNBJets->Write("hPredNBJets");

  TGraphAsymmErrors* gNBJetsFull = MakeNBJetsProjection(hPredFineBins, gFull);
  gNBJetsFull->Write("gNBJetsFull");
  TGraphAsymmErrors* gNBJetsStat = MakeNBJetsProjection(hPredFineBins, gStat);
  gNBJetsStat->Write("gNBJetsStat");
  TGraphAsymmErrors* gNBJetsSyst = MakeNBJetsProjection(hPredFineBins, gSyst);
  gNBJetsSyst->Write("gNBJetsSyst");

  TH1D* hPredNBJets_MHT500 = MakeNBJetsProjection(hPredFineBins, 6, 10);
  hPredNBJets_MHT500->Write("hPredNBJets_MHT500");

  TGraphAsymmErrors* gNBJets_MHT500Full = MakeNBJetsProjection(hPredFineBins, gFull, 6, 10);
  gNBJets_MHT500Full->Write("gNBJets_MHT500Full");
  TGraphAsymmErrors* gNBJets_MHT500Stat = MakeNBJetsProjection(hPredFineBins, gStat, 6, 10);
  gNBJets_MHT500Stat->Write("gNBJets_MHT500Stat");
  TGraphAsymmErrors* gNBJets_MHT500Syst = MakeNBJetsProjection(hPredFineBins, gSyst, 6, 10);
  gNBJets_MHT500Syst->Write("gNBJets_MHT500Syst");

 
  TH1D* hPredMHT = MakeMHTProjection(hPredFineBins);
  hPredMHT->Write("hPredMHT");

  TGraphAsymmErrors* gMHTFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull);
  gMHTFull->Write("gMHTFull");
  TGraphAsymmErrors* gMHTStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat);
  gMHTStat->Write("gMHTStat");
  TGraphAsymmErrors* gMHTSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst);
  gMHTSyst->Write("gMHTSyst");

  TH1D* hPredMHT_3b = MakeMHTProjection(hPredFineBins, 0, 4, 3, 3);
  hPredMHT_3b->Write("hPredMHT_3b");

  TGraphAsymmErrors* gMHT_3bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 0, 4, 3, 3);
  gMHT_3bFull->Write("gMHT_3bFull");
  TGraphAsymmErrors* gMHT_3bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 0, 4, 3, 3);
  gMHT_3bStat->Write("gMHT_3bStat");
  TGraphAsymmErrors* gMHT_3bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 0, 4, 3, 3);
  gMHT_3bSyst->Write("gMHT_3bSyst");

  TH1D* hPredMHT_9j_2b = MakeMHTProjection(hPredFineBins, 4, 4, 2, 3);
  hPredMHT_9j_2b->Write("hPredMHT_9j_2b");

  TGraphAsymmErrors* gMHT_9j_2bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 4, 4, 2, 3);
  gMHT_9j_2bFull->Write("gMHT_9j_2bFull");
  TGraphAsymmErrors* gMHT_9j_2bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 4, 4, 2, 3);
  gMHT_9j_2bStat->Write("gMHT_9j_2bStat");
  TGraphAsymmErrors* gMHT_9j_2bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 4, 4, 2, 3);
  gMHT_9j_2bSyst->Write("gMHT_9j_2bSyst");

  TH1D* hPredMHT_6j_0b = MakeMHTProjection(hPredFineBins, 2, 4, 0, 0);
  hPredMHT_6j_0b->Write("hPredMHT_6j_0b");

  TGraphAsymmErrors* gMHT_6j_0bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 2, 4, 0, 0);
  gMHT_6j_0bFull->Write("gMHT_6j_0bFull");
  TGraphAsymmErrors* gMHT_6j_0bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 2, 4, 0, 0);
  gMHT_6j_0bStat->Write("gMHT_6j_0bStat");
  TGraphAsymmErrors* gMHT_6j_0bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 2, 4, 0, 0);
  gMHT_6j_0bSyst->Write("gMHT_6j_0bSyst");

  TH1D* hPredHT = MakeHTProjection(hPredFineBins);
  hPredHT->Write("hPredHT");

  TGraphAsymmErrors* gHTFull = MakeHTProjection(hPredFineBins, gFull);
  gHTFull->Write("gHTFull");
  TGraphAsymmErrors* gHTStat = MakeHTProjection(hPredFineBins, gStat);
  gHTStat->Write("gHTStat");
  TGraphAsymmErrors* gHTSyst = MakeHTProjection(hPredFineBins, gSyst);
  gHTSyst->Write("gHTSyst");

  TH1D* hPredHT_MHT500 = MakeHTProjection(hPredFineBins, 2, 4, 0, 0, 2, 3);
  hPredHT_MHT500->Write("hPredHT_MHT500");

  TGraphAsymmErrors* gHT_MHT500Full = MakeHTProjection(hPredFineBins, gFull, 2, 4, 0, 0, 2, 3);
  gHT_MHT500Full->Write("gHT_MHT500Full");
  TGraphAsymmErrors* gHT_MHT500Stat = MakeHTProjection(hPredFineBins, gStat, 2, 4, 0, 0, 2, 3);
  gHT_MHT500Stat->Write("gHT_MHT500Stat");
  TGraphAsymmErrors* gHT_MHT500Syst = MakeHTProjection(hPredFineBins, gSyst, 2, 4, 0, 0, 2, 3);
  gHT_MHT500Syst->Write("gHT_MHT500Syst");

  TH1D* hPredHT_MHT500_0b_7j = MakeHTProjection(hPredFineBins, 3, 4, 0, 0, 2, 3);
  hPredHT_MHT500_0b_7j->Write("hPredHT_MHT500_0b_7j");

  TGraphAsymmErrors* gHT_MHT500_0b_7jFull = MakeHTProjection(hPredFineBins, gFull, 3, 4, 0, 0, 2, 3);
  gHT_MHT500_0b_7jFull->Write("gHT_MHT500_0b_7jFull");
  TGraphAsymmErrors* gHT_MHT500_0b_7jStat = MakeHTProjection(hPredFineBins, gStat, 3, 4, 0, 0, 2, 3);
  gHT_MHT500_0b_7jStat->Write("gHT_MHT500_0b_7jStat");
  TGraphAsymmErrors* gHT_MHT500_0b_7jSyst = MakeHTProjection(hPredFineBins, gSyst, 3, 4, 0, 0, 2, 3);
  gHT_MHT500_0b_7jSyst->Write("gHT_MHT500_0b_7jSyst");
  
  TH1D* hPredHT_MHT750 = MakeHTProjection(hPredFineBins, 0, 4, 0, 0, 3, 3);
  hPredHT_MHT750->Write("hPredHT_MHT750");

  TGraphAsymmErrors* gHT_MHT750Full = MakeHTProjection(hPredFineBins, gFull, 0, 4, 0, 0, 3, 3);
  gHT_MHT750Full->Write("gHT_MHT750Full");
  TGraphAsymmErrors* gHT_MHT750Stat = MakeHTProjection(hPredFineBins, gStat, 0, 4, 0, 0, 3, 3);
  gHT_MHT750Stat->Write("gHT_MHT750Stat");
  TGraphAsymmErrors* gHT_MHT750Syst = MakeHTProjection(hPredFineBins, gSyst, 0, 4, 0, 0, 3, 3);
  gHT_MHT750Syst->Write("gHT_MHT750Syst");

  TH1D* hPredHT_3b = MakeHTProjection(hPredFineBins, 0, 4, 3, 3, 0, 3);
  hPredHT_3b->Write("hPredHT_3b");

  TGraphAsymmErrors* gHT_3bFull = MakeHTProjection(hPredFineBins, gFull, 0, 4, 3, 3, 0, 3);
  gHT_3bFull->Write("gHT_3bFull");
  TGraphAsymmErrors* gHT_3bStat = MakeHTProjection(hPredFineBins, gStat, 0, 4, 3, 3, 0, 3);
  gHT_3bStat->Write("gHT_3bStat");
  TGraphAsymmErrors* gHT_3bSyst = MakeHTProjection(hPredFineBins, gSyst, 0, 4, 3, 3, 0, 3);
  gHT_3bSyst->Write("gHT_3bSyst");

  TH1D* hPredMHT_9j_3b = MakeMHTProjection(hPredFineBins, 4, 4, 3, 3);
  hPredMHT_9j_3b->Write("hPredMHT_9j_3b");

  TGraphAsymmErrors* gMHT_9j_3bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 4, 4, 3, 3, true);
  gMHT_9j_3bFull->Write("gMHT_9j_3bFull");
  TGraphAsymmErrors* gMHT_9j_3bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 4, 4, 3, 3);
  gMHT_9j_3bStat->Write("gMHT_9j_3bStat");
  TGraphAsymmErrors* gMHT_9j_3bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 4, 4, 3, 3);
  gMHT_9j_3bSyst->Write("gMHT_9j_3bSyst");

  TH1D* hPredHT_7j_0b = MakeHTProjection(hPredFineBins, 3, 4, 0, 0, 1, 3);
  hPredHT_7j_0b->Write("hPredHT_7j_0b");

  TGraphAsymmErrors* gHT_7j_0bFull = MakeHTProjection(hPredFineBins, gFull, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bFull->Write("gHT_7j_0bFull");
  TGraphAsymmErrors* gHT_7j_0bStat = MakeHTProjection(hPredFineBins, gStat, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bStat->Write("gHT_7j_0bStat");
  TGraphAsymmErrors* gHT_7j_0bSyst = MakeHTProjection(hPredFineBins, gSyst, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bSyst->Write("gHT_7j_0bSyst");

  // for T1qqqqVV
  
  TH1D* hPredNJets_0b_MHT500 = MakeNJetsProjectionV2(hPredFineBins, 0, 0, 6, 10);
  hPredNJets_0b_MHT500->Write("hPredNJets_0b_MHT500");

  TGraphAsymmErrors* gNJets_0b_MHT500Full = MakeNJetsProjectionV2(hPredFineBins, gFull, 0, 0, 6, 10);
  gNJets_0b_MHT500Full->Write("gNJets_0b_MHT500Full");
  TGraphAsymmErrors* gNJets_0b_MHT500Stat = MakeNJetsProjectionV2(hPredFineBins, gStat, 0, 0, 6, 10);
  gNJets_0b_MHT500Stat->Write("gNJets_0b_MHT500Stat");
  TGraphAsymmErrors* gNJets_0b_MHT500Syst = MakeNJetsProjectionV2(hPredFineBins, gSyst, 0, 0, 6, 10);
  gNJets_0b_MHT500Syst->Write("gNJets_0b_MHT500Syst");
  
  outfile->Close();

  
  return;
  
}

