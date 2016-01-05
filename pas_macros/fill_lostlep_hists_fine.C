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
#include "TProfile.h"
#include "TFileCollection.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "inc/jack_style.h"
#include "inc/make_1D_projections.h"

using namespace std;


void fill_lostlep_hists_fine() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile_fine = new TFile("bg_hists/LLPrediction_QCD_HDP_Dez08.root", "read");
  TH1D* hin = (TH1D*)infile_fine->Get("Prediction_data/totalPred_LL");
  TH1D* hsystup = (TH1D*)infile_fine->Get("AdditionalContent/totalPropSysUp_LL");
  TH1D* hsystdown = (TH1D*)infile_fine->Get("AdditionalContent/totalPropSysDown_LL");
  TH1D* hnonclosureup = (TH1D*)infile_fine->Get("Prediction_data/totalPredNonClosureUp_LL");
  TH1D* hnonclosuredown = (TH1D*)infile_fine->Get("Prediction_data/totalPredNonClosureDown_LL");
  
  // TH1D* hsystsup[12];
  // hsystsup[0] = (TH1D*)infile_fine->Get("Prediction_data/totalPredIsoTrackSysUp_LL");
  // hsystsup[1] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMTWSysUp_LL");
  // hsystsup[2] = (TH1D*)infile_fine->Get("Prediction_data/totalPredPuritySysUp_LL");
  // hsystsup[3] = (TH1D*)infile_fine->Get("Prediction_data/totalPredSingleLepPuritySysUp_LL");
  // hsystsup[4] = (TH1D*)infile_fine->Get("Prediction_data/totalPredDiLepFoundSysUp_LL");
  // hsystsup[5] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuIsoSysUp_LL");
  // hsystsup[6] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuRecoSysUp_LL");
  // hsystsup[7] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuAccSysUp_LL");
  // hsystsup[8] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecIsoSysUp_LL");
  // hsystsup[9] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecRecoSysUp_LL");
  // hsystsup[10] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecAccSysUp_LL");
  // hsystsup[11] = (TH1D*)infile_fine->Get("Prediction_data/totalPredNonClosureUp_LL");

  // TH1D* hsystsdown[12];
  // hsystsdown[0] = (TH1D*)infile_fine->Get("Prediction_data/totalPredIsoTrackSysDown_LL");
  // hsystsdown[1] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMTWSysDown_LL");
  // hsystsdown[2] = (TH1D*)infile_fine->Get("Prediction_data/totalPredPuritySysDown_LL");
  // hsystsdown[3] = (TH1D*)infile_fine->Get("Prediction_data/totalPredSingleLepPuritySysDown_LL");
  // hsystsdown[4] = (TH1D*)infile_fine->Get("Prediction_data/totalPredDiLepFoundSysDown_LL");
  // hsystsdown[5] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuIsoSysDown_LL");
  // hsystsdown[6] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuRecoSysDown_LL");
  // hsystsdown[7] = (TH1D*)infile_fine->Get("Prediction_data/totalPredMuAccSysDown_LL");
  // hsystsdown[8] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecIsoSysDown_LL");
  // hsystsdown[9] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecRecoSysDown_LL");
  // hsystsdown[10] = (TH1D*)infile_fine->Get("Prediction_data/totalPredElecAccSysDown_LL");
  // hsystsdown[11] = (TH1D*)infile_fine->Get("Prediction_data/totalPredNonClosureDown_LL");


  // for (unsigned int bin(0); bin<72; bin++) {
  //   double sumw2up(0.), sumw2down(0.);
  //   for (unsigned int ihist(0); ihist<12; ihist++) {
  //     if (hsystsup[ihist]->GetBinContent(bin+1)>0) sumw2up+=(hsystsup[ihist]->GetBinContent(bin+1)-1)*(hsystsup[ihist]->GetBinContent(bin+1)-1);
  //     if (hsystsdown[ihist]->GetBinContent(bin+1)>0) sumw2down+=(hsystsdown[ihist]->GetBinContent(bin+1)-1)*(hsystsdown[ihist]->GetBinContent(bin+1)-1);
  //   }
  //   syst_up[bin] = sumw2up*hin->GetBinContent(bin+1)*hin->GetBinContent(bin+1);
  //   syst_down[bin] = sumw2down*hin->GetBinContent(bin+1)*hin->GetBinContent(bin+1);
  // }
  
  //  TH1D* hNEvts = (TH1D*)infile_fine->Get("Prediction_data/nEvtsCS_LL");
  TProfile* hAvgWeight = (TProfile*)infile_fine->Get("Prediction_MC/avgWeight_LL_MC");

  TFile* infile = new TFile("bg_hists/lostlep_hists.root","read");
  TH1D* hPredAllBins = (TH1D*) infile->Get("hPredAllBins");
  TGraphAsymmErrors* gAllFull = (TGraphAsymmErrors*) infile->Get("Graph;1");
  TGraphAsymmErrors* gAllStat = (TGraphAsymmErrors*) infile->Get("Graph;2");
  TGraphAsymmErrors* gAllSyst = (TGraphAsymmErrors*) infile->Get("Graph;3");
  
  TFile* outfile = new TFile("bg_hists/lostlep_hists_fine.root", "recreate");
  outfile->cd();

  TH1D* hPredFineBins = new TH1D("hPredFineBins", ";Search Bin;Events / Bin", 220, 0.5, 220.5);
  TH1D* hPredFineStatUp = new TH1D("hPredFineStatUp", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineStatDown = new TH1D("hPredFineStatDown", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineSystUp = new TH1D("hPredFineSystUp", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineSystDown = new TH1D("hPredFineSystDown", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hRLL = new TH1D("hRLL", ";Search Bin;RLL", 220, 0.5, 220.5);
  Double_t pred_cv[220];
  Double_t x[220];
  Double_t xl[220];
  Double_t xh[220];
  Double_t stat_up[220];
  Double_t stat_down[220];
  Double_t syst_up[220];
  Double_t syst_down[220];
  Double_t full_err_up[220];
  Double_t full_err_down[220];
  for (unsigned int bin(0); bin<220; bin++) {
    hPredFineBins->SetBinContent(bin+1, hin->GetBinContent(bin+1));
    pred_cv[bin] = hin->GetBinContent(bin+1);
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    stat_up[bin] = pow(hAvgWeight->GetBinContent(bin+1)*1.84102, 2.);
    hRLL->SetBinContent(bin+1, hAvgWeight->GetBinContent(bin+1)*1.84102);
    stat_up[bin] += pow(hin->GetBinError(bin+1), 2.);
    stat_down[bin] = pow(hin->GetBinError(bin+1), 2.);
    syst_up[bin] = 0.;
    syst_down[bin] = 0.;
    if (hsystup->GetBinContent(bin+1)>0) {
      syst_up[bin] += pow((hsystup->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1), 2.);
    }
    if (hsystdown->GetBinContent(bin+1)>0) {
      syst_down[bin] += pow((1.-hsystdown->GetBinContent(bin+1))*hin->GetBinContent(bin+1), 2.);
    }
    syst_up[bin] += pow((hnonclosureup->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1), 2.);
    syst_down[bin] += pow((1.-hnonclosuredown->GetBinContent(bin+1))*hin->GetBinContent(bin+1), 2.);

    // //    double stat_statup2 = hAvgWeight->GetBinContent(bin)*1.84102*hAvgWeight->GetBinContent(bin)*1.84102;
    // // stat_stat2 += (hin->GetBinError(bin)*hin->GetBinError(bin));
    // double stat_syst = total_systs[bin-1];
    // double stat_full = sqrt(stat_stat2 + stat_syst*stat_syst);
    // hPredFineBins->SetBinError(bin, stat_full);
    // printf("Bin %d NCR, 1.84102*MCAvgWeight: %f, %f\n", bin, hNEvts->GetBinContent(bin+1), hAvgWeight->GetBinContent(bin+1)*1.84102);
    // printf("Bin %d pred +/- stat + syst: %f +/- (%f + %f)\n", bin, hin->GetBinContent(bin+1), sqrt(stat_stat2), stat_syst);
    stat_up[bin]=sqrt(stat_up[bin]);
    stat_down[bin]=sqrt(stat_down[bin]);
    syst_up[bin]=sqrt(syst_up[bin]);
    syst_down[bin]=sqrt(syst_down[bin]);
    if (stat_down[bin]==pred_cv[bin]) syst_down[bin]=0;
    if (stat_down[bin]+syst_down[bin]>pred_cv[bin]) syst_down[bin]=pred_cv[bin]-stat_down[bin];

    hPredFineStatUp->SetBinContent(bin+1, stat_up[bin]);
    hPredFineStatDown->SetBinContent(bin+1, stat_down[bin]);
    hPredFineSystUp->SetBinContent(bin+1, syst_up[bin]);
    hPredFineSystDown->SetBinContent(bin+1, syst_down[bin]);    

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
  hPredFineStatUp->Write();
  hPredFineStatDown->Write();
  hPredFineSystUp->Write();
  hPredFineSystDown->Write();
  hRLL->Write();

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

  TH1D* hPredMHT_3b = MakeMHTProjection(hPredFineBins, 0, 4, 3, 3);
  hPredMHT_3b->Write("hPredMHT_3b");

  TGraphAsymmErrors* gMHT_3bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 0, 4, 3, 3);
  gMHT_3bFull->Write("gMHT_3bFull");
  TGraphAsymmErrors* gMHT_3bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 0, 4, 3, 3);
  gMHT_3bStat->Write("gMHT_3bStat");
  TGraphAsymmErrors* gMHT_3bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 0, 4, 3, 3);
  gMHT_3bSyst->Write("gMHT_3bSyst");

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

