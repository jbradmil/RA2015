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

using namespace std;


void fill_lostlep_hists() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile = new TFile("bg_hists/LLPrediction.root", "read");
  TH1D* hin = (TH1D*)infile->Get("Prediction_data/totalPred_LL");
  TH1D* hsystup = (TH1D*)infile->Get("AdditionalContent/totalPropSysUp_LL");
  TH1D* hsystdown = (TH1D*)infile->Get("AdditionalContent/totalPropSysDown_LL");
  TH1D* hnonclosureup = (TH1D*)infile->Get("Prediction_data/totalPredNonClosureUp_LL");
  TH1D* hnonclosuredown = (TH1D*)infile->Get("Prediction_data/totalPredNonClosureDown_LL");
  
  // TH1D* hsystsup[12];
  // hsystsup[0] = (TH1D*)infile->Get("Prediction_data/totalPredIsoTrackSysUp_LL");
  // hsystsup[1] = (TH1D*)infile->Get("Prediction_data/totalPredMTWSysUp_LL");
  // hsystsup[2] = (TH1D*)infile->Get("Prediction_data/totalPredPuritySysUp_LL");
  // hsystsup[3] = (TH1D*)infile->Get("Prediction_data/totalPredSingleLepPuritySysUp_LL");
  // hsystsup[4] = (TH1D*)infile->Get("Prediction_data/totalPredDiLepFoundSysUp_LL");
  // hsystsup[5] = (TH1D*)infile->Get("Prediction_data/totalPredMuIsoSysUp_LL");
  // hsystsup[6] = (TH1D*)infile->Get("Prediction_data/totalPredMuRecoSysUp_LL");
  // hsystsup[7] = (TH1D*)infile->Get("Prediction_data/totalPredMuAccSysUp_LL");
  // hsystsup[8] = (TH1D*)infile->Get("Prediction_data/totalPredElecIsoSysUp_LL");
  // hsystsup[9] = (TH1D*)infile->Get("Prediction_data/totalPredElecRecoSysUp_LL");
  // hsystsup[10] = (TH1D*)infile->Get("Prediction_data/totalPredElecAccSysUp_LL");
  // hsystsup[11] = (TH1D*)infile->Get("Prediction_data/totalPredNonClosureUp_LL");

  // TH1D* hsystsdown[12];
  // hsystsdown[0] = (TH1D*)infile->Get("Prediction_data/totalPredIsoTrackSysDown_LL");
  // hsystsdown[1] = (TH1D*)infile->Get("Prediction_data/totalPredMTWSysDown_LL");
  // hsystsdown[2] = (TH1D*)infile->Get("Prediction_data/totalPredPuritySysDown_LL");
  // hsystsdown[3] = (TH1D*)infile->Get("Prediction_data/totalPredSingleLepPuritySysDown_LL");
  // hsystsdown[4] = (TH1D*)infile->Get("Prediction_data/totalPredDiLepFoundSysDown_LL");
  // hsystsdown[5] = (TH1D*)infile->Get("Prediction_data/totalPredMuIsoSysDown_LL");
  // hsystsdown[6] = (TH1D*)infile->Get("Prediction_data/totalPredMuRecoSysDown_LL");
  // hsystsdown[7] = (TH1D*)infile->Get("Prediction_data/totalPredMuAccSysDown_LL");
  // hsystsdown[8] = (TH1D*)infile->Get("Prediction_data/totalPredElecIsoSysDown_LL");
  // hsystsdown[9] = (TH1D*)infile->Get("Prediction_data/totalPredElecRecoSysDown_LL");
  // hsystsdown[10] = (TH1D*)infile->Get("Prediction_data/totalPredElecAccSysDown_LL");
  // hsystsdown[11] = (TH1D*)infile->Get("Prediction_data/totalPredNonClosureDown_LL");


  // for (unsigned int bin(0); bin<72; bin++) {
  //   double sumw2up(0.), sumw2down(0.);
  //   for (unsigned int ihist(0); ihist<12; ihist++) {
  //     if (hsystsup[ihist]->GetBinContent(bin+1)>0) sumw2up+=(hsystsup[ihist]->GetBinContent(bin+1)-1)*(hsystsup[ihist]->GetBinContent(bin+1)-1);
  //     if (hsystsdown[ihist]->GetBinContent(bin+1)>0) sumw2down+=(hsystsdown[ihist]->GetBinContent(bin+1)-1)*(hsystsdown[ihist]->GetBinContent(bin+1)-1);
  //   }
  //   syst_up[bin] = sumw2up*hin->GetBinContent(bin+1)*hin->GetBinContent(bin+1);
  //   syst_down[bin] = sumw2down*hin->GetBinContent(bin+1)*hin->GetBinContent(bin+1);
  // }
  
  //  TH1D* hNEvts = (TH1D*)infile->Get("Prediction_data/nEvtsCS_LL");
  TProfile* hAvgWeight = (TProfile*)infile->Get("Prediction_MC/avgWeight_LL_MC");
  
  TFile* outfile = new TFile("bg_hists/lostlep_hists.root", "recreate");
  outfile->cd();

  TH1D* hPredAllBins = new TH1D("hPredAllBins", ";Search Bin;Events / Bin", 72, 0.5, 72.5);
  TH1D* hRLL = new TH1D("hRLL", ";Search Bin;RLL", 72, 0.5, 72.5);
  Double_t pred_cv[72];
  Double_t x[72];
  Double_t xl[72];
  Double_t xh[72];
  Double_t stat_up[72];
  Double_t stat_down[72];
  Double_t syst_up[72];
  Double_t syst_down[72];
  Double_t full_err_up[72];
  Double_t full_err_down[72];
  for (unsigned int bin(0); bin<72; bin++) {
    hPredAllBins->SetBinContent(bin+1, hin->GetBinContent(bin+1));
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
    // hPredAllBins->SetBinError(bin, stat_full);
    // printf("Bin %d NCR, 1.84102*MCAvgWeight: %f, %f\n", bin, hNEvts->GetBinContent(bin+1), hAvgWeight->GetBinContent(bin+1)*1.84102);
    // printf("Bin %d pred +/- stat + syst: %f +/- (%f + %f)\n", bin, hin->GetBinContent(bin+1), sqrt(stat_stat2), stat_syst);
    stat_up[bin]=sqrt(stat_up[bin]);
    stat_down[bin]=sqrt(stat_down[bin]);
    syst_up[bin]=sqrt(syst_up[bin]);
    syst_down[bin]=sqrt(syst_down[bin]);
    if (stat_down[bin]==pred_cv[bin]) syst_down[bin]=0;
    if (stat_down[bin]+syst_down[bin]>pred_cv[bin]) syst_down[bin]=pred_cv[bin]-stat_down[bin];

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
  hRLL->Write();
  outfile->Close();

  
  return;
  
}

