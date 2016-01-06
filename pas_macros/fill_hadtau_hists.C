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

void fill_hadtau_hists() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile = new TFile("bg_hists/HadTauEstimation.root", "read");
  TH1D* hin = (TH1D*)infile->Get("searchBin_nominal");
  TH1D* hin_wpoisson = (TH1D*)infile->Get("searchBin_nominal_fullstatuncertainty");
  TH1D* hin_closure = (TH1D*)infile->Get("searchBin_closureUncertainty");
  TH1D* searchBin_BMistagUp = (TH1D*)infile->Get("searchBin_BMistagUp");
  TH1D* searchBin_BMistagDn = (TH1D*)infile->Get("searchBin_BMistagDn");
  TH1D* searchBin_MuRecoIsoUp = (TH1D*)infile->Get("searchBin_MuRecoIsoUp");
  TH1D* searchBin_MuRecoIsoDn = (TH1D*)infile->Get("searchBin_MuRecoIsoDn");
  TH1D* searchBin_MuRecoSysUp = (TH1D*)infile->Get("searchBin_MuRecoSysUp");
  TH1D* searchBin_MuRecoSysDn = (TH1D*)infile->Get("searchBin_MuRecoSysDn");
  TH1D* searchBin_MuIsoSysUp = (TH1D*)infile->Get("searchBin_MuIsoSysUp");
  TH1D* searchBin_MuIsoSysDn = (TH1D*)infile->Get("searchBin_MuIsoSysDn");
  //TH1D* searchBin_StatUncertainties = (TH1D*)infile->Get("");
  TH1D*searchBin_JECSysUp = (TH1D*)infile->Get("searchBin_JECSysUp");
  TH1D*searchBin_JECSysDn = (TH1D*)infile->Get("searchBin_JECSysDn");
  TH1D*searchBin_MTSysUp = (TH1D*)infile->Get("searchBin_MTSysUp");
  TH1D*searchBin_MTSysDn = (TH1D*)infile->Get("searchBin_MTSysDn");
  TH1D*seaerchBin_MtEffStat = (TH1D*)infile->Get("seaerchBin_MtEffStat");
  TH1D*seaerchBin_IsoTrkVetoEffUncertaintyStat = (TH1D*)infile->Get("seaerchBin_IsoTrkVetoEffUncertaintyStat");
  TH1D*seaerchBin_IsoTrkVetoEffUncertaintySys = (TH1D*)infile->Get("seaerchBin_IsoTrkVetoEffUncertaintySys");
  TH1D*seaerchBin_AccStat = (TH1D*)infile->Get("seaerchBin_AccStat");
  TH1D*seaerchBin_AccSysPDFUp = (TH1D*)infile->Get("seaerchBin_AccSysPDFUp");
  TH1D*seaerchBin_AccSysPDFDn = (TH1D*)infile->Get("seaerchBin_AccSysPDFDn");
  TH1D*seaerchBin_AccSysScaleUp = (TH1D*)infile->Get("seaerchBin_AccSysScaleUp");
  TH1D*seaerchBin_AccSysScaleDn = (TH1D*)infile->Get("seaerchBin_AccSysScaleDn");
  TH1D*seaerchBin_MuFromTauStat = (TH1D*)infile->Get("seaerchBin_MuFromTauStat");
  TH1D*searchBin_TrigEffUncertainty = (TH1D*)infile->Get("searchBin_TrigEffUncertainty");

  
  TFile* outfile = new TFile("bg_hists/hadtau_hists.root", "recreate");
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
  TH1D* hRTau = new TH1D("hRTau", ";Search Bin;RTau", 72, 0.5, 72.5);
  for (unsigned int bin(0); bin<72; bin++) {
    hPredAllBins->SetBinContent(bin+1, hin->GetBinContent(bin+1));
    hRTau->SetBinContent(bin+1, 0.275);
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    pred_cv[bin] = hin->GetBinContent(bin+1);
    stat_up[bin]=hin_wpoisson->GetBinError(bin+1)*hin_wpoisson->GetBinError(bin+1);
    stat_down[bin]=hin->GetBinError(bin+1)*hin->GetBinError(bin+1);
    
    double err_closure2 = pow((hin_closure->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1), 2.);
    syst_up[bin]=err_closure2;
    syst_down[bin]=err_closure2;
    syst_up[bin]+=pow((seaerchBin_MtEffStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((seaerchBin_MtEffStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_IsoTrkVetoEffUncertaintyStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((seaerchBin_IsoTrkVetoEffUncertaintyStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_IsoTrkVetoEffUncertaintySys->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((seaerchBin_IsoTrkVetoEffUncertaintySys->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_AccStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((seaerchBin_AccStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_MuFromTauStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((seaerchBin_MuFromTauStat->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);

    syst_up[bin]+=pow((searchBin_MuRecoIsoDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-searchBin_MuRecoIsoUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((searchBin_MuRecoSysDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-searchBin_MuRecoSysUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((searchBin_MuIsoSysDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-searchBin_MuIsoSysUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((searchBin_JECSysDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-searchBin_JECSysUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_AccSysPDFDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-seaerchBin_AccSysPDFUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    syst_up[bin]+=pow((seaerchBin_AccSysScaleDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    syst_down[bin]+=pow((1.-seaerchBin_AccSysScaleUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);

    syst_up[bin]+=pow((1-searchBin_TrigEffUncertainty->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    
    if (searchBin_BMistagUp->GetBinContent(bin+1)>1) syst_up[bin]+=pow((searchBin_BMistagUp->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    else syst_down[bin]+=pow((1.-searchBin_BMistagUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    if (searchBin_BMistagDn->GetBinContent(bin+1)>1) syst_up[bin]+=pow((searchBin_BMistagDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    else syst_down[bin]+=pow((1.-searchBin_BMistagDn->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    if (searchBin_MTSysUp->GetBinContent(bin+1)>1) syst_up[bin]+=pow((searchBin_MTSysUp->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    else syst_down[bin]+=pow((1.-searchBin_MTSysUp->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);
    if (searchBin_MTSysDn->GetBinContent(bin+1)>1) syst_up[bin]+=pow((searchBin_MTSysDn->GetBinContent(bin+1)-1.)*hin->GetBinContent(bin+1),2.);
    else syst_down[bin]+=pow((1.-searchBin_MTSysDn->GetBinContent(bin+1))*hin->GetBinContent(bin+1),2.);

    stat_up[bin]=sqrt(stat_up[bin]);
    stat_down[bin]=sqrt(stat_down[bin]);
    syst_up[bin]=sqrt(syst_up[bin]);
    syst_down[bin]=sqrt(syst_down[bin]);
    //if (stat_down[bin]==pred_cv[bin]) syst_down[bin]=0;
    if (stat_down[bin]+syst_down[bin]>=pred_cv[bin]) syst_down[bin]=pred_cv[bin]-stat_down[bin];
    if (syst_down[bin]<0.0) syst_down[bin]=0.0; // we seem to have a rounding error somewhere--this should fix it
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
  hRTau->Write();
  outfile->Close();

  
  return;
  
}

