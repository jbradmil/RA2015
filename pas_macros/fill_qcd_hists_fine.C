
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
#include "Math/QuantFuncMathCore.h"
#include "inc/jack_style.h"
#include "inc/make_1D_projections.h"

using namespace std;
//const double predSF = 1.63947463;

void fill_qcd_hists_fine() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* infile = new TFile("bg_hists/qcd_hists.root","read");
  TH1D* hPredAllBins = (TH1D*) infile->Get("hPredAllBins");
  TGraphAsymmErrors* gAllFull = (TGraphAsymmErrors*) infile->Get("Graph;1");
  TGraphAsymmErrors* gAllStat = (TGraphAsymmErrors*) infile->Get("Graph;2");
  TGraphAsymmErrors* gAllSyst = (TGraphAsymmErrors*) infile->Get("Graph;3");
  
  TFile* outfile = new TFile("bg_hists/qcd_hists_fine.root", "recreate");
  outfile->cd();

  TH1D* hPredFineBins = new TH1D("hPredFineBins", ";Search Bin;Events / Bin", 220, 0.5, 220.5);
  TH1D* hPredFineStatUp = new TH1D("hPredFineStatUp", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineStatDown = new TH1D("hPredFineStatDown", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineSystUp = new TH1D("hPredFineSystUp", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hPredFineSystDown = new TH1D("hPredFineSystDown", ";Search Bin;", 220, 0.5, 220.5);
  TH1D* hLDPStats = new TH1D("hLDPStats", ";Search Bin;Events / Bin", 220, 0.5, 220.5);
  hLDPStats->SetBinErrorOption(TH1::kPoisson);
  TH1D* hRQCD = new TH1D("hRQCD", ";Search Bin;RQCD", 220, 0.5, 220.5);

  double pred_cv[220] = {93.86, 33.79, 22.10, 3.50, 2.67, 3.34, 0.00, 0.12, 0.26, 0.01, 0.08, 27.48, 10.85, 6.78, 0.62, 0.92, 0.99, 0.01, 0.07, 0.07, 0.00, 0.00, 4.50, 2.25, 0.75, 0.26, 0.14, 0.20, 0.00, 0.03, 0.01, 0.00, 0.00, 0.18, 0.08, 0.06, 0.00, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 80.08, 38.26, 32.93, 1.94, 2.75, 3.97, 0.18, 0.04, 0.40, 0.00, 0.10, 26.90, 13.86, 11.34, 0.69, 0.94, 1.38, 0.00, 0.00, 0.20, 0.00, 0.00, 6.57, 3.74, 1.77, 0.00, 0.09, 0.24, 0.00, 0.00, 0.00, 0.00, 0.03, 0.68, 0.34, 0.31, 0.15, 0.03, 0.07, 0.05, 0.00, 0.00, 0.00, 0.00, 32.25, 19.45, 21.83, 0.42, 1.36, 3.16, 0.03, 0.02, 0.27, 0.00, 0.05, 11.26, 8.76, 9.23, 0.06, 0.45, 1.47, 0.00, 0.00, 0.10, 0.00, 0.02, 3.96, 2.95, 2.20, 0.00, 0.08, 0.32, 0.00, 0.00, 0.06, 0.00, 0.00, 0.21, 0.24, 0.21, 0.04, 0.01, 0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 13.75, 13.95, 19.64, 0.00, 1.08, 2.14, 0.00, 0.06, 0.25, 0.00, 0.19, 7.60, 6.95, 11.53, 0.01, 0.31, 1.31, 0.00, 0.00, 0.12, 0.00, 0.00, 0.20, 1.86, 2.64, 0.00, 0.24, 0.66, 0.00, 0.00, 0.03, 0.00, 0.00, 1.06, 0.56, 0.05, 0.00, 0.03, 0.00, 0.00, 0.04, 0.00, 0.00, 0.00, 0.84, 0.66, 3.08, 0.00, 0.05, 0.38, 0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.73, 1.93, 0.05, 0.01, 0.27, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.21, 1.18, 0.00, 0.02, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.03, 0.00, 0.00, 0.00, 0.00, 0.00};
  double syst[220] = {50.47, 11.55, 5.58, 2.62, 1.63, 1.88, 0.07, 0.14, 0.28, 0.03, 0.09, 14.82, 3.74, 1.75, 0.51, 0.57, 0.57, 0.04, 0.08, 0.08, 0.01, 0.02, 2.48, 0.80, 0.24, 0.22, 0.09, 0.13, 0.02, 0.04, 0.02, 0.01, 0.01, 0.16, 0.05, 0.04, 0.04, 0.02, 0.01, 0.02, 0.01, 0.01, 0.01, 0.01, 47.42, 16.15, 11.64, 1.59, 1.82, 2.45, 0.24, 0.08, 0.44, 0.02, 0.12, 16.02, 5.89, 4.06, 0.62, 0.64, 0.87, 0.05, 0.03, 0.23, 0.02, 0.03, 4.02, 1.64, 0.70, 0.15, 0.11, 0.17, 0.05, 0.02, 0.03, 0.02, 0.04, 0.52, 0.22, 0.17, 0.15, 0.04, 0.07, 0.07, 0.02, 0.03, 0.02, 0.03, 20.07, 8.98, 8.72, 0.45, 0.95, 2.05, 0.08, 0.05, 0.32, 0.03, 0.07, 7.13, 4.09, 3.74, 0.22, 0.34, 0.96, 0.06, 0.04, 0.13, 0.03, 0.04, 2.62, 1.43, 0.97, 0.12, 0.11, 0.24, 0.06, 0.03, 0.08, 0.03, 0.03, 0.27, 0.20, 0.16, 0.07, 0.03, 0.06, 0.06, 0.03, 0.03, 0.03, 0.03, 8.74, 6.58, 8.04, 0.24, 0.78, 1.42, 0.09, 0.09, 0.29, 0.04, 0.22, 4.95, 3.35, 4.77, 0.19, 0.26, 0.90, 0.09, 0.04, 0.16, 0.04, 0.05, 0.88, 1.01, 1.21, 0.14, 0.22, 0.48, 0.09, 0.04, 0.06, 0.04, 0.05, 0.87, 0.37, 0.15, 0.09, 0.06, 0.05, 0.09, 0.06, 0.05, 0.04, 0.05, 0.76, 0.52, 2.07, 0.08, 0.07, 0.34, 0.09, 0.04, 0.06, 0.04, 0.05, 0.23, 0.57, 1.32, 0.09, 0.05, 0.25, 0.08, 0.04, 0.05, 0.04, 0.05, 0.16, 0.24, 0.84, 0.08, 0.04, 0.08, 0.08, 0.04, 0.04, 0.04, 0.04, 0.15, 0.22, 0.10, 0.08, 0.04, 0.05, 0.08, 0.04, 0.04, 0.04, 0.04};
  double Nldp[220] = {2337., 1589., 827., 282., 289., 258., 6., 22., 23., 3., 7., 692., 515., 255., 62., 98., 77., 2., 9., 7., 0., 1., 123., 110., 29., 18., 16., 16., 0., 3., 1., 0., 0., 6., 4., 2., 1., 2., 0., 0., 0., 0., 0., 0., 1055., 943., 637., 95., 165., 162., 6., 8., 19., 0., 5., 398., 357., 223., 40., 61., 58., 0., 0., 9., 0., 0., 112., 103., 37., 7., 13., 11., 0., 0., 0., 0., 1., 12., 12., 6., 4., 2., 3., 1., 0., 0., 0., 0., 347., 387., 339., 23., 65., 102., 1., 2., 11., 0., 2., 147., 191., 146., 11., 26., 49., 0., 1., 5., 0., 1., 58., 65., 37., 3., 10., 12., 0., 0., 2., 0., 0., 4., 7., 4., 1., 1., 2., 0., 0., 0., 0., 0., 101., 193., 211., 5., 37., 49., 0., 2., 6., 0., 4., 67., 114., 126., 4., 14., 35., 0., 0., 4., 0., 0., 18., 38., 34., 2., 12., 17., 0., 0., 1., 0., 0., 10., 10., 2., 0., 2., 0., 0., 1., 0., 0., 0., 6., 12., 38., 0., 2., 10., 0., 0., 1., 0., 0., 2., 13., 24., 1., 2., 8., 0., 0., 0., 0., 0., 1., 6., 16., 0., 1., 2., 0., 0., 0., 0., 0., 0., 4., 1., 0., 0., 1., 0., 0., 0., 0., 0.};
  double Nldp_nonqcd[220] = {352.70, 93.70, 26.10, 133.80, 53.00, 16.10, 5.80, 11.70, 3.80, 1.80, 1.50, 111.00, 34.70, 9.30, 35.70, 16.20, 5.40, 1.70, 3.00, 1.80, 0.50, 1.20, 27.80, 10.40, 1.70, 7.10, 3.90, 1.40, 0.30, 0.30, 0.20, 0.00, 0.00, 2.30, 0.60, 0.00, 1.70, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 175.80, 63.80, 17.50, 52.40, 38.60, 12.70, 2.00, 6.20, 4.10, 0.50, 1.40, 102.70, 38.50, 9.60, 24.90, 18.00, 5.90, 0.20, 1.60, 1.50, 0.10, 0.70, 39.90, 17.10, 3.70, 9.40, 8.90, 1.80, 0.20, 1.00, 0.10, 0.00, 0.00, 4.50, 4.20, 0.10, 0.60, 0.50, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 66.00, 32.30, 13.10, 15.60, 15.40, 7.50, 0.40, 1.40, 2.80, 0.30, 0.50, 48.90, 31.20, 8.10, 9.90, 9.70, 5.20, 0.40, 1.80, 2.00, 0.20, 0.50, 23.50, 11.20, 4.10, 6.10, 6.90, 2.30, 0.00, 0.50, 0.20, 0.00, 0.00, 2.20, 2.70, 0.90, 0.30, 0.70, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 19.40, 19.80, 11.30, 7.20, 10.10, 5.50, 0.20, 0.50, 1.00, 0.20, 0.20, 21.90, 27.70, 8.80, 3.90, 6.40, 8.30, 0.20, 0.20, 1.60, 0.10, 0.10, 16.80, 14.90, 7.20, 2.90, 6.00, 3.50, 0.00, 0.60, 0.40, 0.00, 0.00, 3.70, 3.10, 1.50, 0.40, 1.30, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 2.40, 1.50, 0.60, 0.60, 1.10, 0.60, 0.10, 1.40, 0.60, 0.60, 1.80, 2.40, 1.10, 0.30, 1.60, 1.60, 0.40, 0.00, 0.40, 0.40, 0.40, 1.10, 2.90, 2.00, 0.50, 0.30, 0.80, 0.20, 0.40, 0.20, 0.20, 0.20, 0.30, 0.40, 1.00, 0.00, 0.50, 0.20, 0.10, 0.00, 0.10, 0.10, 0.10};
  double Rqcd[220] = {0.047, 0.023, 0.028, 0.024, 0.011, 0.014, 0.024, 0.011, 0.014, 0.011, 0.014, 0.047, 0.023, 0.028, 0.024, 0.011, 0.014, 0.024, 0.011, 0.014, 0.011, 0.014, 0.047, 0.023, 0.028, 0.024, 0.011, 0.014, 0.024, 0.011, 0.014, 0.011, 0.014, 0.047, 0.023, 0.028, 0.024, 0.011, 0.014, 0.024, 0.011, 0.014, 0.011, 0.014, 0.091, 0.044, 0.053, 0.046, 0.022, 0.027, 0.046, 0.022, 0.027, 0.022, 0.027, 0.091, 0.044, 0.053, 0.046, 0.022, 0.027, 0.046, 0.022, 0.027, 0.022, 0.027, 0.091, 0.044, 0.053, 0.046, 0.022, 0.027, 0.046, 0.022, 0.027, 0.022, 0.027, 0.091, 0.044, 0.053, 0.046, 0.022, 0.027, 0.046, 0.022, 0.027, 0.022, 0.027, 0.115, 0.055, 0.067, 0.057, 0.027, 0.033, 0.057, 0.027, 0.033, 0.027, 0.033, 0.115, 0.055, 0.067, 0.057, 0.027, 0.033, 0.057, 0.027, 0.033, 0.027, 0.033, 0.115, 0.055, 0.067, 0.057, 0.027, 0.033, 0.057, 0.027, 0.033, 0.027, 0.033, 0.115, 0.055, 0.067, 0.057, 0.027, 0.033, 0.057, 0.027, 0.033, 0.027, 0.033, 0.169, 0.081, 0.098, 0.084, 0.040, 0.049, 0.084, 0.040, 0.049, 0.040, 0.049, 0.169, 0.081, 0.098, 0.084, 0.040, 0.049, 0.084, 0.040, 0.049, 0.040, 0.049, 0.169, 0.081, 0.098, 0.084, 0.040, 0.049, 0.084, 0.040, 0.049, 0.040, 0.049, 0.169, 0.081, 0.098, 0.084, 0.040, 0.049, 0.084, 0.040, 0.049, 0.040, 0.049, 0.145, 0.069, 0.084, 0.072, 0.035, 0.042, 0.072, 0.035, 0.042, 0.035, 0.042, 0.145, 0.069, 0.084, 0.072, 0.035, 0.042, 0.072, 0.035, 0.042, 0.035, 0.042, 0.145, 0.069, 0.084, 0.072, 0.035, 0.042, 0.072, 0.035, 0.042, 0.035, 0.042, 0.145, 0.069, 0.084, 0.072, 0.035, 0.042, 0.072, 0.035, 0.042, 0.035, 0.042};

  const double alpha = 1 - 0.6827;
  
  Double_t stat_up[220];
  Double_t stat_down[220];
  Double_t syst_up[220];
  Double_t syst_down[220];
  Double_t full_err_up[220];
  Double_t full_err_down[220];
  Double_t x[220];
  Double_t xl[220];
  Double_t xh[220];

  
  
  for (unsigned int bin(0); bin<220; bin++) {
    double N=std::max(Nldp[bin]-Nldp_nonqcd[bin], 0.);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    
    hPredFineBins->SetBinContent(bin+1, pred_cv[bin]);
    hPredFineStatUp->SetBinContent(bin+1, (U-N)*Rqcd[bin]);
    hPredFineStatDown->SetBinContent(bin+1, (N-L)*Rqcd[bin]);
    hPredFineSystUp->SetBinContent(bin+1, syst[bin]);
    hPredFineSystDown->SetBinContent(bin+1, syst[bin]);
    hLDPStats->SetBinContent(bin+1, std::max(Nldp[bin]-Nldp_nonqcd[bin], 0.));
    hRQCD->SetBinContent(bin+1, Rqcd[bin]);   
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    // if (hLDPStats->GetBinContent(bin+1)>0){
      //      stat_up[bin]=pred_cv[bin]*sqrt(Nldp[bin])/Nldp[bin];
      //      stat_down[bin]=pred_cv[bin]*sqrt(Nldp[bin])/Nldp[bin];

    stat_up[bin]=(U-N)*Rqcd[bin];
    stat_down[bin]=(N-L)*Rqcd[bin];
    //    printf("Bin %d: Nldp=%d, Rqcd=%3.3f, estat=+%3.3f-%3.3f\n", bin+1, (int)hLDPStats->GetBinContent(bin+1), Rqcd[bin], stat_up[bin], stat_down[bin]);
    // }
    // else {
    //   stat_up[bin]=1.1*Rqcd[bin];
    //   stat_down[bin]=0.;
    // }
    syst_up[bin]=syst[bin];
    syst_down[bin]=syst[bin];
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
  hRQCD->Write();
  hLDPStats->Write();

  hPredFineStatUp->Write();
  hPredFineStatDown->Write();
  hPredFineSystUp->Write();
  hPredFineSystDown->Write();

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



  TH1D* hPredHT_7j_0b = MakeHTProjection(hPredFineBins, 3, 4, 0, 0, 1, 3);
  hPredHT_7j_0b->Write("hPredHT_7j_0b");

  TGraphAsymmErrors* gHT_7j_0bFull = MakeHTProjection(hPredFineBins, gFull, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bFull->Write("gHT_7j_0bFull");
  TGraphAsymmErrors* gHT_7j_0bStat = MakeHTProjection(hPredFineBins, gStat, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bStat->Write("gHT_7j_0bStat");
  TGraphAsymmErrors* gHT_7j_0bSyst = MakeHTProjection(hPredFineBins, gSyst, 3, 4, 0, 0, 1, 3);
  gHT_7j_0bSyst->Write("gHT_7j_0bSyst");

  
  
  TH1D* hPredMHT_9j_3b = MakeMHTProjection(hPredFineBins, 4, 4, 3, 3);
  hPredMHT_9j_3b->Write("hPredMHT_9j_3b");

  TGraphAsymmErrors* gMHT_9j_3bFull = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllFull, gFull, 4, 4, 3, 3, true);
  gMHT_9j_3bFull->Write("gMHT_9j_3bFull");
  TGraphAsymmErrors* gMHT_9j_3bStat = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllStat, gStat, 4, 4, 3, 3);
  gMHT_9j_3bStat->Write("gMHT_9j_3bStat");
  TGraphAsymmErrors* gMHT_9j_3bSyst = MakeMHTProjection(hPredAllBins, hPredFineBins, gAllSyst, gSyst, 4, 4, 3, 3);
  gMHT_9j_3bSyst->Write("gMHT_9j_3bSyst");

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

