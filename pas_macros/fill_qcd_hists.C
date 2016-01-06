
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
#include "Math.h"
#include "Math/QuantFuncMathCore.h"
#include "jack_style.h"

using namespace std;
//const double predSF = 1.63947463;

void fill_qcd_hists() {

  TH1::SetDefaultSumw2(1);
  //gROOT->SetBatch(1);

  TFile* outfile = new TFile("qcd_hists.root", "recreate");
  outfile->cd();

  TH1D* hPredAllBins = new TH1D("hPredAllBins", ";Search Bin;Events / Bin", 72, 0.5, 72.5);
  TH1D* hLDPStats = new TH1D("hLDPStats", ";Search Bin;Events / Bin", 72, 0.5, 72.5);
  hLDPStats->SetBinErrorOption(TH1::kPoisson);
  TH1D* hRQCD = new TH1D("hRQCD", ";Search Bin;RQCD", 72, 0.5, 72.5);

  double pred_cv[72] = {212.06, 98.29, 87.33, 0.39, 0.94, 0.24, 67.01, 35.79, 31.20, 0.07, 0.37, 0.02, 15.29, 9.25, 5.50, 0.03, 0.07, 0.03, 1.26, 0.71, 0.72, 0.05, 0.00, 0.00, 13.75, 15.03, 21.78, 0.06, 0.25, 0.19, 7.61, 7.26, 12.84, 0.00, 0.12, 0.00, 0.20, 2.10, 3.30, 0.00, 0.03, 0.00, 1.06, 0.58, 0.05, 0.04, 0.00, 0.00, 0.84, 0.71, 3.46, 0.00, 0.00, 0.00, 0.08, 0.75, 2.20, 0.00, 0.00, 0.00, 0.00, 0.24, 1.23, 0.00, 0.00, 0.00, 0.00, 0.25, 0.03, 0.00, 0.00, 0.00};
  double syst[72] = {116.2, 35.77, 25.47,  0.44,  0.98,  0.26, 39.26, 14.96, 11.20,  0.09,  0.41,  0.05,  9.10,  3.94,  2.08,  0.04,  0.09,  0.04,  0.90,  0.39,  0.37,  0.08,  0.01,  0.02,  8.72,  7.08,  8.94,  0.09,  0.29,  0.22,  5.13,  3.76,  5.91,  0.03,  0.15,  0.02,  0.59,  1.16,  1.61,  0.03,  0.05,  0.02,  0.88,  0.40,  0.11,  0.06,  0.01,  0.02,  0.77,  0.54,  2.31,  0.04,  0.04,  0.04,  0.13,  0.58,  1.56,  0.03,  0.02,  0.02,  0.08,  0.23,  0.90,  0.03,  0.02,  0.02,  0.05,  0.23,  0.06,  0.02,  0.01,  0.02};
  double Nldp[72] = {4139., 3438., 2325., 45., 53., 17., 1350., 1248., 808., 12., 21., 2., 321., 317., 142., 3., 3., 1., 28., 28., 17., 1., 0., 0., 106., 230., 260., 2., 6., 4., 71., 128., 161., 0., 4., 0., 20., 50., 51., 0., 1., 0., 10., 12., 2., 1., 0., 0., 6., 14., 48., 0., 1., 0., 3., 15., 32., 0., 0., 0., 1., 7., 18., 0., 0., 0., 0., 4., 2., 0., 0., 0.};
  double Nldp_nonqcd[72] = {796.20, 296.90, 93.20, 27.30, 10.60,  6.00, 333.10, 148.30, 43.30,  8.70,  5.30,  3.10, 114.00, 58.20, 15.00,  2.30,  0.50,  0.30, 11.80,  9.00,  1.40,  0.00,  0.00,  0.00, 26.60, 29.80, 16.70,  0.70,  1.00,  0.40, 25.80, 34.10, 17.10,  0.40,  1.60,  0.20, 19.70, 20.80, 10.50,  0.80,  0.40,  0.10,  4.20,  4.40,  2.00,  0.00,  0.00,  0.00,  0.90,  3.00,  2.80,  0.70,  1.40,  1.30,  2.10,  3.80,  2.70,  0.40,  0.40,  0.70,  1.70,  3.20,  2.70,  0.60,  0.20,  0.30,  0.40,  0.90,  1.20,  0.10,  0.10,  0.10};
  double Rqcd[72] = {0.063, 0.031, 0.039, 0.022, 0.022, 0.020, 0.066, 0.033, 0.041, 0.012, 0.024, 0.033, 0.072, 0.036, 0.043, 0.011, 0.027, 0.027, 0.074, 0.037, 0.045, 0.046, 0.025, 0.022, 0.169, 0.075, 0.090, 0.040, 0.049, 0.049, 0.168, 0.077, 0.089, 0.062, 0.049, 0.045, 0.169, 0.072, 0.082, 0.062, 0.049, 0.045, 0.169, 0.077, 0.098, 0.040, 0.049, 0.045, 0.145, 0.065, 0.076, 0.053, 0.042, 0.038, 0.088, 0.068, 0.075, 0.053, 0.042, 0.038, 0.109, 0.063, 0.081, 0.053, 0.042, 0.038, 0.109, 0.069, 0.042, 0.053, 0.042, 0.038};

  const double alpha = 1 - 0.6827;
  
  Double_t stat_up[72];
  Double_t stat_down[72];
  Double_t syst_up[72];
  Double_t syst_down[72];
  Double_t full_err_up[72];
  Double_t full_err_down[72];
  Double_t x[72];
  Double_t xl[72];
  Double_t xh[72];

  
  
  for (unsigned int bin(0); bin<72; bin++) {
    hPredAllBins->SetBinContent(bin+1, pred_cv[bin]);
    hLDPStats->SetBinContent(bin+1, std::max(Nldp[bin]-Nldp_nonqcd[bin], 0.));
    hRQCD->SetBinContent(bin+1, Rqcd[bin]);   
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    // if (hLDPStats->GetBinContent(bin+1)>0){
      //      stat_up[bin]=pred_cv[bin]*sqrt(Nldp[bin])/Nldp[bin];
      //      stat_down[bin]=pred_cv[bin]*sqrt(Nldp[bin])/Nldp[bin];
    double N=std::max(Nldp[bin]-Nldp_nonqcd[bin], 0.);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) ;
    stat_up[bin]=(U-N)*Rqcd[bin];
    stat_down[bin]=(N-L)*Rqcd[bin];
    printf("Bin %d: Nldp=%d, Rqcd=%3.3f, estat=+%3.3f-%3.3f\n", bin+1, (int)hLDPStats->GetBinContent(bin+1), Rqcd[bin], stat_up[bin], stat_down[bin]);
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

  TGraphAsymmErrors* gFull = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, full_err_down, full_err_up);
  TGraphAsymmErrors* gStat = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, stat_down, stat_up);
  TGraphAsymmErrors* gSyst = new TGraphAsymmErrors(72, x, pred_cv, xl, xh, syst_down, syst_up);
  gFull->Write();
  gStat->Write();
  gSyst->Write();
  hPredAllBins->Write();
  hRQCD->Write();
  hLDPStats->Write();
  outfile->Close();
  
  return;
  
  
}

