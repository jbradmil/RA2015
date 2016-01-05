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

TH1D* MakeNJetsProjection(TH1D* hPredFineBins, int lowbin=0, int highbin=43) {
 
  Double_t njets_bins[6] = {3.5, 4.5, 5.5, 6.5, 8.5, 12.5};
  TString hname = Form("hPredNJets_%d_%d", lowbin+1, highbin+1);
  TH1D* hPredNJets = new TH1D(hname,";N_{jet} (p_{T} > 30 GeV);Events", 5, njets_bins);

  Double_t pred_cv[5];

  for (int nj(0); nj<5; nj++) {
    pred_cv[nj] = 0.;
    for (int subbin(0); subbin<44; subbin++) {
      if (subbin<lowbin||subbin>highbin) continue;
      pred_cv[nj]+=hPredFineBins->GetBinContent(nj*44+subbin+1);
    }
    //    printf("Total for nj %d: %3.2f\n", nj, pred_cv[nj]);
    hPredNJets->SetBinContent(nj+1, pred_cv[nj]);
    hPredNJets->SetBinError(nj+1, 0);
  }

  return hPredNJets;

}

TGraphAsymmErrors* MakeNJetsProjection(TH1D* hPredFineBins, TGraphAsymmErrors* gInFine, int lowbin=0, int highbin=43) {

  Double_t x[5] = {4.5, 5.5, 6.5, 7, 10.5};
  Double_t xl[5] = {0.5, 0.5, 0.5, 1., 1.5};
  Double_t xh[5] = {0.5, 0.5, 0.5, 1., 1.5};

  Double_t pred_cv[5];
  Double_t err_up[5];
  Double_t err_down[5];

  for (int nj(0); nj<5; nj++) {
    pred_cv[nj] = 0.;
    err_up[nj] = 0.;
    err_down[nj] = 0.;
    for (int subbin(0); subbin<44; subbin++) {
      if (subbin<lowbin||subbin>highbin) continue;
      pred_cv[nj]+=hPredFineBins->GetBinContent(nj*44+subbin+1);
      err_up[nj]+=pow(gInFine->GetErrorYhigh(nj*44+subbin), 2.);
      err_down[nj]+=pow(gInFine->GetErrorYlow(nj*44+subbin), 2.);
    }
    err_up[nj] = sqrt(err_up[nj]);
    err_down[nj] = sqrt(err_down[nj]);
  }

  TGraphAsymmErrors* gOut = new TGraphAsymmErrors(5, x, pred_cv, xl, xh,  err_down, err_up); 

  return gOut;

}

TH1D* MakeNJetsProjectionV2(TH1D* hPredFineBins, int nblow=0, int nbhigh=3, int boxlow=0, int boxhigh=10) {
 
  Double_t njets_bins[6] = {3.5, 4.5, 5.5, 6.5, 8.5, 12.5};
  TString hname = Form("hPredNJets_nb-%d-%d_box-%d-%d", nblow+1, nbhigh+1, boxlow+1, boxhigh+1);
  TH1D* hPredNJets = new TH1D(hname,";N_{jet} (p_{T} > 30 GeV);Events", 5, njets_bins);

  Double_t pred_cv[5];

  for (int nj(0); nj<5; nj++) {
    pred_cv[nj] = 0.;
    for (int nb(0); nb<4; nb++) { // nb bins
      if (nb<nblow||nb>nbhigh) continue;
      for (int ibox(0); ibox<11; ibox++) {
	if (ibox<boxlow||ibox>boxhigh) continue;
	pred_cv[nj]+=hPredFineBins->GetBinContent(nj*44+nb*11+ibox+1);
      }
    }
    // printf("Total for nj %d: %3.2f\n", nj, pred_cv[nj]);
    hPredNJets->SetBinContent(nj+1, pred_cv[nj]);
    hPredNJets->SetBinError(nj+1, 0);
  }

  return hPredNJets;

}

TGraphAsymmErrors* MakeNJetsProjectionV2(TH1D* hPredFineBins, TGraphAsymmErrors* gInFine, int nblow=0, int nbhigh=3, int boxlow=0, int boxhigh=10) {

  Double_t x[5] = {4.5, 5.5, 6.5, 7, 10.5};
  Double_t xl[5] = {0.5, 0.5, 0.5, 1., 1.5};
  Double_t xh[5] = {0.5, 0.5, 0.5, 1., 1.5};

  Double_t pred_cv[5];
  Double_t err_up[5];
  Double_t err_down[5];

  for (int nj(0); nj<5; nj++) {
    pred_cv[nj] = 0.;
    err_up[nj] = 0.;
    err_down[nj] = 0.;
    for (int nb(0); nb<4; nb++) { // nb bins
      if (nb<nblow||nb>nbhigh) continue;
      for (int ibox(0); ibox<11; ibox++) {
	if (ibox<boxlow||ibox>boxhigh) continue;
	pred_cv[nj]+=hPredFineBins->GetBinContent(nj*44+nb*11+ibox+1);
	err_up[nj]+=pow(gInFine->GetErrorYhigh(nj*44+nb*11+ibox), 2.);
	err_down[nj]+=pow(gInFine->GetErrorYlow(nj*44+nb*11+ibox), 2.);
      }
    }
    err_up[nj] = sqrt(err_up[nj]);
    err_down[nj] = sqrt(err_down[nj]);
  }

  TGraphAsymmErrors* gOut = new TGraphAsymmErrors(5, x, pred_cv, xl, xh,  err_down, err_up); 

  return gOut;

}

TH1D* MakeNBJetsProjection(TH1D* hPredFineBins, int lowbin=0, int highbin=10) {
 
  TString hname = Form("hPredNBJets_%d_%d", lowbin+1, highbin+1);
  TH1D* hPredNBJets = new TH1D(hname,";N_{b-jet} (p_{T} > 30 GeV);Events", 4, -0.5, 3.5);

  Double_t pred_cv[4];

  for (int nb(0); nb<4; nb++) {
    pred_cv[nb] = 0.;
    for (int nj(0); nj<5; nj++) {
      for (int subbin(0); subbin<11; subbin++) {
	if (subbin<lowbin||subbin>highbin) continue;
	pred_cv[nb]+=hPredFineBins->GetBinContent(nj*44+nb*11+subbin+1);
	//	printf("nb %d: Bin %d (%3.2f)\n", nb, nj*44+nb*11+subbin+1, hPredFineBins->GetBinContent(nj*44+nb*11+subbin+1));
      }
    }
    //    printf("Total for nb %d: %3.2f\n", nb, pred_cv[nb]);
    hPredNBJets->SetBinContent(nb+1, pred_cv[nb]);
    hPredNBJets->SetBinError(nb+1, 0);
  }
  return hPredNBJets;

}

 
TGraphAsymmErrors* MakeNBJetsProjection(TH1D* hPredFineBins, TGraphAsymmErrors* gInFine, int lowbin=0, int highbin=10) {

  Double_t x[4] = {0, 1, 2, 3};
  Double_t xl[4] = {0.5, 0.5, 0.5, 0.5};
  Double_t xh[4] = {0.5, 0.5, 0.5, 0.5};

  Double_t pred_cv[4];
  Double_t err_up[4];
  Double_t err_down[4];

  for (int nb(0); nb<4; nb++) {
    pred_cv[nb]=0.;
    err_up[nb] = 0.;
    err_down[nb] = 0.;
    for (int nj(0); nj<5; nj++) {
      for (int subbin(0); subbin<11; subbin++) {
	if (subbin<lowbin||subbin>highbin) continue;
	pred_cv[nb]+=hPredFineBins->GetBinContent(nj*44+nb*11+subbin+1);
	err_up[nb]+=pow(gInFine->GetErrorYhigh(nb*44+nb*11+subbin), 2.);
	err_down[nb]+=pow(gInFine->GetErrorYlow(nb*44+nb*11+subbin), 2.);
      }
    }
    err_up[nb] = sqrt(err_up[nb]);
    err_down[nb] = sqrt(err_down[nb]);
  }

  TGraphAsymmErrors* gOut = new TGraphAsymmErrors(4, x, pred_cv, xl, xh,  err_down, err_up);

  return gOut;

}

TH1D* MakeMHTProjection(TH1D* hPredFineBins, int njlow=0, int njhigh=4, int nblow=0, int nbhigh=3) {

  //  printf("MHT projection for NJ [%d,%d], NB [%d,%d]:\n",  njlow, njhigh, nblow, nbhigh);

  
  Double_t mht_bins[5] = {200, 300, 500, 750, 1050};
  TString hname = Form("hPredMHT_nj-%d-%d_nb-%d-%d", njlow+1, njhigh+1, nblow+1, nbhigh+1);
  TH1D* hPredMHT = new TH1D(hname,";H_{T}{}^{miss} [GeV];Events", 4, mht_bins);

  Double_t pred_cv[4];
  
  for (int imht(0); imht<4; imht++) { // mht bins
      pred_cv[imht] = 0.;
      for (int nj(0); nj<5; nj++) { // njets bins
	if (nj<njlow||nj>njhigh) continue;
	for (int nb(0); nb<4; nb++) { // nb bins
	  if (nb<nblow||nb>nbhigh) continue;
	  for (int subbin(0); subbin<3; subbin++) { // ht bins
	    if (imht==3&&subbin==2) continue;
	    pred_cv[imht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+subbin+1);
	    //	    printf("imht %d, nj %d, nb %d (Bin %d): %3.2f\n", imht, nj, nb, nj*44+nb*11+imht*3+subbin+1, hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+subbin+1));
	  } // ht bins
	} // nb bins
      } // njets bins
      //      printf("Total for mht %d: %3.2f\n", imht, pred_cv[imht]);
      hPredMHT->SetBinContent(imht+1, pred_cv[imht]);
      hPredMHT->SetBinError(imht+1, 0);
  }

  return hPredMHT;
  
}

TGraphAsymmErrors* MakeMHTProjection(TH1D* hPredAllBins, TH1D* hPredFineBins, TGraphAsymmErrors* gIn, TGraphAsymmErrors* gInFine, int njlow=0, int njhigh=4, int nblow=0, int nbhigh=3, bool verbose=false) {
 
  Double_t x[4] = {250, 400, 625, 900};
  Double_t xl[4] = {50, 100, 125, 150};
  Double_t xh[4] = {50, 100, 125, 150};

  Double_t pred_cv[4];
  Double_t err_up[4];
  Double_t err_down[4];

  for (int imht(0); imht<4; imht++) { // mht bins
      pred_cv[imht] = 0.;
      for (int nb(0); nb<4; nb++) { // nb bins
	if (nb<nblow||nb>nbhigh) continue;
	  for (int nj(0); nj<5; nj++) { // njets bins
	    if (nj<njlow||nj>njhigh) continue;
	    for (int subbin(0); subbin<3; subbin++) { // ht bins
	      if (imht==3&&subbin==1) continue;
	      pred_cv[imht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+subbin+1);
	      // printf("imht %d: Bin %d\n", imht, nj*44+nb*11+imht*3+subbin+1);
	    } // ht bins
	  } // njets bins
      } // nb bins
  }
  
  for (int imht(0); imht<4; imht++) {
    err_up[imht] = 0.;
    err_down[imht] = 0.;
    for (int nb(0); nb<4; nb++) { // nb bins
      if (nb<nblow||nb>nbhigh) continue;
      if (imht<2) { // fine binning
	for (int nj(0); nj<5; nj++) { // njets bins
	  if (nj<njlow||nj>njhigh) continue;
	  for (int subbin(0); subbin<3; subbin++) { // ht bins
	    err_up[imht]+=pow(gInFine->GetErrorYhigh(nj*44+nb*11+imht*3+subbin), 2.);
	    err_down[imht]+=pow(gInFine->GetErrorYlow(nj*44+nb*11+imht*3+subbin), 2.);
	    //if (verbose) printf("imht %d, Bin %d: %3.3f + %3.3f - %3.3f\n", imht, nj*44+nb*11+imht*3+subbin+1, hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+subbin+1), gInFine->GetErrorYhigh(nj*44+nb*11+imht*3+subbin), gInFine->GetErrorYlow(nj*44+nb*11+imht*3+subbin));
	    // printf("imht %d: Bin %d\n", imht, nj*44+nb*11+imht*3+subbin+1);
	  } // ht bins
	} // njets bins
      } else { // coarse binning
	  for (int nj(0); nj<3; nj++) { // njets bins
	    if (nj<(njlow-2)||nj>(njhigh-2)) continue;
	    for (int subbin(0); subbin<2; subbin++) { // ht bins
	      if (imht==3&&subbin==1) continue;
	      err_up[imht]+=pow(gIn->GetErrorYhigh(nj*24+nb*6+2*imht+subbin-1), 2.);
	      err_down[imht]+=pow(gIn->GetErrorYlow(nj*24+nb*6+2*imht+subbin-1), 2.);
	      // printf("imht %d: Bin %d\n", imht, nj*44+nb*11+imht*3+subbin+1);
	    } // ht bins
	  }
      }
    }
    err_up[imht] = sqrt(err_up[imht]);
    err_down[imht] = sqrt(err_down[imht]);
    //    if (verbose) printf("MHT bin %d: %3.3f + %3.3f - %3.3f\n", imht, pred_cv[imht], err_up[imht], err_down[imht]);
  }

  TGraphAsymmErrors* gOut = new TGraphAsymmErrors(4, x, pred_cv, xl, xh,  err_down, err_up);
  return gOut;

}

TH1D* MakeHTProjection(TH1D* hPredFineBins, int njlow=0, int njhigh=4, int nblow=0, int nbhigh=3, int mhtlow=0, int mhthigh=3) {
 
  Double_t ht_bins[4] = {500, 800, 1200, 1650};
  TString hname = Form("hPredHT_nj-%d-%d_nb-%d-%d_mht-%d-%d", njlow+1, njhigh+1, nblow+1, nbhigh+1, mhtlow+1, mhthigh+1);
  TH1D* hPredHT = new TH1D(hname,";H_{T} [GeV];Events", 3, ht_bins);

  Double_t pred_cv[4];

  for (int iht(0); iht<3; iht++) {
    pred_cv[iht] = 0.;
    for (int nj(0); nj<5; nj++) {
      if (nj<njlow||nj>njhigh) continue;
      for (int nb(0); nb<4; nb++) {
	if (nb<nblow||nb>nbhigh) continue;
	for (int imht(0); imht<4; imht++) {
	  if (imht<mhtlow||imht>mhthigh) continue;
	  if (imht==3) {
	    if (iht==0) continue;
	    pred_cv[iht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+iht);
	    //printf("iht %d: Bin %d (MHT %d)\n", iht, nj*44+nb*11+imht*3+iht, imht);
	  }
	  else {
	    pred_cv[iht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+iht+1);
	    //printf("iht %d: Bin %d (MHT %d)\n", iht, nj*44+nb*11+imht*3+iht+1, imht);
	  }
	}
      }
    }
    //    printf("Total for ht %d: %3.2f\n", iht, pred_cv[iht]);
    hPredHT->SetBinContent(iht+1, pred_cv[iht]);
    hPredHT->SetBinError(iht+1, 0);
  }
  return hPredHT;

}

TGraphAsymmErrors* MakeHTProjection(TH1D* hPredFineBins, TGraphAsymmErrors* gInFine, int njlow=0, int njhigh=4, int nblow=0, int nbhigh=3, int mhtlow=0, int mhthigh=3) {
 
  Double_t x[3] = {650, 1000, 1425};
  Double_t xl[3] = {150, 200, 225};
  Double_t xh[3] = {150, 200, 225};

  Double_t pred_cv[3];
  Double_t err_up[3]; 
  Double_t err_down[3];

  for (int iht(0); iht<3; iht++) {
    pred_cv[iht] = 0.;
    err_up[iht] = 0.; 
    err_down[iht] = 0.;
    for (int nj(0); nj<5; nj++) {
      if (nj<njlow||nj>njhigh) continue;
      for (int nb(0); nb<4; nb++) {
	if (nb<nblow||nb>nbhigh) continue;
	for (int imht(0); imht<4; imht++) {
	  if (imht<mhtlow||imht>mhthigh) continue;
	  if (imht==3) {
	    if (iht==0) continue;
	    pred_cv[iht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+iht);
	    err_up[iht]+=pow(gInFine->GetErrorYhigh(nj*44+nb*11+imht*3+iht-1), 2.); 
	    err_down[iht]+=pow(gInFine->GetErrorYlow(nj*44+nb*11+imht*3+iht-1), 2.);
	    //printf("iht %d: Bin %d (MHT %d)\n", iht, nj*44+nb*11+imht*3+iht, imht);
	  }
	  else {
	    pred_cv[iht]+=hPredFineBins->GetBinContent(nj*44+nb*11+imht*3+iht+1);
	    err_up[iht]+=pow(gInFine->GetErrorYhigh(nj*44+nb*11+imht*3+iht), 2.); 
	    err_down[iht]+=pow(gInFine->GetErrorYlow(nj*44+nb*11+imht*3+iht), 2.);
	    //printf("iht %d: Bin %d (MHT %d)\n", iht, nj*44+nb*11+imht*3+iht+1, imht);
	  }
	}
      }
    }
    err_up[iht] = sqrt(err_up[iht]); 
    err_down[iht] = sqrt(err_down[iht]);
  }
 
  TGraphAsymmErrors* gOut = new TGraphAsymmErrors(3, x, pred_cv, xl, xh,  err_down, err_up); 
  return gOut; 

}
