
#include <iostream>
#include <vector>
#include <string>
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

#ifndef __CINT__
#include "RooStats/NumberCountingUtils.h"
#endif

using namespace std;

TString plotdir = "plots/";
//const double predSF = 1.63947463;

FILE *  texfile;

std::string njets_cuts[3] = {"4-6", "7-8", "9+"};
std::string nbjets_cuts[72] = {"0", "0", "0", "0", "0", "0",
			      "1", "1", "1", "1", "1", "1",
			      "2", "2", "2", "2", "2", "2",
			       "3+", "3+", "3+", "3+", "3+", "3+",
			       "0", "0", "0", "0", "0", "0",
			      "1", "1", "1", "1", "1", "1",
			      "2", "2", "2", "2", "2", "2",
			       "3+", "3+", "3+", "3+", "3+", "3+",
			       "0", "0", "0", "0", "0", "0",
			      "1", "1", "1", "1", "1", "1",
			      "2", "2", "2", "2", "2", "2",
			       "3+", "3+", "3+", "3+", "3+", "3+"};
std::string mht_cuts[6] = {"200-500", "200-500", "200-500", "500-750", "500-750", "750+"};
std::string ht_cuts[6] = {"500-800", "800-1200", "1200+", "500-800", "1200+", "800+"};

void PrintTable1Header() {

  fprintf(texfile, "\\begin{sidewaystable}\n");
  fprintf(texfile, "\\renewcommand{\\arraystretch}{1.25}\n");
  fprintf(texfile, "\\centering\n");
  fprintf(texfile, "\\caption{Observed number of events and pre-fit background predictions in the $4\\leq\\njets\\leq6$ search bins.}\n");
  fprintf(texfile, "\\label{tab:pre-fit-results-nj1}\n");
  fprintf(texfile, "\\begin{tabular}{ |c|c|c|c|c||c|c|c|c||c|c| }\n");
  fprintf(texfile, "\\hline\n");
  fprintf(texfile, "Bin & \\MHT [GeV] & \\HT [GeV] & \\njets & \\nbjets & Lost-$e/\\mu$ & $\\tau\\rightarrow\\mathrm{had}$ & $Z\\rightarrow\\nu\\bar{\\nu}$ & QCD & Total Pred. & Obs. \\\\ \\hline\n");

}

void PrintTable2Header() {

  fprintf(texfile, "\\begin{sidewaystable}\n");
  fprintf(texfile, "\\renewcommand{\\arraystretch}{1.25}\n");
  fprintf(texfile, "\\centering\n");
  fprintf(texfile, "\\caption{Observed number of events and pre-fit background predictions in the $7\\leq\\njets\\leq8$ search bins.}\n");
  fprintf(texfile, "\\label{tab:pre-fit-results-nj2}\n");
  fprintf(texfile, "\\begin{tabular}{ |c|c|c|c|c||c|c|c|c||c|c| }\n");
  fprintf(texfile, "\\hline\n");
  fprintf(texfile, "Bin & \\MHT [GeV] & \\HT [GeV] & \\njets & \\nbjets & Lost-$e/\\mu$ & $\\tau\\rightarrow\\mathrm{had}$ & $Z\\rightarrow\\nu\\bar{\\nu}$ & QCD & Total Pred. & Obs. \\\\ \\hline\n");

}

void PrintTable3Header() {

  fprintf(texfile, "\\begin{sidewaystable}\n");
  fprintf(texfile, "\\renewcommand{\\arraystretch}{1.25}\n");
  fprintf(texfile, "\\centering\n");
  fprintf(texfile, "\\caption{Observed number of events and pre-fit background predictions in the $\\njets\\qeq9$ search bins.}\n");
  fprintf(texfile, "\\label{tab:pre-fit-results-nj2}\n");
  fprintf(texfile, "\\begin{tabular}{ |c|c|c|c|c||c|c|c|c||c|c| }\n");
  fprintf(texfile, "\\hline\n");
  fprintf(texfile, "Bin & \\MHT [GeV] & \\HT [GeV] & \\njets & \\nbjets & Lost-$e/\\mu$ & $\\tau\\rightarrow\\mathrm{had}$ & $Z\\rightarrow\\nu\\bar{\\nu}$ & QCD & Total Pred. & Obs. \\\\ \\hline\n");

}

void PrintTableTrailer() {

  fprintf(texfile, "\\end{tabular}\n");
  fprintf(texfile, "\\end{sidewaystable}\n");
  
}


void PrintResults() {

  texfile = fopen ("pre-fit-tables.tex","w");

  TFile* f_lostlep = new TFile("bg_hists/lostlep_hists.root", "read");
  TFile* f_hadtau = new TFile("bg_hists/hadtau_hists.root", "read");
  TFile* f_qcd = new TFile("bg_hists/qcd_hists.root", "read");
  TFile* f_znn = new TFile("bg_hists/znn_hists.root", "read");
  TFile* f_data_obs = new TFile("data_hists/data_hists.root", "read");
  TH1D* hdata_obs = (TH1D*) f_data_obs->Get("hObsAllBins");

  TGraphAsymmErrors* glostlepstat = (TGraphAsymmErrors*) f_lostlep->Get("Graph;2");
  TGraphAsymmErrors* ghadtaustat = (TGraphAsymmErrors*) f_hadtau->Get("Graph;2");
  TGraphAsymmErrors* glostlepsyst = (TGraphAsymmErrors*) f_lostlep->Get("Graph;3");
  TGraphAsymmErrors* ghadtausyst = (TGraphAsymmErrors*) f_hadtau->Get("Graph;3");
  TGraphAsymmErrors* gqcdstat = (TGraphAsymmErrors*) f_qcd->Get("Graph;2");
  TGraphAsymmErrors* gznnstat = (TGraphAsymmErrors*) f_znn->Get("Graph;2");
  TGraphAsymmErrors* gqcdsyst = (TGraphAsymmErrors*) f_qcd->Get("Graph;3");
  TGraphAsymmErrors* gznnsyst = (TGraphAsymmErrors*) f_znn->Get("Graph;3");


  
  // special values for bins with 0 expected events
  // TH1D* hLL0 = (TH1D*) f_lostlep->Get("hRLL");
  // TH1D* hTau0 = (TH1D*) f_lostlep->Get("hRTau");
  // TH1D* hZnn0 = (TH1D*) f_lostlep->Get("hRZ0");
  // TH1D* hQCD0 = (TH1D*) f_lostlep->Get("hRQCD");

  
  Double_t pred_cv[72];
  Double_t full_stat_up[72];
  Double_t full_stat_down[72];
  Double_t full_syst_up[72];
  Double_t full_syst_down[72];
  Double_t x[72];
  Double_t xl[72];
  Double_t xh[72];

  PrintTable1Header();
  
  for (unsigned int bin(0); bin<24; bin++) {
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    pred_cv[bin]=(glostlepstat->Eval(bin+1)+ghadtaustat->Eval(bin+1)+gqcdstat->Eval(bin+1)+gznnstat->Eval(bin+1));
    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    full_stat_up[bin] = (sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.)));
    full_stat_down[bin] = (sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.)));
    // printf("Bin %d, LL stat err: + %3.3f - %3.3f\n",bin+1, glostlepstat->GetErrorYhigh(bin), glostlepstat->GetErrorYlow(bin));
    // printf("Bin %d, LL syst err: + %3.3f - %3.3f\n",bin+1, glostlepsyst->GetErrorYhigh(bin), glostlepsyst->GetErrorYlow(bin));
    // printf("Bin %d, wtop err: + %3.3f - %3.3f\n",bin+1,wtop_err_up, wtop_err_down);
    full_syst_up[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.)));
    full_syst_down[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.)));
    double pull(0.);
    // if (pred_cv[bin]-hdata_obs->GetBinContent(bin+1)!=0) {
    if (pred_cv[bin] < hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_up[bin]*full_stat_up[bin] + full_syst_up[bin]*full_syst_up[bin] + pred_cv[bin]);
      // printf("Obs: %d, Pred: %3.3f, Err: %3.3f, Pull: %3.2f\n", (int)hdata_obs->GetBinContent(bin+1), pred_cv[bin], sqrt(full_stat_up[bin]*full_stat_up[bin] + full_syst_up[bin]*full_syst_up[bin]), pull);
    }
    else if (pred_cv[bin] > hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin] + pred_cv[bin]);
      //  printf("Obs: %d, Pred: %3.3f, Err: %3.3f, Pull: %3.2f\n", (int)hdata_obs->GetBinContent(bin+1), pred_cv[bin], sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin]), pull);
    }
      // }

    // double ZBi(0.);
    // if (hdata_obs->GetBinContent(bin+1)>pred_cv[bin]) {
    //   if (pred_cv[bin]>0) {
    // 	ZBi = RooStats::NumberCountingUtils::BinomialObsZ(hdata_obs->GetBinContent(bin+1), pred_cv[bin], sqrt(full_stat_up[bin]*full_stat_up[bin]+full_syst_up[bin]*full_syst_up[bin])/pred_cv[bin]);
    //   } else {
    // 	double nom_0_val = hLL0->GetBinContent(bin+1)+hTau0->GetBinContent(bin+1)+hZnn0->GetBinContent(bin+1)+hQCD0->GetBinContent(bin+1);
    // 	ZBi = RooStats::NumberCountingUtils::BinomialObsZ(hdata_obs->GetBinContent(bin+1), nom_0_val, sqrt(full_stat_up[bin]*full_stat_up[bin]+full_syst_up[bin]*full_syst_up[bin])/nom_0_val);
    //   }
    // }
    fprintf(texfile, "%d & %s & %s & %s & %s & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & %d \\\\ \\hline\n", bin+1,
	    mht_cuts[bin%6].c_str(), ht_cuts[bin%6].c_str(), njets_cuts[0].c_str(), nbjets_cuts[bin].c_str(),
	   glostlepstat->Eval(bin+1), glostlepstat->GetErrorYhigh(bin), glostlepsyst->GetErrorYhigh(bin), glostlepstat->GetErrorYlow(bin), glostlepsyst->GetErrorYlow(bin),
	   ghadtaustat->Eval(bin+1), ghadtaustat->GetErrorYhigh(bin), ghadtausyst->GetErrorYhigh(bin), ghadtaustat->GetErrorYlow(bin), ghadtausyst->GetErrorYlow(bin),
	   gznnstat->Eval(bin+1), gznnstat->GetErrorYhigh(bin), gznnsyst->GetErrorYhigh(bin), gznnstat->GetErrorYlow(bin), gznnsyst->GetErrorYlow(bin),
	   gqcdstat->Eval(bin+1), gqcdstat->GetErrorYhigh(bin), gqcdsyst->GetErrorYhigh(bin), gqcdstat->GetErrorYlow(bin), gqcdsyst->GetErrorYlow(bin),
	    pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin], (int)hdata_obs->GetBinContent(bin+1)	   );
    //if (pred_cv[bin]==0) printf("%d & $%f^{+%f+%f}_{-%f-%f}$ & %d & %f \\\\ \\hline\n", bin+1, pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin], (int)hdata_obs->GetBinContent(bin+1), pull);
   printf("%d | %f + %f + %f - %f - %f\n", bin+1,
	 pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin]);
 }
  PrintTableTrailer();


  PrintTable2Header();
  
  for (unsigned int bin(24); bin<48; bin++) {
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    pred_cv[bin]=(glostlepstat->Eval(bin+1)+ghadtaustat->Eval(bin+1)+gqcdstat->Eval(bin+1)+gznnstat->Eval(bin+1));
    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    full_stat_up[bin] = (sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.)));
    full_stat_down[bin] = (sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.)));
    full_syst_up[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.)));
    full_syst_down[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.)));
    double pull(0.);
    if (pred_cv[bin] < hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_up[bin]*full_stat_up[bin] + full_syst_up[bin]*full_syst_up[bin] + pred_cv[bin]);
    }
    else if (pred_cv[bin] > hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin] + pred_cv[bin]);
      // printf("Obs: %d, Pred: %3.3f, Err: %3.3f, Pull: %3.2f\n", (int)hdata_obs->GetBinContent(bin+1), pred_cv[bin], sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin]), pull);
    }

    fprintf(texfile, "%d & %s & %s & %s & %s & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & %d \\\\ \\hline\n", bin+1,
	    mht_cuts[(bin-24)%6].c_str(), ht_cuts[(bin-24)%6].c_str(), njets_cuts[1].c_str(), nbjets_cuts[bin].c_str(),
	   glostlepstat->Eval(bin+1), glostlepstat->GetErrorYhigh(bin), glostlepsyst->GetErrorYhigh(bin), glostlepstat->GetErrorYlow(bin), glostlepsyst->GetErrorYlow(bin),
	   ghadtaustat->Eval(bin+1), ghadtaustat->GetErrorYhigh(bin), ghadtausyst->GetErrorYhigh(bin), ghadtaustat->GetErrorYlow(bin), ghadtausyst->GetErrorYlow(bin),
	   gznnstat->Eval(bin+1), gznnstat->GetErrorYhigh(bin), gznnsyst->GetErrorYhigh(bin), gznnstat->GetErrorYlow(bin), gznnsyst->GetErrorYlow(bin),
	   gqcdstat->Eval(bin+1), gqcdstat->GetErrorYhigh(bin), gqcdsyst->GetErrorYhigh(bin), gqcdstat->GetErrorYlow(bin), gqcdsyst->GetErrorYlow(bin),
	    pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin], (int)hdata_obs->GetBinContent(bin+1));
  printf("%d | %f + %f + %f - %f - %f\n", bin+1,
	 pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin]);
  }
  PrintTableTrailer();

  PrintTable3Header();
  
  for (unsigned int bin(48); bin<72; bin++) {
    x[bin] = bin+1;
    xl[bin]=0.5;
    xh[bin]=0.5;
    pred_cv[bin]=(glostlepstat->Eval(bin+1)+ghadtaustat->Eval(bin+1)+gqcdstat->Eval(bin+1)+gznnstat->Eval(bin+1));
    double wtop_stat_up = sqrt(pow(glostlepstat->GetErrorYhigh(bin)+ghadtaustat->GetErrorYhigh(bin),2.));
    double wtop_stat_down = sqrt(pow(glostlepstat->GetErrorYlow(bin)+ghadtaustat->GetErrorYlow(bin),2.));
    full_stat_up[bin] = (sqrt(pow(wtop_stat_up,2.)+pow(gqcdstat->GetErrorYhigh(bin),2.)+pow(gznnstat->GetErrorYhigh(bin),2.)));
    full_stat_down[bin] = (sqrt(pow(wtop_stat_down,2.)+pow(gqcdstat->GetErrorYlow(bin),2.)+pow(gznnstat->GetErrorYlow(bin),2.)));
    full_syst_up[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYhigh(bin),2.)+pow(gznnsyst->GetErrorYhigh(bin),2.)));
    full_syst_down[bin] = (sqrt(pow(glostlepsyst->GetErrorYlow(bin),2.)+pow(ghadtausyst->GetErrorYlow(bin),2.)+pow(gqcdsyst->GetErrorYlow(bin),2.)+pow(gznnsyst->GetErrorYlow(bin),2.)));
    double pull(0.);
    if (pred_cv[bin] < hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_up[bin]*full_stat_up[bin] + full_syst_up[bin]*full_syst_up[bin] + pred_cv[bin]);
    }
    else if (pred_cv[bin] > hdata_obs->GetBinContent(bin+1)) {
      pull = (hdata_obs->GetBinContent(bin+1)-pred_cv[bin])/sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin] + pred_cv[bin]);
      // printf("Obs: %d, Pred: %3.3f, Err: %3.3f, Pull: %3.2f\n", (int)hdata_obs->GetBinContent(bin+1), pred_cv[bin], sqrt(full_stat_down[bin]*full_stat_down[bin] + full_syst_down[bin]*full_syst_down[bin]), pull);
    }

    fprintf(texfile, "%d & %s & %s & %s & %s & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & $%3.2f^{+%3.2f+%3.2f}_{-%3.2f-%3.2f}$ & %d \\\\ \\hline\n", bin+1,
	    mht_cuts[(bin-48)%6].c_str(), ht_cuts[(bin-48)%6].c_str(), njets_cuts[2].c_str(), nbjets_cuts[bin].c_str(),
	   glostlepstat->Eval(bin+1), glostlepstat->GetErrorYhigh(bin), glostlepsyst->GetErrorYhigh(bin), glostlepstat->GetErrorYlow(bin), glostlepsyst->GetErrorYlow(bin),
	   ghadtaustat->Eval(bin+1), ghadtaustat->GetErrorYhigh(bin), ghadtausyst->GetErrorYhigh(bin), ghadtaustat->GetErrorYlow(bin), ghadtausyst->GetErrorYlow(bin),
	   gznnstat->Eval(bin+1), gznnstat->GetErrorYhigh(bin), gznnsyst->GetErrorYhigh(bin), gznnstat->GetErrorYlow(bin), gznnsyst->GetErrorYlow(bin),
	   gqcdstat->Eval(bin+1), gqcdstat->GetErrorYhigh(bin), gqcdsyst->GetErrorYhigh(bin), gqcdstat->GetErrorYlow(bin), gqcdsyst->GetErrorYlow(bin),
	    pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin], (int)hdata_obs->GetBinContent(bin+1));
  printf("%d | %f + %f + %f - %f - %f\n", bin+1,
	 pred_cv[bin], full_stat_up[bin], full_syst_up[bin], full_stat_down[bin], full_syst_down[bin]);
  }


  PrintTableTrailer();
  
  return;
  
}

