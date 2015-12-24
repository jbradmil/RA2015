#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPave.h"
#include "TLatex.h"
#include "TString.h"
#include "TH1.h"
#include "THStack.h"
#include "TLine.h"
#include "TColor.h"

/* TColor uscb_blue(1000, 1/255.,57/255.,166/255.); */
/* TColor uscb_gold(1001, 255/255.,200/255.,47/255); */
/* TColor penn_red(1002, 149/255.,0/255.,26/255.); */
/* TColor ut_orange(1003, 191/255.,87/255.,0/255.); */
/* //TColor pantone(1003, 255/255.,88/255.,0/255.); */
/* TColor uo_green(1004, 0/255.,79/255.,39/255.); */
/* TColor tcu_purple(1005, 52/255.,42/255.,123/255.); */
/* TColor tar_heel_blue(1006, 86/255.,160/255.,211/255.); */
/* TColor cool_teal(1007, 2/255.,132/255.,130/255.); */
/* TColor signal_gold(1008, 204/255.,145/255.,38/255.); */

TColor uscb_blue(5000, 1/255.,57/255.,166/255.); 
 TColor penn_red(5002, 149/255.,0/255.,26/255.); 
 TColor uo_green(5004, 0/255.,79/255.,39/255.);


TColor ucsb_blue_lite(2000, 13/255.,67/255.,255/255.);
TColor carrot(2003, 230/255., 91/255., 20/255.);

TColor uo_green_lite(2004, 23/255.,160/255.,74/255.);
TColor tcu_purple_lite(2005, 86/255.,69/255.,203/255.);

TColor qcd_uscb_gold(2001, 255/255.,200/255.,47/255);
TColor znn_penn_red(2002, 255/255.,0/255.,43/255.);
TColor lost_lep_dusk_blue(2006, 105/255.,166/255., 202/255.);
TColor had_tau_grayed_jade(2007, 133/255.,189/255., 164/255.);
TColor heather(2008, 133/255.,99/255., 163/255.);
TColor c_iso(3002, 201/255., 191/255., 138/255.);

TColor uw_gold(3003, 145/255., 123/255., 76/255.);
TColor pretty_purple(3004, 185/255., 141/255., 226/255.);
TColor rose_quartz(3005, 240/255., 187/255., 183/255.);


void set_style(TH1 * h, const TString& p, const bool dashes=false) {


  if (p.Contains("t1bbbb")) {
    h->SetFillColor (2);
    h->SetMarkerColor(2);
  } else if(p.Contains("t1tttt")) {
    h->SetFillColor(6);
    h->SetMarkerColor(6);
  }  else if (p == "data_obs") {
    h->SetMarkerSize(1);
    h->SetLineWidth(1);
    h->SetMarkerStyle(20);
    h->SetLineColor(1);
  } else {
    if (dashes) h->SetFillStyle(3135);
    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    if (p == "znn") {
      h->SetFillColor(1002);
      h->SetMarkerColor(1002);
    } else if (p == "qcd") {
      h->SetFillColor (2001);
      h->SetMarkerColor(2001);
    } else if (p == "ttbar") {
      h->SetFillColor(1000);
      h->SetMarkerColor(1000);
    } else if (p == "wjets") {
      h->SetFillColor (1004);
      h->SetMarkerColor(1004);
    } else if (p == "single_top") {
      h->SetFillColor (1005);
      h->SetMarkerColor(1005);
    } else if (p == "other") {
      h->SetFillColor (1003);
      h->SetMarkerColor(1003);
    } else if (p== "bg") {
      h->SetFillColor(kAzure+1);
      h->SetMarkerColor(kAzure+1);
    } else if (p== "signal") {
      h->SetFillColor(1008);
      h->SetMarkerColor(1008);
    } else if (p== "pred") {
      h->SetFillColor(1007);
      h->SetMarkerColor(1007);
    }
  }
  return;
}

void set_style_lite(TH1 * h, const TString& p, const bool dashes=false) {


  if (p.Contains("t1bbbb")) {
    h->SetFillColor (2);
    h->SetMarkerColor(2);
  } else if(p.Contains("t1tttt")) {
    h->SetFillColor(6);
    h->SetMarkerColor(6);
  }  else if (p == "data_obs") {
    h->SetMarkerSize(1);
    h->SetLineWidth(1);
    h->SetMarkerStyle(20);
    h->SetLineColor(1);
  } else {
    if (dashes) h->SetFillStyle(3135);
    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    if (p == "znn") {
      h->SetFillColor(2002);
      h->SetMarkerColor(2002);
    } else if (p == "qcd") {
      h->SetFillColor (2001);
      h->SetMarkerColor(2001);
    } else if (p == "ttbar") {
      h->SetFillColor(2000);
      h->SetMarkerColor(2000);
    } else if (p == "wjets") {
      h->SetFillColor (2004);
      h->SetMarkerColor(2004);
    } else if (p == "single_top") {
      h->SetFillColor (2008);
      h->SetMarkerColor(2008);
    } else if (p == "other") {
      h->SetFillColor (2003);
      h->SetMarkerColor(2003);
    } else if (p== "bg") {
      h->SetFillColor(kAzure+1);
      h->SetMarkerColor(kAzure+1);
    } else if (p== "signal") {
      h->SetFillColor(1008);
      h->SetMarkerColor(1008);
    } else if (p== "pred") {
      h->SetFillColor(1007);
      h->SetMarkerColor(1007);
    }
    else if (p== "lost_lep") {
      h->SetFillColor(2006);
      h->SetMarkerColor(2007);
    }
    else if (p== "had_tau") {
      h->SetFillColor(2007);
      h->SetMarkerColor(2007);
    }
  }
  return;
}

void set_style(TLegend* leg, float text_size=0.045, bool transparent=false) {
  if (transparent) leg->SetFillStyle(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(text_size);
  leg->SetBorderSize(0);
  return;
}

void set_style(TLine* line) {
  line->SetLineColor(kGray+3);
  line->SetLineWidth(4);
  line->SetLineStyle(7);
}

void setTDRStyle(TStyle* tdrStyle) {
  //TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

/// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  
//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

}


void set_palette_style()
{
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void set_deep_sea_palette()
{
  const Int_t nRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 1.00, 0.09, 0.18, 0.09, 0.00 };
  // Double_t green[nRGBs] = { 0.97, 0.68, 0.39, 0.02, 0.01 };
  // Double_t blue[nRGBs]  = { 0.97, 0.79, 0.62, 0.39, 0.17 };
  Double_t green[nRGBs] = { 1.00, 0.97, 0.57, 0.20, 0.10 };
  Double_t blue[nRGBs]  = { 1.00, 0.97, 0.80, 0.57, 0.47 };
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void set_avocado_palette()
{
  const Int_t nRGBs = 9;
  const Int_t NCont = 255;
  Double_t stops[nRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[nRGBs]   = { 0./255.,  4./255., 12./255.,  30./255.,  52./255., 101./255., 142./255., 190./255., 237./255.};
  Double_t green[nRGBs] = { 0./255., 40./255., 86./255., 121./255., 140./255., 172./255., 187./255., 213./255., 240./255.};
  Double_t blue[nRGBs]  = { 0./255.,  9./255., 14./255.,  18./255.,  21./255.,  23./255.,  27./255.,  35./255., 101./255.};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void set_temperature_palette()
{
  const Int_t nRGBs = 10;
  const Int_t NCont = 255;
  Double_t stops[nRGBs] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5, 0.5625, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[nRGBs]   = {  30./255.,  55./255., 103./255., 147./255., 255./255., 255./255., 203./255., 188./255., 151./255., 105./255.};
  Double_t green[nRGBs] = {   0./255.,  65./255., 138./255., 182./255., 255./255., 255./255., 175./255., 121./255.,  53./255.,   9./255.};
  Double_t blue[nRGBs]  = { 191./255., 202./255., 212./255., 208./255., 255./255., 255./255., 140./255.,  97./255.,  57./255.,  30./255.};
  TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void Pal1()
{
  gStyle->SetPalette(51);
}

void Pal2()
{
  gStyle->SetPalette(53);
}
