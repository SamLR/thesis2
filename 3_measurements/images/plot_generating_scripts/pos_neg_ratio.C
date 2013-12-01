#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"


TH1D* makeHist(int id, float* f, float* f_er, float* c, float* c_er);
float getRelativeError(float f, float f_er, float c, float c_er);

const int n_ch = 5;
const int n_files = 6;
const TString ch_names [n_ch] = {
  TString("D1"), TString("D2"), TString("D3"), TString("D4"), TString("D5")
};
const TString file_names [n_files] = {
  TString("448 (0)"), TString("451 (0.5)"), TString("452 (0.5)"), TString("455 (1)"), TString("458 (5)"), TString("459 (5)")
};

void pos_neg_ratio(){
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  
  // Data take from output_txt/analysis_g4bl_sin_exec_d4_d5_tight/rates_and_integrals.txt
  float c_ints [n_ch][n_files] = {
    {5.57e+03, 2.99e+02, 2.11e+03, 2.69e+03, 2.22e+03, 1.17e+03},
    {1.55e+04, 1.55e+03, 7.17e+03, 9.10e+03, 8.12e+03, 3.75e+03},
    {6.01e+03, 8.69e+02, 3.40e+03, 4.45e+03, 2.40e+03, 1.12e+03},
    {1.98e+03, 2.63e+02, 1.30e+03, 1.98e+03, 7.37e+02, 3.98e+02},
    {2.55e+02, 7.48e+01, 2.66e+02, 4.31e+02, 1.96e+02, 8.21e+01}
  };
  
  float c_er_ints [n_ch][n_files] = {
    {3.2e+02, 8.2e+01, 1.7e+02, 1.9e+02, 1.5e+02, 1.0e+02},
    {3.9e+02, 1.1e+02, 2.4e+02, 2.6e+02, 2.0e+02, 1.4e+02},
    {2.1e+02, 6.7e+01, 1.4e+02, 1.5e+02, 1.1e+02, 7.5e+01},
    {1.2e+02, 3.5e+01, 7.7e+01, 8.9e+01, 6.0e+01, 5.0e+01},
    {8.0e+01, 2.6e+01, 5.3e+01, 6.0e+01, 4.3e+01, 2.8e+01}
  };
  
  float f_ints [n_ch][n_files] = {
    {1.89e+05, 1.39e+04, 5.98e+04, 7.40e+04, 4.81e+04, 2.13e+04},
    {2.66e+05, 2.45e+04, 1.08e+05, 1.33e+05, 7.72e+04, 3.49e+04},
    {8.52e+04, 7.93e+03, 3.66e+04, 4.48e+04, 2.45e+04, 1.11e+04},    
    {2.37e+04, 2.28e+03, 1.05e+04, 1.20e+04, 6.76e+03, 3.02e+03},
    {1.25e+04, 1.20e+03, 5.31e+03, 6.45e+03, 3.86e+03, 1.75e+03}
  };
  
  float f_er_ints [n_ch][n_files] = {
    {7.8e+02, 2.1e+02, 4.2e+02, 4.6e+02, 3.4e+02, 2.2e+02},
    {9.5e+02, 2.7e+02, 5.5e+02, 6.0e+02, 4.2e+02, 2.8e+02},
    {4.7e+02, 1.4e+02, 2.8e+02, 3.1e+02, 2.2e+02, 1.5e+02},
    {2.3e+02, 6.8e+01, 1.4e+02, 1.5e+02, 1.1e+02, 7.5e+01},
    {1.9e+02, 5.5e+01, 1.1e+02, 1.2e+02, 8.9e+01, 5.9e+01}
  };
  
  TH1D* hists [n_ch];
  TCanvas* can = new TCanvas("c", "c", 1436,856);
  TLegend* leg = new TLegend(0.90, 0.35, 0.99, 0.65, "Channel");
  leg->SetFillStyle(0);
  for(int ch = 0; ch < n_ch; ++ch) {
    hists[ch] = makeHist(ch, f_ints[ch], f_er_ints[ch], c_ints[ch], c_er_ints[ch]);
    leg->AddEntry(hists[ch], ch_names[ch]);
    if (ch == 0) {
      hists[ch]->SetTitle("Ratio of muonic-copper decays to free decays");
      hists[ch]->Draw("E1 P 9");
    } else {
      hists[ch]->Draw("E1 P 9 SAME");
    }
  }
  leg->Draw();
  can->SaveAs("Ratio.svg");
  can->SaveAs("Ratio.eps");
}

TH1D* makeHist(int id, float* f, float* f_er, float* c, float* c_er) {

  float offset [n_ch] = {0.00, -0.05, -0.10, 0.05, 0.10,};
  int non_horrific_colors [n_ch] = {2, 3, 4, 6, 7};
  float x_max = static_cast<float>(n_files);
  
  TH1D* hist = NULL;
  
  if (id == 0) {
    hist = new TH1D(ch_names[id].Data(), ch_names[id].Data(), n_files, 0, x_max);
  } else {
    hist = new TH1D(ch_names[id].Data(), ch_names[id].Data(), 100*n_files, 0, x_max);
  }
  
  
  hist->SetLineColor(non_horrific_colors[id]);
  hist->SetLineWidth(2);
  hist->SetMarkerColor(non_horrific_colors[id]);
  hist->SetMarkerStyle(20);
  hist->GetYaxis()->SetTitle("Ratio copper:free");
  hist->GetYaxis()->SetRangeUser(0.00, 0.2);
  hist->GetXaxis()->SetTicks("-");
  hist->GetXaxis()->SetTitle("Run (Degrader (mm))");
  hist->GetXaxis()->SetTitleOffset(-1.25);
  
  for(int i = 0; i < n_files; ++i) {
    float y    = c[i]/f[i];
    float y_er = y*getRelativeError(f[i], f_er[i], c[i], c_er[i]);
    float x = static_cast<float>(i) + offset[id] + 0.5;
    hist->Fill(x, y);
    if (id == 0) {
      int bin = i + 1;
      hist->GetXaxis()->SetBinLabel(bin, file_names[i].Data());
      hist->SetBinError(bin, y_er);
    } else {
      int bin = hist->FindBin(x);
      hist->SetBinError(bin, y_er);
    }
  }
  return hist;
}

float getRelativeError(float f, float f_er, float c, float c_er) {
  float fer_f2 = (f_er/f)*(f_er/f);
  float cer_c2 = (c_er/c)*(c_er/c);
  return TMath::Sqrt(fer_f2 + cer_c2);
}
