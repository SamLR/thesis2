#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"

TH1D* makeDatPlot(TString name, float* y, float* y_er, int i);
// TGraphErrors* makeDatPlot(TString name, float* y, float* y_er, int id);
TGraphErrors* makeSimPlot(TString name, float* y, float* y_er);

const int n_ch = 5;
const int n_files = 6;
const int n_sim_files = 4;

// All data taken from fits with the other tau parameter held constant
// copper data taken from output_txt/analysis_g4bl_sin_phase_loose_cu
// free data taken from output_txt/analysis_g4bl_sin_phase_loose_f

void fitted_mu_lifetimes(){
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);
  
  TString names[n_ch] = {TString("D1"), TString("D2"), TString("D3"), 
    TString("D4"), TString("D5")};
    
  float f [n_ch][n_files] = {
    {2159.65, 2106.15, 2345.72, 2285.74, 2287.15, 2251.11}, 
    {2061.01, 1936.76, 2116.00, 2099.17, 2222.17, 2224.48},
    {2268.50, 2335.78, 2293.76, 2261.60, 2336.52, 2348.64},
    {2230.74, 2255.66, 2254.51, 2227.35, 2323.17, 2367.91},
    {2210.65, 2117.36, 2226.03, 2234.24, 2237.82, 2271.48}
  };
  float f_er [n_ch][n_files] = {
    {  35.61,   99.55,   52.52,   47.37,   56.52,   83.19},
    {  56.00,  139.94,   79.02,   68.96,   85.81,  123.57},
    {  14.69,   45.14,   20.39,   17.95,   21.59,   31.56},
    {  21.22,   64.73,   28.71,   25.39,   33.12,   49.13},
    {  16.39,   55.55,   27.96,   24.89,   26.46,   40.21}
  };
  
  TH1D* f_hists[n_ch];
  TCanvas* f_can = new TCanvas("f","f", 1436,856);
  TLegend* f_l = new TLegend(0.90, 0.35, 0.99, 0.65, "Channel");
  for(int ch_id = 0; ch_id < n_ch; ++ch_id) {
    f_hists[ch_id] = makeDatPlot(names[ch_id], f[ch_id], f_er[ch_id], ch_id);
    f_hists[ch_id]->GetYaxis()->SetRangeUser(1750, 2450);
    f_l->AddEntry(f_hists[ch_id], names[ch_id].Data(), "LPFE");
    if (ch_id == 0) {
      f_hists[ch_id]->SetTitle("Measured free muon lifetime");
      f_hists[ch_id]->Draw("E1 P 9");
    } else {
      f_hists[ch_id]->Draw("E1 P 9 SAME");
    }
  }
  TF1* f_func = new TF1("f1","[0]", 0.0, static_cast<float>(n_ch+1));
  f_func->FixParameter(0, 2196.9811);
  f_func->SetLineStyle(9);
  f_func->SetLineWidth(2);
  f_func->SetLineColor(1);
  f_l->SetFillStyle(0); // no fill
  f_l->Draw();
  f_func->Draw("SAME");
  f_can->SaveAs("per_ch_free_lifetime.svg");
  f_can->SaveAs("per_ch_free_lifetime.eps");
  
  
    
  float c [n_ch][n_files] = {
    {197.88, 104.64, 112.30, 125.01, 129.64, 162.16},
    {  2.83,   3.57,   0.21,   2.08, 308.06, 267.55},
    { 45.03,  40.90,  39.82,  39.98,  43.15,  44.57},
    { 51.61,  55.66,  49.29,  44.86,  45.31,  47.69},
    {  0.02,  25.28,  63.99,  49.22,  29.37,  31.31}
  };
  float c_er [n_ch][n_files] = {
    {   17.99,    20.48,     9.35,    8.71, 14.83,  30.02},
    {    6.92, 11286.80, 15711.30, 4586.57, 82.28, 122.26},
    {    2.14,     5.39,     0.68,    2.51,  0.56,   0.87},
    {    3.41,     8.17,     3.76,    0.87,  1.16,   1.84},
    { 4596.97,     5.20,    11.81,    7.16,  0.74,   6.12}
  };
  
  TH1D* c_hists[n_ch];
  TCanvas* c_can = new TCanvas("c","c", 1436,856);
  TLegend* c_l = new TLegend(0.90, 0.35, 0.99, 0.65, "Channel");
  for(int ch_id = 0; ch_id < n_ch; ++ch_id) {
    TString n = TString(" ")+names[ch_id];
    c_hists[ch_id] = makeDatPlot(n, c[ch_id], c_er[ch_id], ch_id);
    c_hists[ch_id]->GetYaxis()->SetRangeUser(0, 350);
    c_l->AddEntry(c_hists[ch_id], names[ch_id].Data(), "LPFE");
    if (ch_id == 0) {
      c_hists[ch_id]->SetTitle("Measured muonic copper lifetime");
      c_hists[ch_id]->Draw("E1 P 9");
    } else {
      c_hists[ch_id]->Draw("E1 P 9 SAME");
    }
  }
  TF1* c_func = new TF1("f2","[0]", 0.0, static_cast<float>(n_ch+1));
  c_func->FixParameter(0, 163.5);
  c_func->SetLineStyle(9);
  c_func->SetLineWidth(2);
  c_func->SetLineColor(1);
  c_l->SetFillStyle(0); // no fill
  c_l->Draw();
  c_func->Draw("SAME");
  c_can->SaveAs("per_ch_copper_lifetime.svg");
  c_can->SaveAs("per_ch_copper_lifetime.eps");
 
  // Sim free                         0mm    0.5mm      1mm      5mm
  float f_sim    [n_sim_files] = {2188.98, 2169.02, 2194.00, 2161.80};
  float f_er_sim [n_sim_files] = {   9.21,   22.36,   22.80,   39.38};
  TString f_s_name = "Simulated free muon lifetime";
  TGraphErrors* g_f_sim = makeSimPlot(f_s_name,f_sim, f_er_sim);
  g_f_sim->GetYaxis()->SetRangeUser(2150, 2200); 
  TCanvas* f_s_can = new TCanvas("f_s", "f_s", 1436,856);
  g_f_sim->Draw("AP");
  f_s_can->SaveAs("simulated_free_lifetime.svg");
  f_s_can->SaveAs("simulated_free_lifetime.eps");
  
  // Sim copper                       0mm    0.5mm       1mm       5mm
  float c_sim    [n_sim_files] = {1282.23,    4.65,   492.55,  1455.40};
  float c_er_sim [n_sim_files] = {2404.62, 1635.78, 12319.21, 10219.90};
  TString c_s_name = "Simulated copper muon lifetime";
  TGraphErrors* g_c_sim = makeSimPlot(c_s_name,c_sim, c_er_sim);
  TCanvas* c_s_can = new TCanvas("c_s", "c_s", 1436,856);
  g_c_sim->GetYaxis()->SetRangeUser(-50, 1500); 
  g_c_sim->Draw("AP");
  c_s_can->SaveAs("simulated_copper_lifetime.svg");
  c_s_can->SaveAs("simulated_copper_lifetime.eps");
}


TGraphErrors* makeSimPlot(TString name, float* y, float* y_er) {
  float x    [n_sim_files] = {0.0,  0.5,  1.0,  5.0};
  float x_er [n_sim_files] = {0.05, 0.05, 0.05, 0.05};
  TGraphErrors* g = new TGraphErrors(n_sim_files, x, y, x_er, y_er);
  g->SetTitle(name.Data());
  g->GetXaxis()->SetTitle("Degrader thickness (mm)");
  g->GetXaxis()->SetRangeUser(-0.5, 5.5); 
  g->GetYaxis()->SetTitle("Lifetime (ns)");
  return g;
}


// 
// 
TH1D* makeDatPlot(TString name, float* y, float* y_er, int id) {
  
  // one for each data set
  int non_horrific_colors [n_ch] = {2, 3, 4, 6, 7};
  float offset [n_ch] = {0.00, -0.10, -0.05, 0.05, 0.10};
  
  TH1D* hist;
  if (id == 0) {
     hist = new TH1D(name.Data(), name.Data(), n_files, 0, static_cast<float>(n_files));
  } else {
    hist = new TH1D(name.Data(), name.Data(), 100*n_files, 0, static_cast<float>(n_files));
  }
  
  TString names [n_files] = {
    TString("448(0)"), TString("451(0.5)"), TString("452(0.5)"), TString("455(1)"), TString("458(5)"), TString("459(5)")
  };
  
  for(int i = 0; i < n_files; ++i) {
    // Find the centre of the bin & apply the per ch offset
    float x = static_cast<float>(i) + offset[id] + 0.5;
    hist->Fill(x, y[i]);
    if (id == 0) {
      int bin = i + 1;
      hist->GetXaxis()->SetBinLabel(bin, names[i].Data());
      hist->SetBinError(bin, y_er[i]);
    } else {
      int bin = hist->FindBin(x);
      hist->SetBinError(bin, y_er[i]);
    }
  }
  
  hist->SetLineColor(non_horrific_colors[id]);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(20); // circle
  hist->GetYaxis()->SetTitle("Lifetime (ns)");
  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetXaxis()->SetTicks("-");
  hist->GetXaxis()->SetTitle("Run (Degrader (mm))");
  hist->GetXaxis()->SetTitleOffset(-1.25);
  return hist;
}