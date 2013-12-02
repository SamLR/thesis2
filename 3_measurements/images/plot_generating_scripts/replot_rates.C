#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"

float get_efficiency (const bool return_val);
float add_in_quad(const float a, const float a_er, const float b, const float b_er);
float sqrt_sum_errors2(const float a_er, const float b_er);

void replot_rates() {
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  const int n_runs = 6;
  const int n_bins = 4;
  // From fit systematic.ipnb taken the largest difference from 50ns bin/50 ns l_bound (for 75ns l_bound)
  // const float max_systematic_err = 381.26;
  TString filename("/Users/scook/code/MuSIC/offline_analysis_music5/MuSIC5_offline_analysis/MuSIC5_offline_analysis/scripts/images/analysis_g4bl_sin_phase_exec_d4_d5_l_bound50bin_16tight/out.root");
  
  TFile* in_file = new TFile(filename.Data(), "READ");
  TH1D* in_hist = (TH1D*) in_file->Get("Rate of muons decaying in f");
  TH1D* sim_hist = (TH1D*) in_file->Get("Simulated rate of muons decaying in f");
  
  // to scale by:
  // 1. Acceptance DONT as this adjusts back to total number of muons
  // 2. Efficiency 
  
  // Acceptances as calculated in section 4.3.1 measurements part
  // each has an error of +/- 0.001
  // const double acceptance [n_bins] = {0.047, 0.042, 0.039, 0.026};
  // const double acceptance_err      =  0.001;
  
  // Efficiency (get efficiency returns either the value or its error (the variance))
  // eff**2 as it is per MPPC and we have 2
  const double efficiency     = get_efficiency(true)*get_efficiency(true);
  // x +/- y, f = x**2 +/- xy 2**0.5
  const double efficiency_err = TMath::Sqrt(2)*get_efficiency(true)*get_efficiency(false);
  printf("Efficiency per MPPC: %.3f +/- %.3f\n", get_efficiency(true), get_efficiency(false));
  printf("Efficiency (2 MPPCs): %.3f +/- %.3f\n", efficiency, efficiency_err);
  // photon acceptance is the number of decays / number of times the threshold & 50ns cut are passed
  // From calc acceptance
  const double sim_mu_decay    = 7687.0;
  const double sim_mu_decay_er = TMath::Sqrt(sim_mu_decay);
  printf("Number of muon decays in detector: %.0f +/- %.0f\n", sim_mu_decay, sim_mu_decay_er);
  const double sim_mu_decay_detected    = 5974.0;
  const double sim_mu_decay_detected_er = TMath::Sqrt(sim_mu_decay_detected);
  printf("Number photon producing muon decays in detector: %.0f +/- %.0f\n", sim_mu_decay_detected, sim_mu_decay_detected_er);
  
  const double photon_acceptance    = sim_mu_decay_detected/sim_mu_decay;
  const double photon_acceptance_er = photon_acceptance * add_in_quad(sim_mu_decay, sim_mu_decay_er, sim_mu_decay_detected, sim_mu_decay_detected_er);
  printf("Photon acceptance: %.3f +/- %.3f\n", photon_acceptance, photon_acceptance_er);
  
  // Combine efficiency & acceptance
  double scale = 1.0/(efficiency*photon_acceptance);
  double scale_er = scale * add_in_quad(efficiency, efficiency_err, photon_acceptance, photon_acceptance_er);
  
  // original bin values
  float o_bins    [n_runs];
  float o_bin_ers [n_runs];

  for(int i = 0; i < n_runs; ++i) {
    o_bins   [i] = in_hist->GetBinContent(i+1);
    o_bin_ers[i] = in_hist->GetBinError(i+1);
  }
  
  // for the simulation values
  double d_bins_s    [n_bins];
  double d_bin_s_ers [n_bins];
  for(int i = 0; i < n_bins; ++i) {
    d_bins_s    [i] = sim_hist->GetBinContent(i+1);
    d_bin_s_ers [i] = sim_hist->GetBinError(i+1);
  }
  
  // destination bins
  double d_bins    [n_bins];
  double d_bin_ers [n_bins];
  // Set of errors excluding the affect of the efficiency
  double d_bin_ers_exec_eff [n_bins];
  
  // Average things so we have one value / bin
  // Air 
  // origin index = 0, destination index = 0
  d_bins   [0] = scale * o_bins[0];
  d_bin_ers[0] = d_bins[0] * add_in_quad(scale, scale_er, o_bins[0], o_bin_ers[0]);
  d_bin_ers_exec_eff[0] = d_bins[0] * add_in_quad(photon_acceptance, photon_acceptance_er, o_bins[0], o_bin_ers[0]);
  // 2x0.5 mm Al values
  // origin index = 1 & 2, destination index = 1
  const double bin_sum1    = (o_bins[1] + o_bins[2])/2.0;
  const double bin_sum_er1 = sqrt_sum_errors2(o_bin_ers[1], o_bin_ers[2])/2.0;
  d_bins   [1] = scale * bin_sum1;
  d_bin_ers[1] = d_bins[1] * add_in_quad(scale, scale_er, bin_sum1, bin_sum_er1);
  d_bin_ers_exec_eff[1] = d_bins[1] * add_in_quad(photon_acceptance, photon_acceptance_er, bin_sum1, bin_sum_er1);
  // 1mm Al
  // origin index = 3, destination index = 2
  d_bins   [2] = scale * o_bins[3];
  d_bin_ers[2] = d_bins[2] * add_in_quad(scale, scale_er, o_bins[3], o_bin_ers[3]);
  d_bin_ers_exec_eff[2] = d_bins[2] * add_in_quad(photon_acceptance, photon_acceptance_er, o_bins[3], o_bin_ers[3]);
  // 2x5mm Al
  // origin indexs = 4 & 5, destination index = 3
  const double bin_sum3    = (o_bins[4] + o_bins[5])/2.0;
  const double bin_sum_er3 = sqrt_sum_errors2(o_bin_ers[4], o_bin_ers[5])/2.0;
  d_bins   [3] = scale * bin_sum3;
  d_bin_ers[3] = d_bins[3] * add_in_quad(scale, scale_er, bin_sum3, bin_sum_er3);
  d_bin_ers_exec_eff[3] = d_bins[3] * add_in_quad(photon_acceptance, photon_acceptance_er, bin_sum3, bin_sum_er3);
  
  double bin_means_r[n_bins] = {41.13, 46.78, 50.39, 65.76}; // run values round & take 
  double bin_means_s[n_bins] = {40.93, 46.58, 50.19, 65.56}; // offsets at +/- 0.25
  double bin_sigmas [n_bins] = { 0.00,  0.00,  0.00,  0.00};
  
  char name [] = "Adjusted rate of freely decaying muons";
  TGraphErrors* out_f_hist_exec_eff = new TGraphErrors(n_bins, bin_means_r, d_bins, bin_sigmas, d_bin_ers_exec_eff);
  out_f_hist_exec_eff->SetTitle(name);
  out_f_hist_exec_eff->GetXaxis()->SetTitle("Momentum (MeV/c)");
  out_f_hist_exec_eff->GetYaxis()->SetTitle("Muon rate (nA^{-1})");
  out_f_hist_exec_eff->GetYaxis()->SetTitleOffset(1.45);
  out_f_hist_exec_eff->GetYaxis()->SetRangeUser(0, 50000);
  
  TCanvas* can = new TCanvas("c1", "c1", 1436, 856);
  out_f_hist_exec_eff->SetLineColor(4);
  out_f_hist_exec_eff->SetFillColor(0);
  out_f_hist_exec_eff->SetLineWidth(1);
  out_f_hist_exec_eff->SetMarkerSize(0.75);
  out_f_hist_exec_eff->SetMarkerColor(4);
  out_f_hist_exec_eff->SetMarkerStyle(8);
  out_f_hist_exec_eff->GetYaxis()->SetRangeUser(0,55000);
  out_f_hist_exec_eff->Draw("A P");
  
  char name2 [] = "Simulated rate of freely decaying muons";
  TGraphErrors* out_s_hist = new TGraphErrors(n_bins, bin_means_s, d_bins_s, bin_sigmas, d_bin_s_ers);
  out_s_hist->SetTitle(name2);
  out_s_hist->SetLineColor(2);
  out_s_hist->SetFillColor(0);
  out_s_hist->SetLineWidth(1);
  out_s_hist->SetMarkerSize(0.75);
  out_s_hist->SetMarkerColor(2);
  out_s_hist->SetMarkerStyle(8);
  out_s_hist->Draw("SAME P");
  
  TLegend* leg = can->BuildLegend(0.6, 0.7, 0.9, 0.9, "Muon momentum distribution");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  
  char name99 [] = "Adjusted rate of freely decaying muons (errors due to MPPC efficiency)";
  TGraphErrors* out_f_hist = new TGraphErrors(n_bins, bin_means_r, d_bins, bin_sigmas, d_bin_ers);
  out_f_hist->SetTitle(name99);
  out_f_hist->SetLineColor(4);
  out_f_hist->SetFillColor(0);
  out_f_hist->SetLineWidth(0);
  out_f_hist->SetMarkerSize(0.75);
  out_f_hist->SetMarkerColor(4);
  out_f_hist->SetMarkerStyle(8);
  out_f_hist->Draw("SAME P");
  
  can->SaveAs("adjusted_muon_rates.eps");
  can->SaveAs("adjusted_muon_rates.svg");
  
  
  float means [n_bins] =    {41.03,      46.68, 50.29,      65.66};
  char* run_bins [n_bins] = {"448", "451, 452", "455", "458, 459"};
  for(int i = 0; i < n_bins; ++i) {
    printf("%.1f  &  %8s  &  %.0f  &  %.0f  &  %.0f  \\\\\n", means[i], run_bins[i], d_bins[i], d_bin_ers[i], d_bin_ers_exec_eff[i]);
  }
}


float get_efficiency (const bool return_val) {
  // Values based on MuSIC 2 measurements
  const float vals [3] = {0.5265,  0.3786,  0.3889};
  float sum = 0;
  for(int i = 0; i < 3; ++i) {
    sum += vals[i];
  }
  const float ave = sum/3.0;
  if (return_val) {
    return ave;
  }
  
  // calculate the variance 
  float var = 0.0;
  for(int i = 0; i < 3; ++i) {
    var += (vals[i] - ave)*(vals[i] - ave);
  }
  return TMath::Sqrt(var/3.0);
}

float add_in_quad(const float a, const float a_er, const float b, const float b_er) {
  float aer_a2 = (a_er/a) * (a_er/a);
  float ber_b2 = (b_er/b) * (b_er/b);
  return TMath::Sqrt(aer_a2 + ber_b2);
}

float sqrt_sum_errors2(const float a_er, const float b_er) {
  return TMath::Sqrt(a_er*a_er + b_er*b_er);
}
