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

void replot_rates_cu() {
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  const int n_runs = 6;
  const int n_bins = 4;
  // From fit systematic.ipnb taken the largest difference from 50ns bin/50 ns l_bound (for 75ns l_bound)
  // const float max_systematic_err = 381.26;
  const float max_systematic = 0.5585; // fraction of value
  TString filename("/Users/scook/code/MuSIC/offline_analysis_music5/MuSIC5_offline_analysis/MuSIC5_offline_analysis/scripts/images/analysis_g4bl_sin_phase_exec_d4_d5_l_bound50bin_16tight/out.root");
  
  TFile* in_file = new TFile(filename.Data(), "READ");
  TH1D* in_hist = (TH1D*) in_file->Get("Rate of muons decaying in cu");
  // TH1D* sim_hist = (TH1D*) in_file->Get("Simulated rate of muons decaying in f");
  
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
  // const double sim_mu_decay    = 7687.0;
  const double sim_mu_decay    = 10500.0;
  const double sim_mu_decay_er = TMath::Sqrt(sim_mu_decay);
  printf("Number of muon decays in detector: %.0f +/- %.0f\n", sim_mu_decay, sim_mu_decay_er);
  // const double sim_mu_decay_detected    = 5974.0;
  const double sim_mu_decay_detected    = 8362.0;
  const double sim_mu_decay_detected_er = TMath::Sqrt(sim_mu_decay_detected);
  printf("Number photon producing muon decays in detector: %.0f +/- %.0f\n", sim_mu_decay_detected, sim_mu_decay_detected_er);
  
  const double photon_acceptance    = sim_mu_decay_detected/sim_mu_decay;
  const double photon_acceptance_er = photon_acceptance * add_in_quad(sim_mu_decay, sim_mu_decay_er, sim_mu_decay_detected, sim_mu_decay_detected_er);
  printf("Photon acceptance: %.3f +/- %.3f\n", photon_acceptance, photon_acceptance_er);
  
  // Combine efficiency & acceptance
  const double scale = 1.0/(efficiency*photon_acceptance);
  const double scale_er = scale * add_in_quad(efficiency, efficiency_err, photon_acceptance, photon_acceptance_er);
  const double scale_sys_er = (scale_er/max_systematic);
  
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
    // d_bins_s    [i] = sim_hist->GetBinContent(i+1);
    // d_bin_s_ers [i] = sim_hist->GetBinError(i+1);
  }
  
  // destination bins
  double d_bins    [n_bins];
  double d_bin_ers [n_bins];
  // Errors with the added systematic uncertainty
  double d_bin_ers_incl_sys [n_bins];
  // Set of errors excluding the affect of the efficiency
  double d_bin_ers_exec_eff [n_bins];
  
  // Average things so we have one value / bin
  // Air 
  // origin index = 0, destination index = 0
  d_bins   [0] = scale * o_bins[0];
  d_bin_ers[0] = d_bins[0] * add_in_quad(scale, scale_er, o_bins[0], o_bin_ers[0]);
  d_bin_ers_incl_sys[0] = d_bins[0] * add_in_quad(scale, scale_sys_er, o_bins[0], o_bin_ers[0]);
  d_bin_ers_exec_eff[0] = d_bins[0] * add_in_quad(photon_acceptance, photon_acceptance_er, o_bins[0], o_bin_ers[0]);
  // 2x0.5 mm Al values
  // origin index = 1 & 2, destination index = 1
  const double bin_sum1    = (o_bins[1] + o_bins[2])/2.0;
  const double bin_sum_er1 = sqrt_sum_errors2(o_bin_ers[1], o_bin_ers[2])/2.0;
  d_bins   [1] = scale * bin_sum1;
  d_bin_ers[1] = d_bins[1] * add_in_quad(scale, scale_er, bin_sum1, bin_sum_er1);
  d_bin_ers_incl_sys[1] = d_bins[1] * add_in_quad(scale, scale_sys_er, o_bins[1], o_bin_ers[1]);
  d_bin_ers_exec_eff[1] = d_bins[1] * add_in_quad(photon_acceptance, photon_acceptance_er, bin_sum1, bin_sum_er1);
  // 1mm Al
  // origin index = 3, destination index = 2
  d_bins   [2] = scale * o_bins[3];
  d_bin_ers[2] = d_bins[2] * add_in_quad(scale, scale_er, o_bins[3], o_bin_ers[3]);
  d_bin_ers_incl_sys[2] = d_bins[2] * add_in_quad(scale, scale_sys_er, o_bins[2], o_bin_ers[2]);
  d_bin_ers_exec_eff[2] = d_bins[2] * add_in_quad(photon_acceptance, photon_acceptance_er, o_bins[3], o_bin_ers[3]);
  // 2x5mm Al
  // origin indexs = 4 & 5, destination index = 3
  const double bin_sum3    = (o_bins[4] + o_bins[5])/2.0;
  const double bin_sum_er3 = sqrt_sum_errors2(o_bin_ers[4], o_bin_ers[5])/2.0;
  d_bins   [3] = scale * bin_sum3;
  d_bin_ers[3] = d_bins[3] * add_in_quad(scale, scale_er, bin_sum3, bin_sum_er3);
  d_bin_ers_incl_sys[3] = d_bins[3] * add_in_quad(scale, scale_sys_er, bin_sum3, bin_sum_er3);
  d_bin_ers_exec_eff[3] = d_bins[3] * add_in_quad(photon_acceptance, photon_acceptance_er, bin_sum3, bin_sum_er3);
  
  const char bin_label_fmt [] = "%.0f #pm %.0f";
  // Actual values are 45, 50, 52 and 66 
  double bin_means_r[n_bins]  = {44.95, 49.95, 51.95, 65.95};
  double bin_means_s[n_bins] =  {45.05, 50.05, 52.05, 66.05}; 
  double bin_sigmas [n_bins] = { 0.00,  0.00,  0.00,  0.00};
  
  char name [] = "Corrected rate of muons decaying in copper";
  TGraphErrors* out_f_hist_exec_eff = new TGraphErrors(n_bins, bin_means_r, d_bins, bin_sigmas, d_bin_ers_exec_eff);
  out_f_hist_exec_eff->SetTitle(name);
  out_f_hist_exec_eff->GetXaxis()->SetTitle("Momentum (MeV/c)");
  out_f_hist_exec_eff->GetYaxis()->SetTitle("Muon rate (nA^{-1})");
  out_f_hist_exec_eff->GetYaxis()->SetTitleOffset(1.45);
  out_f_hist_exec_eff->GetYaxis()->SetRangeUser(0, 50000);
  
  TCanvas* can = new TCanvas("c1", "c1", 1436, 856);
  out_f_hist_exec_eff->SetLineColor(1);
  out_f_hist_exec_eff->SetFillColor(0);
  out_f_hist_exec_eff->SetLineWidth(1);
  out_f_hist_exec_eff->SetMarkerSize(0.75);
  out_f_hist_exec_eff->SetMarkerColor(1);
  out_f_hist_exec_eff->SetMarkerStyle(5);
  out_f_hist_exec_eff->GetYaxis()->SetRangeUser(0,2500);
  out_f_hist_exec_eff->Draw("A P");
  
  char name99 [] = "Errors due to MPPC efficiency";
  TGraphErrors* out_f_hist = new TGraphErrors(n_bins, bin_means_r, d_bins, bin_sigmas, d_bin_ers);
  out_f_hist->SetTitle(name99);
  out_f_hist->SetLineColor(2);
  out_f_hist->SetFillColor(0);
  out_f_hist->SetLineWidth(1);
  out_f_hist->SetMarkerSize(0.75);
  out_f_hist->SetMarkerColor(2);
  out_f_hist->SetMarkerStyle(1);
  out_f_hist->Draw("SAME P");

  char name999 [] = "Systematic errors";
  TGraphErrors* out_f_hist_incl_sys = new TGraphErrors(n_bins, bin_means_r, d_bins, bin_sigmas, d_bin_ers_incl_sys);
  out_f_hist_incl_sys->SetTitle(name999);
  out_f_hist_incl_sys->SetLineColor(4);
  out_f_hist_incl_sys->SetFillColor(0);
  out_f_hist_incl_sys->SetLineWidth(0);
  out_f_hist_incl_sys->SetMarkerSize(0.75);
  out_f_hist_incl_sys->SetMarkerColor(4);
  out_f_hist_incl_sys->SetMarkerStyle(1);
  out_f_hist_incl_sys->Draw("SAME P");
  
  TLegend* leg = can->BuildLegend(0.6, 0.8, 0.9, 0.9, "Muon momentum distribution");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  
  out_f_hist->Draw("SAME P");
  out_f_hist_exec_eff->Draw("SAME P");
  
  can->SaveAs("adjusted_muon_rates_cu.eps");
  can->SaveAs("adjusted_muon_rates_cu.svg");
  
  
  float means [n_bins] =    {41.03,      46.68, 50.29,      65.66};
  char* run_bins [n_bins] = {"448", "451, 452", "455", "458, 459"};
  printf("means[i], run_bins[i], d_bins[i], d_bin_ers_incl_sys[i], d_bin_ers[i], d_bin_ers_exec_eff[i]\n");
  for(int i = 0; i < n_bins; ++i) {
    printf("%.1f  &  %8s  &  %.0f  &  %.0f  &  %.0f  &  %.0f  \\\\\n", means[i], run_bins[i], d_bins[i], d_bin_ers_incl_sys[i], d_bin_ers[i], d_bin_ers_exec_eff[i]);
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
