#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TStyle.h"

#define MAX_HITS 5000

struct truth_branch {
  int nhits;
  double in_px;
  double in_py;
  double in_pz;
  
  int counter  [MAX_HITS];
  int pdgid    [MAX_HITS];
  int parentid [MAX_HITS];
  int trkid    [MAX_HITS];
  bool first_step [MAX_HITS];
  
  double px  [MAX_HITS];
  double py  [MAX_HITS];
  double pz  [MAX_HITS];
  double tof [MAX_HITS];
};

struct mppc_branch {
  int mppc_hits;
  double mppc_x    [MAX_HITS];
  double mppc_y    [MAX_HITS];
  double mppc_z    [MAX_HITS];
  double mppc_time [MAX_HITS];
};

double get_momentum(truth_branch &branch);
void set_truth_branches(TTree* tree, truth_branch &branch);
void set_mppc_branches(TTree* tree, mppc_branch &branch);
int get_mppc(double y, double z);
void move_stats_box(TCanvas* can, TH1D* hist, float x1, float y1, float x2, float y2);

void calc_acceptance() {
  // Want to find how many muons stop that we don't detect photons for
  // requirement:
  // 1. muon in upstream detector
  // 2. electron in downstream
  // 3. muon produces > X photons (for both U and D)
  // 4. dt(D-U) > 50
  
  gStyle->SetOptStat(2211);
  
  // This is a quick check that the get_mppc calculation works
  printf("point (832.83839, 75.019284, 3923.2438) is MPPC %d\n",get_mppc(75.019284, 3923.2438));
  
  // const int n_ch = 5;
  const int n_d_ch = 5;
  const int n_u_ch = 8;
  
  // Thresholds from Akira's slides
  const int u_thrs = 8;
  const int d_thrs = 10;
  
  TFile* in_file = new TFile("/Users/scook/code/MuSIC/simulation/MuSIC_5_detector_sim/MuSIC5/output/output_from_hep_batch/optical_processes.root", "READ");
  
  TTree* truth = (TTree*) in_file->Get("truth");
  truth_branch t_branch;
  set_truth_branches(truth, t_branch);
  
  TTree* mppc  = (TTree*) in_file->Get("mppc");
  mppc_branch m_branch;
  set_mppc_branches(mppc, m_branch);
  
  char name [] = "Decay muon momentum spectrum";
  TH1D* momentum = new TH1D(name, name, 100, 0, 300);

  char name2 [] = "Decay muon momentum spectrum with #Delta t > 50ns";
  TH1D* momentum50ns = new TH1D(name2, name2, 100, 0, 300);
  
  char name3 [] = "Detected photon counts at all scintillators for muon decays";
  TH1D* photon_counts = new TH1D(name3, name3, 1500, 0, 1500);
  
  char name4 [] = "Detected photon counts at all scintillators (no cuts)";
  TH1D* gross_photon_counts = new TH1D(name4, name4, 3000, 0, 3000);
  
  char name5 [] = "Detected photon counts at all scintillators for muon decays with #Delta t > 50ns";
  TH1D* photon_counts50ns = new TH1D(name5, name5, 1500, 0, 1500);
  
  char name6 [] = "Upstream photon counts passing trigger for muon decays with #Delta t > 50ns";
  TH1D* photon_counts50ns_uthrs = new TH1D(name6, name6, 1500, 0, 1500);
  
  char name7 [] = "Downstream photon counts passing trigger for muon decays with #Delta t > 50ns";
  TH1D* photon_counts50ns_dthrs = new TH1D(name7, name7, 1500, 0, 1500);
  
  int decay_muons      = 0; // number of muons we see at U with a daughter electron at D
  int decay_muons_gt50 = 0; // decay muons with dt > 50 ns
  int muon_photon_events      = 0;  // muon decays with photons @ U & D
  int muon_photon_events_thrs = 0;  // muon decays photons > thrs @ U & D
  
  const int n_entries = truth->GetEntries();
  for(int entry = 0; entry < n_entries; ++entry) {
    
    truth->GetEntry(entry);
    mppc->GetEntry(entry);
    
    int  mu_trkid       = 0;     // Store the mu id to check against e's parent
    bool u_stream_mu    = false; // Have we seen an upstream mu?
    bool d_stream_e     = false; // Have we seen a downstream e?
    bool mu_decay_gt_50 = false; // Is the dt between U&D greater than 50 ns?
    double u_stream_t   = 0.0;   // When did we see the mu?
    double d_stream_t   = 0.0;   // When did we see the e?
    double mu_momentum  = 0.0;   // The muon's initial momentum
    
    // Calculate the number of true decays
    for(int hit = 0; hit < t_branch.nhits; ++hit) {
      if (t_branch.first_step[hit]) {
        if (t_branch.counter[hit]==1 &&
            TMath::Abs(t_branch.pdgid[hit])==13) {
              
          // U-stream muon
          u_stream_mu = true;
          u_stream_t = t_branch.tof[hit];
          mu_momentum = get_momentum(t_branch);
          mu_trkid = t_branch.trkid[hit];
          
        } else if (t_branch.counter[hit]==3 &&
                   TMath::Abs(t_branch.pdgid[hit])==11) {// &&
                   // t_branch.parentid[hit] == mu_trkid) {
                     
          // D-stream electron
          d_stream_e = true;
          d_stream_t = t_branch.tof[hit];
          
          if (u_stream_mu && d_stream_e) {
            
            // Decay!
            ++decay_muons;
            momentum->Fill(mu_momentum);
            
            if ((t_branch.tof[hit] - u_stream_t) > 50) {
              
              // Decay we'd detect!
              ++decay_muons_gt50;
              momentum50ns->Fill(mu_momentum);
              mu_decay_gt_50 = true;
              
            }
          }
        }
      }
    }
    
    
    if (m_branch.mppc_hits>0) {
      // Curiosity plot
      gross_photon_counts->Fill(m_branch.mppc_hits);
    }
    
    // Calculate the number of decays that produce MPPC detected photons
    
    int n_u_hits [n_u_ch] = {0, 0, 0, 0, 0, 0, 0, 0}; // counts in U scints
    int n_d_hits [n_d_ch] = {0, 0, 0, 0, 0};          // counts in D scints
    bool gt_u_thrs = false;   // Have we exceeded the U detection threshold in one counter?
    bool gt_d_thrs = false;   // Have we exceeded the D detection threshold in one counter?
    
    for(int hit = 0; hit < m_branch.mppc_hits; ++hit) {
      
      // Which MPPC is it
      const int mppc_id = get_mppc(m_branch.mppc_y[hit], m_branch.mppc_z[hit]);
      
      if (mppc_id > 9) {
        
        // Upstream scint
        const int index = (mppc_id/10) - 1;
        ++n_u_hits[index];
        if (n_u_hits[index]>u_thrs) {
          gt_u_thrs = true;
        }
      } else {
        
        // Downstream scint
        const int index = mppc_id - 1;
        ++n_d_hits[index];
        if (n_d_hits[index]>d_thrs) {
          gt_d_thrs = true;
        }
      }
    }
    
    if (u_stream_mu && d_stream_e) {
      // We know these photons are from a muon decay 
      photon_counts->Fill(m_branch.mppc_hits);
      ++muon_photon_events;
      if (mu_decay_gt_50){
        photon_counts50ns->Fill(m_branch.mppc_hits);
        ++muon_photon_events_thrs;
      }
    }
  }
  
  TCanvas* can = new TCanvas("c1", "c1", 1436, 856);
  momentum->GetXaxis()->SetTitle("Initial momentum (MeV/c)");
  momentum->GetYaxis()->SetTitle("Count");
  momentum->Draw();
  momentum50ns->SetLineColor(2);
  momentum50ns->Draw("SAMES");
  TLegend* leg = can->BuildLegend(  0.6, 0.8, 0.9, 0.90);
  move_stats_box(can, momentum,     0.6, 0.6, 0.9, 0.79);
  move_stats_box(can, momentum50ns, 0.6, 0.4, 0.9, 0.59);
  can->Update();
  leg->SetFillStyle(0);
  
  can->SaveAs("stopped_mom.eps");
  can->SaveAs("stopped_mom.svg");
  
  
  // Disable the name portion of the stats boxes
  gStyle->SetOptStat(2210);
  
  TCanvas* can_ph = new TCanvas("photon_counts","photon_counts");
  can_ph->SetLogy();
  photon_counts->Draw();
  photon_counts->GetXaxis()->SetTitle("Photons per event");
  photon_counts->GetYaxis()->SetTitle("Count");
  move_stats_box(can_ph, photon_counts, 0.6, 0.7, 0.9, 0.9);
  can_ph->SaveAs("photon_counts.eps");
  can_ph->SaveAs("photon_counts.svg");
  
  TCanvas* can_gr = new TCanvas("gross_photon_counts","gross_photon_counts");
  can_gr->SetLogy();
  gross_photon_counts->Draw();
  gross_photon_counts->GetXaxis()->SetTitle("Photons per event");
  gross_photon_counts->GetYaxis()->SetTitle("Count");
  move_stats_box(can_gr, gross_photon_counts, 0.6, 0.7, 0.9, 0.9);
  can_gr->SaveAs("gross_photon_counts.eps");
  can_gr->SaveAs("gross_photon_counts.svg");
  
  TCanvas* can_p50 = new TCanvas("photon_counts50ns","photon_counts50ns");
  can_p50->SetLogy();
  photon_counts50ns->Draw();
  photon_counts50ns->GetXaxis()->SetTitle("Photons per event");
  photon_counts50ns->GetYaxis()->SetTitle("Count");
  move_stats_box(can_p50, photon_counts50ns, 0.6, 0.7, 0.9, 0.9);
  can_p50->SaveAs("photon_counts50ns.eps");
  can_p50->SaveAs("photon_counts50ns.svg");
  
  printf("decay_muons = %d\n",             decay_muons);
  printf("decay_muons_gt50 = %d\n",        decay_muons_gt50);
  printf("muon_photon_events = %d\n",      muon_photon_events);
  printf("muon_photon_events_thrs = %d\n", muon_photon_events_thrs);
}

void move_stats_box(TCanvas* can, TH1D* hist, float x1, float y1, float x2, float y2) {
  can->Update();
  TPaveStats* st = (TPaveStats*) hist->FindObject("stats");
  st->SetX1NDC(x1);
  st->SetX2NDC(x2);
  st->SetY1NDC(y1);
  st->SetY2NDC(y2);
}

double get_momentum(truth_branch &branch) {
  double px2 = branch.in_px*branch.in_px;
  double py2 = branch.in_py*branch.in_py;
  double pz2 = branch.in_pz*branch.in_pz;
  return TMath::Sqrt(px2 + py2 + pz2);
}

void set_mppc_branches(TTree* tree, mppc_branch &branch) {
   tree->SetBranchAddress("mppc_hits",  &branch.mppc_hits);
   tree->SetBranchAddress("mppc_x",      branch.mppc_x);
   tree->SetBranchAddress("mppc_y",      branch.mppc_y);
   tree->SetBranchAddress("mppc_z",      branch.mppc_z);
   tree->SetBranchAddress("mppc_time",   branch.mppc_time);
}

void set_truth_branches(TTree* tree, truth_branch &branch) {
  tree->SetBranchAddress("nhit",  &branch.nhits);
  tree->SetBranchAddress("in_Px", &branch.in_px);
  tree->SetBranchAddress("in_Py", &branch.in_py);
  tree->SetBranchAddress("in_Pz", &branch.in_pz);
  
  tree->SetBranchAddress("counter",     branch.counter);
  tree->SetBranchAddress("pdgid",       branch.pdgid);
  tree->SetBranchAddress("parentid",    branch.parentid);
  tree->SetBranchAddress("trkid",       branch.trkid);
  tree->SetBranchAddress("first_step",  branch.first_step);
  tree->SetBranchAddress("px",          branch.px);
  tree->SetBranchAddress("py",          branch.py);
  tree->SetBranchAddress("pz",          branch.pz);
  tree->SetBranchAddress("tof",         branch.tof);
}

int get_mppc(double y, double z) {
  // We can ignore the X values as they only differentiate left/right MPPC
  // and for the experiment they were multiplex anyway (also optical 
  // symmetry is a reasonable assumption).
  
  if ((3635.5 < z && z < 3638.5) || (3930.5 < z && z < 3933.5)) {
    // Downstream scintillator
    if        (   99.0<y && y<101.0 ) {
      return 1;
    } else if (   49.0<y && y< 51.0 ) {
      return 2;
    } else if (   -1.0<y && y<  1.0 ) {
      return 3;
    } else if (  -51.0<y && y<-49.0 ) {
      return 4;
    } else if ( -101.0<y && y<-99.0 ) {
      return 5;
    } 
  } else if ((3627.5 < z && z < 3630.5) || (3922.5 < z && z < 3925.5)){
    // Upstream
    if        (   104.0<y && y< 106.0 ) {
      return 11;
    } else if (    74.0<y && y<  76.0 ) {
      return 12;
    } else if (    44.0<y && y<  46.0 ) {
      return 13;
    } else if (    14.0<y && y<  16.0 ) {
      return 14;
    } else if (   -16.0<y && y< -14.0 ) {
      return 15;
    } else if (   -46.0<y && y< -44.0 ) {
      return 16;
    } else if (   -76.0<y && y< -74.0 ) {
      return 17;
    } else if (  -106.0<y && y<-104.0 ) {
      return 18;
    } 
  } 
  printf("Balls up: y:%.1f z:%.1f", y, z);
  return -1;
}

/*
Simulated downstream MPPC positions, all +/- 0.5
X =  832.5, 838.5
Y = -100, -50, 0, 50, 100
Z = 3637, 3932 

Upstream positions, all +/- 0.5
X = 1238, 832.7
Y = -105, -75, -45, -15, 15, 45, 75, 105
Z = 3629, 3924
*/