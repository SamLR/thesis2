#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"

struct n_photons_t
{
  Double_t x[3000];
  Double_t y[3000];
  Double_t z[3000];
  Int_t nhits;
};

void n_photons(){
  gStyle->SetOptStat(2210);
  TString filename = TString("/Users/scook/code/MuSIC/simulation/") +
    TString("MuSIC_5_detector_sim/MuSIC5/output/output_from_hep_batch/") +
    TString("optical_processes.root");
  TFile* file = new TFile(filename.Data(), "READ");
  TTree* mppc = (TTree*) file->Get("mppc");
  TH1F* hist = new TH1F("h", "Number of detected photons", 1000, 0, 1000);
  hist->GetXaxis()->SetTitle("Number of photons");
  hist->GetYaxis()->SetTitle("Count");
  n_photons_t branch;
  mppc->SetBranchAddress("mppc_x",     branch.x);
  mppc->SetBranchAddress("mppc_y",     branch.y);
  mppc->SetBranchAddress("mppc_z",     branch.z);
  mppc->SetBranchAddress("mppc_hits", &branch.nhits);
  
  for(int entry = 0; entry < mppc->GetEntries(); ++entry) {
    mppc->GetEntry(entry);
    
    int count = 0;
    for(Int_t hit = 0; hit < branch.nhits; ++hit) {
      bool correct_z = 3635.0 < branch.z[hit] && branch.z[hit] < 3640.0;
      bool correct_y = 97.0   < branch.y[hit] && branch.y[hit] < 102.0;
      
      if(correct_z && correct_y) {
           ++count;
      }
    }
    if (count > 0) hist->Fill(count);
  }
  TCanvas* can = new TCanvas("c1","c1",1436,856);
  can->SetLogx();
  can->SetLogy();
  
  hist->Draw();
  can->Update();
  TPaveStats* stats = (TPaveStats*) hist->FindObject("stats");
  stats->SetX1NDC(0.65);
  stats->SetX2NDC(0.90);
  stats->SetY1NDC(0.70);
  stats->SetY2NDC(0.90);
  can->Update();
  
  can->SaveAs("n_photons.png");
  can->SaveAs("n_photons.svg");
  
}