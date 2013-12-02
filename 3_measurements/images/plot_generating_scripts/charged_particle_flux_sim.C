#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLegend.h"

void charged_particle_flux_sim (){
    gStyle->SetOptStat(2220);
    // Scale from per 9e8 protons to nA, 6.242e9e == 1nA 
    Double_t scale_factor = (9e8/6.242e9);
    TFile* file = new TFile("/Users/scook/code/MuSIC/simulation/MuSIC_5_detector_sim/MuSIC5/g4blout/from_hep_1Bn/out_36_rotate.root");
    TTree* tree = (TTree*) file->Get("t");
    Double_t low_edge [10] = {-165, -135, -65, -35, -15, 15, 35, 65, 135, 165};
    TH1F* hist_1d = new TH1F("hist_1D", "Simulated charged particle flux", 9, low_edge);
    hist_1d->GetXaxis()->SetTitle("Vertical position (mm)");
    hist_1d->GetYaxis()->SetTitle("Charged particle rate (nA^{-1})");
    hist_1d->GetYaxis()->SetTitleOffset(1.4);
    TString selection = TString("(abs(PDGid)==211||abs(PDGid)==13||abs(PDGid)==11||PDGid==2212)&&")+     // charged particles only
                        TString("(abs(y)<15||(35<abs(y)&&abs(y)<65)||(135<abs(y)&&abs(y)<165))&&")+ // vertical positions
                        TString("(abs(x)<190)");                                                         // horizontal range
    TCanvas* can1 = new TCanvas("c1", "c1", 1436,856);
    tree->Draw("-y>>hist_1D", selection.Data(), "e");
    can1->Update();
    hist_1d->Scale(scale_factor);
    TPaveStats* stats_1d = (TPaveStats*) hist_1d->FindObject("stats");
    stats_1d->SetX1NDC(0.1);
    stats_1d->SetX2NDC(0.3);
    stats_1d->SetY1NDC(0.7);
    stats_1d->SetY2NDC(0.9);
    can1->Update();
    
    TCanvas* can2 = new TCanvas("c2", "c2", 1436,856);
    Double_t low_edge_x[6] = {-205, -135, -35, 35, 135, 205};
    Double_t low_edge_y[7] = {-195, -125, -35, 35, 165, 225, 285};
    TH2F* hist_2d = new TH2F("hist_2D", "Simulated charged particle flux", 5, low_edge_x, 6, low_edge_y);
    hist_2d->GetXaxis()->SetTitle("Horizontal position (mm)");
    hist_2d->GetXaxis()->SetTitleOffset(1.6);
    hist_2d->GetYaxis()->SetTitle("Vertical position (mm)");
    hist_2d->GetYaxis()->SetTitleOffset(1.85);
    hist_2d->GetZaxis()->SetTitle("Charged particle rate (nA^{-1})");
    hist_2d->GetZaxis()->SetTitleOffset(1.5);
    
    TString selection2 = TString("(abs(PDGid)==211||abs(PDGid)==13||abs(PDGid)==11||PDGid==2212)&&")+ // charged particles only
                         TString("(abs(x)<35&&abs(y)<35)||")+                  //   0,   0
                         TString("(-205<x&&x<-135&&165<(-y)&&(-y)<225)||")+    // -17,  20
                         TString("(-205<x&&x<-135&&-195<(-y)&&(-y)<-125)||")+  // -17, -16
                         TString("(-205<x&&x<-135&&abs(-y)<35)||")+            // -17,   0
                         TString("(135<x&&x<205&&abs(-y)<35)||")+              //  17,   0
                         TString("(135<x&&x<205&&165<(-y)&&(-y)<225)||")+      //  17,  20
                         TString("(135<x&&x<205&&-195<(-y)&&(-y)<-125)||")+    //  17, -16
                         TString("(abs(x)<35&&-195<(-y)&&(-y)<-125)||")+       //   0, -16
                         TString("(abs(x)<35&&165<(-y)&&(-y)<225)||")+         //   0,  20
                         TString("(abs(x)<35&&225<(-y)&&(-y)<285)||")+         //   0,  20
                         TString("(-205<x&&x<-135&&225<(-y)&&(-y)<285)");      // -17,  20
    tree->Draw("-y:x>>hist_2D", selection2.Data(), "LEGO2FB0");
    
    // Values for the TGraphErrors
    const Int_t n_points = 9;
    Double_t pos     [n_points] = {  -00.5,   -50.0,  -149.5,  -150.5,    50.5,    49.5,   150.5,   149.5,    00.5};
    Double_t pos_er  [n_points] = {   15.0,    15.0,    15.0,    15.0,    15.0,    15.0,    15.0,    15.0,    15.0};
    Double_t flux    [n_points] = {72200.0, 47000.0, 15800.0, 17800.0, 89300.0, 92200.0, 55600.0, 58800.0, 70600.0};
    Double_t flux_er [n_points] = { 3600.0,  2300.0,   810.0,   920.0,  4600.0,  5300.0,  2800.0,  2900.0,  3600.0};
    TCanvas* can3 = new TCanvas("c3", "c3", 1436,856);
    TGraphErrors* measured_flux_1d = new TGraphErrors(n_points,pos,flux,pos_er, flux_er);
    measured_flux_1d->SetTitle("Measured charged particle rate (1D)");
    measured_flux_1d->GetXaxis()->SetTitle("Vertical position (mm)");
    measured_flux_1d->GetYaxis()->SetTitle("Charged particle rate (nA^{-1})");
    measured_flux_1d->GetYaxis()->SetRangeUser(0.0,1e5);
    measured_flux_1d->Draw("AP");
    
    
    // for the 2D graph error
    const Int_t n_points_2d = 13;
    Double_t p_er [n_points_2d]  = {  35.0,   35.0,   35.0,   35.0,   35.0,   35.0,   35.0,   35.0,   35.0,   35.0,   35.0,    35.0,   35.0};
    Double_t x    [n_points_2d]  = {  00.0, -170.0, -170.0, -170.0,  170.0,  170.0,  170.0,   00.0,   00.0,   00.0,   00.0,  -170.0,   00.0};
    Double_t y    [n_points_2d]  = {  00.0,  200.0, -160.0,   00.0,   00.0,  200.0, -160.0, -160.0,  200.0,   00.0,  250.0,   250.0,  200.0};
    Double_t z    [n_points_2d]  = {5400.0, 9300.0, 2200.0, 5010.0, 3110.0, 5590.0, 1540.0, 2400.0, 8600.0, 4360.0, 7300.0, 10800.0, 9200.0};
    Double_t z_er [n_points_2d]  = {1000.0, 1700.0,  330.0,  710.0,  430.0,  810.0,  230.0,  350.0, 1300.0,  660.0, 1500.0,  2100.0, 1500.0};
    // Double_t z    [n_points_2d]  = {54.0,   93.0,   22.1,   50.0,  31.0,  56.0,   15.3,   24.0,  86.0, 44.0,  73.0,  108.0,  92.0};
    // Double_t z_er [n_points_2d]  = {22.0,   37.0,    8.3,   19.0,  11.0,  21.0,    5.8,    9.0,  32.0, 16.0,  38.0,   56.0,  46.0};
    TCanvas* can4 = new TCanvas("c4", "c4", 1436,856);
    TGraph2DErrors* measured_flux_2d = new TGraph2DErrors(n_points_2d, x, y, z, p_er, p_er, z_er);
    // TGraph2D* measured_flux_2d = new TGraph2D(n_points_2d, x, y, z);
    measured_flux_2d->SetTitle("Measured Charged particle rate (2D)");
    measured_flux_2d->GetXaxis()->SetTitle("Horizontal position (mm)");
    measured_flux_2d->GetXaxis()->SetTitleOffset(1.65);
    measured_flux_2d->GetYaxis()->SetTitle("Vertical position (mm)");
    measured_flux_2d->GetYaxis()->SetTitleOffset(1.85);
    measured_flux_2d->GetZaxis()->SetTitle("Charged particle rate (nA^{-1})");
    measured_flux_2d->GetZaxis()->SetTitleOffset(1.5);
    
    measured_flux_2d->SetMarkerStyle(20);
    measured_flux_2d->SetMarkerSize(3);
    measured_flux_2d->Draw("PCOL ERR FB");

    can1->SaveAs("sim_1d_charged_flux.eps");
    can1->SaveAs("sim_1d_charged_flux.svg");
    can2->SaveAs("sim_2d_charged_flux.eps");
    can2->SaveAs("sim_2d_charged_flux.svg");
    
    can3->SaveAs("measured_1d_charged_flux.eps");
    can3->SaveAs("measured_1d_charged_flux.svg");
    can4->SaveAs("measured_2d_charged_flux.eps");
    can4->SaveAs("measured_2d_charged_flux.svg");
    
    
    TCanvas* can5 = new TCanvas("c5", "c5", 1436,856);
    measured_flux_1d->SetLineColor(4);
    measured_flux_1d->SetTitle("Measured rate");
    measured_flux_1d->SetFillColor(0);
    measured_flux_1d->Draw("A P");
    hist_1d->SetLineColor(2);
    hist_1d->SetTitle("Simulated rate");
    hist_1d->SetStats(false);
    hist_1d->SetFillColor(0);
    hist_1d->Draw("SAME E");
    TLegend* leg = can5->BuildLegend(0.1, 0.7, 0.4, 0.9);
    leg->SetFillColor(0);
    leg->Draw();

    measured_flux_1d->SetTitle("1D Charged particle rate");
    can5->Update();
    can5->SaveAs("1D_charged_particle_flux.eps");
    can5->SaveAs("1D_charged_particle_flux.svg");
    
    
}
