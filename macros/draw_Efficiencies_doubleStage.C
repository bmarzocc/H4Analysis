#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

Double_t fitfunc(Double_t *v, Double_t *par) {
  Double_t arg = 0;
  if (par[2] != 0) arg = v[0]/par[2];
  Double_t fitval;
  Double_t effEmission=par[1]*(1 - 1./(par[3]*log(arg)+1));
  if (arg>1)
    fitval = par[0]*(1-effEmission)+effEmission;
  else
    fitval = par[0];
  
  return fitval;
}

void draw_Efficiencies_doubleStage() {

  // inputs
  TFile *inputBinps = new TFile("myOutputFileBinps.root");  
  TFile *inputCamer = new TFile("myOutputFileCameroniBtf.root");  
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *binp1 = (TGraphAsymmErrors*)inputBinps->Get("binp1");
  binp1->RemovePoint(0);
  binp1->RemovePoint(0);
  
  TGraphAsymmErrors *binp3 = (TGraphAsymmErrors*)inputBinps->Get("binp3");
  binp3->RemovePoint(0);
  binp3->RemovePoint(0);
  binp3->RemovePoint(0);

  TGraphAsymmErrors *binp4 = (TGraphAsymmErrors*)inputBinps->Get("binp4");
  binp4->RemovePoint(0);
  binp4->RemovePoint(0);
  binp4->RemovePoint(0);
  
  TGraphAsymmErrors *effMib25Corr_1200 = (TGraphAsymmErrors*)inputCamer->Get("effMib25Corr_1200");
  effMib25Corr_1200->RemovePoint(0);
  effMib25Corr_1200->RemovePoint(0);
  effMib25Corr_1200->RemovePoint(0);

  TGraphAsymmErrors *effRm8Corr_1200   = (TGraphAsymmErrors*)inputCamer->Get("effRm8Corr_1200");
  TGraphAsymmErrors *effRm5Corr_1200   = (TGraphAsymmErrors*)inputCamer->Get("effRm5Corr_1200");

  TGraphAsymmErrors *effZS2Corr  = (TGraphAsymmErrors*)inputH4->Get("effZS2Corr");
  TGraphAsymmErrors *effMib3Corr = (TGraphAsymmErrors*)inputH4->Get("effMib3Corr");

  // cosmetics
  binp1->SetMarkerStyle(24);
  binp1->SetMarkerColor(1);
  binp1->SetLineColor(1);
  binp1->SetLineStyle(6);
  binp1->SetLineWidth(2);
  binp1->SetMarkerSize(1.1);
  binp3->SetMarkerStyle(25);
  binp3->SetMarkerColor(kRed+2);
  binp3->SetLineColor(kRed+2);
  binp3->SetLineStyle(2);
  binp3->SetLineWidth(2);
  binp4->SetMarkerStyle(26);
  binp4->SetMarkerColor(kGreen+3);
  binp4->SetLineColor(kGreen+3);
  binp4->SetLineStyle(3);
  binp4->SetLineWidth(2);
  binp4->SetMarkerSize(1.4);

  effMib25Corr_1200->SetMarkerStyle(32);
  effMib25Corr_1200->SetMarkerColor(kRed+3);
  effMib25Corr_1200->SetLineColor(kRed+3);
  effMib25Corr_1200->SetLineStyle(4);
  effMib25Corr_1200->SetLineWidth(2);
  effMib25Corr_1200->SetMarkerSize(1.4);
  effRm8Corr_1200->SetMarkerStyle(24);
  effRm8Corr_1200->SetMarkerColor(5);
  effRm8Corr_1200->SetLineColor(5);
  effRm8Corr_1200->SetLineStyle(1);
  effRm8Corr_1200->SetLineWidth(2);
  effRm8Corr_1200->SetMarkerSize(1);

  effRm5Corr_1200->SetMarkerStyle(27);
  effRm5Corr_1200->SetMarkerColor(kCyan+2);
  effRm5Corr_1200->SetLineColor(kCyan+2);
  effRm5Corr_1200->SetLineStyle(5);
  effRm5Corr_1200->SetLineWidth(2);
  effRm5Corr_1200->SetMarkerSize(1.7);

  effZS2Corr->SetMarkerStyle(24);
  effZS2Corr->SetMarkerColor(kMagenta+3);
  effZS2Corr->SetLineColor(kMagenta+3);
  effZS2Corr->SetLineStyle(6);
  effZS2Corr->SetLineWidth(2);
  effZS2Corr->SetMarkerSize(1.1);
  effMib3Corr->SetMarkerStyle(25);
  effMib3Corr->SetMarkerColor(kBlue-1);
  effMib3Corr->SetLineColor(kBlue-1);
  effMib3Corr->SetLineStyle(7);
  effMib3Corr->SetLineWidth(2);

  TF1* f1 = new TF1("f1",fitfunc,0.,2700.,4);
  f1->SetParameters(0.2,1.,600.,2.);  
  f1->SetParLimits(0,0.,1.);
  f1->SetParLimits(1,0.,1.5);
  f1->SetParLimits(2,0.,1000.);
  f1->SetParLimits(3,0.,30.);
  f1->SetLineColor(kRed+2);
  f1->SetLineStyle(2);
  f1->SetLineWidth(2);
  binp3->Fit("f1","BR","",650,1500);
 
  TF1* f2 = new TF1("f2",fitfunc,0.,2700.,4);
  f2->SetParameters(0.,1.,1400.,2.);  
  f2->FixParameter(0,0.);
  f2->FixParameter(1,1.);
  //  f2->SetParLimits(0,0.,1.);
  //  f2->SetParLimits(1,0.,2.);
  f2->SetParLimits(2,0.,2000.);
  f2->SetParLimits(3,0.,40.);
  f2->SetLineColor(kRed+3);
  f2->SetLineStyle(4);
  f2->SetLineWidth(2);
  effMib25Corr_1200->Fit("f2","BR","",1490,2500);

  TF1* f3 = new TF1("f3",fitfunc,0.,2700.,4);
  f3->SetParameters(0.,1.,1400.,2.);  
  f3->FixParameter(0,0.);
  f3->FixParameter(1,1.);
  //  f3->SetParLimits(0,0.,1.);
  //  f3->SetParLimits(1,0.,2.);
  f3->SetParLimits(2,0.,2000.);
  f3->SetParLimits(3,0.,40.);
  f3->SetLineColor(kCyan+2);
  f3->SetLineStyle(5);
  f3->SetLineWidth(2);
  effRm5Corr_1200->Fit("f3","BR","",880,1700);

  TF1* f4 = new TF1("f4",fitfunc,0.,2700.,4);
  f4->SetParameters(0.,1.,1480.,2.);  
  //f4->FixParameter(0,0.);
  f4->FixParameter(1,1.);
  f4->SetParLimits(0,0.,1.);
  //f4->SetParLimits(1,0.,2.);
  f4->SetParLimits(2,0.,2000.);
  f4->SetParLimits(3,0.,40.);
  f4->SetLineColor(1);
  f4->SetLineStyle(6);
  f4->SetLineWidth(2);
  binp1->Fit("f4","BR","",1480,2500);

  TF1* f5 = new TF1("f5",fitfunc,0.,2700.,4);
  f5->SetParameters(0.,1.,1400.,2.);  
  //  f5->FixParameter(0,0.);
  f5->FixParameter(1,1.);
  f5->SetParLimits(0,0.,1.);
  //  f5->SetParLimits(1,0.,2.);
  f5->SetParLimits(2,0.,2000.);
  f5->SetParLimits(3,0.,40.);
  f5->SetLineColor(kGreen+3);
  f5->SetLineStyle(3);
  f5->SetLineWidth(2);
  binp4->Fit("f5","BR","",1480,2500);

  // plotting
  setStyle();

  //TH2F* H2 = new TH2F("H2","",3200,0.,3200.,100,0.,1.4);
  TH2F* H2 = new TH2F("H2","",2200,500.,2700.,100,0.,1.4);
  H2->GetXaxis()->SetTitle("MCP-emission stage bias (V)");
  H2->GetYaxis()->SetTitle("Efficiency");
  H2->GetXaxis()->SetNdivisions(505);

  TLine *line = new TLine(500,1,2700,1);

  TLegend *legA;
  legA = new TLegend(0.18,0.70,0.38,0.90);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.03);
  legA->SetFillColor(0);
  legA->AddEntry(binp1, "40x2+40x2, Amp bias=2100V", "pl");
  legA->AddEntry(binp3, "90x1+40x2, Amp bias=2100V", "pl");
  legA->AddEntry(binp4, "90x2+40x2, Amp bias=2100V", "pl");
  //legA->AddEntry(effZS2Corr,  "40x3, electrons at 50 GeV", "pl");

  TLegend *legB;
  legB = new TLegend(0.58,0.70,0.78,0.90);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.03);
  legB->SetFillColor(0);
  legB->AddEntry(effMib25Corr_1200, "80x2+40x1, Amp bias=1200V", "pl");
  //legB->AddEntry(effRm8Corr_1200,   "80x1+40x1, Amp bias=1200V", "pl");
  legB->AddEntry(effRm5Corr_1200,   "80x1+40x1, Amp bias=1200V", "pl");
  //legB->AddEntry(effMib3Corr, "40x2, electrons at 50 GeV", "pl");

  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  binp1->Draw("Psame");
  binp3->Draw("Psame");
  binp4->Draw("Psame");
  effMib25Corr_1200->Draw("Psame");
  //effRm8Corr_1200->Draw("PCsame");
  effRm5Corr_1200->Draw("Psame");
  //effZS2Corr->Draw("PCsame");
  //effMib3Corr->Draw("PCsame");
  legA->Draw();
  legB->Draw();
  TLatex latex2(0.69, 0.94,"#bf{#bf{491 MeV electrons}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw(); 
  c1->SaveAs("summaryEff_doubleStage.png");
  c1->SaveAs("summaryEff_doubleStage.pdf");


}
