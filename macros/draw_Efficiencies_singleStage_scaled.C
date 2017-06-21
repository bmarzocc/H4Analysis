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

Double_t fitfunc(Double_t *v, Double_t *par);
Double_t fitfunc2(Double_t *v, Double_t *par);
std::pair<float,std::pair<float,float> > computeStat( std::vector<double>, std::vector<double>, std::vector<double>); 
void ScaleGraph(TGraphAsymmErrors*,float);

void draw_Efficiencies_singleStage_scaled() {

  gStyle->SetOptTitle(0); 
  //gStyle->SetOptStat(1110); 
  gStyle->SetOptStat(0000); 
  //gStyle->SetOptFit(1); 
  gStyle->SetOptFit(0); 
  //gStyle->SetErrorX(0); 

  // inputs
  TFile *input    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *effSee = (TGraphAsymmErrors*)input->Get("effSee");
  TGraphAsymmErrors *effMgO_PadA = (TGraphAsymmErrors*)input->Get("effMgO_PadA");
  TGraphAsymmErrors *effMgO_PadB = (TGraphAsymmErrors*)input->Get("effMgO_PadB");
  TGraphAsymmErrors *effMgO_PadC = (TGraphAsymmErrors*)input->Get("effMgO_PadC");
  TGraphAsymmErrors *effMgO_PadD = (TGraphAsymmErrors*)input->Get("effMgO_PadD");

  TGraphAsymmErrors *effMgO = new TGraphAsymmErrors();
  
  std::vector<std::vector<double> > mean;
  mean.resize(effMgO_PadA->GetN());
  std::vector<std::vector<double> > error_up;
  error_up.resize(effMgO_PadA->GetN());
  std::vector<std::vector<double> > error_down;
  error_down.resize(effMgO_PadA->GetN());

  for(int ii = 0; ii<effMgO_PadA->GetN(); ii++)
  {
      double x,y;
      double errorUp, errorDown;

      effMgO_PadA->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadA->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadA->GetErrorYlow(ii));

      effMgO_PadB->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadB->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadB->GetErrorYlow(ii));

      effMgO_PadC->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadC->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadC->GetErrorYlow(ii));

      effMgO_PadD->GetPoint(ii,x,y);
      mean.at(ii).push_back(y);
      error_up.at(ii).push_back(effMgO_PadD->GetErrorYhigh(ii));
      error_down.at(ii).push_back(effMgO_PadD->GetErrorYlow(ii));
  }

  for(int ii = 0; ii<effMgO_PadA->GetN(); ii++)
  {
      float Mean = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).first;
      float Error_up = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).second.first;
      float Error_down = computeStat( mean.at(ii), error_up.at(ii), error_down.at(ii)).second.second;
        
      double x,y;
      effMgO_PadA->GetPoint(ii,x,y);
      effMgO->SetPoint(ii,x,Mean);
  }

  // inputsH4
  TFile *inputH4    = new TFile("myOutputFileH4.root");  

  TGraphAsymmErrors *effZS2Corr  = (TGraphAsymmErrors*)inputH4->Get("effZS2Corr");
  TGraphAsymmErrors *effMib3Corr = (TGraphAsymmErrors*)inputH4->Get("effMib3Corr");
  TGraphAsymmErrors *effZS2PMTCorr = (TGraphAsymmErrors*)inputH4->Get("effZS2PMTCorr");
 
  TF1* f1 = new TF1("f1",fitfunc,1000.,3200.,4);
  f1->SetParameters(0.2,0.1,1.,1450.);  
  f1->SetParLimits(1,0.,1.);
  f1->SetParLimits(3,1400.,1450.);
  //f1->FixParameter(0,0.);
  effSee->Fit("f1","B,0");

  TF1* f1_post = new TF1("f1_post",fitfunc,0.8,2.,4);
  f1_post->FixParameter(0,f1->GetParameter(0));
  f1_post->FixParameter(1,f1->GetParameter(1));
  f1_post->FixParameter(2,f1->GetParameter(2));
  f1_post->FixParameter(3,1.);
  f1_post->SetLineColor(kCyan+2);
  f1_post->SetLineStyle(7);
  f1_post->SetLineWidth(2);

  ScaleGraph(effSee,f1->GetParameter(3));

  TF1* f2 = new TF1("f2",fitfunc,1000.,3200.,4);
  f2->SetParameters(0.2,0.1,1.,2100.);  
  f2->SetParLimits(1,0.,1.);
  f2->SetParLimits(3,2000.,2100.);
  f2->FixParameter(0,0.);
  effMgO->Fit("f2","B,0");

  TF1* f2_post = new TF1("f2_post",fitfunc,0.8,2.,4);
  f2_post->FixParameter(0,f2->GetParameter(0));
  f2_post->FixParameter(1,f2->GetParameter(1));
  f2_post->FixParameter(2,f2->GetParameter(2));
  f2_post->FixParameter(3,1.);
  f2_post->SetLineColor(kGreen+2);
  f2_post->SetLineStyle(8);
  f2_post->SetLineWidth(2);
   
  ScaleGraph(effMgO,f2->GetParameter(3));

  TF1* f3 = new TF1("f3",fitfunc,1300.,3200.,4);
  f3->SetParameters(0.2,1.,1.,2250.);  
  f3->SetParLimits(1,0.,1.);
  f3->SetParLimits(3,1800.,2250.);
  //f3->FixParameter(0,0.);
  effZS2Corr->Fit("f3","B,0");

  TF1* f3_post = new TF1("f3_post",fitfunc,0.8,2.,4);
  f3_post->FixParameter(0,f3->GetParameter(0));
  f3_post->FixParameter(1,f3->GetParameter(1));
  f3_post->FixParameter(2,f3->GetParameter(2));
  f3_post->FixParameter(3,1.);
  f3_post->SetLineColor(kViolet+1);
  f3_post->SetLineStyle(7);
  f3_post->SetLineWidth(2);

  ScaleGraph(effZS2Corr,f3->GetParameter(3));

  TF1* f4 = new TF1("f4",fitfunc,1300.,3200.,4);
  f4->SetParameters(0.2,1.,1.,2000.);  
  f4->SetParLimits(1,0.,1.);
  f4->SetParLimits(3,1800.,2050.);
  //f4->FixParameter(0,0.);
  effMib3Corr->Fit("f4","B,0");

  TF1* f4_post = new TF1("f4_post",fitfunc,0.8,2.,4);
  f4_post->FixParameter(0,f4->GetParameter(0));
  f4_post->FixParameter(1,f4->GetParameter(1));
  f4_post->FixParameter(2,f4->GetParameter(2));
  f4_post->FixParameter(3,1.);
  f4_post->SetLineColor(kRed+2);
  f4_post->SetLineStyle(8);
  f4_post->SetLineWidth(2);

  ScaleGraph(effMib3Corr,f4->GetParameter(3));

  ScaleGraph(effZS2PMTCorr,2030.);

  TF1* f5 = new TF1("f5",fitfunc2,0.8,2.,3);
  f5->SetLineColor(kOrange+1);  
  f5->SetLineStyle(6);
  f5->SetLineWidth(2);
  f5->SetParameters(1.,1.,1.,1.); 
  f5->SetParLimits(0,0.,1.);  
  f5->SetParLimits(2,-10.,0.); 
  effZS2PMTCorr->Fit("f5","B");

  // cosmetics
  effSee->SetMarkerStyle(24);
  effSee->SetMarkerColor(kCyan+2);
  effSee->SetLineColor(kCyan+2);  
  effSee->SetLineStyle(4);
  effSee->SetLineWidth(2);
  effSee->SetMarkerSize(1.1);

  effMgO->SetMarkerStyle(25);
  effMgO->SetMarkerColor(kGreen+2);
  effMgO->SetLineColor(kGreen+2);  
  effMgO->SetLineStyle(5);
  effMgO->SetLineWidth(2);
  effMgO->SetMarkerSize(1.1);

  effZS2Corr->SetMarkerStyle(26);
  effZS2Corr->SetMarkerColor(kViolet+1);
  effZS2Corr->SetLineColor(kViolet+1);
  effZS2Corr->SetLineStyle(6);
  effZS2Corr->SetLineWidth(2);
  effZS2Corr->SetMarkerSize(1.3);

  effMib3Corr->SetMarkerStyle(32);
  effMib3Corr->SetMarkerColor(kRed+2);
  effMib3Corr->SetLineColor(kRed+2);
  effMib3Corr->SetLineStyle(7);
  effMib3Corr->SetLineWidth(2);
  effMib3Corr->SetMarkerSize(1.3);

  effZS2PMTCorr->SetMarkerStyle(27);
  effZS2PMTCorr->SetMarkerColor(kOrange+1);
  effZS2PMTCorr->SetLineColor(kOrange+1);  
  effZS2PMTCorr->SetLineStyle(8);
  effZS2PMTCorr->SetLineWidth(2);
  effZS2PMTCorr->SetMarkerSize(1.5);

  // plotting
  setStyle();

  TH2F* H2 = new TH2F("H2","",220,0.8,2.,100,0.,1.25);
  H2->GetXaxis()->SetTitle("MCP-stack bias/threshold bias");
  H2->GetYaxis()->SetTitle("Efficiency");
  H2->GetXaxis()->SetNdivisions(505);

  TLine *line = new TLine(0.8,1,2.,1);

  TLegend *legA;
  legA = new TLegend(0.23,0.77,0.43,0.90);
  legA->SetFillStyle(0);
  legA->SetBorderSize(0);
  legA->SetTextSize(0.03);
  legA->SetFillColor(0);
  //legA->AddEntry(effMib3Corr,  "40x2 iMCP", "pl");
  legA->AddEntry(effSee,  "40x2 SEE-iMCP", "pl");
  legA->AddEntry(effMgO,  "40x2 MgO-iMCP, 150 GeV muons", "pl");

  TLegend *legB;
  legB = new TLegend(0.68,0.77,0.88,0.90);
  legB->SetFillStyle(0);
  legB->SetBorderSize(0);
  legB->SetTextSize(0.03);
  legB->SetFillColor(0);
  legB->AddEntry(effMib3Corr,  "40x2 iMCP", "pl");
  legB->AddEntry(effZS2Corr,  "40x3 iMCP", "pl");
  //legB->AddEntry(effZS2PMTCorr,  "40x3 PMT-MCP", "pl");

  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  line->Draw();
  effSee->Draw("Psame");
  effMgO->Draw("Psame");
  effZS2Corr->Draw("Psame");
  effMib3Corr->Draw("Psame");
  //effZS2PMTCorr->Draw("Psame");
  //c1->RedrawAxis("sameaxis");
  f1_post->Draw("Lsame");
  f2_post->Draw("Lsame");
  f3_post->Draw("Lsame");
  f4_post->Draw("Lsame");
  legA->Draw("same");
  legB->Draw("same");
  TLatex latex2(0.705, 0.94,"#bf{#bf{50 GeV electrons}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw(); 
  c1->SaveAs("summaryEff_singleStage_scaled.png");
  c1->SaveAs("summaryEff_singleStage_scaled.pdf");
}

Double_t fitfunc(Double_t *v, Double_t *par) {
  Double_t arg = 0;
  if (par[3] != 0) arg = v[0]/par[3];
  Double_t fitval;
  if (arg>1)
    fitval = par[0] + (par[1]*(1.- (1./(par[2]*log(arg)+1))));
  else
    fitval = 0.;
  
  return fitval;
}

Double_t fitfunc2(Double_t *v, Double_t *par) {

  Double_t fitval;
  fitval = 1./(1. + par[0]*exp(par[1]*(v[0]+par[2])));
 
  return fitval;
}
std::pair<float,std::pair<float,float> > computeStat( std::vector<double> mean, std::vector<double> error_up, std::vector<double> error_down)
{
  float Mean = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Mean = Mean +mean.at(ii);

  Mean = Mean/mean.size();
  
  float Error_up = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Error_up = Error_up + (mean.at(ii)-Mean)*(mean.at(ii)-Mean);

  Error_up = Error_up/(mean.size()-1);
  Error_up = sqrt(Error_up);

  float Error_down = 0.;
  for(unsigned int ii = 0; ii<mean.size(); ii++)
      Error_down = Error_down + (mean.at(ii)-Mean)*(mean.at(ii)-Mean);

  Error_down = Error_down/(mean.size()-1);
  Error_down = sqrt(Error_down);

  std::pair<float,std::pair<float,float> > outpair = std::make_pair(Mean,std::make_pair(Error_up,Error_down));
  return outpair ;
}

void ScaleGraph(TGraphAsymmErrors* g,float scale)
{
   for(int ii = 0; ii<g->GetN(); ii++)
   {
      double x,y;
      g->GetPoint(ii,x,y);
      g->SetPoint(ii,x/scale,y);
      g->SetPointEXhigh(ii,g->GetErrorXhigh(ii)/scale);
      g->SetPointEXlow(ii,g->GetErrorXlow(ii)/scale);
   } 
} 
