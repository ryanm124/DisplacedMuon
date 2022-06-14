// ----------------------------------------------------------------------------------------------------------------
// Feasibility study of using L1 Tracks to identify Displaced Vertex
//
// By Bharadwaj Harikrishnan, May 2021
// Edited by Ryan McCarthy, Sept 2021
// ----------------------------------------------------------------------------------------------------------------

#define Analyzer_DisplacedMuon_cxx
#include "Analyzer_DisplacedMuon.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <valarray>
#include <deque>
using namespace std;

bool VERBOSE[]={false,false,false,false};

//float d0_res = 0.0105; //cm
//float z0_res = 0.3477; //cm
float d0_res = 0.0554; //cm
float z0_res = 0.7045; //cm
float CUTOFF = 1.0; // Common z cut off in cm for Displaced Vertex and verification

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char *text);

void displayProgress(long current, long max)
{
   using std::cerr;
   if (max < 2500)
      return;
   if (current % (max / 2500) != 0 && current < max - 1)
      return;

   int width = 52; // Hope the terminal is at least that wide.
   int barWidth = width - 2;
   cerr << "\x1B[2K";    // Clear line
   cerr << "\x1B[2000D"; // Cursor left
   cerr << '[';
   for (int i = 0; i < barWidth; ++i)
   {
      if (i < barWidth * current / max)
      {
         cerr << '=';
      }
      else
      {
         cerr << ' ';
      }
   }
   cerr << ']';
   cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
   cerr.flush();
}


bool ComparePtTrack(Track_Parameters a, Track_Parameters b) { return a.pt > b.pt; }
bool CompareZ0Track(Track_Parameters a, Track_Parameters b) { return a.z0 > b.z0; }
bool CompareD0Track(Track_Parameters a, Track_Parameters b) { return a.d0 > b.d0; }

std::valarray<float> calcPVec(Track_Parameters a, double_t v_x, double_t v_y)
{
  std::valarray<float> r_vec = {v_x-a.x0,v_y-a.y0};
  std::valarray<float> p_vec = {-r_vec[1],r_vec[0]};
  if(a.charge>0){
    p_vec *= -1;
  }
  p_vec /= TMath::Sqrt(pow(p_vec[0],2)+pow(p_vec[1],2));
  p_vec *= a.pt;
  return p_vec;
}

template<typename T>
std::vector<T> linspace(T start, T end, int num){
  std::vector<T> out;
  T delta = (end - start) / (num-1);
  for(int i=0; i<num-1; i++){
    out.push_back(start+delta*i);
  }
  out.push_back(end);
  return out;
}

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
   Double_t dPhi = phi1 - phi2;
   if (dPhi > TMath::Pi())
      dPhi -= 2. * TMath::Pi();
   if (dPhi < -TMath::Pi())
      dPhi += 2. * TMath::Pi();
   return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
   Double_t dEta, dPhi;
   dEta = eta1 - eta2;
   dPhi = deltaPhi(phi1, phi2);
   return sqrt(dEta * dEta + dPhi * dPhi);
}

Double_t dist(Double_t x1, Double_t y1 , Double_t x2=0, Double_t y2=0){ // Distance between 2 points
   return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
}

Double_t dist_Vertex(Double_t x_vtx, Double_t y_vtx, Track_Parameters *a){ // Distance between track and displaced vertex
   float R = dist(x_vtx,y_vtx,a->x0,a->y0);
   return (fabs(R-(a->rho)));
}

Double_t dist_TPs(Track_Parameters a, Track_Parameters b); // Closest distance between 2 tracks
bool CompareDeltaXY(Vertex_Parameters v1, Vertex_Parameters v2) {return dist_TPs(v1.a,v1.b) < dist_TPs(v2.a,v2.b); }


Int_t Vertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx); 
// Identify the displaced vertex (x_vtx,y_vtx,z_vtx) and return the status 
//-2 = Circles with same center. No Intersection
//-1 = Circles don't Intersect. A point on the line connecting the centers is chosen.
// 0 = Only 1 Intersection not satisfying Z cutoff
// 1 = Only 1 Intersection satisfying Z cutoff
// 2 = Only 1 Intersection detectable dist(x,y)<20
// 3 = 2 Intersections 

void Analyzer_DisplacedMuon::Loop(TString type,
                            TString type_dir = "",
                            float TP_maxD0 = 10.0,
                            float TP_minD0 = 0.0004196,
                            int TP_select_injet = 0)
{
   gROOT->SetBatch();
   gErrorIgnoreLevel = kWarning;

   SetPlotStyle();

   float TP_minPt = 3.0;
   float TP_maxPt = 10000.0;
   float TP_maxEta = 2.4;

   TH1F *h_trk_d0 = new TH1F("h_trk_d0","h_trk_d0; Track d_{0} Distribution (cm) ; Events / 0.1 cm",200,-10,10);
   TH1F *h_trk_pt = new TH1F("h_trk_pt","h_trk_pt; Track p_{T} Distribution (GeV); Events / 1.0 GeV", 100, 0, 100.0);
   TH1F *h_trk_eta = new TH1F("h_trk_eta","h_trk_eta; Track #eta Distribution; Events / 0.1", 50, -2.5, 2.5);
   TH1F *h_tp_pt = new TH1F("h_tp_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F *h_tp_eta = new TH1F("h_tp_eta", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F *h_tp_d0 = new TH1F("h_tp_d0", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   TH1F *h_trueVertex_errorCode = new TH1F("h_trueVertex_errorCode","h_trueVertex_errorCode; Error Code; Events / 1.0",11,0,11);
   TH1F *h_trueVertex_numAllCuts = new TH1F("h_trueVertex_numAllCuts","h_trueVertex_numAllCuts; TP Vertices; Events / 1.0",40,0,40);
   TH1F *h_trueVertex_numNoCuts = new TH1F("h_trueVertex_numNoCuts","h_trueVertex_numNoCuts; TP Vertices; Events / 1.0",40,0,40);
   TH1F *h_trueVertex_x = new TH1F("h_trueVertex_x","h_trueVertex_x; TP Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trueVertex_y = new TH1F("h_trueVertex_y","h_trueVertex_y; TP Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trueVertex_z = new TH1F("h_trueVertex_z","h_trueVertex_z; TP Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
   TH1F *h_trueVertex_sumPt = new TH1F("h_trueVertex_sumPt","h_trueVertex_sumPt; TP Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trueVertex_highPt = new TH1F("h_trueVertex_highPt","h_trueVertex_highPt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trueVertex_lowPt = new TH1F("h_trueVertex_lowPt","h_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trueVertexCuts_x = new TH1F("h_trueVertexCuts_x","h_trueVertexCuts_x; TP Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trueVertexCuts_y = new TH1F("h_trueVertexCuts_y","h_trueVertexCuts_y; TP Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trueVertexCuts_z = new TH1F("h_trueVertexCuts_z","h_trueVertexCuts_z; TP Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
   TH1F *h_trueVertexCuts_sumPt = new TH1F("h_trueVertexCuts_sumPt","h_trueVertexCuts_sumPt; TP Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trueVertexCuts_indexPt = new TH1F("h_trueVertexCuts_indexPt","h_trueVertexCuts_indexPt; Pt Ranking of TP; Events / 1.0",15,0,15.0);
   TH1F *h_trackVertex_numAllCuts = new TH1F("h_trackVertex_numAllCuts","h_trackVertex_numAllCuts; Track Vertices; Events / 1.0",40,0,40);
   TH1F *h_trackVertex_numNoCuts = new TH1F("h_trackVertex_numNoCuts","h_trackVertex_numNoCuts; Track Vertices; Events / 1.0",40,0,40);
   TH1F *h_trackVertex_x = new TH1F("h_trackVertex_x","h_trackVertex_x; Track Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trackVertex_y = new TH1F("h_trackVertex_y","h_trackVertex_y; Track Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trackVertex_z = new TH1F("h_trackVertex_z","h_trackVertex_z; Track Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
   TH1F *h_trackVertex_sumPt = new TH1F("h_trackVertex_sumPt","h_trackVertex_sumPt; Track Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trackVertexCuts_x = new TH1F("h_trackVertexCuts_x","h_trackVertexCuts_x; Track Vertex X Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trackVertexCuts_y = new TH1F("h_trackVertexCuts_y","h_trackVertexCuts_y; Track Vertex Y Position (cm); Events / 0.1 cm",100,-5.0,5.0);
   TH1F *h_trackVertexCuts_z = new TH1F("h_trackVertexCuts_z","h_trackVertexCuts_z; Track Vertex Z Position (cm); Events / 0.025 cm",4000,-50,50);
   TH1F *h_trackVertexCuts_sumPt = new TH1F("h_trackVertexCuts_sumPt","h_trackVertexCuts_sumPt; Track Vertex Momentum (GeV); Events / 1.0 GeV",200,0,200.0);
   TH1F *h_trackVertexCuts_indexPt = new TH1F("h_trackVertexCuts_indexPt","h_trackVertexCuts_indexPt; Pt Ranking of Track; Events / 1.0",50,0,50.0);

   // Chi2 plots
   TH1F *h_trk_chi2rphidof = new TH1F("h_trk_chi2rphidof","h_trk_chi2rphidof; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rphidof_H = new TH1F("h_trk_chi2rphidof_H","h_trk_chi2rphidof_H; Track #chi^{2}_{r#phi}/d.o.f (p_{T}>8 GeV); Events / 0.1",50,0,5);   
   TH1F *h_trk_chi2rphidof_L = new TH1F("h_trk_chi2rphidof_L","h_trk_chi2rphidof_L; Track #chi^{2}_{r#phi}/d.o.f (p_{T}<8 GeV); Events / 0.1",50,0,5);   
   TH1F *h_trk_chi2rphidof_C = new TH1F("h_trk_chi2rphidof_C","h_trk_chi2rphidof_C; Track #chi^{2}_{r#phi}/d.o.f (|#eta|<0.8); Events / 0.1",50,0,5);   
   TH1F *h_trk_chi2rphidof_I = new TH1F("h_trk_chi2rphidof_I","h_trk_chi2rphidof_I; Track #chi^{2}_{r#phi}/d.o.f (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);   
   TH1F *h_trk_chi2rphidof_F = new TH1F("h_trk_chi2rphidof_F","h_trk_chi2rphidof_F; Track #chi^{2}_{r#phi}/d.o.f (|#eta|>1.6); Events / 0.1",50,0,5);   
   TH1F *h_trk_chi2rphidof_P = new TH1F("h_trk_chi2rphidof_P","h_trk_chi2rphidof_P; Track #chi^{2}_{r#phi}/d.o.f (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rphidof_D = new TH1F("h_trk_chi2rphidof_D","h_trk_chi2rphidof_D; Track #chi^{2}_{r#phi}/d.o.f (d_{0}>1 cm); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof = new TH1F("h_trk_chi2rzdof","h_trk_chi2rzdof; Track #chi^{2}_{rz}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_H = new TH1F("h_trk_chi2rzdof_H","h_trk_chi2rzdof_H; Track #chi^{2}_{rz}/d.o.f (p_{T}>8 GeV); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_L = new TH1F("h_trk_chi2rzdof_L","h_trk_chi2rzdof_L; Track #chi^{2}_{rz}/d.o.f (p_{T}<8 GeV); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_C = new TH1F("h_trk_chi2rzdof_C","h_trk_chi2rzdof_C; Track #chi^{2}_{rz}/d.o.f (|#eta|<0.8); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_I = new TH1F("h_trk_chi2rzdof_I","h_trk_chi2rzdof_I; Track #chi^{2}_{rz}/d.o.f (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_F = new TH1F("h_trk_chi2rzdof_F","h_trk_chi2rzdof_F; Track #chi^{2}_{rz}/d.o.f (|#eta|>1.6); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_P = new TH1F("h_trk_chi2rzdof_P","h_trk_chi2rzdof_P; Track #chi^{2}_{rz}/d.o.f (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
   TH1F *h_trk_chi2rzdof_D = new TH1F("h_trk_chi2rzdof_D","h_trk_chi2rzdof_D; Track #chi^{2}_{rz}/d.o.f (d_{0}>1 cm); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2 = new TH1F("h_trk_bendchi2","h_trk_bendchi2; Track bend #chi^{2} ; Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_H = new TH1F("h_trk_bendchi2_H","h_trk_bendchi2_H; Track bend #chi^{2} (p_{T}>8 GeV); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_L = new TH1F("h_trk_bendchi2_L","h_trk_bendchi2_L; Track bend #chi^{2} (p_{T}<8 GeV); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_C = new TH1F("h_trk_bendchi2_C","h_trk_bendchi2_C; Track bend #chi^{2} (|#eta|<0.8); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_I = new TH1F("h_trk_bendchi2_I","h_trk_bendchi2_I; Track bend #chi^{2} (0.8#leq|#eta|#leq1.6); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_F = new TH1F("h_trk_bendchi2_F","h_trk_bendchi2_F; Track bend #chi^{2} (|#eta|>1.6); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_P = new TH1F("h_trk_bendchi2_P","h_trk_bendchi2_P; Track bend #chi^{2} (d_{0}#leq 1 cm); Events / 0.1",50,0,5);
   TH1F *h_trk_bendchi2_D = new TH1F("h_trk_bendchi2_D","h_trk_bendchi2_D; Track bend #chi^{2} (d_{0}>1 cm); Events / 0.1",50,0,5);

   // Efficiency of Identifying Tracks Plots
   TH1F* h_tp_pt_noCuts = new TH1F("h_tp_pt_noCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_noCuts = new TH1F("h_match_tp_pt_noCuts", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_noCuts = new TH1F("h_tp_eta_noCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_noCuts = new TH1F("h_match_tp_eta_noCuts", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_noCuts = new TH1F("h_tp_d0_noCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_noCuts = new TH1F("h_match_tp_d0_noCuts", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   TH1F* h_tp_pt_maxD0Cut = new TH1F("h_tp_pt_maxD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_maxD0Cut = new TH1F("h_match_tp_pt_maxD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_maxD0Cut = new TH1F("h_tp_eta_maxD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_maxD0Cut = new TH1F("h_match_tp_eta_maxD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_maxD0Cut = new TH1F("h_tp_d0_maxD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_maxD0Cut = new TH1F("h_match_tp_d0_maxD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   
   TH1F* h_tp_pt_minD0Cut = new TH1F("h_tp_pt_minD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_minD0Cut = new TH1F("h_match_tp_pt_minD0Cut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_minD0Cut = new TH1F("h_tp_eta_minD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_minD0Cut = new TH1F("h_match_tp_eta_minD0Cut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_minD0Cut = new TH1F("h_tp_d0_minD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_minD0Cut = new TH1F("h_match_tp_d0_minD0Cut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   TH1F* h_tp_pt_minPtCut = new TH1F("h_tp_pt_minPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_minPtCut = new TH1F("h_match_tp_pt_minPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_minPtCut = new TH1F("h_tp_eta_minPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_minPtCut = new TH1F("h_match_tp_eta_minPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_minPtCut = new TH1F("h_tp_d0_minPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_minPtCut = new TH1F("h_match_tp_d0_minPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   TH1F* h_tp_pt_maxPtCut = new TH1F("h_tp_pt_maxPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_maxPtCut = new TH1F("h_match_tp_pt_maxPtCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_maxPtCut = new TH1F("h_tp_eta_maxPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_maxPtCut = new TH1F("h_match_tp_eta_maxPtCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_maxPtCut = new TH1F("h_tp_d0_maxPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_maxPtCut = new TH1F("h_match_tp_d0_maxPtCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   TH1F* h_tp_pt_maxEtaCut = new TH1F("h_tp_pt_maxEtaCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt_maxEtaCut = new TH1F("h_match_tp_pt_maxEtaCut", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta_maxEtaCut = new TH1F("h_tp_eta_maxEtaCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta_maxEtaCut = new TH1F("h_match_tp_eta_maxEtaCut", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0_maxEtaCut = new TH1F("h_tp_d0_maxEtaCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0_maxEtaCut = new TH1F("h_match_tp_d0_maxEtaCut", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);

   // Displaced Vertex Plots
   TH1F *h_delta_dist_xy = new TH1F("h_delta_dist_xy","h_delta_dist_xy; (a) Distance between chosen TPs in x-y (cm) ; Events / 0.02 cm",50,0,1);
   TH1F *h_error_delta_x = new TH1F("h_error_delta_x","h_error_delta_x; (b) #Delta x_{displaced} with chosen TP true x (cm) ; Events / 0.02 cm",100,0,2);
   TH1F *h_delta_dist_z = new TH1F("h_delta_dist_z","h_delta_dist_z; (c) Distance between chosen TPs in z (cm) ; Events / 0.05 cm",40,0,2);
   TH1F *h_error_delta_z = new TH1F("h_error_delta_z","h_error_delta_z; (d) Chosen TP error in z (cm) ; Events / 0.1 cm",100,0,10);
   TH1F *h_delta_x = new TH1F("h_delta_x","h_delta_x; #Delta x between chosen TPs true x (cm) ; Events / 0.05 cm",40,0,2);

   TH1I *h_trk_Counter_TPcombination = new TH1I("h_trk_Counter_TPcombination","h_trk_Counter_TPcombination; Track combination chosen; Events / 1.0",6,0,6);
   TH1F *h_trk_delta_dist_xy = new TH1F("h_trk_delta_dist_xy","h_trk_delta_dist_xy; Distance between chosen Tracks in x-y (cm) ; Events / 2E-8 cm",50,0,0.000001);
   TH1F *h_trk_delta_dist_z = new TH1F("h_trk_delta_dist_z","h_trk_delta_dist_z; Distance between chosen Tracks in z (cm) ; Events / 0.05 cm",40,0,2);
   TH1F *h_trackVertex_cos_T = new TH1F("h_trackVertex_cos_T","h_trackVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.05",40,-1,1);
   TH1F *h_trackVertex_alpha_T = new TH1F("h_trackVertex_alpha_T","h_trackVertex_alpha_T; angle bwtn parent momentum and vertex position; Events / 0.315",40,-6.3,6.3);
   TH1F *h_trackVertex_d_T = new TH1F("h_trackVertex_d_T","h_trackVertex_d_T; Impact parameter of parent (cm) ; Events / 0.025 cm",40,0,1);
   TH1F *h_trackVertex_R_T = new TH1F("h_trackVertex_R_T","h_trackVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",40,0,10);
   TH1F *h_trueVertex_cos_T = new TH1F("h_trueVertex_cos_T","h_trueVertex_cos_T; Cos(angle): parent momentum and vertex position; Events / 0.05",40,-1,1);
   TH1F *h_trueVertex_alpha_T = new TH1F("h_trueVertex_alpha_T","h_trueVertex_alpha_T; angle bwtn parent momentum and vertex position; Events / 0.315",40,-6.3,6.3);
   TH1F *h_trueVertex_d_T = new TH1F("h_trueVertex_d_T","h_trueVertex_d_T; Impact parameter of parent (cm) ; Events / 0.005 cm",40,0,0.2);
   TH1F *h_trueVertex_R_T = new TH1F("h_trueVertex_R_T","h_trueVertex_R_T; Tranverse distance of vertex (cm) ; Events / 0.25 cm",40,0,10);

   TH1F *h_res_tp_trk_x = new TH1F("h_res_tp_trk_x","h_res_tp_trk_x; x residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
   TH1F *h_res_tp_trk_y = new TH1F("h_res_tp_trk_y","h_res_tp_trk_y; y residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
   TH1F *h_res_tp_trk_z = new TH1F("h_res_tp_trk_z","h_res_tp_trk_z; z residual of vertex (cm) ; Events / 0.02 cm",1000,-10,10);
   TH1F *h_res_tp_trk_r = new TH1F("h_res_tp_trk_r","h_res_tp_trk_r; r residual of vertex (cm) ; Events / 0.02 cm",100,-1,1);
   TH1F *h_res_tp_trk_phi = new TH1F("h_res_tp_trk_phi","h_res_tp_trk_phi; phi residual of vertex ; Events / 0.02",100,-1,1);

   // Efficiency of Identifying Displaced Vertex Plots
   TH1F *h_all_trueVertex_pt = new TH1F("h_all_trueVertex_pt","h_all_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_trueVertex_minD0 = new TH1F("h_all_trueVertex_minD0","h_all_trueVertex_minD0; minimum d_{0} of TPs from vertex (cm); Events / 0.05 cm",80,-2,2);
   TH1F *h_all_trueVertex_lowPt = new TH1F("h_all_trueVertex_lowPt","h_all_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_oneMatch_trueVertex_pt = new TH1F("h_all_oneMatch_trueVertex_pt","h_all_oneMatch_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_oneMatch_trueVertex_lowPt = new TH1F("h_all_oneMatch_trueVertex_lowPt","h_all_oneMatch_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_oneMatchAlt_trueVertex_pt = new TH1F("h_all_oneMatchAlt_trueVertex_pt","h_all_oneMatchAlt_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_trackVertex_pt = new TH1F("h_all_trackVertex_pt","h_all_trackVertex_pt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_all_trackVertex_minD0 = new TH1F("h_all_trackVertex_minD0","h_all_trackVertex_minD0; minimum d_{0} of trks from vertex (cm); Events / 0.05 cm",80,-2,2);
   TH1F *h_all_trackVertex_lowPt = new TH1F("h_all_trackVertex_lowPt","h_all_trackVertex_lowPt; p_{T} of Lower p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_correct_trueVertex_pt = new TH1F("h_correct_trueVertex_pt","h_correct_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_correct_trueVertex_lowPt = new TH1F("h_correct_trueVertex_lowPt","h_correct_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_correct_oneMatch_trueVertex_pt = new TH1F("h_correct_oneMatch_trueVertex_pt","h_correct_oneMatch_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_correct_oneMatch_trueVertex_lowPt = new TH1F("h_correct_oneMatch_trueVertex_lowPt","h_correct_oneMatch_trueVertex_lowPt; p_{T} of Lower p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_correct_oneMatchAlt_trueVertex_pt = new TH1F("h_correct_oneMatchAlt_trueVertex_pt","h_correct_oneMatchAlt_trueVertex_pt; p_{T} of Leading p_{T} TP from vertex (GeV); Events / 1.0 GeV",100,0,100);
   TH1F *h_false_trackVertex_pt = new TH1F("h_false_trackVertex_pt","h_false_trackVertex_pt; p_{T} of Leading p_{T} trk from vertex (GeV); Events / 1.0 GeV",100,0,100);

   TH1F *h_false_trackVertex_delta_dist_z = new TH1F("h_false_trackVertex_delta_dist_z","h_false_trackVertex_delta_dist_z; Distance in z btwn trks from vertex (cm); Events / 0.5 cm",100,0,5);
   TH1F *h_correct_trackVertex_delta_dist_z = new TH1F("h_correct_trackVertex_delta_dist_z","h_correct_trackVertex_delta_dist_z; Distance in z btwn trks from vertex (cm); Events / 0.5 cm",100,0,5);
   TH1F *h_all_trackVertex_numStubs = new TH1F("h_all_trackVertex_numStubs","h_all_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
   TH1F *h_correct_trackVertex_numStubs = new TH1F("h_correct_trackVertex_numStubs","h_correct_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
   TH1F *h_correct_trackVertex_numStubsSum = new TH1F("h_correct_trackVertex_numStubsSum","h_correct_trackVertex_numStubsSum; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
   TH1F *h_false_trackVertex_numStubs = new TH1F("h_false_trackVertex_numStubs","h_false_trackVertex_numStubs; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
   TH1F *h_false_trackVertex_numStubsSum = new TH1F("h_false_trackVertex_numStubsSum","h_false_trackVertex_numStubsSum; Number of stubs of trks from vertex; Events / 1.0",12,0,12);
   TH1F *h_false_trackVertex_d0 = new TH1F("h_false_trackVertex_d0","h_false_trackVertex_d0; d_{0} of trks from vertex; Events / 0.05 cm",100,0,5);
   TH1F *h_false_trackVertex_d_T = new TH1F("h_false_trackVertex_d_T","h_false_trackVertex_d_T; Impact parameter of parent (cm) ; Events / 0.025 cm",40,0,1);
   TH1F *h_false_trackVertex_chi2rphidofSum = new TH1F("h_false_trackVertex_chi2rphidofSum","h_false_trackVertex_chi2rphidofSum; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_false_trackVertex_chi2rzdofSum = new TH1F("h_false_trackVertex_chi2rzdofSum","h_false_trackVertex_chi2rzdofSum; Track #chi^{2}_{rz}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_false_trackVertex_bendchi2Sum = new TH1F("h_false_trackVertex_bendchi2Sum","h_false_trackVertex_bendchi2Sum; Track bend #chi^{2} ; Events / 0.1",50,0,5);
   TH1F *h_correct_trackVertex_chi2rphidofSum = new TH1F("h_correct_trackVertex_chi2rphidofSum","h_correct_trackVertex_chi2rphidofSum; Track #chi^{2}_{r#phi}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_correct_trackVertex_chi2rzdofSum = new TH1F("h_correct_trackVertex_chi2rzdofSum","h_correct_trackVertex_chi2rzdofSum; Track #chi^{2}_{rz}/d.o.f ; Events / 0.1",50,0,5);
   TH1F *h_correct_trackVertex_bendchi2Sum = new TH1F("h_correct_trackVertex_bendchi2Sum","h_correct_trackVertex_bendchi2Sum; Track bend #chi^{2} ; Events / 0.1",50,0,5);

   TH1F *h_all_trueVertex_eta = new TH1F("h_all_trueVertex_eta","h_all_trueVertex_eta; #eta of Leading p_{T} TP from vertex; Events / 0.096",50,-2.4,2.4);
   TH1F *h_all_oneMatch_trueVertex_eta = new TH1F("h_all_oneMatch_trueVertex_eta","h_all_oneMatch_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_all_oneMatchAlt_trueVertex_eta = new TH1F("h_all_oneMatchAlt_trueVertex_eta","h_all_oneMatchAlt_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_all_trackVertex_eta = new TH1F("h_all_trackVertex_eta","h_all_trackVertex_eta; #eta of Leading p_{T} trk from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_correct_trueVertex_eta = new TH1F("h_correct_trueVertex_eta","h_correct_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_correct_oneMatch_trueVertex_eta = new TH1F("h_correct_oneMatch_trueVertex_eta","h_correct_oneMatch_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_correct_oneMatchAlt_trueVertex_eta = new TH1F("h_correct_oneMatchAlt_trueVertex_eta","h_correct_oneMatchAlt_trueVertex_eta; #eta of Leading p_{T} TP from vertex ; Events / 0.096",50,-2.4,2.4);
   TH1F *h_false_trackVertex_eta = new TH1F("h_false_trackVertex_eta","h_false_trackVertex_eta; #eta of Leading p_{T} trk from vertex ; Events / 0.096",50,-2.4,2.4);
   
   TH1F *h_all_trueVertex_dxy = new TH1F("h_all_trueVertex_dxy","h_all_trueVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
   TH1F *h_all_oneMatch_trueVertex_dxy = new TH1F("h_all_oneMatch_trueVertex_dxy","h_all_oneMatch_trueVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
   TH1F *h_all_oneMatchAlt_trueVertex_dxy = new TH1F("h_all_oneMatchAlt_trueVertex_dxy","h_all_oneMatchAlt_trueVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
   TH1F *h_all_trackVertex_dxy = new TH1F("h_all_trackVertex_dxy","h_all_trackVertex_dxy; dxy of vertex (cm); Events / 0.1 cm",20,0,2);
   TH1F *h_correct_trueVertex_dxy = new TH1F("h_correct_trueVertex_dxy","h_correct_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);
   TH1F *h_correct_oneMatch_trueVertex_dxy = new TH1F("h_correct_oneMatch_trueVertex_dxy","h_correct_oneMatch_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);
   TH1F *h_correct_oneMatchAlt_trueVertex_dxy = new TH1F("h_correct_oneMatchAlt_trueVertex_dxy","h_correct_oneMatchAlt_trueVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);
   TH1F *h_false_trackVertex_dxy = new TH1F("h_false_trackVertex_dxy","h_false_trackVertex_dxy; dxy of vertex (cm) ; Events / 0.1 cm",20,0,2);

   // Trigger Rates Study 
   TH2F *h_Count_trk_pt_d0 = new TH2F("h_Count_trk_pt_d0","h_Count_trk_pt_d0; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks
   TH2F *h_Count_trk_pt_d0_dv = new TH2F("h_Count_trk_pt_d0_dv","h_Count_trk_pt_d0_dv;(DV selection) Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks including the displaced vertex selection
   

   float pt_cuts[5] = {2.0,4.0,6.0,8.0,10.0};     // Cuts to control event rate
   float d0_cuts[5] = {0.1,0.2,0.3,0.4,0.5}; 
   std::vector<double> dxy_cuts = linspace(0.0,0.0000001,100); //cut <dxy
   std::vector<double> dz_cuts = linspace(0.0,3.0,100); //cut <dz
   std::vector<double> cos_T_cuts = linspace(.8,.95,100); //cut >cos
   std::vector<double> cos_T_cuts_2 = linspace(.95,1.0,900);
   cos_T_cuts.pop_back();
   cos_T_cuts.insert(cos_T_cuts.end(),cos_T_cuts_2.begin(),cos_T_cuts_2.end());
   std::vector<double> d_T_cuts = linspace(0.0,0.1,100); //cut <d_T
   std::vector<double> R_T_cuts = linspace(0.0,2.0,100); //cut >R_T
   double true_vertices = 0.0;

   std::vector<double> correct_vert_dxy_cut(dxy_cuts.size()); 
   std::vector<double> false_vert_dxy_cut(dxy_cuts.size());
   std::vector<double> all_vert_dxy_cut(dxy_cuts.size());

   std::vector<double> correct_vert_dz_cut(dz_cuts.size());
   std::vector<double> false_vert_dz_cut(dz_cuts.size());
   std::vector<double> all_vert_dz_cut(dz_cuts.size());

   std::vector<double> correct_vert_cos_T_cut(cos_T_cuts.size());
   std::vector<double> false_vert_cos_T_cut(cos_T_cuts.size());
   std::vector<double> all_vert_cos_T_cut(cos_T_cuts.size());
   
   std::vector<double> correct_vert_d_T_cut(d_T_cuts.size());
   std::vector<double> false_vert_d_T_cut(d_T_cuts.size());
   std::vector<double> all_vert_d_T_cut(d_T_cuts.size());

   std::vector<double> correct_vert_R_T_cut(R_T_cuts.size());
   std::vector<double> false_vert_R_T_cut(R_T_cuts.size());
   std::vector<double> all_vert_R_T_cut(R_T_cuts.size());

   if (fChain == 0) return;

   
   Long64_t nevt = fChain->GetEntries();
   //Long64_t nevt = 10;

   for (Long64_t i_evnt=0; i_evnt<nevt; i_evnt++) {
     //std::cout<<"event number: "<<i_evnt<<std::endl;
       fChain->GetEntry(i_evnt);
       displayProgress(i_evnt, nevt);
       
       std::deque<Track_Parameters> selectedTracks;      // Tracks 
       std::vector<Track_Parameters> selectedTPs;         // Tracking particles
       std::vector<Vertex_Parameters> trueVertices;
       std::vector<Vertex_Parameters> trackVertices;
       std::vector<Vertex_Parameters> tempVertices;
       int maxPT_i = 0;
       int firstMatch_j = 0;
       bool oneMatch = false;
       // ----------------------------------------------------------------------------------------------------------------
       // track loop
       //std::cout<<"test1"<<std::endl;
       for (int it = 0; it < (int)trk_pt->size(); it++){
	 h_trk_eta->Fill(trk_eta->at(it));
	 if (fabs(trk_eta->at(it)) > TP_maxEta)
	   continue;
	 h_trk_pt->Fill(trk_pt->at(it));
	 if (trk_pt->at(it) < TP_minPt)
	   continue;
	 if (trk_pt->at(it) > TP_maxPt)
	   continue;
	 h_trk_d0->Fill(trk_d0->at(it));
	 if (std::fabs(trk_d0->at(it)) > TP_maxD0)
	   continue;
	 if (std::fabs(trk_d0->at(it)) < TP_minD0)
	   continue;
	 
	 int ndof = 2 * trk_nstub->at(it) - 5;
	 float chi2rphidof = (float)trk_chi2rphi->at(it) / ndof;
	 float chi2rzdof = (float)trk_chi2rz->at(it) / ndof;
	     
	 
	 h_trk_chi2rphidof->Fill(chi2rphidof);
	 h_trk_chi2rzdof->Fill(chi2rzdof);
	 h_trk_bendchi2->Fill(trk_bendchi2 ->at(it));
	 
	 if(trk_pt->at(it) > 8){   
	   h_trk_chi2rphidof_H->Fill(chi2rphidof);
	   h_trk_chi2rzdof_H->Fill(chi2rzdof);
	   h_trk_bendchi2_H->Fill(trk_bendchi2 ->at(it));
	 }
	 else{
	   h_trk_chi2rphidof_L->Fill(chi2rphidof);
	   h_trk_chi2rzdof_L->Fill(chi2rzdof);
	   h_trk_bendchi2_L->Fill(trk_bendchi2 ->at(it));
	 }
	 if(fabs(trk_eta->at(it)) < 0.8){
	   h_trk_chi2rphidof_C->Fill(chi2rphidof);
	   h_trk_chi2rzdof_C->Fill(chi2rzdof);
	   h_trk_bendchi2_C->Fill(trk_bendchi2 ->at(it));
	 }
	 else if(fabs(trk_eta->at(it)) <= 1.6){
	   h_trk_chi2rphidof_I->Fill(chi2rphidof);
	   h_trk_chi2rzdof_I->Fill(chi2rzdof);
	   h_trk_bendchi2_I->Fill(trk_bendchi2 ->at(it));
	 }
	 else{
	   h_trk_chi2rphidof_F->Fill(chi2rphidof);
	   h_trk_chi2rzdof_F->Fill(chi2rzdof);
	   h_trk_bendchi2_F->Fill(trk_bendchi2 ->at(it));
	 }
	 if(trk_d0->at(it) <= 1){   
	   h_trk_chi2rphidof_P->Fill(chi2rphidof);
	   h_trk_chi2rzdof_P->Fill(chi2rzdof);
	   h_trk_bendchi2_P->Fill(trk_bendchi2 ->at(it));
	 }
	 else{   
	   h_trk_chi2rphidof_D->Fill(chi2rphidof);
	   h_trk_chi2rzdof_D->Fill(chi2rzdof);
	   h_trk_bendchi2_D->Fill(trk_bendchi2 ->at(it));
	 }
#if 0 
			if(chi2rphidof > 2)
				continue;
			if(chi2rzdof   > 2)
				continue;
			if(trk_bendchi2 ->at(it) > 5)
				continue;
#endif	   
			
			selectedTracks.push_back(Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999, trk_matchtp_pt->at(it) ));
       }

       int pt_check[5] = {0,0,0,0,0};
       int d0_check[5] = {0,0,0,0,0};
       int pt_d0_check[5][5];
       for (int k=0; k < 5; k++){
         for (int l=0; l < 5; l++){
	   pt_d0_check[k][l] = 0;
         }
       }
       
       if ((selectedTracks.size() >= 2)){ 
	 for (int j = 0; j < selectedTracks.size(); j++){
	   // for (int k = j+1; k < selectedTracks->size(); k++){
	   // if((fabs((*selectedTracks)[j]->z(x_dv_trk,y_dv_trk,z_dv_trk) - z_dv_trk)) < 1.0){
	       if (selectedTracks[j].pt   >   pt_cuts[0]){ pt_check[0]++; }
               if (selectedTracks[j].pt   >   pt_cuts[1]){ pt_check[1]++; }
               if (selectedTracks[j].pt   >   pt_cuts[2]){ pt_check[2]++; }
               if (selectedTracks[j].pt   >   pt_cuts[3]){ pt_check[3]++; }
               if (selectedTracks[j].pt   >   pt_cuts[4]){ pt_check[4]++; }

               if (fabs(selectedTracks[j].d0)>d0_cuts[0]){ d0_check[0]++;}
               if (fabs(selectedTracks[j].d0)>d0_cuts[1]){ d0_check[1]++;}
               if (fabs(selectedTracks[j].d0)>d0_cuts[2]){ d0_check[2]++;}
               if (fabs(selectedTracks[j].d0)>d0_cuts[3]){ d0_check[3]++;}
               if (fabs(selectedTracks[j].d0)>d0_cuts[4]){ d0_check[4]++;}
               for (int l=0; l < 5; l++){
                  if (selectedTracks[j].pt   >   pt_cuts[0] && fabs(selectedTracks[j].d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
                  if (selectedTracks[j].pt   >   pt_cuts[1] && fabs(selectedTracks[j].d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
                  if (selectedTracks[j].pt   >   pt_cuts[2] && fabs(selectedTracks[j].d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
                  if (selectedTracks[j].pt   >   pt_cuts[3] && fabs(selectedTracks[j].d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
                  if (selectedTracks[j].pt   >   pt_cuts[4] && fabs(selectedTracks[j].d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
               }
            // }
        // }
      }

      for(int k=0;k<5;k++){
	pt_check[k] = 0;
	d0_check[k] = 0;
	for (int l=0; l < 5; l++){
	  if(pt_d0_check[k][l]>=2){
	    h_Count_trk_pt_d0->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0->GetBinContent(l+1,k+1) + 1));
	  }
	  pt_d0_check[k][l] = 0;
	}
      }
      }
      
      // ----------------------------------------------------------------------------------------------------------------
      // tracking particle loop
       //std::cout<<"test2"<<std::endl;
      for (int it = 0; it < (int)tp_pt->size(); it++){
      
	float tmp_d0 = -tp_d0->at(it);	// Sign difference in the NTupleMaker
	float tmp_z0 = tp_z0->at(it);
	
	h_tp_pt_noCuts->Fill(tp_pt->at(it));
	h_tp_eta_noCuts->Fill(tp_eta->at(it));
	h_tp_d0_noCuts->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_noCuts->Fill(tp_pt->at(it));
	  h_match_tp_eta_noCuts->Fill(tp_eta->at(it));
	  h_match_tp_d0_noCuts->Fill(tp_d0->at(it));
	}
	h_tp_d0->Fill(tp_d0->at(it));
	if (std::fabs(tmp_d0) > TP_maxD0)
	  continue;

	h_tp_pt_maxD0Cut->Fill(tp_pt->at(it));
	h_tp_eta_maxD0Cut->Fill(tp_eta->at(it));
	h_tp_d0_maxD0Cut->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_maxD0Cut->Fill(tp_pt->at(it));
	  h_match_tp_eta_maxD0Cut->Fill(tp_eta->at(it));
	  h_match_tp_d0_maxD0Cut->Fill(tp_d0->at(it));
	}

	if (std::fabs(tmp_d0) < TP_minD0)
	  continue;

	h_tp_pt_minD0Cut->Fill(tp_pt->at(it));
	h_tp_eta_minD0Cut->Fill(tp_eta->at(it));
	h_tp_d0_minD0Cut->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_minD0Cut->Fill(tp_pt->at(it));
	  h_match_tp_eta_minD0Cut->Fill(tp_eta->at(it));
	  h_match_tp_d0_minD0Cut->Fill(tp_d0->at(it));
	}
	h_tp_pt->Fill(tp_pt->at(it));
	if (tp_pt->at(it) < TP_minPt)
	  continue;

	h_tp_pt_minPtCut->Fill(tp_pt->at(it));
	h_tp_eta_minPtCut->Fill(tp_eta->at(it));
	h_tp_d0_minPtCut->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_minPtCut->Fill(tp_pt->at(it));
	  h_match_tp_eta_minPtCut->Fill(tp_eta->at(it));
	  h_match_tp_d0_minPtCut->Fill(tp_d0->at(it));
	}

	if (tp_pt->at(it) > TP_maxPt)
	  continue;

	h_tp_pt_maxPtCut->Fill(tp_pt->at(it));
	h_tp_eta_maxPtCut->Fill(tp_eta->at(it));
	h_tp_d0_maxPtCut->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_maxPtCut->Fill(tp_pt->at(it));
	  h_match_tp_eta_maxPtCut->Fill(tp_eta->at(it));
	  h_match_tp_d0_maxPtCut->Fill(tp_d0->at(it));
	}
	h_tp_eta->Fill(tp_eta->at(it));
	if (std::fabs(tp_eta->at(it)) > TP_maxEta)
	  continue;

	h_tp_pt_maxEtaCut->Fill(tp_pt->at(it));
	h_tp_eta_maxEtaCut->Fill(tp_eta->at(it));
	h_tp_d0_maxEtaCut->Fill(tp_d0->at(it));
	
	if (tp_nmatch->at(it) >= 1){
	  h_match_tp_pt_maxEtaCut->Fill(tp_pt->at(it));
	  h_match_tp_eta_maxEtaCut->Fill(tp_eta->at(it));
	  h_match_tp_d0_maxEtaCut->Fill(tp_d0->at(it));
	}
	//std::cout<<"TP pos: "<<tp_x->at(it)<<" "<<tp_y->at(it)<<" "<<tp_z->at(it)<<" TP mom: "<<tp_pt->at(it)<<std::endl;
	selectedTPs.push_back(Track_Parameters(tp_pt->at(it), tmp_d0, tmp_z0, tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, abs(tp_pdgid->at(it))));
      }
    
      // --------------------------------------------------------------------------------------------
      //         Vertex finding in Tracking Particles
      // --------------------------------------------------------------------------------------------
      if (!(selectedTracks.size() >= 2)) continue;
      bool true_DV = false;
      //std::cout<<"test3"<<std::endl;
      if(VERBOSE[0])
	std::cout<<"End of z-vertex Finding"<<endl;
      //std::cout<<"event"<<std::endl;
      //sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);
      
      Double_t x_dv = -9999.0;// (tp_x->at((*selectedTPs)[0]->index));//+tp_x->at((*selectedTPs)[1]->index))/2.0;
      Double_t y_dv = -9999.0;// (tp_y->at((*selectedTPs)[0]->index));//+tp_y->at((*selectedTPs)[1]->index))/2.0;
      Double_t z_dv = -9999.0;// (tp_z->at((*selectedTPs)[0]->index));//+tp_z->at((*selectedTPs)[1]->index))/2.0;
            
      int Vertex_check = -1;
      int selected_track_j0 = -1;
      int selected_track_j1 = -1;
      std::vector<int> selected_tps_j;

      Double_t x_tmp = x_dv;
      Double_t y_tmp = y_dv;
      Double_t z_tmp = z_dv;
      Int_t Vertex_check_tmp = -1;
      //std::cout<<"selectedTPs size: "<<selectedTPs.size()<<std::endl;
      
      int noCuts_trueVertices = 0;

      if(selectedTPs.size()>=2){
	
	//sort(selectedTPs.begin(), selectedTPs.end(), *this);
	sort(selectedTPs.begin(), selectedTPs.end(), [&](Track_Parameters a, Track_Parameters b) {
	    return tp_x->at(a.index) > tp_x->at(b.index);
	  });
	
	for( uint i=0; i<(selectedTPs.size() - 1); i++ ){
	   int index0 = selectedTPs[i].index;
	   int index1 = selectedTPs[i+1].index;
	   if( fabs(tp_x->at(index0)-tp_x->at(index1))<d0_res && fabs(tp_y->at(index0)-tp_y->at(index1))<d0_res && fabs(tp_z->at(index0)-tp_z->at(index1))<z0_res ){
	     x_dv = (tp_x->at(index0) + tp_x->at(index1))/2;
	     y_dv = (tp_y->at(index0) + tp_y->at(index1))/2;
	     z_dv = (tp_z->at(index0) + tp_z->at(index1))/2;
	     h_trueVertex_x->Fill(x_dv);
	     h_trueVertex_y->Fill(y_dv);
	     h_trueVertex_z->Fill(z_dv);
	     h_trueVertex_sumPt->Fill(selectedTPs[i].pt+selectedTPs[i+1].pt);
	     if(selectedTPs[i].pt>selectedTPs[i+1].pt){
	       h_trueVertex_lowPt->Fill(selectedTPs[i+1].pt);
	       h_trueVertex_highPt->Fill(selectedTPs[i].pt);
	     }
	     else{
	       h_trueVertex_lowPt->Fill(selectedTPs[i].pt);
	       h_trueVertex_highPt->Fill(selectedTPs[i+1].pt);
	     }
	     noCuts_trueVertices++;
	     if(dist(x_dv,y_dv)>d0_res && dist(x_dv,y_dv)<20){
	       //std::cout<<"true vertex: "<<x_dv<<" "<<y_dv<<" "<<z_dv<<" tp_pt: "<<selectedTPs[i].pt<<" "<<selectedTPs[i+1].pt<<std::endl;
	       if(selectedTPs[i].pt>selectedTPs[i+1].pt){
		 trueVertices.push_back(Vertex_Parameters(x_dv, y_dv, z_dv, selectedTPs[i], selectedTPs[i+1]) );
	       }
	       else{
		 trueVertices.push_back(Vertex_Parameters(x_dv, y_dv, z_dv, selectedTPs[i+1], selectedTPs[i]) );
	       }
	       selected_tps_j.push_back(i);
	       i++;
	     }
	   }
	 }
#if 0
	 uint i_dupe = 1;
	 while( i_dupe<trueVertices.size() ){
	   Vertex_Parameters currentVert = trueVertices[i_dupe];
	   Vertex_Parameters lastVert = trueVertices[i_dupe-1];
	   if( (currentVert.x_dv == lastVert.x_dv) && (currentVert.y_dv == lastVert.y_dv) && (currentVert.z_dv == lastVert.z_dv) ){
	     trueVertices.erase(trueVertices.begin()+i_dupe);
	   }
	   else{
	     i_dupe++;
	   }
	 }
#endif
         // Selection of leading 2 p_T tracks from true common vertex
         // if(!(!((fabs(tp_x->at((*selectedTPs)[0]->index) - tp_x->at((*selectedTPs)[1]->index))<0.01) && (fabs(tp_y->at((*selectedTPs)[0]->index) - tp_y->at((*selectedTPs)[1]->index))<0.01)))){nCount_leadingpt++;}// continue;
#if 0
         Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[1],x_dv,y_dv,z_dv); // 0+1
         Vertex_check_tmp = Vertex_check;

         if (Vertex_check>0){
            selected_track_j0 = 0;
            selected_track_j1 = 1;
         } 
         else if(selectedTPs->size()>=3){
            Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[2],x_dv,y_dv,z_dv);
            if(Vertex_check>0){ //0+2
               selected_track_j0 = 0;
               selected_track_j1 = 2;
            }
            else if((Vertex_check=Vertex((*selectedTPs)[1],(*selectedTPs)[2],x_dv,y_dv,z_dv))>0){ //1+2
               selected_track_j0 = 1;
               selected_track_j1 = 2;
            }
            else if(selectedTPs->size()>=4){
               Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[3],x_dv,y_dv,z_dv);
               if(Vertex_check>0){ //0+3
                  selected_track_j0 = 0;
                  selected_track_j1 = 3;
               }
               else if((Vertex_check=Vertex((*selectedTPs)[1],(*selectedTPs)[3],x_dv,y_dv,z_dv))>0){ //1+3
                  selected_track_j0 = 1;
                  selected_track_j1 = 3;
               }
               else if((Vertex_check=Vertex((*selectedTPs)[2],(*selectedTPs)[3],x_dv,y_dv,z_dv))>0){ //2+3
                  selected_track_j0 = 2;
                  selected_track_j1 = 3;
               }
               else{               
                  x_dv = x_tmp;
                  y_dv = y_tmp;
                  z_dv = z_tmp;//(*selectedTPs)[0]->z0;
                  Vertex_check = Vertex_check_tmp;
                  
                  selected_track_j0 = 0;
                  selected_track_j1 = 1;
               }
            }
            else{
               
               x_dv = x_tmp;
               y_dv = y_tmp;
               z_dv = z_tmp;//(*selectedTPs)[0]->z0;
               Vertex_check = Vertex_check_tmp;
               
               selected_track_j0 = 0;
               selected_track_j1 = 1;
            }
         } 
         else{
            x_dv = x_tmp;
            y_dv = y_tmp;
            z_dv = z_tmp;//(*selectedTPs)[0]->z0;
            Vertex_check = Vertex_check_tmp;
            
            selected_track_j0 = 0;
            selected_track_j1 = 1;
         }
         if(VERBOSE[0])
            std::cout<<"Selected Tracks = "<<selected_track_j0<<"\t and  "<<selected_track_j1<<endl;      
#endif

	 //std::cout<<"trueVertices size: "<<trueVertices.size()<<std::endl;
	 h_trueVertex_numNoCuts->Fill(noCuts_trueVertices);
	 h_trueVertex_numAllCuts->Fill(trueVertices.size());
	 true_vertices += trueVertices.size();
	 float maxPT = 0.0;
	 sort(selectedTPs.begin(), selectedTPs.end(), ComparePtTrack);
	 for(uint i=0; i<trueVertices.size(); i++){
	   if(trueVertices[i].a.pt>maxPT){
	     maxPT = trueVertices[i].a.pt;
	     maxPT_i = i;
	   }

	   std::valarray<float> p_trk_1 = calcPVec(trueVertices[i].a,trueVertices[i].x_dv,trueVertices[i].y_dv);
	   std::valarray<float> p_trk_2 = calcPVec(trueVertices[i].b,trueVertices[i].x_dv,trueVertices[i].y_dv);
	   std::valarray<float> p_tot = p_trk_1+p_trk_2;
	   float R_T = TMath::Sqrt(pow(trueVertices[i].x_dv,2)+pow(trueVertices[i].y_dv,2));
	   float cos_T = (p_tot[0]*trueVertices[i].x_dv+p_tot[1]*trueVertices[i].y_dv)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
	   float alpha_T = acos(cos_T);
	   float theta = atan2(p_tot[1],p_tot[0]);
	   float d_T = std::abs(cos(theta)*trueVertices[i].y_dv-sin(theta)*trueVertices[i].x_dv);
	
	   h_trueVertex_d_T->Fill(d_T);
	   h_trueVertex_cos_T->Fill(cos_T);
	   h_trueVertex_alpha_T->Fill(alpha_T);
	   h_trueVertex_R_T->Fill(R_T);
	   

	   uint i_pt = 9999;
	   uint j_pt = 9999;
	   for( uint j=0; j<selectedTPs.size();j++ ){
	     if(trueVertices[i].a.index == selectedTPs[j].index){
	       i_pt = j;
	     }
	     if(trueVertices[i].b.index == selectedTPs[j].index){
	       j_pt = j;
	     }
	   }
	   h_trueVertexCuts_indexPt->Fill(i_pt);
	   h_trueVertexCuts_indexPt->Fill(j_pt);
	   //h_delta_dist_xy->Fill(dist_TPs(selectedTPs[selected_tps_j[i]],selectedTPs[selected_tps_j[i]+1]));
	   h_delta_dist_xy->Fill(dist(tp_x->at(trueVertices[i].a.index), tp_y->at(trueVertices[i].a.index), tp_x->at(trueVertices[i].b.index), tp_y->at(trueVertices[i].b.index) ));
	   //h_error_delta_x->Fill(fabs((tp_x->at((*selectedTPs)[selected_track_j0]->index) - x_dv)));
	   
	   //float z_dv_1 = selectedTPs[selected_tps_j[i]].z(trueVertices[i].x_dv,trueVertices[i].y_dv);
	   //float z_dv_2 = selectedTPs[selected_tps_j[i]+1].z(trueVertices[i].x_dv,trueVertices[i].y_dv);

	   h_delta_dist_z->Fill(fabs(tp_z->at(trueVertices[i].a.index) - tp_z->at(trueVertices[i].b.index) ));
	   //h_error_delta_z->Fill(fabs((tp_z->at((*selectedTPs)[selected_track_j0]->index) - z_dv_1)));
	   
	   if(trueVertices[i].x_dv!=-9999.0){
	     h_trueVertexCuts_x->Fill(trueVertices[i].x_dv);
	     h_trueVertexCuts_y->Fill(trueVertices[i].y_dv);
	     h_trueVertexCuts_z->Fill(trueVertices[i].z_dv);
	     h_trueVertexCuts_sumPt->Fill(trueVertices[i].a.pt+trueVertices[i].b.pt);
	   }

	   h_delta_x->Fill(fabs((tp_x->at(trueVertices[i].a.index) - (tp_x->at(trueVertices[i].b.index)))));
	   true_DV = true;
	   
	   h_all_trueVertex_pt->Fill(trueVertices[i].a.pt);
	   if(fabs(trueVertices[i].a.d0) < fabs(trueVertices[i].b.d0)){
	     h_all_trueVertex_minD0->Fill(trueVertices[i].a.d0);
	   }
	   else{
	     h_all_trueVertex_minD0->Fill(trueVertices[i].b.d0);
	   }
	   h_all_trueVertex_lowPt->Fill(trueVertices[i].b.pt);
	   h_all_trueVertex_eta->Fill(trueVertices[i].a.eta);
	   h_all_trueVertex_dxy->Fill(dist(trueVertices[i].x_dv,trueVertices[i].y_dv));
	   
	   

	 }
	 if(trueVertices.size()>0){
	   h_all_oneMatch_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
	   h_all_oneMatch_trueVertex_lowPt->Fill(trueVertices[maxPT_i].b.pt);
	   h_all_oneMatch_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
	   h_all_oneMatch_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
	 }
      }
      
      // --------------------------------------------------------------------------------------------
      //                Vertex finding in Tracks
      // --------------------------------------------------------------------------------------------
      //std::cout<<"test4"<<std::endl;
      sort(selectedTracks.begin(), selectedTracks.end(), ComparePtTrack);
      std::vector<Track_Parameters> copyTracks;
      for (uint i=0; i<selectedTracks.size(); i++){
	copyTracks.push_back(selectedTracks[i]);
      }
      Double_t x_dv_trk = -9999.0;// (tp_x->at((*selectedTracks)[0]->index));//+tp_x->at((*selectedTracks)[1]->index))/2.0;
      Double_t y_dv_trk = -9999.0;// (tp_y->at((*selectedTracks)[0]->index));//+tp_y->at((*selectedTracks)[1]->index))/2.0;
      Double_t z_dv_trk = -9999.0;// (tp_z->at((*selectedTracks)[0]->index));//+tp_z->at((*selectedTracks)[1]->index))/2.0;
      std::vector<int> selected_tracks_j;
      //std::cout<<"test 1"<<std::endl;
      //std::cout<<"selectedTracks size: "<<selectedTracks.size()<<std::endl;
      int i_sel = 0;
      int noCuts_trackVertices = 0;
      
      if(selectedTracks.size()>4){
	selectedTracks.erase(selectedTracks.begin()+4,selectedTracks.end());
      }
      while ( i_sel<(int(selectedTracks.size())-1) ){
	for( uint j=i_sel+1; j<selectedTracks.size(); j++){
	  int k_match = -1;
	  for( uint k=0; k<trueVertices.size(); k++ ){
	    if( (selectedTracks[i_sel].tp_pt == trueVertices[k].a.pt) && (selectedTracks[j].tp_pt == trueVertices[k].b.pt) ){
	      k_match = k;
	      trueVertices[k_match].errorCode = 2;
	    }
	  }
	  if( dist_TPs( selectedTracks[i_sel], selectedTracks[j] ) == 0 ){
	    if(k_match != -1){
	      trueVertices[k_match].errorCode = 1;
	    }
	    //std::cout<<"i j: "<<i_sel<<" "<<j<<std::endl;
	    
	    Vertex_check = Vertex(selectedTracks[i_sel],selectedTracks[j],x_dv_trk,y_dv_trk,z_dv_trk);
	    //std::cout<<"no cut track vertex: "<<x_dv_trk<<" "<<y_dv_trk<<" "<<z_dv_trk<<" indices: "<<selectedTracks[i_sel].index<<" "<<selectedTracks[j].index<<std::endl;
	    h_trackVertex_x->Fill(x_dv_trk);
	    h_trackVertex_y->Fill(y_dv_trk);
	    h_trackVertex_z->Fill(z_dv_trk);
	    h_trackVertex_sumPt->Fill(selectedTracks[i_sel].pt+selectedTracks[j].pt);
	    noCuts_trackVertices++;

	    std::valarray<float> p_trk_1 = calcPVec(selectedTracks[i_sel],x_dv_trk,y_dv_trk);
	    std::valarray<float> p_trk_2 = calcPVec(selectedTracks[j],x_dv_trk,y_dv_trk);
	    std::valarray<float> p_tot = p_trk_1+p_trk_2;
	    float R_T = TMath::Sqrt(pow(x_dv_trk,2)+pow(y_dv_trk,2));
	    float cos_T = (p_tot[0]*x_dv_trk+p_tot[1]*y_dv_trk)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
	    float alpha_T = acos(cos_T);
	    float theta = atan2(p_tot[1],p_tot[0]);
	    float d_T = std::abs(cos(theta)*y_dv_trk-sin(theta)*x_dv_trk);

	    int ndof_1 = 2 * trk_nstub->at(selectedTracks[i_sel].index) - 5;
	    float chi2rphidof_1 = (float)trk_chi2rphi->at(selectedTracks[i_sel].index) / ndof_1;
	    float chi2rzdof_1 = (float)trk_chi2rz->at(selectedTracks[i_sel].index) / ndof_1;
	    float bendchi2_1 = trk_bendchi2->at(selectedTracks[i_sel].index);
	    int ndof_2 = 2 * trk_nstub->at(selectedTracks[j].index) - 5;
	    float chi2rphidof_2 = (float)trk_chi2rphi->at(selectedTracks[j].index) / ndof_2;
	    float chi2rzdof_2 = (float)trk_chi2rz->at(selectedTracks[j].index) / ndof_2;
	    float bendchi2_2 = trk_bendchi2->at(selectedTracks[j].index);
	    float chi2rphidofSum = chi2rphidof_1 + chi2rphidof_2;
	    float chi2rzdofSum = chi2rzdof_1 + chi2rzdof_2;
	    float bendchi2Sum = bendchi2_1 + bendchi2_2;
	    int numStubsSum = trk_nstub->at(selectedTracks[i_sel].index) + trk_nstub->at(selectedTracks[j].index);

	    if(fabs(selectedTracks[i_sel].z(x_dv_trk,y_dv_trk)-selectedTracks[j].z(x_dv_trk,y_dv_trk))>= 1.25 && k_match != -1){
	      //std::cout<<"delta z too large: "<<fabs(selectedTracks[i_sel].z(x_dv_trk,y_dv_trk)-selectedTracks[j].z(x_dv_trk,y_dv_trk))<<std::endl;
	      trueVertices[k_match].errorCode = 3;
	    }
	    else if(dist(x_dv_trk,y_dv_trk)<=d0_res && k_match != -1){
	      trueVertices[k_match].errorCode = 4;
	    }
	    else if(dist(x_dv_trk,y_dv_trk)>=20 && k_match != -1){
	      trueVertices[k_match].errorCode = 5;
	    }
	    else if(trackVertices.size()>=2 && k_match != -1){
	      trueVertices[k_match].errorCode = 6;
	    }
	    else if(cos_T<=0.995 && k_match != -1){
	      trueVertices[k_match].errorCode = 7;
	    }
	    else if(d_T>=.05 && k_match != -1){
	      trueVertices[k_match].errorCode = 8;
	    }
	    else if((selectedTracks[i_sel].charge+selectedTracks[j].charge)!=0 && k_match != -1){
	      trueVertices[k_match].errorCode = 9;
	    }
	    else if(chi2rzdofSum>=5 && k_match != -1){
	      trueVertices[k_match].errorCode = 10;
	    }
	    else if(numStubsSum<=8 && k_match != -1){
	      trueVertices[k_match].errorCode = 11;
	    }
	    else if((fabs(selectedTracks[i_sel].d0)<=0.05 && fabs(selectedTracks[j].d0)<=0.05) && k_match != -1){
	      trueVertices[k_match].errorCode = 12;
	    }

	    //if(fabs(selectedTracks[i_sel].z(x_dv_trk,y_dv_trk)-selectedTracks[j].z(x_dv_trk,y_dv_trk))< 2 && dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackVertices.size()<2 && R_T>0.25 && cos_T > 0.95 && d_T < .05 && (selectedTracks[i_sel].charge+selectedTracks[j].charge)==0 ){ //vertex cuts go here
	    if(fabs(selectedTracks[i_sel].z(x_dv_trk,y_dv_trk)-selectedTracks[j].z(x_dv_trk,y_dv_trk))< 1.25 && dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20 && trackVertices.size()<2 && cos_T > 0.995 && d_T < .05 && (selectedTracks[i_sel].charge+selectedTracks[j].charge)==0 && chi2rzdofSum<5 && numStubsSum>8 && (fabs(selectedTracks[i_sel].d0)>0.05 || fabs(selectedTracks[j].d0)>0.05) ){ //vertex cuts go here
	      //std::cout<<"z pos at vertex: "<<selectedTracks[i_sel].z(x_dv_trk,y_dv_trk)<<" "<<selectedTracks[j].z(x_dv_trk,y_dv_trk)<<std::endl;
	      trackVertices.push_back(Vertex_Parameters(x_dv_trk, y_dv_trk, z_dv_trk, selectedTracks[i_sel], selectedTracks[j] ) );
	      //std::cout<<"track vertex: "<<x_dv_trk<<" "<<y_dv_trk<<" "<<z_dv_trk<<" tp_pt: "<<selectedTracks[i_sel].tp_pt<<" "<<selectedTracks[j].tp_pt<<std::endl;
	      if(j==(selectedTracks.size()-1)){
		//selectedTracks.pop_back();
	      }
	      else{
		//selectedTracks.erase(selectedTracks.begin()+j);
	      }
	      
	      //break;
	    }
	  }
	}	
	selectedTracks.pop_front();
      }

      //std::cout<<"test 2"<<std::endl;
      //std::cout<<"selectedTracks size: "<<selectedTracks.size()<<std::endl;
#if 0
      for( int i=0; i<(int(selectedTracks.size())-1); i++ ){
	for( uint j=i+1; j<(selectedTracks.size()); j++){
	  Vertex_check = Vertex(selectedTracks[i],selectedTracks[j],x_dv_trk,y_dv_trk,z_dv_trk);
	  if(fabs(selectedTracks[i].z(x_dv_trk,y_dv_trk)-selectedTracks[j].z(x_dv_trk,y_dv_trk))< z0_res && dist(x_dv_trk,y_dv_trk)>d0_res && dist(x_dv_trk,y_dv_trk)<20){ //same vertex cuts as above
	    tempVertices.push_back(Vertex_Parameters(x_dv_trk, y_dv_trk, z_dv_trk, selectedTracks[i], selectedTracks[j]) );
	  }
	}
      }
      //std::cout<<"test 3"<<std::endl;
      sort(tempVertices.begin(), tempVertices.end(), CompareDeltaXY);
      while( tempVertices.size() > 0 ){
	trackVertices.push_back( tempVertices[0] );
	int index0 = tempVertices[0].a.index;
	int index1 = tempVertices[0].b.index;
	int i_temp = 1;
	while( i_temp<int(tempVertices.size()) ){
	  int temp_index0 = tempVertices[i_temp].a.index; 
	  int temp_index1 = tempVertices[i_temp].b.index;
	  if( (index0==temp_index0) || (index0==temp_index1) || (index1==temp_index0) || (index1==temp_index1) ){
	    tempVertices.erase(tempVertices.begin()+i_temp);
	  }
	  else{
	    i_temp++;
	  }
	}
	tempVertices.erase(tempVertices.begin());
      }
      //std::cout<<"test 4"<<std::endl;
#endif
#if 0           
      Vertex_check = -1;
      selected_track_j0 = -1;
      selected_track_j1 = -1;
      Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[1],x_dv_trk,y_dv_trk,z_dv_trk); // 0+1
      x_tmp = x_dv_trk;
      y_tmp = y_dv_trk;
      z_tmp = z_dv_trk;
      Vertex_check_tmp = Vertex_check;

         if (Vertex_check>0){
            selected_track_j0 = 0;
            selected_track_j1 = 1;
            h_trk_Counter_TPcombination->AddBinContent(1); 
         } 
         else if(selectedTracks->size()>=3){
            Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[2],x_dv_trk,y_dv_trk,z_dv_trk);
            if(Vertex_check>0){ //0+2
               selected_track_j0 = 0;
               selected_track_j1 = 2;
               h_trk_Counter_TPcombination->AddBinContent(2); 
            }
            else if((Vertex_check=Vertex((*selectedTracks)[1],(*selectedTracks)[2],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //1+2
               selected_track_j0 = 1;
               selected_track_j1 = 2;
               h_trk_Counterr_TPcombination->AddBinContent(3); 
            }
            else if(selectedTracks->size()>=4){
               Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk);
               if(Vertex_check>0){ //0+3
                  selected_track_j0 = 0;
                  selected_track_j1 = 3;
                  h_trk_Counter_TPcombination->AddBinContent(4); 
               }
               else if((Vertex_check=Vertex((*selectedTracks)[1],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //1+3
                  selected_track_j0 = 1;
                  selected_track_j1 = 3;
                  h_trk_Counter_TPcombination->AddBinContent(5); 
               }
               else if((Vertex_check=Vertex((*selectedTracks)[2],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //2+3
                  selected_track_j0 = 2;
                  selected_track_j1 = 3;
                  h_trk_Counter_TPcombination->AddBinContent(6); 
               }
               else{
                  
                  x_dv_trk = x_tmp;
                  y_dv_trk = y_tmp;
                  z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
                  Vertex_check = Vertex_check_tmp;
                  
                  selected_track_j0 = 0;
                  selected_track_j1 = 1;
                  h_trk_Counter_TPcombination->AddBinContent(0); 
               }
            }
            else{
               
               x_dv_trk = x_tmp;
               y_dv_trk = y_tmp;
               z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
               Vertex_check = Vertex_check_tmp;
               
               selected_track_j0 = 0;
               selected_track_j1 = 1;
               h_trk_Counter_TPcombination->AddBinContent(0); 
            }
         } 
         else{
            x_dv_trk = x_tmp;
            y_dv_trk = y_tmp;
            z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
            Vertex_check = Vertex_check_tmp;
            
            selected_track_j0 = 0;
            selected_track_j1 = 1;
            h_trk_Counter_TPcombination->AddBinContent(0); 
         }
#endif
   
	 //std::cout<<"trackVertices size: "<<trackVertices.size()<<std::endl;
	 h_trackVertex_numNoCuts->Fill(noCuts_trackVertices);
	 h_trackVertex_numAllCuts->Fill(trackVertices.size());
      for( uint i=0; i<trackVertices.size(); i++ ){
	uint i_pt = 9999;
	uint j_pt = 9999;
	for( uint j=0; j<copyTracks.size();j++ ){
	  if(trackVertices[i].a.index == copyTracks[j].index){
	    i_pt = j;
	  }
	  if(trackVertices[i].b.index == copyTracks[j].index){
	    j_pt = j;
	  }
	}
	h_trackVertexCuts_indexPt->Fill(i_pt);
	h_trackVertexCuts_indexPt->Fill(j_pt);

	h_trk_delta_dist_xy->Fill(dist_TPs(trackVertices[i].a,trackVertices[i].b));
      
	float z_dv_trk_1 = trackVertices[i].a.z(trackVertices[i].x_dv,trackVertices[i].y_dv);
	float z_dv_trk_2 = trackVertices[i].b.z(trackVertices[i].x_dv,trackVertices[i].y_dv);
      
	std::valarray<float> p_trk_1 = calcPVec(trackVertices[i].a,trackVertices[i].x_dv,trackVertices[i].y_dv);
	std::valarray<float> p_trk_2 = calcPVec(trackVertices[i].b,trackVertices[i].x_dv,trackVertices[i].y_dv);
	std::valarray<float> p_tot = p_trk_1+p_trk_2;
	float R_T = TMath::Sqrt(pow(trackVertices[i].x_dv,2)+pow(trackVertices[i].y_dv,2));
	float cos_T = (p_tot[0]*trackVertices[i].x_dv+p_tot[1]*trackVertices[i].y_dv)/(R_T*TMath::Sqrt(pow(p_tot[0],2)+pow(p_tot[1],2)));
	float alpha_T = acos(cos_T);
	float theta = atan2(p_tot[1],p_tot[0]);
	float d_T = std::abs(cos(theta)*trackVertices[i].y_dv-sin(theta)*trackVertices[i].x_dv);
	
	int ndof_1 = 2 * trk_nstub->at(trackVertices[i].a.index) - 5;
	float chi2rphidof_1 = (float)trk_chi2rphi->at(trackVertices[i].a.index) / ndof_1;
	float chi2rzdof_1 = (float)trk_chi2rz->at(trackVertices[i].a.index) / ndof_1;
	float bendchi2_1 = trk_bendchi2->at(trackVertices[i].a.index);
	int ndof_2 = 2 * trk_nstub->at(trackVertices[i].b.index) - 5;
	float chi2rphidof_2 = (float)trk_chi2rphi->at(trackVertices[i].b.index) / ndof_2;
	float chi2rzdof_2 = (float)trk_chi2rz->at(trackVertices[i].b.index) / ndof_2;
	float bendchi2_2 = trk_bendchi2->at(trackVertices[i].b.index);
	float chi2rphidofSum = chi2rphidof_1 + chi2rphidof_2;
	float chi2rzdofSum = chi2rzdof_1 + chi2rzdof_2;
	float bendchi2Sum = bendchi2_1 + bendchi2_2;

	h_trk_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2)); //* Use this condition
	h_trackVertex_d_T->Fill(d_T);
	h_trackVertex_cos_T->Fill(cos_T);
	h_trackVertex_alpha_T->Fill(alpha_T);
	h_trackVertex_R_T->Fill(R_T);
	
	if(trackVertices[i].x_dv!=-9999.0){
	  h_trackVertexCuts_x->Fill(trackVertices[i].x_dv);
	  h_trackVertexCuts_y->Fill(trackVertices[i].y_dv);
	  h_trackVertexCuts_z->Fill(trackVertices[i].z_dv);
	  h_trackVertexCuts_sumPt->Fill(trackVertices[i].a.pt+trackVertices[i].b.pt);


	  h_all_trackVertex_pt->Fill(trackVertices[i].a.pt);
	  if(fabs(trackVertices[i].a.d0) < fabs(trackVertices[i].b.d0)){
	     h_all_trackVertex_minD0->Fill(trackVertices[i].a.d0);
	   }
	   else{
	     h_all_trackVertex_minD0->Fill(trackVertices[i].b.d0);
	   }
	  h_all_trackVertex_lowPt->Fill(trackVertices[i].b.pt);
	  h_all_trackVertex_eta->Fill(trackVertices[i].a.eta);
	  h_all_trackVertex_dxy->Fill(dist(trackVertices[i].x_dv,trackVertices[i].y_dv));
	  h_all_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	  h_all_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));

	  for(uint k=0;k<dxy_cuts.size();k++){
	    if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
	      all_vert_dxy_cut[k]++;
	    }
	  }
	  for(uint k=0;k<dz_cuts.size();k++){
	    if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k]){
	      all_vert_dz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<cos_T_cuts.size();k++){
	    if(cos_T > cos_T_cuts[k]){
	      all_vert_cos_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<d_T_cuts.size();k++){
	    if(d_T < d_T_cuts[k]){
	      all_vert_d_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<R_T_cuts.size();k++){
	    if(R_T > R_T_cuts[k]){
	      all_vert_R_T_cut[k]++;
	    }
	  }
	   
	}
	
      
	uint matched_j = 0;
	bool foundMatch = false;
	for(uint j=0; j<trueVertices.size(); j++){
	  if( (trackVertices[i].a.tp_pt == trueVertices[j].a.pt) && (trackVertices[i].b.tp_pt == trueVertices[j].b.pt) && trueVertices[j].matched==false){
	    trueVertices[j].matched = true;
	    matched_j = j;
	    foundMatch = true;
	    if(oneMatch==false){
	      oneMatch = true;
	      firstMatch_j = j;
	    }
	    
	    break;
	  }
	}
	//std::cout<<"track vertex: "<<trackVertices[i].x_dv<<" "<<trackVertices[i].y_dv<<" "<<trackVertices[i].z_dv<<" tp_pt: "<<trackVertices[i].a.tp_pt<<" "<<trackVertices[i].b.tp_pt<<" matched: "<<foundMatch<<std::endl;
	if(true_DV && foundMatch){
	  //h_correct_trk_pt->Fill(trackVertices[i].a.pt);
	  //h_correct_trk_eta->Fill(trackVertices[i].a.eta);
	  //h_correct_trk_dxy->Fill(dist(trackVertices[i].x_dv,trackVertices[i].y_dv));
	  h_correct_trueVertex_pt->Fill(trueVertices[matched_j].a.pt);
	  h_correct_trueVertex_lowPt->Fill(trueVertices[matched_j].b.pt);
	  h_correct_trueVertex_eta->Fill(trueVertices[matched_j].a.eta);
	  h_correct_trueVertex_dxy->Fill(dist(trueVertices[matched_j].x_dv,trueVertices[matched_j].y_dv));
	  h_correct_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	  h_correct_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));
	  h_correct_trackVertex_numStubsSum->Fill(trk_nstub->at(trackVertices[i].a.index)+trk_nstub->at(trackVertices[i].b.index));
	  h_correct_trackVertex_chi2rphidofSum->Fill(chi2rphidofSum);
	  h_correct_trackVertex_chi2rzdofSum->Fill(chi2rzdofSum);
	  h_correct_trackVertex_bendchi2Sum->Fill(bendchi2Sum);
	  h_correct_trackVertex_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2));

	  for(uint k=0;k<dxy_cuts.size();k++){
	    if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
	      correct_vert_dxy_cut[k]++;
	    }
	  }
	  for(uint k=0;k<dz_cuts.size();k++){
	    if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k]){
	      correct_vert_dz_cut[k]++;
	    }
	  }
	  for(uint k=0;k<cos_T_cuts.size();k++){
	    if(cos_T > cos_T_cuts[k]){
	      correct_vert_cos_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<d_T_cuts.size();k++){
	    if(d_T < d_T_cuts[k]){
	      correct_vert_d_T_cut[k]++;
	    }
	  }
	  for(uint k=0;k<R_T_cuts.size();k++){
	    if(R_T > R_T_cuts[k]){
	      correct_vert_R_T_cut[k]++;
	    }
	  }
	  
	  h_res_tp_trk_x->Fill(trueVertices[matched_j].x_dv-trackVertices[i].x_dv);
	  h_res_tp_trk_y->Fill(trueVertices[matched_j].y_dv-trackVertices[i].y_dv);
	  h_res_tp_trk_z->Fill(trueVertices[matched_j].z_dv-trackVertices[i].z_dv);
	  float r_j = dist(trueVertices[matched_j].x_dv,trueVertices[matched_j].y_dv);
	  float r_i = dist(trackVertices[i].x_dv,trackVertices[i].y_dv);
	  float phi_j = atan2(trueVertices[matched_j].y_dv,trueVertices[matched_j].x_dv);
	  float phi_i = atan2(trackVertices[i].y_dv,trackVertices[i].x_dv);
          h_res_tp_trk_r->Fill(r_j-r_i);
	  h_res_tp_trk_phi->Fill(phi_j-phi_i);
       
	  if (trackVertices[i].a.pt   >   pt_cuts[0]){ pt_check[0]++; }
	  if (trackVertices[i].a.pt   >   pt_cuts[1]){ pt_check[1]++; }
	  if (trackVertices[i].a.pt   >   pt_cuts[2]){ pt_check[2]++; }
	  if (trackVertices[i].a.pt   >   pt_cuts[3]){ pt_check[3]++; }
	  if (trackVertices[i].a.pt   >   pt_cuts[4]){ pt_check[4]++; }
	  
	  if (fabs(trackVertices[i].a.d0)>d0_cuts[0]){ d0_check[0]++;}
	  if (fabs(trackVertices[i].a.d0)>d0_cuts[1]){ d0_check[1]++;}
	  if (fabs(trackVertices[i].a.d0)>d0_cuts[2]){ d0_check[2]++;}
	  if (fabs(trackVertices[i].a.d0)>d0_cuts[3]){ d0_check[3]++;}
	  if (fabs(trackVertices[i].a.d0)>d0_cuts[4]){ d0_check[4]++;}
	  for (int l=0; l < 5; l++){
	    if (trackVertices[i].a.pt   >   pt_cuts[0] && fabs(trackVertices[i].a.d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
	    if (trackVertices[i].a.pt   >   pt_cuts[1] && fabs(trackVertices[i].a.d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
	    if (trackVertices[i].a.pt   >   pt_cuts[2] && fabs(trackVertices[i].a.d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
	    if (trackVertices[i].a.pt   >   pt_cuts[3] && fabs(trackVertices[i].a.d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
	    if (trackVertices[i].a.pt   >   pt_cuts[4] && fabs(trackVertices[i].a.d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
	  }

	  if (trackVertices[i].b.pt   >   pt_cuts[0]){ pt_check[0]++; }
	  if (trackVertices[i].b.pt   >   pt_cuts[1]){ pt_check[1]++; }
	  if (trackVertices[i].b.pt   >   pt_cuts[2]){ pt_check[2]++; }
	  if (trackVertices[i].b.pt   >   pt_cuts[3]){ pt_check[3]++; }
	  if (trackVertices[i].b.pt   >   pt_cuts[4]){ pt_check[4]++; }
	  
	  if (fabs(trackVertices[i].b.d0)>d0_cuts[0]){ d0_check[0]++;}
	  if (fabs(trackVertices[i].b.d0)>d0_cuts[1]){ d0_check[1]++;}
	  if (fabs(trackVertices[i].b.d0)>d0_cuts[2]){ d0_check[2]++;}
	  if (fabs(trackVertices[i].b.d0)>d0_cuts[3]){ d0_check[3]++;}
	  if (fabs(trackVertices[i].b.d0)>d0_cuts[4]){ d0_check[4]++;}
	  for (int l=0; l < 5; l++){
	    if (trackVertices[i].b.pt   >   pt_cuts[0] && fabs(trackVertices[i].b.d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
	    if (trackVertices[i].b.pt   >   pt_cuts[1] && fabs(trackVertices[i].b.d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
	    if (trackVertices[i].b.pt   >   pt_cuts[2] && fabs(trackVertices[i].b.d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
	    if (trackVertices[i].b.pt   >   pt_cuts[3] && fabs(trackVertices[i].b.d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
	    if (trackVertices[i].b.pt   >   pt_cuts[4] && fabs(trackVertices[i].b.d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
	  }
		    
	  for(int k=0;k<5;k++){
            pt_check[k] = 0;
            d0_check[k] = 0;
            for (int l=0; l < 5; l++){
	      if(pt_d0_check[k][l]>=2){
		h_Count_trk_pt_d0_dv->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0_dv->GetBinContent(l+1,k+1) + 1));
	      }
	      pt_d0_check[k][l] = 0;
            }
	  }
	}
	else{
	  if(trackVertices[i].x_dv!=-9999.0){
	    h_false_trackVertex_pt->Fill(trackVertices[i].a.pt);
	    h_false_trackVertex_eta->Fill(trackVertices[i].a.eta);
	    h_false_trackVertex_dxy->Fill(dist(trackVertices[i].x_dv,trackVertices[i].y_dv));
	    h_false_trackVertex_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2));
	    h_false_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].a.index));
	    h_false_trackVertex_numStubs->Fill(trk_nstub->at(trackVertices[i].b.index));
	    h_false_trackVertex_numStubsSum->Fill(trk_nstub->at(trackVertices[i].a.index)+trk_nstub->at(trackVertices[i].b.index));
	    h_false_trackVertex_chi2rphidofSum->Fill(chi2rphidofSum);
	    h_false_trackVertex_chi2rzdofSum->Fill(chi2rzdofSum);
	    h_false_trackVertex_bendchi2Sum->Fill(bendchi2Sum);
	    h_false_trackVertex_d0->Fill(trackVertices[i].a.d0);
	    h_false_trackVertex_d0->Fill(trackVertices[i].b.d0);
	    h_false_trackVertex_d_T->Fill(d_T);

	    for(uint k=0;k<dxy_cuts.size();k++){
	      if(dist_TPs(trackVertices[i].a,trackVertices[i].b) < dxy_cuts[k]){
		false_vert_dxy_cut[k]++;
	      }
	    }
	    for(uint k=0;k<dz_cuts.size();k++){
	      if(fabs(z_dv_trk_1-z_dv_trk_2) < dz_cuts[k]){
		false_vert_dz_cut[k]++;
	      }
	    }
	    for(uint k=0;k<cos_T_cuts.size();k++){
	      if(cos_T > cos_T_cuts[k]){
		false_vert_cos_T_cut[k]++;
	      }
	    }
	    for(uint k=0;k<d_T_cuts.size();k++){
	      if(d_T < d_T_cuts[k]){
		false_vert_d_T_cut[k]++;
	      }
	    }
	    for(uint k=0;k<R_T_cuts.size();k++){
	      if(R_T > R_T_cuts[k]){
		false_vert_R_T_cut[k]++;
	      }
	    }
	    
	  }
	}
      } // End of Loop of TrackVertices
      if(oneMatch){
	h_correct_oneMatch_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
	h_correct_oneMatch_trueVertex_lowPt->Fill(trueVertices[maxPT_i].b.pt);
	h_correct_oneMatch_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
	h_correct_oneMatch_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
	h_correct_oneMatchAlt_trueVertex_pt->Fill(trueVertices[firstMatch_j].a.pt);
	h_correct_oneMatchAlt_trueVertex_eta->Fill(trueVertices[firstMatch_j].a.eta);
	h_correct_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[firstMatch_j].x_dv,trueVertices[firstMatch_j].y_dv));
	h_all_oneMatchAlt_trueVertex_pt->Fill(trueVertices[firstMatch_j].a.pt);
	h_all_oneMatchAlt_trueVertex_eta->Fill(trueVertices[firstMatch_j].a.eta);
	h_all_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[firstMatch_j].x_dv,trueVertices[firstMatch_j].y_dv));
      }
      else{
	if(trueVertices.size()>0){
	  h_all_oneMatchAlt_trueVertex_pt->Fill(trueVertices[maxPT_i].a.pt);
	  h_all_oneMatchAlt_trueVertex_eta->Fill(trueVertices[maxPT_i].a.eta);
	  h_all_oneMatchAlt_trueVertex_dxy->Fill(dist(trueVertices[maxPT_i].x_dv,trueVertices[maxPT_i].y_dv));
	}
      }
      for(uint i=0; i<trueVertices.size(); i++){
	int ai = trueVertices[i].a.index;
	int bi = trueVertices[i].b.index;
	h_trueVertex_errorCode->Fill(trueVertices[i].errorCode);
	//std::cout<<"true vertex: "<<trueVertices[i].x_dv<<" "<<trueVertices[i].y_dv<<" "<<trueVertices[i].z_dv<<" tp_pt: "<<trueVertices[i].a.pt<<" "<<trueVertices[i].b.pt<<" matched: "<<trueVertices[i].matched<<" tp pos 1: "<<tp_x->at(ai)<<" "<<tp_y->at(ai)<<" "<<tp_z->at(ai)<<" tp pos 2: "<<tp_x->at(bi)<<" "<<tp_y->at(bi)<<" "<<tp_z->at(bi)<<std::endl;
      }
   } // End of Event Loop

   // ---------------------------------------------------------------------------------------------------------
   //some Histograms

   char ctxt[500];
   if(type.Contains("cT0")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=0mm");
   }
   else if(type.Contains("cT10000")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=10000mm");
   }
   else if(type.Contains("cT5000")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=5000mm");
   }   
   else if(type.Contains("cT100")){
      sprintf(ctxt, "Dark Photon, PU=0, #tau=100mm");
   }
   else if(type.Contains("cT10")){
      if(type.Contains("PU200")){
         sprintf(ctxt, "Dark Photon, PU=200, #tau=10mm");
      }
      else{
         sprintf(ctxt, "Dark Photon, PU=0, #tau=10mm");
      }
   }
   else if(type.Contains("NeutrinoGun")){
      sprintf(ctxt, "Neutrino Gun, PU=200");
   }
   else if(type.Contains("DisplacedMu")){
      if(type.Contains("PU200")){
         sprintf(ctxt, "Displaced Mu, PU=200");
      }
      else{
         sprintf(ctxt, "Displaced Mu, PU=0");
      }
   }
   else{
      sprintf(ctxt, "");
   }
   TCanvas c;

   TString DIR = type_dir + "AnalyzerTrkPlots/";
   TString makedir = "mkdir -p " + DIR;
   const char *mkDIR = makedir.Data();
   gSystem->Exec(mkDIR);

   TFile *fout;
   fout = new TFile(type_dir + "output_" + type + ".root", "recreate");

   h_trk_d0->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_d0->GetXaxis()->SetRange(1, h_trk_d0->GetNbinsX() + 2);
   h_trk_d0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_d0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_d0->GetName() + ".jpeg");
   delete h_trk_d0;

   h_tp_d0->GetYaxis()->SetNoExponent(kTRUE);
   h_tp_d0->GetXaxis()->SetRange(1, h_tp_d0->GetNbinsX() + 2);
   h_tp_d0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_tp_d0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_tp_d0->GetName() + ".jpeg");
   delete h_tp_d0;

   h_trk_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_pt->GetXaxis()->SetRange(1, h_trk_pt->GetNbinsX() + 2);
   h_trk_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_pt->GetName() + ".jpeg");
   delete h_trk_pt;

   h_tp_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_tp_pt->GetXaxis()->SetRange(1, h_tp_pt->GetNbinsX() + 2);
   h_tp_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_tp_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_tp_pt->GetName() + ".jpeg");
   delete h_tp_pt;

   h_trk_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_eta->GetXaxis()->SetRange(1, h_trk_eta->GetNbinsX() + 2);
   h_trk_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_eta->GetName() + ".jpeg");
   delete h_trk_eta;

   h_tp_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_tp_eta->GetXaxis()->SetRange(1, h_tp_eta->GetNbinsX() + 2);
   h_tp_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_tp_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_tp_eta->GetName() + ".jpeg");
   delete h_tp_eta;

   h_trueVertex_errorCode->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_errorCode->GetXaxis()->SetRange(1, h_trueVertex_errorCode->GetNbinsX() + 2);
   h_trueVertex_errorCode->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_errorCode->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_errorCode->GetName() + ".jpeg");
   delete h_trueVertex_errorCode;

   h_trueVertex_numAllCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_numAllCuts->GetXaxis()->SetRange(1, h_trueVertex_numAllCuts->GetNbinsX() + 2);
   h_trueVertex_numAllCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_numAllCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_numAllCuts->GetName() + ".jpeg");
   delete h_trueVertex_numAllCuts;

   h_trueVertex_numNoCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_numNoCuts->GetXaxis()->SetRange(1, h_trueVertex_numNoCuts->GetNbinsX() + 2);
   h_trueVertex_numNoCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_numNoCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_numNoCuts->GetName() + ".jpeg");
   delete h_trueVertex_numNoCuts;

   h_trackVertex_numAllCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_numAllCuts->GetXaxis()->SetRange(1, h_trackVertex_numAllCuts->GetNbinsX() + 2);
   //c.SetLogy();
   h_trackVertex_numAllCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_numAllCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_numAllCuts->GetName() + ".jpeg");
   delete h_trackVertex_numAllCuts;

   h_trackVertex_numNoCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_numNoCuts->GetXaxis()->SetRange(1, h_trackVertex_numNoCuts->GetNbinsX() + 2);
   h_trackVertex_numNoCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_numNoCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_numNoCuts->GetName() + ".jpeg");
   delete h_trackVertex_numNoCuts;

   h_trueVertex_x->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_x->GetXaxis()->SetRange(1, h_trueVertex_x->GetNbinsX() + 2);
   h_trueVertex_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_x->GetName() + ".jpeg");
   delete h_trueVertex_x;

   h_trueVertex_y->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_y->GetXaxis()->SetRange(1, h_trueVertex_y->GetNbinsX() + 2);
   h_trueVertex_y->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_y->GetName() + ".jpeg");
   delete h_trueVertex_y;

   h_trueVertex_z->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_z->GetXaxis()->SetRange(1, h_trueVertex_z->GetNbinsX() + 2);
   h_trueVertex_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_z->GetName() + ".jpeg");
   delete h_trueVertex_z; 

   h_trueVertex_sumPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_sumPt->GetXaxis()->SetRange(1, h_trueVertex_sumPt->GetNbinsX() + 2);
   h_trueVertex_sumPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_sumPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_sumPt->GetName() + ".jpeg");
   delete h_trueVertex_sumPt;

   h_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_lowPt->GetXaxis()->SetRange(1, h_trueVertex_lowPt->GetNbinsX() + 2);
   h_trueVertex_lowPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_lowPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_lowPt->GetName() + ".jpeg");
   delete h_trueVertex_lowPt;

   h_trueVertex_highPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_highPt->GetXaxis()->SetRange(1, h_trueVertex_highPt->GetNbinsX() + 2);
   h_trueVertex_highPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_highPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_highPt->GetName() + ".jpeg");
   delete h_trueVertex_highPt;

   h_false_trackVertex_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_delta_dist_z->GetXaxis()->SetRange(1, h_false_trackVertex_delta_dist_z->GetNbinsX() + 2);
   h_false_trackVertex_delta_dist_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_delta_dist_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_delta_dist_z->GetName() + ".jpeg");
   delete h_false_trackVertex_delta_dist_z;

   h_correct_trackVertex_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_delta_dist_z->GetXaxis()->SetRange(1, h_correct_trackVertex_delta_dist_z->GetNbinsX() + 2);
   h_correct_trackVertex_delta_dist_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_delta_dist_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_delta_dist_z->GetName() + ".jpeg");
   delete h_correct_trackVertex_delta_dist_z;

   h_false_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_numStubs->GetXaxis()->SetRange(1, h_false_trackVertex_numStubs->GetNbinsX() + 2);
   h_false_trackVertex_numStubs->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_numStubs->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_numStubs->GetName() + ".jpeg");
   delete h_false_trackVertex_numStubs;

   h_false_trackVertex_numStubsSum->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_numStubsSum->GetXaxis()->SetRange(1, h_false_trackVertex_numStubsSum->GetNbinsX() + 2);
   h_false_trackVertex_numStubsSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_numStubsSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_numStubsSum->GetName() + ".jpeg");
   delete h_false_trackVertex_numStubsSum;

   h_false_trackVertex_chi2rphidofSum->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_chi2rphidofSum->GetXaxis()->SetRange(1, h_false_trackVertex_chi2rphidofSum->GetNbinsX() + 2);
   h_false_trackVertex_chi2rphidofSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_chi2rphidofSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_chi2rphidofSum->GetName() + ".jpeg");
   delete h_false_trackVertex_chi2rphidofSum;

   h_false_trackVertex_chi2rzdofSum->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_chi2rzdofSum->GetXaxis()->SetRange(1, h_false_trackVertex_chi2rzdofSum->GetNbinsX() + 2);
   h_false_trackVertex_chi2rzdofSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_chi2rzdofSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_chi2rzdofSum->GetName() + ".jpeg");
   delete h_false_trackVertex_chi2rzdofSum;

   h_false_trackVertex_bendchi2Sum->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_bendchi2Sum->GetXaxis()->SetRange(1, h_false_trackVertex_bendchi2Sum->GetNbinsX() + 2);
   h_false_trackVertex_bendchi2Sum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_bendchi2Sum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_bendchi2Sum->GetName() + ".jpeg");
   delete h_false_trackVertex_bendchi2Sum;

   h_correct_trackVertex_chi2rphidofSum->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_chi2rphidofSum->GetXaxis()->SetRange(1, h_correct_trackVertex_chi2rphidofSum->GetNbinsX() + 2);
   h_correct_trackVertex_chi2rphidofSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_chi2rphidofSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_chi2rphidofSum->GetName() + ".jpeg");
   delete h_correct_trackVertex_chi2rphidofSum;

   h_correct_trackVertex_chi2rzdofSum->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_chi2rzdofSum->GetXaxis()->SetRange(1, h_correct_trackVertex_chi2rzdofSum->GetNbinsX() + 2);
   h_correct_trackVertex_chi2rzdofSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_chi2rzdofSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_chi2rzdofSum->GetName() + ".jpeg");
   delete h_correct_trackVertex_chi2rzdofSum;

   h_correct_trackVertex_bendchi2Sum->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_bendchi2Sum->GetXaxis()->SetRange(1, h_correct_trackVertex_bendchi2Sum->GetNbinsX() + 2);
   h_correct_trackVertex_bendchi2Sum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_bendchi2Sum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_bendchi2Sum->GetName() + ".jpeg");
   delete h_correct_trackVertex_bendchi2Sum;

   h_all_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
   h_all_trackVertex_numStubs->GetXaxis()->SetRange(1, h_all_trackVertex_numStubs->GetNbinsX() + 2);
   h_all_trackVertex_numStubs->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trackVertex_numStubs->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trackVertex_numStubs->GetName() + ".jpeg");
   delete h_all_trackVertex_numStubs;

   h_correct_trackVertex_numStubs->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_numStubs->GetXaxis()->SetRange(1, h_correct_trackVertex_numStubs->GetNbinsX() + 2);
   h_correct_trackVertex_numStubs->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_numStubs->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_numStubs->GetName() + ".jpeg");
   delete h_correct_trackVertex_numStubs;

   h_correct_trackVertex_numStubsSum->GetYaxis()->SetNoExponent(kTRUE);
   h_correct_trackVertex_numStubsSum->GetXaxis()->SetRange(1, h_correct_trackVertex_numStubsSum->GetNbinsX() + 2);
   h_correct_trackVertex_numStubsSum->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trackVertex_numStubsSum->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trackVertex_numStubsSum->GetName() + ".jpeg");
   delete h_correct_trackVertex_numStubsSum;

   h_false_trackVertex_d0->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_d0->GetXaxis()->SetRange(1, h_false_trackVertex_d0->GetNbinsX() + 2);
   h_false_trackVertex_d0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_d0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_d0->GetName() + ".jpeg");
   delete h_false_trackVertex_d0;

   h_false_trackVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
   h_false_trackVertex_d_T->GetXaxis()->SetRange(1, h_false_trackVertex_d_T->GetNbinsX() + 2);
   h_false_trackVertex_d_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_false_trackVertex_d_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_false_trackVertex_d_T->GetName() + ".jpeg");
   delete h_false_trackVertex_d_T;

   h_trackVertex_x->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_x->GetXaxis()->SetRange(1, h_trackVertex_x->GetNbinsX() + 2);
   h_trackVertex_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_x->GetName() + ".jpeg");
   delete h_trackVertex_x;

   h_trackVertex_y->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_y->GetXaxis()->SetRange(1, h_trackVertex_y->GetNbinsX() + 2);
   h_trackVertex_y->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_y->GetName() + ".jpeg");
   delete h_trackVertex_y;

   h_trackVertex_z->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_z->GetXaxis()->SetRange(1, h_trackVertex_z->GetNbinsX() + 2);
   h_trackVertex_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_z->GetName() + ".jpeg");
   delete h_trackVertex_z;

   h_trackVertex_sumPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_sumPt->GetXaxis()->SetRange(1, h_trackVertex_sumPt->GetNbinsX() + 2);
   h_trackVertex_sumPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_sumPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_sumPt->GetName() + ".jpeg");
   delete h_trackVertex_sumPt;

   h_trackVertexCuts_x->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertexCuts_x->GetXaxis()->SetRange(1, h_trackVertexCuts_x->GetNbinsX() + 2);
   h_trackVertexCuts_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertexCuts_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertexCuts_x->GetName() + ".jpeg");
   delete h_trackVertexCuts_x;

   h_trackVertexCuts_y->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertexCuts_y->GetXaxis()->SetRange(1, h_trackVertexCuts_y->GetNbinsX() + 2);
   h_trackVertexCuts_y->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertexCuts_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertexCuts_y->GetName() + ".jpeg");
   delete h_trackVertexCuts_y;

   h_trackVertexCuts_z->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertexCuts_z->GetXaxis()->SetRange(1, h_trackVertexCuts_z->GetNbinsX() + 2);
   h_trackVertexCuts_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertexCuts_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertexCuts_z->GetName() + ".jpeg");
   delete h_trackVertexCuts_z;

   h_trackVertexCuts_sumPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertexCuts_sumPt->GetXaxis()->SetRange(1, h_trackVertexCuts_sumPt->GetNbinsX() + 2);
   h_trackVertexCuts_sumPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertexCuts_sumPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertexCuts_sumPt->GetName() + ".jpeg");
   delete h_trackVertexCuts_sumPt;

   h_trackVertexCuts_indexPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertexCuts_indexPt->GetXaxis()->SetRange(1, h_trackVertexCuts_indexPt->GetNbinsX() + 2);
   h_trackVertexCuts_indexPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertexCuts_indexPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertexCuts_indexPt->GetName() + ".jpeg");
   delete h_trackVertexCuts_indexPt;

   h_trueVertexCuts_x->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertexCuts_x->GetXaxis()->SetRange(1, h_trueVertexCuts_x->GetNbinsX() + 2);
   h_trueVertexCuts_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertexCuts_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertexCuts_x->GetName() + ".jpeg");
   delete h_trueVertexCuts_x;

   h_trueVertexCuts_y->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertexCuts_y->GetXaxis()->SetRange(1, h_trueVertexCuts_y->GetNbinsX() + 2);
   h_trueVertexCuts_y->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertexCuts_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertexCuts_y->GetName() + ".jpeg");
   delete h_trueVertexCuts_y;

   h_trueVertexCuts_z->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertexCuts_z->GetXaxis()->SetRange(1, h_trueVertexCuts_z->GetNbinsX() + 2);
   h_trueVertexCuts_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertexCuts_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertexCuts_z->GetName() + ".jpeg");
   delete h_trueVertexCuts_z;

   h_trueVertexCuts_sumPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertexCuts_sumPt->GetXaxis()->SetRange(1, h_trueVertexCuts_sumPt->GetNbinsX() + 2);
   h_trueVertexCuts_sumPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertexCuts_sumPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertexCuts_sumPt->GetName() + ".jpeg");
   delete h_trueVertexCuts_sumPt;

   h_trueVertexCuts_indexPt->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertexCuts_indexPt->GetXaxis()->SetRange(1, h_trueVertexCuts_indexPt->GetNbinsX() + 2);
   h_trueVertexCuts_indexPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertexCuts_indexPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertexCuts_indexPt->GetName() + ".jpeg");
   delete h_trueVertexCuts_indexPt;

   h_trk_chi2rphidof->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof->GetXaxis()->SetRange(1, h_trk_chi2rphidof->GetNbinsX() + 2);
   h_trk_chi2rphidof->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof->GetName() + ".jpeg");
   delete h_trk_chi2rphidof;

   h_trk_chi2rzdof->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof->GetXaxis()->SetRange(1, h_trk_chi2rzdof->GetNbinsX() + 2);
   h_trk_chi2rzdof->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof->GetName() + ".jpeg");
   delete h_trk_chi2rzdof;
   
   h_trk_bendchi2->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2->GetXaxis()->SetRange(1, h_trk_bendchi2->GetNbinsX() + 2);
   h_trk_bendchi2->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2->GetName() + ".jpeg");
   delete h_trk_bendchi2;

   h_trk_chi2rphidof_H->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_H->GetXaxis()->SetRange(1, h_trk_chi2rphidof_H->GetNbinsX() + 2);
   h_trk_chi2rphidof_H->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_H->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_H->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_H;

   h_trk_chi2rzdof_H->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_H->GetXaxis()->SetRange(1, h_trk_chi2rzdof_H->GetNbinsX() + 2);
   h_trk_chi2rzdof_H->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_H->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_H->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_H;

   h_trk_bendchi2_H->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_H->GetXaxis()->SetRange(1, h_trk_bendchi2_H->GetNbinsX() + 2);
   h_trk_bendchi2_H->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_H->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_H->GetName() + ".jpeg");
   delete h_trk_bendchi2_H;

   h_trk_chi2rphidof_L->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_L->GetXaxis()->SetRange(1, h_trk_chi2rphidof_L->GetNbinsX() + 2);
   h_trk_chi2rphidof_L->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_L->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_L->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_L;

   h_trk_chi2rzdof_L->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_L->GetXaxis()->SetRange(1, h_trk_chi2rzdof_L->GetNbinsX() + 2);
   h_trk_chi2rzdof_L->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_L->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_L->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_L;

   h_trk_bendchi2_L->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_L->GetXaxis()->SetRange(1, h_trk_bendchi2_L->GetNbinsX() + 2);
   h_trk_bendchi2_L->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_L->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_L->GetName() + ".jpeg");
   delete h_trk_bendchi2_L;

   h_trk_chi2rphidof_C->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_C->GetXaxis()->SetRange(1, h_trk_chi2rphidof_C->GetNbinsX() + 2);
   h_trk_chi2rphidof_C->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_C->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_C->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_C;

   h_trk_chi2rzdof_C->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_C->GetXaxis()->SetRange(1, h_trk_chi2rzdof_C->GetNbinsX() + 2);
   h_trk_chi2rzdof_C->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_C->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_C->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_C;

   h_trk_bendchi2_C->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_C->GetXaxis()->SetRange(1, h_trk_bendchi2_C->GetNbinsX() + 2);
   h_trk_bendchi2_C->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_C->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_C->GetName() + ".jpeg");
   delete h_trk_bendchi2_C;

   h_trk_chi2rphidof_I->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_I->GetXaxis()->SetRange(1, h_trk_chi2rphidof_I->GetNbinsX() + 2);
   h_trk_chi2rphidof_I->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_I->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_I->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_I;

   h_trk_chi2rzdof_I->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_I->GetXaxis()->SetRange(1, h_trk_chi2rzdof_I->GetNbinsX() + 2);
   h_trk_chi2rzdof_I->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_I->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_I->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_I;

   h_trk_bendchi2_I->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_I->GetXaxis()->SetRange(1, h_trk_bendchi2_I->GetNbinsX() + 2);
   h_trk_bendchi2_I->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_I->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_I->GetName() + ".jpeg");
   delete h_trk_bendchi2_I;

   h_trk_chi2rphidof_F->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_F->GetXaxis()->SetRange(1, h_trk_chi2rphidof_F->GetNbinsX() + 2);
   h_trk_chi2rphidof_F->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_F->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_F->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_F;

   h_trk_chi2rzdof_F->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_F->GetXaxis()->SetRange(1, h_trk_chi2rzdof_F->GetNbinsX() + 2);
   h_trk_chi2rzdof_F->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_F->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_F->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_F;

   h_trk_bendchi2_F->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_F->GetXaxis()->SetRange(1, h_trk_bendchi2_F->GetNbinsX() + 2);
   h_trk_bendchi2_F->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_F->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_F->GetName() + ".jpeg");
   delete h_trk_bendchi2_F;

   h_trk_chi2rphidof_P->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_P->GetXaxis()->SetRange(1, h_trk_chi2rphidof_P->GetNbinsX() + 2);
   h_trk_chi2rphidof_P->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_P->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_P->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_P;

   h_trk_chi2rzdof_P->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_P->GetXaxis()->SetRange(1, h_trk_chi2rzdof_P->GetNbinsX() + 2);
   h_trk_chi2rzdof_P->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_P->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_P->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_P;

   h_trk_bendchi2_P->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_P->GetXaxis()->SetRange(1, h_trk_bendchi2_P->GetNbinsX() + 2);
   h_trk_bendchi2_P->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_P->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_P->GetName() + ".jpeg");
   delete h_trk_bendchi2_P;

   h_trk_chi2rphidof_D->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_D->GetXaxis()->SetRange(1, h_trk_chi2rphidof_D->GetNbinsX() + 2);
   h_trk_chi2rphidof_D->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_D->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_D->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_D;

   h_trk_chi2rzdof_D->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_D->GetXaxis()->SetRange(1, h_trk_chi2rzdof_D->GetNbinsX() + 2);
   h_trk_chi2rzdof_D->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_D->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_D->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_D;

   h_trk_bendchi2_D->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_D->GetXaxis()->SetRange(1, h_trk_bendchi2_D->GetNbinsX() + 2);
   h_trk_bendchi2_D->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_D->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_D->GetName() + ".jpeg");
   delete h_trk_bendchi2_D;

   h_match_tp_pt_noCuts->Sumw2();
   h_tp_pt_noCuts->Sumw2();
   TH1F* h_eff_pt_noCuts = (TH1F*)h_match_tp_pt_noCuts->Clone();
   h_eff_pt_noCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_noCuts->GetXaxis()->SetRange(1, h_eff_pt_noCuts->GetNbinsX() + 2);
   h_eff_pt_noCuts->SetName("eff_pt_noCuts");
   h_eff_pt_noCuts->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_noCuts->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_noCuts->Divide(h_match_tp_pt_noCuts, h_tp_pt_noCuts, 1.0, 1.0, "B");
   h_eff_pt_noCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_noCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_noCuts->GetName() + ".jpeg");
   delete h_eff_pt_noCuts;
   delete h_match_tp_pt_noCuts;
   delete h_tp_pt_noCuts;

   h_match_tp_eta_noCuts->Sumw2();
   h_tp_eta_noCuts->Sumw2();
   TH1F* h_eff_eta_noCuts = (TH1F*)h_match_tp_eta_noCuts->Clone();
   h_eff_eta_noCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_noCuts->GetXaxis()->SetRange(1, h_eff_eta_noCuts->GetNbinsX() + 2);
   h_eff_eta_noCuts->SetName("eff_eta_noCuts");
   h_eff_eta_noCuts->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_noCuts->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_noCuts->Divide(h_match_tp_eta_noCuts, h_tp_eta_noCuts, 1.0, 1.0, "B");
   h_eff_eta_noCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_noCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_noCuts->GetName() + ".jpeg");
   delete h_eff_eta_noCuts;
   delete h_match_tp_eta_noCuts;
   delete h_tp_eta_noCuts;

   h_match_tp_d0_noCuts->Sumw2();
   h_tp_d0_noCuts->Sumw2();
   TH1F* h_eff_d0_noCuts = (TH1F*)h_match_tp_d0_noCuts->Clone();
   h_eff_d0_noCuts->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_noCuts->GetXaxis()->SetRange(1, h_eff_d0_noCuts->GetNbinsX() + 2);
   h_eff_d0_noCuts->SetName("eff_d0_noCuts");
   h_eff_d0_noCuts->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_noCuts->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_noCuts->Divide(h_match_tp_d0_noCuts, h_tp_d0_noCuts, 1.0, 1.0, "B");
   h_eff_d0_noCuts->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_noCuts->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_noCuts->GetName() + ".jpeg");
   delete h_eff_d0_noCuts;
   delete h_match_tp_d0_noCuts;
   delete h_tp_d0_noCuts;

   h_match_tp_pt_maxD0Cut->Sumw2();
   h_tp_pt_maxD0Cut->Sumw2();
   TH1F* h_eff_pt_maxD0Cut = (TH1F*)h_match_tp_pt_maxD0Cut->Clone();
   h_eff_pt_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_maxD0Cut->GetXaxis()->SetRange(1, h_eff_pt_maxD0Cut->GetNbinsX() + 2);
   h_eff_pt_maxD0Cut->SetName("eff_pt_maxD0Cut");
   h_eff_pt_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_maxD0Cut->Divide(h_match_tp_pt_maxD0Cut, h_tp_pt_maxD0Cut, 1.0, 1.0, "B");
   h_eff_pt_maxD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_maxD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_maxD0Cut->GetName() + ".jpeg");
   delete h_eff_pt_maxD0Cut;
   delete h_match_tp_pt_maxD0Cut;
   delete h_tp_pt_maxD0Cut;

   h_match_tp_eta_maxD0Cut->Sumw2();
   h_tp_eta_maxD0Cut->Sumw2();
   TH1F* h_eff_eta_maxD0Cut = (TH1F*)h_match_tp_eta_maxD0Cut->Clone();
   h_eff_eta_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_maxD0Cut->GetXaxis()->SetRange(1, h_eff_eta_maxD0Cut->GetNbinsX() + 2);
   h_eff_eta_maxD0Cut->SetName("eff_eta_maxD0Cut");
   h_eff_eta_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_maxD0Cut->Divide(h_match_tp_eta_maxD0Cut, h_tp_eta_maxD0Cut, 1.0, 1.0, "B");
   h_eff_eta_maxD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_maxD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_maxD0Cut->GetName() + ".jpeg");
   delete h_eff_eta_maxD0Cut;
   delete h_match_tp_eta_maxD0Cut;
   delete h_tp_eta_maxD0Cut;

   h_match_tp_d0_maxD0Cut->Sumw2();
   h_tp_d0_maxD0Cut->Sumw2();
   TH1F* h_eff_d0_maxD0Cut = (TH1F*)h_match_tp_d0_maxD0Cut->Clone();
   h_eff_d0_maxD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_maxD0Cut->GetXaxis()->SetRange(1, h_eff_d0_maxD0Cut->GetNbinsX() + 2);
   h_eff_d0_maxD0Cut->SetName("eff_d0_maxD0Cut");
   h_eff_d0_maxD0Cut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_maxD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_maxD0Cut->Divide(h_match_tp_d0_maxD0Cut, h_tp_d0_maxD0Cut, 1.0, 1.0, "B");
   h_eff_d0_maxD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_maxD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_maxD0Cut->GetName() + ".jpeg");
   delete h_eff_d0_maxD0Cut;
   delete h_match_tp_d0_maxD0Cut;
   delete h_tp_d0_maxD0Cut;

   h_match_tp_pt_minD0Cut->Sumw2();
   h_tp_pt_minD0Cut->Sumw2();
   TH1F* h_eff_pt_minD0Cut = (TH1F*)h_match_tp_pt_minD0Cut->Clone();
   h_eff_pt_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_minD0Cut->GetXaxis()->SetRange(1, h_eff_pt_minD0Cut->GetNbinsX() + 2);
   h_eff_pt_minD0Cut->SetName("eff_pt_minD0Cut");
   h_eff_pt_minD0Cut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_minD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_minD0Cut->Divide(h_match_tp_pt_minD0Cut, h_tp_pt_minD0Cut, 1.0, 1.0, "B");
   h_eff_pt_minD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_minD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_minD0Cut->GetName() + ".jpeg");
   delete h_eff_pt_minD0Cut;
   delete h_match_tp_pt_minD0Cut;
   delete h_tp_pt_minD0Cut;

   h_match_tp_eta_minD0Cut->Sumw2();
   h_tp_eta_minD0Cut->Sumw2();
   TH1F* h_eff_eta_minD0Cut = (TH1F*)h_match_tp_eta_minD0Cut->Clone();
   h_eff_eta_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_minD0Cut->GetXaxis()->SetRange(1, h_eff_eta_minD0Cut->GetNbinsX() + 2);
   h_eff_eta_minD0Cut->SetName("eff_eta_minD0Cut");
   h_eff_eta_minD0Cut->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_minD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_minD0Cut->Divide(h_match_tp_eta_minD0Cut, h_tp_eta_minD0Cut, 1.0, 1.0, "B");
   h_eff_eta_minD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_minD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_minD0Cut->GetName() + ".jpeg");
   delete h_eff_eta_minD0Cut;
   delete h_match_tp_eta_minD0Cut;
   delete h_tp_eta_minD0Cut;

   h_match_tp_d0_minD0Cut->Sumw2();
   h_tp_d0_minD0Cut->Sumw2();
   TH1F* h_eff_d0_minD0Cut = (TH1F*)h_match_tp_d0_minD0Cut->Clone();
   h_eff_d0_minD0Cut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_minD0Cut->GetXaxis()->SetRange(1, h_eff_d0_minD0Cut->GetNbinsX() + 2);
   h_eff_d0_minD0Cut->SetName("eff_d0_minD0Cut");
   h_eff_d0_minD0Cut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_minD0Cut->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_minD0Cut->Divide(h_match_tp_d0_minD0Cut, h_tp_d0_minD0Cut, 1.0, 1.0, "B");
   h_eff_d0_minD0Cut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_minD0Cut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_minD0Cut->GetName() + ".jpeg");
   delete h_eff_d0_minD0Cut;
   delete h_match_tp_d0_minD0Cut;
   delete h_tp_d0_minD0Cut;

   h_match_tp_pt_minPtCut->Sumw2();
   h_tp_pt_minPtCut->Sumw2();
   TH1F* h_eff_pt_minPtCut = (TH1F*)h_match_tp_pt_minPtCut->Clone();
   h_eff_pt_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_minPtCut->GetXaxis()->SetRange(1, h_eff_pt_minPtCut->GetNbinsX() + 2);
   h_eff_pt_minPtCut->SetName("eff_pt_minPtCut");
   h_eff_pt_minPtCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_minPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_minPtCut->Divide(h_match_tp_pt_minPtCut, h_tp_pt_minPtCut, 1.0, 1.0, "B");
   h_eff_pt_minPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_minPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_minPtCut->GetName() + ".jpeg");
   delete h_eff_pt_minPtCut;
   delete h_match_tp_pt_minPtCut;
   delete h_tp_pt_minPtCut;

   h_match_tp_eta_minPtCut->Sumw2();
   h_tp_eta_minPtCut->Sumw2();
   TH1F* h_eff_eta_minPtCut = (TH1F*)h_match_tp_eta_minPtCut->Clone();
   h_eff_eta_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_minPtCut->GetXaxis()->SetRange(1, h_eff_eta_minPtCut->GetNbinsX() + 2);
   h_eff_eta_minPtCut->SetName("eff_eta_minPtCut");
   h_eff_eta_minPtCut->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_minPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_minPtCut->Divide(h_match_tp_eta_minPtCut, h_tp_eta_minPtCut, 1.0, 1.0, "B");
   h_eff_eta_minPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_minPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_minPtCut->GetName() + ".jpeg");
   delete h_eff_eta_minPtCut;
   delete h_match_tp_eta_minPtCut;
   delete h_tp_eta_minPtCut;

   h_match_tp_d0_minPtCut->Sumw2();
   h_tp_d0_minPtCut->Sumw2();
   TH1F* h_eff_d0_minPtCut = (TH1F*)h_match_tp_d0_minPtCut->Clone();
   h_eff_d0_minPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_minPtCut->GetXaxis()->SetRange(1, h_eff_d0_minPtCut->GetNbinsX() + 2);
   h_eff_d0_minPtCut->SetName("eff_d0_minPtCut");
   h_eff_d0_minPtCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_minPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_minPtCut->Divide(h_match_tp_d0_minPtCut, h_tp_d0_minPtCut, 1.0, 1.0, "B");
   h_eff_d0_minPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_minPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_minPtCut->GetName() + ".jpeg");
   delete h_eff_d0_minPtCut;
   delete h_match_tp_d0_minPtCut;
   delete h_tp_d0_minPtCut;

   h_match_tp_pt_maxPtCut->Sumw2();
   h_tp_pt_maxPtCut->Sumw2();
   TH1F* h_eff_pt_maxPtCut = (TH1F*)h_match_tp_pt_maxPtCut->Clone();
   h_eff_pt_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_maxPtCut->GetXaxis()->SetRange(1, h_eff_pt_maxPtCut->GetNbinsX() + 2);
   h_eff_pt_maxPtCut->SetName("eff_pt_maxPtCut");
   h_eff_pt_maxPtCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_maxPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_maxPtCut->Divide(h_match_tp_pt_maxPtCut, h_tp_pt_maxPtCut, 1.0, 1.0, "B");
   h_eff_pt_maxPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_maxPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_maxPtCut->GetName() + ".jpeg");
   delete h_eff_pt_maxPtCut;
   delete h_match_tp_pt_maxPtCut;
   delete h_tp_pt_maxPtCut;

   h_match_tp_eta_maxPtCut->Sumw2();
   h_tp_eta_maxPtCut->Sumw2();
   TH1F* h_eff_eta_maxPtCut = (TH1F*)h_match_tp_eta_maxPtCut->Clone();
   h_eff_eta_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_maxPtCut->GetXaxis()->SetRange(1, h_eff_eta_maxPtCut->GetNbinsX() + 2);
   h_eff_eta_maxPtCut->SetName("eff_eta_maxPtCut");
   h_eff_eta_maxPtCut->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_maxPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_maxPtCut->Divide(h_match_tp_eta_maxPtCut, h_tp_eta_maxPtCut, 1.0, 1.0, "B");
   h_eff_eta_maxPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_maxPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_maxPtCut->GetName() + ".jpeg");
   delete h_eff_eta_maxPtCut;
   delete h_match_tp_eta_maxPtCut;
   delete h_tp_eta_maxPtCut;

   h_match_tp_d0_maxPtCut->Sumw2();
   h_tp_d0_maxPtCut->Sumw2();
   TH1F* h_eff_d0_maxPtCut = (TH1F*)h_match_tp_d0_maxPtCut->Clone();
   h_eff_d0_maxPtCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_maxPtCut->GetXaxis()->SetRange(1, h_eff_d0_maxPtCut->GetNbinsX() + 2);
   h_eff_d0_maxPtCut->SetName("eff_d0_maxPtCut");
   h_eff_d0_maxPtCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_maxPtCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_maxPtCut->Divide(h_match_tp_d0_maxPtCut, h_tp_d0_maxPtCut, 1.0, 1.0, "B");
   h_eff_d0_maxPtCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_maxPtCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_maxPtCut->GetName() + ".jpeg");
   delete h_eff_d0_maxPtCut;
   delete h_match_tp_d0_maxPtCut;
   delete h_tp_d0_maxPtCut;

   h_match_tp_pt_maxEtaCut->Sumw2();
   h_tp_pt_maxEtaCut->Sumw2();
   TH1F* h_eff_pt_maxEtaCut = (TH1F*)h_match_tp_pt_maxEtaCut->Clone();
   h_eff_pt_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt_maxEtaCut->GetXaxis()->SetRange(1, h_eff_pt_maxEtaCut->GetNbinsX() + 2);
   h_eff_pt_maxEtaCut->SetName("eff_pt_maxEtaCut");
   h_eff_pt_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt_maxEtaCut->Divide(h_match_tp_pt_maxEtaCut, h_tp_pt_maxEtaCut, 1.0, 1.0, "B");
   h_eff_pt_maxEtaCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt_maxEtaCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt_maxEtaCut->GetName() + ".jpeg");
   delete h_eff_pt_maxEtaCut;
   delete h_match_tp_pt_maxEtaCut;
   delete h_tp_pt_maxEtaCut;

   h_match_tp_eta_maxEtaCut->Sumw2();
   h_tp_eta_maxEtaCut->Sumw2();
   TH1F* h_eff_eta_maxEtaCut = (TH1F*)h_match_tp_eta_maxEtaCut->Clone();
   h_eff_eta_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta_maxEtaCut->GetXaxis()->SetRange(1, h_eff_eta_maxEtaCut->GetNbinsX() + 2);
   h_eff_eta_maxEtaCut->SetName("eff_eta_maxEtaCut");
   h_eff_eta_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta_maxEtaCut->Divide(h_match_tp_eta_maxEtaCut, h_tp_eta_maxEtaCut, 1.0, 1.0, "B");
   h_eff_eta_maxEtaCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta_maxEtaCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta_maxEtaCut->GetName() + ".jpeg");
   delete h_eff_eta_maxEtaCut;
   delete h_match_tp_eta_maxEtaCut;
   delete h_tp_eta_maxEtaCut;

   h_match_tp_d0_maxEtaCut->Sumw2();
   h_tp_d0_maxEtaCut->Sumw2();
   TH1F* h_eff_d0_maxEtaCut = (TH1F*)h_match_tp_d0_maxEtaCut->Clone();
   h_eff_d0_maxEtaCut->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0_maxEtaCut->GetXaxis()->SetRange(1, h_eff_d0_maxEtaCut->GetNbinsX() + 2);
   h_eff_d0_maxEtaCut->SetName("eff_d0_maxEtaCut");
   h_eff_d0_maxEtaCut->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0_maxEtaCut->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0_maxEtaCut->Divide(h_match_tp_d0_maxEtaCut, h_tp_d0_maxEtaCut, 1.0, 1.0, "B");
   h_eff_d0_maxEtaCut->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0_maxEtaCut->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0_maxEtaCut->GetName() + ".jpeg");
   delete h_eff_d0_maxEtaCut;
   delete h_match_tp_d0_maxEtaCut;
   delete h_tp_d0_maxEtaCut;

   h_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
   h_delta_dist_xy->GetXaxis()->SetRange(1, h_delta_dist_xy->GetNbinsX() + 2);
   h_delta_dist_xy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_delta_dist_xy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_delta_dist_xy->GetName() + ".jpeg");
   delete h_delta_dist_xy;
   
   h_error_delta_x->GetYaxis()->SetNoExponent(kTRUE);
   h_error_delta_x->GetXaxis()->SetRange(1, h_error_delta_x->GetNbinsX() + 2);
   h_error_delta_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_error_delta_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_error_delta_x->GetName() + ".jpeg");
   delete h_error_delta_x;
   
   h_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
   h_delta_dist_z->GetXaxis()->SetRange(1, h_delta_dist_z->GetNbinsX() + 2);
   h_delta_dist_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_delta_dist_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_delta_dist_z->GetName() + ".jpeg");
   delete h_delta_dist_z;

   h_error_delta_z->GetYaxis()->SetNoExponent(kTRUE);
   h_error_delta_z->GetXaxis()->SetRange(1, h_error_delta_z->GetNbinsX() + 2);
   h_error_delta_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_error_delta_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_error_delta_z->GetName() + ".jpeg");
   delete h_error_delta_z;
   
   h_delta_x->GetYaxis()->SetNoExponent(kTRUE);
   h_delta_x->GetXaxis()->SetRange(1, h_delta_x->GetNbinsX() + 2);
   h_delta_x->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_delta_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_delta_x->GetName() + ".jpeg");
   delete h_delta_x;

   char binlabel[1000];
   sprintf(binlabel, "0 + 1,   0 + 2,   1 + 2,   0 + 3,   1 + 3,   2 + 3");

   h_trk_Counter_TPcombination->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_Counter_TPcombination->GetXaxis()->SetRange(0, h_trk_Counter_TPcombination->GetNbinsX() + 1);
   h_trk_Counter_TPcombination->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   mySmallText(0.25, 0.95, 1, binlabel);
   h_trk_Counter_TPcombination->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_Counter_TPcombination->GetName() + ".jpeg");
   delete h_trk_Counter_TPcombination;

   h_trk_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_delta_dist_xy->GetXaxis()->SetRange(1, h_trk_delta_dist_xy->GetNbinsX() + 2);
   h_trk_delta_dist_xy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_delta_dist_xy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_delta_dist_xy->GetName() + ".jpeg");
   delete h_trk_delta_dist_xy;

   h_trackVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_cos_T->GetXaxis()->SetRange(1, h_trackVertex_cos_T->GetNbinsX() + 2);
   h_trackVertex_cos_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_cos_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_cos_T->GetName() + ".jpeg");
   delete h_trackVertex_cos_T;

   h_trackVertex_alpha_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_alpha_T->GetXaxis()->SetRange(1, h_trackVertex_alpha_T->GetNbinsX() + 2);
   h_trackVertex_alpha_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_alpha_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_alpha_T->GetName() + ".jpeg");
   delete h_trackVertex_alpha_T;

   h_trackVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_d_T->GetXaxis()->SetRange(1, h_trackVertex_d_T->GetNbinsX() + 2);
   h_trackVertex_d_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_d_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_d_T->GetName() + ".jpeg");
   delete h_trackVertex_d_T;

   h_trackVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trackVertex_R_T->GetXaxis()->SetRange(1, h_trackVertex_R_T->GetNbinsX() + 2);
   h_trackVertex_R_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trackVertex_R_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trackVertex_R_T->GetName() + ".jpeg");
   delete h_trackVertex_R_T;

   h_trueVertex_cos_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_cos_T->GetXaxis()->SetRange(1, h_trueVertex_cos_T->GetNbinsX() + 2);
   h_trueVertex_cos_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_cos_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_cos_T->GetName() + ".jpeg");
   delete h_trueVertex_cos_T;

   h_trueVertex_alpha_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_alpha_T->GetXaxis()->SetRange(1, h_trueVertex_alpha_T->GetNbinsX() + 2);
   h_trueVertex_alpha_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_alpha_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_alpha_T->GetName() + ".jpeg");
   delete h_trueVertex_alpha_T;

   h_trueVertex_d_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_d_T->GetXaxis()->SetRange(1, h_trueVertex_d_T->GetNbinsX() + 2);
   h_trueVertex_d_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_d_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_d_T->GetName() + ".jpeg");
   delete h_trueVertex_d_T;

   h_trueVertex_R_T->GetYaxis()->SetNoExponent(kTRUE);
   h_trueVertex_R_T->GetXaxis()->SetRange(1, h_trueVertex_R_T->GetNbinsX() + 2);
   h_trueVertex_R_T->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trueVertex_R_T->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trueVertex_R_T->GetName() + ".jpeg");
   delete h_trueVertex_R_T;

   h_trk_delta_dist_z->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_delta_dist_z->GetXaxis()->SetRange(1, h_trk_delta_dist_z->GetNbinsX() + 2);
   h_trk_delta_dist_z->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_delta_dist_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_delta_dist_z->GetName() + ".jpeg");
   delete h_trk_delta_dist_z;

   char res[1000];
   float rms = 0;
   TF1* fit;
   fit = new TF1("fit", "gaus", -1, 1);
   h_res_tp_trk_x->GetXaxis()->SetRange(1, h_res_tp_trk_x->GetNbinsX() + 2);
   h_res_tp_trk_x->Fit("fit");
   h_res_tp_trk_x->Draw();
   rms = fit->GetParameter(2);
   sprintf(res, "RMS = %.4f", rms);
   mySmallText(0.22, 0.82, 1, res);
   mySmallText(0.4, 0.42, 1, ctxt);
   h_res_tp_trk_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_res_tp_trk_x->GetName() + ".jpeg");
   delete h_res_tp_trk_x;
   delete fit;

   fit = new TF1("fit", "gaus", -1, 1);
   h_res_tp_trk_y->GetXaxis()->SetRange(1, h_res_tp_trk_y->GetNbinsX() + 2);
   h_res_tp_trk_y->Fit("fit");
   h_res_tp_trk_y->Draw();
   rms = fit->GetParameter(2);
   sprintf(res, "RMS = %.4f", rms);
   mySmallText(0.22, 0.82, 1, res);
   mySmallText(0.4, 0.42, 1, ctxt);
   h_res_tp_trk_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_res_tp_trk_y->GetName() + ".jpeg");
   delete h_res_tp_trk_y;
   delete fit;

   fit = new TF1("fit", "gaus", -1, 1);
   h_res_tp_trk_r->GetXaxis()->SetRange(1, h_res_tp_trk_r->GetNbinsX() + 2);
   h_res_tp_trk_r->Fit("fit");
   h_res_tp_trk_r->Draw();
   rms = fit->GetParameter(2);
   sprintf(res, "RMS = %.4f", rms);
   mySmallText(0.22, 0.82, 1, res);
   mySmallText(0.4, 0.42, 1, ctxt);
   h_res_tp_trk_r->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_res_tp_trk_r->GetName() + ".jpeg");
   delete h_res_tp_trk_r;
   delete fit;

   fit = new TF1("fit", "gaus", -1, 1);
   h_res_tp_trk_phi->GetXaxis()->SetRange(1, h_res_tp_trk_phi->GetNbinsX() + 2);
   h_res_tp_trk_phi->Fit("fit");
   h_res_tp_trk_phi->Draw();
   rms = fit->GetParameter(2);
   sprintf(res, "RMS = %.4f", rms);
   mySmallText(0.22, 0.82, 1, res);
   mySmallText(0.4, 0.42, 1, ctxt);
   h_res_tp_trk_phi->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_res_tp_trk_phi->GetName() + ".jpeg");
   delete h_res_tp_trk_phi;
   delete fit;

   fit = new TF1("fit", "gaus", -10, 10);
   h_res_tp_trk_z->GetXaxis()->SetRange(1, h_res_tp_trk_z->GetNbinsX() + 2);
   h_res_tp_trk_z->Fit("fit");
   h_res_tp_trk_z->Draw();
   rms = fit->GetParameter(2);
   sprintf(res, "RMS = %.4f", rms);
   mySmallText(0.22, 0.82, 1, res);
   mySmallText(0.4, 0.42, 1, ctxt);
   h_res_tp_trk_z->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_res_tp_trk_z->GetName() + ".jpeg");
   delete h_res_tp_trk_z;
   delete fit;

   h_all_trueVertex_pt->GetXaxis()->SetRange(1, h_all_trueVertex_pt->GetNbinsX() + 2);
   h_all_trueVertex_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trueVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trueVertex_pt->GetName() + ".jpeg");

   h_all_trueVertex_minD0->GetXaxis()->SetRange(1, h_all_trueVertex_minD0->GetNbinsX() + 2);
   h_all_trueVertex_minD0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trueVertex_minD0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trueVertex_minD0->GetName() + ".jpeg");

   h_all_trueVertex_lowPt->GetXaxis()->SetRange(1, h_all_trueVertex_lowPt->GetNbinsX() + 2);
   h_all_trueVertex_lowPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trueVertex_lowPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trueVertex_lowPt->GetName() + ".jpeg");

   h_all_trackVertex_pt->GetXaxis()->SetRange(1, h_all_trackVertex_pt->GetNbinsX() + 2);
   h_all_trackVertex_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trackVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trackVertex_pt->GetName() + ".jpeg");

   h_all_trackVertex_minD0->GetXaxis()->SetRange(1, h_all_trackVertex_minD0->GetNbinsX() + 2);
   h_all_trackVertex_minD0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trackVertex_minD0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trackVertex_minD0->GetName() + ".jpeg");

   h_all_trackVertex_lowPt->GetXaxis()->SetRange(1, h_all_trackVertex_lowPt->GetNbinsX() + 2);
   h_all_trackVertex_lowPt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trackVertex_lowPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trackVertex_lowPt->GetName() + ".jpeg");

   h_correct_trueVertex_pt->GetXaxis()->SetRange(1, h_correct_trueVertex_pt->GetNbinsX() + 2);
   h_correct_trueVertex_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trueVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trueVertex_pt->GetName() + ".jpeg");
   
   h_all_trueVertex_pt->Sumw2();
   h_correct_trueVertex_pt->Sumw2();
   TH1F* h_eff_trueVertex_pt = (TH1F*)h_correct_trueVertex_pt->Clone();
   h_eff_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trueVertex_pt->GetXaxis()->SetRange(1, h_eff_trueVertex_pt->GetNbinsX() + 2);
   h_eff_trueVertex_pt->SetName("eff_trueVertex_pt");
   h_eff_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
   h_eff_trueVertex_pt->Divide(h_correct_trueVertex_pt,h_all_trueVertex_pt, 1.0, 1.0, "B");
   h_eff_trueVertex_pt->Draw();
   h_eff_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trueVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trueVertex_pt->GetName() + ".jpeg");
   delete h_eff_trueVertex_pt;
   delete h_all_trueVertex_pt;
   delete h_correct_trueVertex_pt;

   h_all_trueVertex_lowPt->Sumw2();
   h_correct_trueVertex_lowPt->Sumw2();
   TH1F* h_eff_trueVertex_lowPt = (TH1F*)h_correct_trueVertex_lowPt->Clone();
   h_eff_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trueVertex_lowPt->GetXaxis()->SetRange(1, h_eff_trueVertex_lowPt->GetNbinsX() + 2);
   h_eff_trueVertex_lowPt->SetName("eff_trueVertex_lowPt");
   h_eff_trueVertex_lowPt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_trueVertex_lowPt->GetYaxis()->SetTitle("Efficiency");
   h_eff_trueVertex_lowPt->Divide(h_correct_trueVertex_lowPt,h_all_trueVertex_lowPt, 1.0, 1.0, "B");
   h_eff_trueVertex_lowPt->Draw();
   h_eff_trueVertex_lowPt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trueVertex_lowPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trueVertex_lowPt->GetName() + ".jpeg");
   delete h_eff_trueVertex_lowPt;
   delete h_all_trueVertex_lowPt;
   delete h_correct_trueVertex_lowPt;

   h_all_oneMatch_trueVertex_pt->Sumw2();
   h_correct_oneMatch_trueVertex_pt->Sumw2();
   TH1F* h_eff_oneMatch_trueVertex_pt = (TH1F*)h_correct_oneMatch_trueVertex_pt->Clone();
   h_eff_oneMatch_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatch_trueVertex_pt->GetXaxis()->SetRange(1, h_eff_oneMatch_trueVertex_pt->GetNbinsX() + 2);
   h_eff_oneMatch_trueVertex_pt->SetName("eff_oneMatch_trueVertex_pt");
   h_eff_oneMatch_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_oneMatch_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatch_trueVertex_pt->Divide(h_correct_oneMatch_trueVertex_pt,h_all_oneMatch_trueVertex_pt, 1.0, 1.0, "B");
   h_eff_oneMatch_trueVertex_pt->Draw();
   h_eff_oneMatch_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatch_trueVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_pt->GetName() + ".jpeg");
   delete h_eff_oneMatch_trueVertex_pt;
   delete h_all_oneMatch_trueVertex_pt;
   delete h_correct_oneMatch_trueVertex_pt;

   h_all_oneMatch_trueVertex_lowPt->Sumw2();
   h_correct_oneMatch_trueVertex_lowPt->Sumw2();
   TH1F* h_eff_oneMatch_trueVertex_lowPt = (TH1F*)h_correct_oneMatch_trueVertex_lowPt->Clone();
   h_eff_oneMatch_trueVertex_lowPt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatch_trueVertex_lowPt->GetXaxis()->SetRange(1, h_eff_oneMatch_trueVertex_lowPt->GetNbinsX() + 2);
   h_eff_oneMatch_trueVertex_lowPt->SetName("eff_oneMatch_trueVertex_lowPt");
   h_eff_oneMatch_trueVertex_lowPt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_oneMatch_trueVertex_lowPt->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatch_trueVertex_lowPt->Divide(h_correct_oneMatch_trueVertex_lowPt,h_all_oneMatch_trueVertex_lowPt, 1.0, 1.0, "B");
   h_eff_oneMatch_trueVertex_lowPt->Draw();
   h_eff_oneMatch_trueVertex_lowPt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatch_trueVertex_lowPt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_lowPt->GetName() + ".jpeg");
   delete h_eff_oneMatch_trueVertex_lowPt;
   delete h_all_oneMatch_trueVertex_lowPt;
   delete h_correct_oneMatch_trueVertex_lowPt;

   h_all_oneMatchAlt_trueVertex_pt->Sumw2();
   h_correct_oneMatchAlt_trueVertex_pt->Sumw2();
   TH1F* h_eff_oneMatchAlt_trueVertex_pt = (TH1F*)h_correct_oneMatchAlt_trueVertex_pt->Clone();
   h_eff_oneMatchAlt_trueVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatchAlt_trueVertex_pt->GetXaxis()->SetRange(1, h_eff_oneMatchAlt_trueVertex_pt->GetNbinsX() + 2);
   h_eff_oneMatchAlt_trueVertex_pt->SetName("eff_oneMatchAlt_trueVertex_pt");
   h_eff_oneMatchAlt_trueVertex_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_oneMatchAlt_trueVertex_pt->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatchAlt_trueVertex_pt->Divide(h_correct_oneMatchAlt_trueVertex_pt,h_all_oneMatchAlt_trueVertex_pt, 1.0, 1.0, "B");
   h_eff_oneMatchAlt_trueVertex_pt->Draw();
   h_eff_oneMatchAlt_trueVertex_pt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatchAlt_trueVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_pt->GetName() + ".jpeg");
   delete h_eff_oneMatchAlt_trueVertex_pt;
   delete h_all_oneMatchAlt_trueVertex_pt;
   delete h_correct_oneMatchAlt_trueVertex_pt;

   std::vector<double> eff_dxy_cuts;
   for(uint i=0;i<dxy_cuts.size();i++){
     eff_dxy_cuts.push_back(correct_vert_dxy_cut[i] / true_vertices);
   }
   TGraph* gr_eff_dxy_cuts = new TGraph(int(dxy_cuts.size()),dxy_cuts.data(),eff_dxy_cuts.data());
   gr_eff_dxy_cuts->SetTitle("Eff vs Dxy Cut; Dxy Cut (cm); Efficiency");
   gr_eff_dxy_cuts->SetName("gr_eff_dxy_cuts");
   gr_eff_dxy_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_dxy_cuts->GetName() + ".jpeg");
   delete gr_eff_dxy_cuts;

   std::vector<double> eff_dz_cuts;
   for(uint i=0;i<dz_cuts.size();i++){
     eff_dz_cuts.push_back(correct_vert_dz_cut[i] / true_vertices);
   }
   TGraph* gr_eff_dz_cuts = new TGraph(int(dz_cuts.size()),dz_cuts.data(),eff_dz_cuts.data());
   gr_eff_dz_cuts->SetTitle("Eff vs Dz Cut; Dz Cut (cm); Efficiency");
   gr_eff_dz_cuts->SetName("gr_eff_dz_cuts");
   gr_eff_dz_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_dz_cuts->GetName() + ".jpeg");
   delete gr_eff_dz_cuts;
   
   std::vector<double> eff_cos_T_cuts;
   for(uint i=0;i<cos_T_cuts.size();i++){
     eff_cos_T_cuts.push_back(correct_vert_cos_T_cut[i] / true_vertices);
   }
   TGraph* gr_eff_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),cos_T_cuts.data(),eff_cos_T_cuts.data());
   gr_eff_cos_T_cuts->SetTitle("Eff vs Cos(alpha) Cut; Cos(alpha) Cut; Efficiency");
   gr_eff_cos_T_cuts->SetName("gr_eff_cos_T_cuts");
   gr_eff_cos_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_cos_T_cuts->GetName() + ".jpeg");
   delete gr_eff_cos_T_cuts;

   std::vector<double> eff_d_T_cuts;
   for(uint i=0;i<d_T_cuts.size();i++){
     eff_d_T_cuts.push_back(correct_vert_d_T_cut[i] / true_vertices);
   }
   TGraph* gr_eff_d_T_cuts = new TGraph(int(d_T_cuts.size()),d_T_cuts.data(),eff_d_T_cuts.data());
   gr_eff_d_T_cuts->SetTitle("Eff vs d_T Cut; d_T Cut (cm); Efficiency");
   gr_eff_d_T_cuts->SetName("gr_eff_d_T_cuts");
   gr_eff_d_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_d_T_cuts->GetName() + ".jpeg");
   delete gr_eff_d_T_cuts;
   
   std::vector<double> eff_R_T_cuts;
   for(uint i=0;i<R_T_cuts.size();i++){
     eff_R_T_cuts.push_back(correct_vert_R_T_cut[i] / true_vertices);
   }
   TGraph* gr_eff_R_T_cuts = new TGraph(int(R_T_cuts.size()),R_T_cuts.data(),eff_R_T_cuts.data());
   gr_eff_R_T_cuts->SetTitle("Eff vs R_T Cut; R_T Cut (cm); Efficiency");
   gr_eff_R_T_cuts->SetName("gr_eff_R_T_cuts");
   gr_eff_R_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_R_T_cuts->GetName() + ".jpeg");
   delete gr_eff_R_T_cuts;

   std::vector<double> false_dxy_cuts;
   for(uint i=0;i<dxy_cuts.size();i++){
     if(all_vert_dxy_cut[i]!=0){
       false_dxy_cuts.push_back(false_vert_dxy_cut[i] / all_vert_dxy_cut[i]);
     }
     else{
       false_dxy_cuts.push_back(0);
     }
   }
   TGraph* gr_false_dxy_cuts = new TGraph(int(dxy_cuts.size()),dxy_cuts.data(),false_dxy_cuts.data());
   gr_false_dxy_cuts->SetTitle("False Rate vs Dxy Cut; Dxy Cut (cm); False Rate");
   gr_false_dxy_cuts->SetName("gr_false_dxy_cuts");
   gr_false_dxy_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_false_dxy_cuts->GetName() + ".jpeg");
   delete gr_false_dxy_cuts;

   std::vector<double> false_dz_cuts;
   for(uint i=0;i<dz_cuts.size();i++){
     if(all_vert_dz_cut[i]!=0){
       false_dz_cuts.push_back(false_vert_dz_cut[i] / all_vert_dz_cut[i]);
     }
     else{
       false_dz_cuts.push_back(0);
     }
   }
   TGraph* gr_false_dz_cuts = new TGraph(int(dz_cuts.size()),dz_cuts.data(),false_dz_cuts.data());
   gr_false_dz_cuts->SetTitle("False Rate vs Dz Cut; Dz Cut (cm); False Rate");
   gr_false_dz_cuts->SetName("gr_false_dz_cuts");
   gr_false_dz_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_false_dz_cuts->GetName() + ".jpeg");
   delete gr_false_dz_cuts;

   TGraph* gr_eff_vs_false_dz_cuts = new TGraph(int(dz_cuts.size()),eff_dz_cuts.data(),false_dz_cuts.data());
   gr_eff_vs_false_dz_cuts->SetTitle("False Rate vs Efficiency for Dz Cut; Efficiency; False Rate");
   gr_eff_vs_false_dz_cuts->SetName("gr_eff_vs_false_dz_cuts");
   gr_eff_vs_false_dz_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_vs_false_dz_cuts->GetName() + ".jpeg");
   delete gr_eff_vs_false_dz_cuts;
   
   std::vector<double> false_cos_T_cuts;
   for(uint i=0;i<cos_T_cuts.size();i++){
     if(all_vert_cos_T_cut[i]!=0){
       false_cos_T_cuts.push_back(false_vert_cos_T_cut[i] / all_vert_cos_T_cut[i]);
     }
     else{
       false_cos_T_cuts.push_back(0);
     }
   }
   TGraph* gr_false_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),cos_T_cuts.data(),false_cos_T_cuts.data());
   gr_false_cos_T_cuts->SetTitle("False Rate vs Cos(alpha) Cut; Cos(alpha) Cut; False Rate");
   gr_false_cos_T_cuts->SetName("gr_false_cos_T_cuts");
   gr_false_cos_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_false_cos_T_cuts->GetName() + ".jpeg");
   delete gr_false_cos_T_cuts;

   TGraph* gr_eff_vs_false_cos_T_cuts = new TGraph(int(cos_T_cuts.size()),eff_cos_T_cuts.data(),false_cos_T_cuts.data());
   gr_eff_vs_false_cos_T_cuts->SetTitle("False Rate vs Efficiency of Cos(alpha) Cut; Efficiency; False Rate");
   gr_eff_vs_false_cos_T_cuts->SetName("gr_eff_vs_false_cos_T_cuts");
   gr_eff_vs_false_cos_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_vs_false_cos_T_cuts->GetName() + ".jpeg");
   delete gr_eff_vs_false_cos_T_cuts;

   std::vector<double> false_d_T_cuts;
   for(uint i=0;i<d_T_cuts.size();i++){
     if(all_vert_d_T_cut[i]!=0){
       false_d_T_cuts.push_back(false_vert_d_T_cut[i] / all_vert_d_T_cut[i]);
     }
     else{
       false_d_T_cuts.push_back(0);
     }
   }
   TGraph* gr_false_d_T_cuts = new TGraph(int(d_T_cuts.size()),d_T_cuts.data(),false_d_T_cuts.data());
   gr_false_d_T_cuts->SetTitle("False Rate vs d_T Cut; d_T Cut (cm); False Rate");
   gr_false_d_T_cuts->SetName("gr_false_d_T_cuts");
   gr_false_d_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_false_d_T_cuts->GetName() + ".jpeg");
   delete gr_false_d_T_cuts;

   TGraph* gr_eff_vs_false_d_T_cuts = new TGraph(int(d_T_cuts.size()),eff_d_T_cuts.data(),false_d_T_cuts.data());
   gr_eff_vs_false_d_T_cuts->SetTitle("False Rate vs Efficiency of d_T Cut; Efficiency; False Rate");
   gr_eff_vs_false_d_T_cuts->SetName("gr_eff_vs_false_d_T_cuts");
   gr_eff_vs_false_d_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_vs_false_d_T_cuts->GetName() + ".jpeg");
   delete gr_eff_vs_false_d_T_cuts;

   std::vector<double> false_R_T_cuts;
   for(uint i=0;i<R_T_cuts.size();i++){
     if(all_vert_R_T_cut[i]!=0){
       false_R_T_cuts.push_back(false_vert_R_T_cut[i] / all_vert_R_T_cut[i]);
     }
     else{
       false_R_T_cuts.push_back(0);
     }
   }
   TGraph* gr_false_R_T_cuts = new TGraph(int(R_T_cuts.size()),R_T_cuts.data(),false_R_T_cuts.data());
   gr_false_R_T_cuts->SetTitle("False Rate vs R_T Cut; R_T Cut (cm); False Rate");
   gr_false_R_T_cuts->SetName("gr_false_R_T_cuts");
   gr_false_R_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_false_R_T_cuts->GetName() + ".jpeg");
   delete gr_false_R_T_cuts;

   TGraph* gr_eff_vs_false_R_T_cuts = new TGraph(int(R_T_cuts.size()),eff_R_T_cuts.data(),false_R_T_cuts.data());
   gr_eff_vs_false_R_T_cuts->SetTitle("False Rate vs Efficiency of R_T Cut; Efficiency; False Rate");
   gr_eff_vs_false_R_T_cuts->SetName("gr_eff_vs_false_R_T_cuts");
   gr_eff_vs_false_R_T_cuts->Draw("AC*");
   c.SaveAs(DIR + "/"+ gr_eff_vs_false_R_T_cuts->GetName() + ".jpeg");
   delete gr_eff_vs_false_R_T_cuts;

   h_all_trackVertex_pt->Sumw2();
   h_false_trackVertex_pt->Sumw2();
   TH1F* h_fake_trackVertex_pt = (TH1F*)h_false_trackVertex_pt->Clone();
   h_fake_trackVertex_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_fake_trackVertex_pt->GetXaxis()->SetRange(1, h_fake_trackVertex_pt->GetNbinsX() + 2);
   h_fake_trackVertex_pt->SetName("fake_trackVertex_pt");
   h_fake_trackVertex_pt->GetXaxis()->SetTitle("Track p_{T} (GeV)");
   h_fake_trackVertex_pt->GetYaxis()->SetTitle("Fake Rate");
   h_fake_trackVertex_pt->Divide(h_false_trackVertex_pt,h_all_trackVertex_pt, 1.0, 1.0, "B");
   h_fake_trackVertex_pt->Draw();
   h_fake_trackVertex_pt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_fake_trackVertex_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_fake_trackVertex_pt->GetName() + ".jpeg");
   delete h_fake_trackVertex_pt;
   delete h_all_trackVertex_pt;
   delete h_false_trackVertex_pt;

   h_all_trueVertex_eta->GetXaxis()->SetRange(1, h_all_trueVertex_eta->GetNbinsX() + 2);
   h_all_trueVertex_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trueVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trueVertex_eta->GetName() + ".jpeg");

   h_correct_trueVertex_eta->GetXaxis()->SetRange(1, h_correct_trueVertex_eta->GetNbinsX() + 2);
   h_correct_trueVertex_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trueVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trueVertex_eta->GetName() + ".jpeg");

   h_all_trueVertex_eta->Sumw2();
   h_correct_trueVertex_eta->Sumw2();
   TH1F* h_eff_trueVertex_eta = (TH1F*)h_correct_trueVertex_eta->Clone();
   h_eff_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trueVertex_eta->GetXaxis()->SetRange(1, h_eff_trueVertex_eta->GetNbinsX() + 2);
   h_eff_trueVertex_eta->SetName("eff_trueVertex_eta");
   h_eff_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
   h_eff_trueVertex_eta->Divide(h_correct_trueVertex_eta,h_all_trueVertex_eta, 1.0, 1.0, "B");
   h_eff_trueVertex_eta->Draw();
   h_eff_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trueVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trueVertex_eta->GetName() + ".jpeg");
   delete h_eff_trueVertex_eta;
   delete h_all_trueVertex_eta;
   delete h_correct_trueVertex_eta;

   h_all_oneMatch_trueVertex_eta->Sumw2();
   h_correct_oneMatch_trueVertex_eta->Sumw2();
   TH1F* h_eff_oneMatch_trueVertex_eta = (TH1F*)h_correct_oneMatch_trueVertex_eta->Clone();
   h_eff_oneMatch_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatch_trueVertex_eta->GetXaxis()->SetRange(1, h_eff_oneMatch_trueVertex_eta->GetNbinsX() + 2);
   h_eff_oneMatch_trueVertex_eta->SetName("eff_oneMatch_trueVertex_eta");
   h_eff_oneMatch_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_oneMatch_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatch_trueVertex_eta->Divide(h_correct_oneMatch_trueVertex_eta,h_all_oneMatch_trueVertex_eta, 1.0, 1.0, "B");
   h_eff_oneMatch_trueVertex_eta->Draw();
   h_eff_oneMatch_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatch_trueVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_eta->GetName() + ".jpeg");
   delete h_eff_oneMatch_trueVertex_eta;
   delete h_all_oneMatch_trueVertex_eta;
   delete h_correct_oneMatch_trueVertex_eta;

   h_all_oneMatchAlt_trueVertex_eta->Sumw2();
   h_correct_oneMatchAlt_trueVertex_eta->Sumw2();
   TH1F* h_eff_oneMatchAlt_trueVertex_eta = (TH1F*)h_correct_oneMatchAlt_trueVertex_eta->Clone();
   h_eff_oneMatchAlt_trueVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatchAlt_trueVertex_eta->GetXaxis()->SetRange(1, h_eff_oneMatchAlt_trueVertex_eta->GetNbinsX() + 2);
   h_eff_oneMatchAlt_trueVertex_eta->SetName("eff_oneMatchAlt_trueVertex_eta");
   h_eff_oneMatchAlt_trueVertex_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_oneMatchAlt_trueVertex_eta->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatchAlt_trueVertex_eta->Divide(h_correct_oneMatchAlt_trueVertex_eta,h_all_oneMatchAlt_trueVertex_eta, 1.0, 1.0, "B");
   h_eff_oneMatchAlt_trueVertex_eta->Draw();
   h_eff_oneMatchAlt_trueVertex_eta->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatchAlt_trueVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_eta->GetName() + ".jpeg");
   delete h_eff_oneMatchAlt_trueVertex_eta;
   delete h_all_oneMatchAlt_trueVertex_eta;
   delete h_correct_oneMatchAlt_trueVertex_eta;

   h_all_trackVertex_eta->Sumw2();
   h_false_trackVertex_eta->Sumw2();
   TH1F* h_fake_trackVertex_eta = (TH1F*)h_false_trackVertex_eta->Clone();
   h_fake_trackVertex_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_fake_trackVertex_eta->GetXaxis()->SetRange(1, h_fake_trackVertex_eta->GetNbinsX() + 2);
   h_fake_trackVertex_eta->SetName("fake_trackVertex_eta");
   h_fake_trackVertex_eta->GetXaxis()->SetTitle("Track #eta");
   h_fake_trackVertex_eta->GetYaxis()->SetTitle("Fake Rate");
   h_fake_trackVertex_eta->Divide(h_false_trackVertex_eta,h_all_trackVertex_eta, 1.0, 1.0, "B");
   h_fake_trackVertex_eta->Draw();
   h_fake_trackVertex_eta->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_fake_trackVertex_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_fake_trackVertex_eta->GetName() + ".jpeg");
   delete h_fake_trackVertex_eta;
   delete h_all_trackVertex_eta;
   delete h_false_trackVertex_eta;

   h_all_trueVertex_dxy->GetXaxis()->SetRange(1, h_all_trueVertex_dxy->GetNbinsX() + 2);
   h_all_trueVertex_dxy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trueVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trueVertex_dxy->GetName() + ".jpeg");

   h_correct_trueVertex_dxy->GetXaxis()->SetRange(1, h_correct_trueVertex_dxy->GetNbinsX() + 2);
   h_correct_trueVertex_dxy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trueVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trueVertex_dxy->GetName() + ".jpeg");
   
   h_all_trueVertex_dxy->Sumw2();
   h_correct_trueVertex_dxy->Sumw2();
   TH1F* h_eff_trueVertex_dxy = (TH1F*)h_correct_trueVertex_dxy->Clone();
   h_eff_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trueVertex_dxy->GetXaxis()->SetRange(1, h_eff_trueVertex_dxy->GetNbinsX() + 2);
   h_eff_trueVertex_dxy->SetName("eff_trueVertex_dxy");
   h_eff_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
   h_eff_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
   h_eff_trueVertex_dxy->Divide(h_correct_trueVertex_dxy,h_all_trueVertex_dxy, 1.0, 1.0, "B");
   h_eff_trueVertex_dxy->Draw();
   h_eff_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trueVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trueVertex_dxy->GetName() + ".jpeg");
   delete h_eff_trueVertex_dxy;
   delete h_all_trueVertex_dxy;
   delete h_correct_trueVertex_dxy;

   h_all_oneMatch_trueVertex_dxy->Sumw2();
   h_correct_oneMatch_trueVertex_dxy->Sumw2();
   TH1F* h_eff_oneMatch_trueVertex_dxy = (TH1F*)h_correct_oneMatch_trueVertex_dxy->Clone();
   h_eff_oneMatch_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatch_trueVertex_dxy->GetXaxis()->SetRange(1, h_eff_oneMatch_trueVertex_dxy->GetNbinsX() + 2);
   h_eff_oneMatch_trueVertex_dxy->SetName("eff_oneMatch_trueVertex_dxy");
   h_eff_oneMatch_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
   h_eff_oneMatch_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatch_trueVertex_dxy->Divide(h_correct_oneMatch_trueVertex_dxy,h_all_oneMatch_trueVertex_dxy, 1.0, 1.0, "B");
   h_eff_oneMatch_trueVertex_dxy->Draw();
   h_eff_oneMatch_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatch_trueVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatch_trueVertex_dxy->GetName() + ".jpeg");
   delete h_eff_oneMatch_trueVertex_dxy;
   delete h_all_oneMatch_trueVertex_dxy;
   delete h_correct_oneMatch_trueVertex_dxy;

   h_all_oneMatchAlt_trueVertex_dxy->Sumw2();
   h_correct_oneMatchAlt_trueVertex_dxy->Sumw2();
   TH1F* h_eff_oneMatchAlt_trueVertex_dxy = (TH1F*)h_correct_oneMatchAlt_trueVertex_dxy->Clone();
   h_eff_oneMatchAlt_trueVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_oneMatchAlt_trueVertex_dxy->GetXaxis()->SetRange(1, h_eff_oneMatchAlt_trueVertex_dxy->GetNbinsX() + 2);
   h_eff_oneMatchAlt_trueVertex_dxy->SetName("eff_oneMatchAlt_trueVertex_dxy");
   h_eff_oneMatchAlt_trueVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
   h_eff_oneMatchAlt_trueVertex_dxy->GetYaxis()->SetTitle("Efficiency");
   h_eff_oneMatchAlt_trueVertex_dxy->Divide(h_correct_oneMatchAlt_trueVertex_dxy,h_all_oneMatchAlt_trueVertex_dxy, 1.0, 1.0, "B");
   h_eff_oneMatchAlt_trueVertex_dxy->Draw();
   h_eff_oneMatchAlt_trueVertex_dxy->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_oneMatchAlt_trueVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_oneMatchAlt_trueVertex_dxy->GetName() + ".jpeg");
   delete h_eff_oneMatchAlt_trueVertex_dxy;
   delete h_all_oneMatchAlt_trueVertex_dxy;
   delete h_correct_oneMatchAlt_trueVertex_dxy;
   
   h_all_trackVertex_dxy->Sumw2();
   h_false_trackVertex_dxy->Sumw2();
   TH1F* h_fake_trackVertex_dxy = (TH1F*)h_false_trackVertex_dxy->Clone();
   h_fake_trackVertex_dxy->GetYaxis()->SetNoExponent(kTRUE);
   h_fake_trackVertex_dxy->GetXaxis()->SetRange(1, h_fake_trackVertex_dxy->GetNbinsX() + 2);
   h_fake_trackVertex_dxy->SetName("fake_trackVertex_dxy");
   h_fake_trackVertex_dxy->GetXaxis()->SetTitle("Vertex d_{xy} (cm)");
   h_fake_trackVertex_dxy->GetYaxis()->SetTitle("Fake Rate");
   h_fake_trackVertex_dxy->Divide(h_false_trackVertex_dxy,h_all_trackVertex_dxy, 1.0, 1.0, "B");
   h_fake_trackVertex_dxy->Draw();
   h_fake_trackVertex_dxy->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_fake_trackVertex_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_fake_trackVertex_dxy->GetName() + ".jpeg");
   delete h_fake_trackVertex_dxy;
   delete h_all_trackVertex_dxy;
   delete h_false_trackVertex_dxy;

   std::stringstream txt;
   float rate;
   if(type.Contains("NeutrinoGun")){
      rate = 40000.0;
      txt << " kHz ";
   }
   else{
      rate = 100.0;
      txt << " % ";
   }

   for(int k=0;k<5;k++){
      for (int l=0; l < 5; l++){
         h_Count_trk_pt_d0->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0->GetBinContent(l+1,k+1) * rate/nevt));
         h_Count_trk_pt_d0_dv->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0_dv->GetBinContent(l+1,k+1) * rate/nevt));
      }
   }
   
   h_Count_trk_pt_d0->Draw("colz");
   h_Count_trk_pt_d0->SetMarkerSize(2);
   h_Count_trk_pt_d0->Draw("textsame");
   h_Count_trk_pt_d0->GetXaxis()->SetNdivisions(5);
   h_Count_trk_pt_d0->GetYaxis()->SetNdivisions(5);
   h_Count_trk_pt_d0->GetXaxis()->CenterLabels();
   h_Count_trk_pt_d0->GetYaxis()->CenterLabels();
   h_Count_trk_pt_d0->GetYaxis()->SetTitle("Transverse Momentum p_T (GeV)");
   h_Count_trk_pt_d0->GetXaxis()->SetTitle("Transverse Impact Parameter d_0 (cm)");
   h_Count_trk_pt_d0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_Count_trk_pt_d0->GetName() + ".pdf");
   delete h_Count_trk_pt_d0;

   h_Count_trk_pt_d0_dv->Draw("colz");
   h_Count_trk_pt_d0_dv->SetMarkerSize(2);
   h_Count_trk_pt_d0_dv->Draw("textsame");
   h_Count_trk_pt_d0_dv->GetXaxis()->SetNdivisions(5);
   h_Count_trk_pt_d0_dv->GetYaxis()->SetNdivisions(5);
   h_Count_trk_pt_d0_dv->GetXaxis()->CenterLabels();
   h_Count_trk_pt_d0_dv->GetYaxis()->CenterLabels();
   h_Count_trk_pt_d0_dv->GetYaxis()->SetTitle("Transverse Momentum p_T (GeV)");
   h_Count_trk_pt_d0_dv->GetXaxis()->SetTitle("(DV selection) Transverse Impact Parameter d_0 (cm)");
   h_Count_trk_pt_d0_dv->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_Count_trk_pt_d0_dv->GetName() + ".pdf");
   delete h_Count_trk_pt_d0_dv;

   fout->Close();
}


Double_t dist_TPs(Track_Parameters a, Track_Parameters b){  
   float x1 = a.x0; //   Centers of the circles
   float y1 = a.y0; // 
   float x2 = b.x0; // 
   float y2 = b.y0; // 
   float R1 = a.rho;   // Radii of the circles
   float R2 = b.rho;
   float R = dist(x1,y1,x2,y2); // Distance between centers
   if((R>=(R1-R2)) && (R<=(R1+R2))){
      return (0);
   }
   else if(R==0){
      return (-99999.0);
   }
   else{

      return(R-R1-R2);
   }
}

Int_t Vertex(Track_Parameters a, Track_Parameters b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx){
   float x1 = a.x0; //   Centers of the circles
   float y1 = a.y0; // 
   float x2 = b.x0; // 
   float y2 = b.y0; // 
   float R1 = a.rho;   // Radii of the circles
   float R2 = b.rho;
   float del1, del2, x_11, x_12, x_21, x_22, y_11, y_12, y_21, y_22;
   float R = dist(x1,y1,x2,y2); // Distance between centers
   float centerdx = x1 - x2;
   float centerdy = y1 - y2;
   int retval = -1;

   if((R>=(R1-R2)) && (R<=(R1+R2))){
      // Circles Intersect
      float R4 = R*R*R*R;
      float A = (R1*R1 - R2*R2) / (2 * R*R);
      float r2r2 = (R1*R1 - R2*R2);
      float C = TMath::Sqrt(2 * (R1*R1 + R2*R2) / (R*R) - (r2r2 * r2r2) / R4 - 1);
      
      float fx = (x1+x2) / 2 + A * (x2 - x1);
      float gx = C * (y2 - y1) / 2;
      float ix1 = fx + gx;
      float ix2 = fx - gx;

      float fy = (y1+y2) / 2 + A * (y2 - y1);
      float gy = C * (x1 - x2) / 2;
      float iy1 = fy + gy;
      float iy2 = fy - gy;
      
      
      float ap1 = a.phi_T(ix1,iy1);
      float bp1 = b.phi_T(ix1,iy1);
      float ap2 = a.phi_T(ix2,iy2);
      float bp2 = b.phi_T(ix2,iy2);/*
      float z11 = a->z(ap1);
      float z21 = b->z(bp1);
      float z12 = a->z(ap2);
      float z22 = b->z(bp2);
      */
      float z11 = a.z(ix1,iy1);
      float z21 = b.z(ix1,iy1);
      float z12 = a.z(ix2,iy2);
      float z22 = b.z(ix2,iy2);
      
      float delz1 = fabs(z11 - z21);//fabs(fabs(ap1 - bp1)-TMath::Pi());// fabs(dist3(ix1,iy1,z11)-dist3(ix1,iy1,z21));
      float delz2 = fabs(z12 - z22);//fabs(fabs(ap2 - bp2)-TMath::Pi());// fabs(dist3(ix2,iy2,z12)-dist3(ix2,iy2,z22));
      
      if(VERBOSE[3]){// &&(fabs(z11-z_vtx)>0.1 && fabs(z12-z_vtx)>0.1 && fabs(z21-z_vtx)>0.1 && fabs(z22-z_vtx)>0.1)){
         std::cout<<Form("ix1 = %5.2f    |   iy1 = %5.2f    |  ix2 = %5.2f    |   iy2 = %5.2f",ix1,iy1,ix2,iy2)<<endl;
         // std::cout<<Form("ap1 = %5.2f    |   bp1 = %5.2f    |  ap2 = %5.2f    |   bp2 = %5.2f",ap1,bp1,ap2,bp2)<<endl;
         std::cout<<Form("z11 = %5.2f    |   z21 = %5.2f    |  z12 = %5.2f    |   z22 = %5.2f    |   delz1 = %5.2f    |   delz2 = %5.2f",z11,z21,z12,z22,delz1,delz2)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
      }
      retval = 0;
      if(gx==0 && gy==0){ // Only 1 intersection
         x_vtx = ix1;
         y_vtx = iy1;
         if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
            retval = 1;
         }
         z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
         // cout<<"----1 intersection ----";
         return retval;//true;
      }
      
      if(dist(ix1,iy1)>20){
         x_vtx = ix2;
         y_vtx = iy2;
         if(delz2<CUTOFF){//(fabs(z12-tp_z)<1.0 || fabs(z22-tp_z)<1.0){//
            retval = 2;
         }
         z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;//fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;// dist3(ix2,iy2,z12)<dist3(ix2,iy2,z22)?z12:z22;//fabs(a->z0-z12)<fabs(b->z0-z22)?z12:z22;//z12;
         return retval;//true;
      }
      if(dist(ix2,iy2)>20){
         x_vtx = ix1;
         y_vtx = iy1;
         if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
            retval = 2;
         }
         z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;//dist3(ix1,iy1,z11)<dist3(ix1,iy1,z21)?z11:z21;//fabs(a->z0-z11)<fabs(b->z0-z21)?z11:z21;//z11;
         return retval;//true;
      }
      
      // Finding minimum z distance to decide between the 2 intersection vertex

      // if((fabs(z11-tp_z)<1.0 || fabs(z12-tp_z)<1.0 || fabs(z21-tp_z)<1.0 || fabs(z22-tp_z)>1.0)){
          retval = 3;
      // }
 
      if (delz1<=delz2){
         x_vtx = ix1;
         y_vtx = iy1;
         z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;// fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
         // cout<<"----2 intersection ----";
         return retval;//true;         
      }
      else{
         x_vtx = ix2;
         y_vtx = iy2;
         z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;// fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;
         return retval;//true;
      }
      

   }
   else if(R==0){
      // Circles with Same center. 
      x_vtx = -999.0;//a->x();  // Considering First track coordinates as vertex
      y_vtx = -999.0;//a->y();
      return -2;
   }/*
   x_vtx = -99999.0;
   y_vtx = -99999.0;
   z_vtx = -99999.0;
   return false;
   else{
      x_vtx = -9999.0;
      y_vtx = -9999.0;
      return;
   }
   */
   //* Circles don't intersect. 

   if(x1==x2){
      if(y1==y2){
         // Circles with Same center. 
         //x_vtx = a->x();  // Considering First track coordinates as vertex
         // y_vtx = a->y();
         return -2;
      }
      x_11 = x1;
      x_12 = x1;

      x_21 = x2;
      x_22 = x2;

      y_11 = y1 + R1;
      y_12 = y1 - R1;

      y_21 = y2 + R2;
      y_22 = y2 - R2;

   }
   else{

      del1 = R1 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));
      del2 = R2 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));

      x_11 = x1 + del1;
      x_12 = x1 - del1;

      x_21 = x2 + del2;
      x_22 = x2 - del2;

      y_11 = y1 + (x_11-x1) * ((y1-y2)/(x1-x2));
      y_12 = y1 + (x_12-x1) * ((y1-y2)/(x1-x2));

      y_21 = y2 + (x_21-x2) * ((y1-y2)/(x1-x2));
      y_22 = y2 + (x_22-x2) * ((y1-y2)/(x1-x2));
   }

   if(dist(x_11,y_11,x2,y2) <= dist(x_12,y_12,x2,y2)){
      x_vtx = x_11;
      y_vtx = y_11;
   }
   else{
      x_vtx = x_12;
      y_vtx = y_12;
   }

   if(dist(x_21,y_21,x1,y1) <= dist(x_22,y_22,x1,y1)){
      x_vtx = (x_vtx + x_21)/2;
      y_vtx = (y_vtx + y_21)/2;
   }
   else{
      x_vtx = (x_vtx + x_22)/2;
      y_vtx = (y_vtx + y_22)/2;
   }
   //! Does this make sense for no intersection case???
   float az = a.z(x_vtx,y_vtx);
   float bz = b.z(x_vtx,y_vtx);
   z_vtx = (az+bz)/2.0;
   if(VERBOSE[3]){
      std::cout<<Form("z_vtx = %5.1f  |  az = %5.1f  |  bz = %5.1f",z_vtx,az,bz)<<endl;
   }
   return -1;
}

void SetPlotStyle()
{
   // from ATLAS plot style macro

   // use plain black on white colors
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameFillColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetStatColor(0);
   gStyle->SetHistLineColor(1);

   gStyle->SetPalette(1);

   // set the paper & margin sizes
   gStyle->SetPaperSize(20, 26);
   gStyle->SetPadTopMargin(0.05);
   gStyle->SetPadRightMargin(0.15);
   gStyle->SetPadBottomMargin(0.16);
   gStyle->SetPadLeftMargin(0.16);

   // set title offsets (for axis label)
   gStyle->SetTitleXOffset(1.4);
   gStyle->SetTitleYOffset(1.4);

   // use large fonts
   gStyle->SetTextFont(42);
   gStyle->SetTextSize(0.05);
   gStyle->SetLabelFont(42, "x");
   gStyle->SetTitleFont(42, "x");
   gStyle->SetLabelFont(42, "y");
   gStyle->SetTitleFont(42, "y");
   gStyle->SetLabelFont(42, "z");
   gStyle->SetTitleFont(42, "z");
   gStyle->SetLabelSize(0.05, "x");
   gStyle->SetTitleSize(0.05, "x");
   gStyle->SetLabelSize(0.05, "y");
   gStyle->SetTitleSize(0.05, "y");
   gStyle->SetLabelSize(0.05, "z");
   gStyle->SetTitleSize(0.05, "z");

   // use bold lines and markers
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerSize(1.2);
   gStyle->SetHistLineWidth(2.);
   gStyle->SetLineStyleString(2, "[12 12]");

   // get rid of error bar caps
   gStyle->SetEndErrorSize(0.);

   // do not display any of the standard histogram decorations
   gStyle->SetOptTitle(0);
   //gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);

   // put tick marks on top and RHS of plots
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
}

void mySmallText(Double_t x, Double_t y, Color_t color, char *text)
{
   Double_t tsize = 0.044;
   TLatex l;
   l.SetTextSize(tsize);
   l.SetNDC();
   l.SetTextColor(color);
   l.DrawLatex(x, y, text);
}
