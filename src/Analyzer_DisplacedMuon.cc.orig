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
#include "TMath.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x, Double_t y, Color_t color, char *text);

bool ComparePtTrack(Track_Parameters *a, Track_Parameters *b) { return a->pt > b->pt; }
bool CompareZ0Track(Track_Parameters *a, Track_Parameters *b) { return a->z0 > b->z0; }
bool CompareD0Track(Track_Parameters *a, Track_Parameters *b) { return a->d0 > b->d0; }

Double_t dist(Double_t x1, Double_t y1 , Double_t x2, Double_t y2){
   return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
}

Double_t dist_Vertex(Double_t x_vtx, Double_t y_vtx, Track_Parameters *a){
   float R = dist(x_vtx,y_vtx,a->x0,a->y0);
   return (fabs(R-(a->rho)));
}

void Vertex(Track_Parameters *a, Track_Parameters*b, Double_t &x_vtx, Double_t &y_vtx){
   float x1 = a->x0; //   Centers of the circles
   float y1 = a->y0; // 
   float x2 = b->x0; // 
   float y2 = b->y0; // 
   float R1 = a->rho;   // Radii of the circles
   float R2 = b->rho;
   float del1, del2, x_11, x_12, x_21, x_22, y_11, y_12, y_21, y_22;
   float R = dist(x1,y1,x2,y2); // Distance between centers
   float centerdx = x1 - x2;
   float centerdy = y1 - y2;

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

      if(gx==0 && gy==0){ // Only 1 intersection
         x_vtx = ix1;
         y_vtx = iy1;
         // cout<<"----1 intersection ----";
         return;
      }
      // Finding minimum z distance to decide between the 2 intersection vertex
      float z11 = a->phi_T(ix1,iy1);
      float z12 = a->phi_T(ix2,iy2);
      float z21 = b->phi_T(ix1,iy1);
      float z22 = b->phi_T(ix2,iy2);

      // cout<<Form("z11 = %5.2f    |   z12 = %5.2f    |  z21 = %5.2f    |   z22 = %5.2f    |   z22 = %5.2f",z11,z12,z21,z22)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;


      float delz1 = fabs(z11 - z21);
      float delz2 = fabs(z12 - z22);
      if (delz1 <= delz2){
         x_vtx = ix1;
         y_vtx = iy1;
         // cout<<"----2 intersection ----";
         return;         
      }
      else{
         x_vtx = ix2;
         y_vtx = iy2;
         return;
      }

   }
   else if(R==0){
      // Circles with Same center. 
      // x_vtx = -999.0;//a->x();  // Considering First track coordinates as vertex
      // y_vtx = -999.0;//a->y();
      return;
   }
   /*
   else{
      x_vtx = -9999.0;
      y_vtx = -9999.0;
      return;
   }
   */
   // Circles don't intersect. 

   if(x1==x2){
      if(y1==y2){
         // Circles with Same center. 
         //x_vtx = a->x();  // Considering First track coordinates as vertex
         // y_vtx = a->y();
         return;
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

void Analyzer_DisplacedMuon::Loop(TString type,
                            TString type_dir = "",
                            float TP_maxD0 = 10.0,
                            float TP_minD0 = 0.0,
                            int TP_select_injet = 0)
{

   gROOT->SetBatch();
   gErrorIgnoreLevel = kWarning;

   SetPlotStyle();

   float TP_minPt = 2.0;
   float TP_maxPt = 10000.0;
   float TP_maxEta = 2.4;
   int TP_select_pdgid = 0;

   // ----------------------------------------------------------------------------------------------------------------
   // histograms
   // ----------------------------------------------------------------------------------------------------------------

   std::vector<TString> regions{"tpgeq2", "tpgeq2osMu"};
   std::vector<TString> vars{"ntp","nMu",
                             "tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_d0", "tp_charge", "tp_pdgid","tp_dist",
                             "tp1pt", "tp1eta", "tp1phi", "tp1z0", "tp1d0", "tp1charge", "tp1pdgid","tp1_dist",
                             "tp2pt", "tp2eta", "tp2phi", "tp2z0", "tp2d0", "tp2charge", "tp2pdgid","tp2_dist",
                             "tp_tp_eta", "tp_tp_phi", "tp_tp_z0", "tp_tp_d0", "tp_tp_dR",
                             "tpz1pt", "tpz1eta", "tpz1phi", "tpz1z0", "tpz1d0", "tpz1charge", "tpz1pdgid","tpz1dist",
                             "tpz2pt", "tpz2eta", "tpz2phi", "tpz2z0", "tpz2d0", "tpz2charge", "tpz2pdgid","tpz2dist",
                             "tp_z0_min_delta_eta", "tp_z0_min_delta_phi", "tp_z0_min_delta_z0", "tp_z0_min_delta_d0", "tp_z0_min_delta_dR",
                             "tpd1pt", "tpd1eta", "tpd1phi", "tpd1z0", "tpd1d0", "tpd1charge", "tpd1pdgid","tpd1dist",
                             "tpd2pt", "tpd2eta", "tpd2phi", "tpd2z0", "tpd2d0", "tpd2charge", "tpd2pdgid","tpd2dist",
                             "tp_d0_min_delta_eta", "tp_d0_min_delta_phi", "tp_d0_min_delta_z0", "tp_d0_min_delta_d0", "tp_d0_min_delta_dR",
                             "muz1pt", "muz1eta", "muz1phi", "muz1z0", "muz1d0", "muz1charge", "muz1pdgid","muz1dist",
                             "muz2pt", "muz2eta", "muz2phi", "muz2z0", "muz2d0", "muz2charge", "muz2pdgid","muz2dist",
                             "mu_z0_min_delta_eta", "mu_z0_min_delta_phi", "mu_z0_min_delta_z0", "mu_z0_min_delta_d0", "mu_z0_min_delta_dR",
                             "mud1pt", "mud1eta", "mud1phi", "mud1z0", "mud1d0", "mud1charge", "mud1pdgid","mud1dist",
                             "mud2pt", "mud2eta", "mud2phi", "mud2z0", "mud2d0", "mud2charge", "mud2pdgid","mud2dist",
                             "mu_d0_min_delta_eta", "mu_d0_min_delta_phi", "mu_d0_min_delta_z0", "mu_d0_min_delta_d0", "mu_d0_min_delta_dR"};
   std::vector<int> nbins{100, 20,
                          20, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100, 100, 400, 100,
                          100, 100, 100, 100, 100};
   std::vector<float> lowEdge{0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, 0, 0, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, 0, 0, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, 0, 0, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, 0, 0, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, -2.4, -3.2, -20, -1, -2, 0, 0,
                              0, 0, 0, 0, 0};
   std::vector<float> highEdge{100, 20,
                               20,  2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               5, 4, 0.5, 0.5, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               5, 4, 0.5, 0.5, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               5, 4, 0.5, 0.5, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               5, 4, 0.5, 0.5, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               100, 2.4, 3.2, 20, 1, 2, 400, 5,
                               5, 4, 0.5, 0.5, 5};

   typedef vector<TH1F *> Dim1;
   typedef vector<Dim1> Dim2;
   typedef vector<Dim2> Dim3;
   typedef vector<Dim3> Dim4;
   TH2F *h_d0_Mu_n_Mu_p = new TH2F("d0_Mu_n_Mu_p","d0_Mu_n_Mu_p",50,-5,5,50,-5,5);
   TH2F *h_z0_Mu_n_Mu_p = new TH2F("z0_Mu_n_Mu_p","z0_Mu_n_Mu_p",100,-20,20,100,-20,20);
   TH2F *h_displaced_vertex_x = new TH2F("displaced_vertex_x","displaced_vertex_x",50,-5,5,50,-5,5);
   TH2F *h_displaced_vertex_y = new TH2F("displaced_vertex_y","displaced_vertex_y",50,-5,5,50,-5,5);
Dim2 Hists(regions.size(), Dim1(vars.size()));
   std::stringstream name;
   TH1F *h_test;
   Float_t **Output_vars = new Float_t *[regions.size()];
   for (int k = 0; k < regions.size(); ++k)
   {
      Output_vars[k] = new Float_t[vars.size()];
      for (int l = 0; l < vars.size(); ++l)
      {
         name << regions[k] << "_" << vars[l];
         h_test = new TH1F((name.str()).c_str(), (name.str()).c_str(), nbins[l], lowEdge[l], highEdge[l]);
         h_test->StatOverflows(kTRUE);
         h_test->Sumw2(kTRUE);
         Hists[k][l] = h_test;
         name.str("");
         Output_vars[k][l] = 0;
      }
   }

   if (fChain == 0) return;

   int nAccept = 0;
   std::vector<Track_Parameters *> *selectedTracks;
   std::vector<Track_Parameters *> *selectedTPs;
   std::vector<Track_Parameters *> *selectedTPs_zmin;
   std::vector<Track_Parameters *> *selectedTPs_dmin;
   std::vector<Track_Parameters *> *selectedMu;
   std::vector<Track_Parameters *> *selectedMu_zmin;
   std::vector<Track_Parameters *> *selectedMu_dmin;
   int TP_index = -1;
   int TP_mind0_nindex = -1;
   int TP_minz0_nindex = -1;
   int TP_mind0_pindex = -1;
   int TP_minz0_pindex = -1;

   int nCount_2trk = 0;
   float pt_cuts[5] = {5.0,7.0,9.0,11.0,13.0};
   float d0_cuts[5] = {0.0,0.2,0.5,0.8,1.0};
   // float z0_cuts[5] = {0.2,0.4,0.6,0.8,1.0};
   int nCount_trk_pt[5] = {0,0,0,0,0};
   int nCount_trk_d0[5] = {0,0,0,0,0};
   int nCount_trk_pt_d0[5][5] = {0};
   // int nCount_trk_z0[5] = {0,0,0,0,0};
   int nCount_2tp = 0;
   int nCount_tp_pt[5] = {0,0,0,0,0};
   int nCount_tp_d0[5] = {0,0,0,0,0};
   int nCount_tp_pt_d0[5][5] = {0};
   for (int k=0; k < 5; k++){
      for (int l=0; l < 5; l++){
         nCount_tp_pt_d0[k][l] = 0;
         nCount_trk_pt_d0[k][l] = 0;
      }
   }
   // int nCount_tp_z0[5] = {0,0,0,0,0};
   int nCount_2mu=0;
   int nCount_2muOS = 0;
   Long64_t nevt = fChain->GetEntries();

   for (Long64_t i=0; i<nevt; i++) { // event = 55461 or 41263 (cT100)  or 68163 (cT5000)
      fChain->GetEntry(i);
      // if (Cut(ientry) < 0) continue;

      displayProgress(i, nevt);

      // ----------------------------------------------------------------------------------------------------------------
      // track particle loop

      int ntrk = 0;
      int h_index;
      selectedTracks = new std::vector<Track_Parameters *>();
      // selectedTracks_zmin = new std::vector<Track_Parameters *>();
      // selectedTracks_dmin = new std::vector<Track_Parameters *>();

      selectedTPs = new std::vector<Track_Parameters *>();
      selectedTPs_zmin = new std::vector<Track_Parameters *>();
      selectedTPs_dmin = new std::vector<Track_Parameters *>();
      selectedMu = new std::vector<Track_Parameters *>();
      selectedMu_zmin = new std::vector<Track_Parameters *>();
      selectedMu_dmin = new std::vector<Track_Parameters *>();
      

      for (int it = 0; it < (int)trk_pt->size(); it++)
      {
         
         if (fabs(trk_eta->at(it)) > TP_maxEta)
            continue;
         if (trk_pt->at(it) < TP_minPt)
            continue;
         if (trk_pt->at(it) > TP_maxPt)
            continue;
         if (std::fabs(trk_d0->at(it)) > TP_maxD0)
            continue;
         if (std::fabs(trk_d0->at(it)) < TP_minD0)
            continue;

         // only look at tracks in (ttbar) jets ?
         if (TP_select_injet > 0)
         {
            if (TP_select_injet == 1 && trk_injet->at(it) == 0)
               continue;
            if (TP_select_injet == 2 && trk_injet_highpt->at(it) == 0)
               continue;
            if (TP_select_injet == 3 && trk_injet_vhighpt->at(it) == 0)
               continue;
         }
         int ndof = 2 * trk_nstub->at(it) - 5;
         float chi2rphidof = (float)trk_chi2rphi->at(it) / ndof;
         float chi2rzdof = (float)trk_chi2rz->at(it) / ndof;

         if(chi2rphidof > 10)
            continue;
         if(chi2rzdof   > 10)
            continue;
         if(trk_bendchi2 ->at(it) > 10)
            continue;

         float K = 0.01*0.5696*(-trk_rinv->at(it))/trk_pt->at(it);
         float rp = ((1/(2*K)) - trk_d0->at(it))/(-trk_rinv->at(it));
         float dxy = sqrt(rp*rp + (trk_d0->at(it)/K) - 1/(4*K*K));

         selectedTracks->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
<<<<<<< HEAD
         // (*selectedTracks)[selectedTracks->size()-1]->Propagate_Transverse(1.0/trk_rinv->at(it));
=======
>>>>>>> origin/main
         // selectedTPs->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
      }

      // ----------------------------------------------------------------------------------------------------------------
      // tracking particle loop
      TH1F htmp_tp_z("htmp_tp_z", ";z (cm); Tracks", 601 , -30, 30);
      TH1F htmp_tp_d("htmp_tp_d", ";d (cm); Tracks", 401 , -20, 20);

      TH1F htmp_tp_x("htmp_tp_x", ";x (cm); Tracks", 401 , -20, 20);
      TH1F htmp_tp_y("htmp_tp_y", ";y (cm); Tracks", 401 , -20, 20);

      // cout<<htmp_tp_d.GetNbinsX()<<endl;

      for (int it = 0; it < (int)tp_pt->size(); it++)
      {
         float tmp_d0 = tp_d0->at(it);
         float tmp_z0 = tp_z0->at(it);
         if (std::fabs(tmp_d0) > TP_maxD0)
            continue;
         if (std::fabs(tmp_d0) < TP_minD0)
            continue;
         if (tp_pt->at(it) < TP_minPt)
            continue;
         if (tp_pt->at(it) > TP_maxPt)
            continue;
         if (std::fabs(tp_eta->at(it)) > TP_maxEta)
            continue;
         if (TMath::IsNaN(tmp_d0))
            continue;
         if (TMath::IsNaN(tmp_z0))
            continue;

         if (TP_select_injet > 0)
         {
            if (TP_select_injet == 1 && tp_injet->at(it) == 0)
               continue;
            if (TP_select_injet == 2 && tp_injet_highpt->at(it) == 0)
               continue;
            if (TP_select_injet == 3 && tp_injet_vhighpt->at(it) == 0)
               continue;
         }
         // cut on PDG ID at plot stage?
         if (TP_select_pdgid != 0)
         {
            if (fabs(tp_pdgid->at(it)) != fabs(TP_select_pdgid))
               continue;
         }

         // if(fabs(tp_d0->at(it))<0.5)
         //    continue;

         float K = 0.01*0.5696*(tp_charge->at(it))/tp_pt->at(it);
         float rp = ((1/(2*K)) - tp_d0->at(it))/(tp_charge->at(it));
         float dxy = sqrt(rp*rp - (tp_d0_prod->at(it)/K) - 1/(4*K*K));
         // if(fabs(dxy-tp_dxy->at(it))>0.1)
            // cout<<i<<"\t dxy = "<<dxy<<"  \t tp_dxy = "<<tp_dxy->at(it)<<endl;
         // cout << i << "\t" << tp_z0->at(it) << "\t" << tp_d0->at(it) << "\t" << tp_pt->at(it) << "\t" << tp_pdgid->at(it)<<"\n";
<<<<<<< HEAD
         selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tmp_d0, tmp_z0, tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, tp_pdgid->at(it)));
         int j = selectedTPs->size()-1;
         // (*selectedTPs)[j]->x0 = tp_x->at(it);
         // (*selectedTPs)[j]->y0 = tp_y->at(it);
         // (*selectedTPs)[j]->Propagate_Transverse();
         // cout<<Form("phi_T =   %7.2f     |  ",
         // (*selectedTPs)[j]->Propagate_Transverse();
         // if(fabs(tp_pdgid->at(it))==13)
         // cout<<Form("dxy_calc = %5.2f    |   dxy = %5.2f    |  dist = %5.2f    |   eta = %5.2f    |   phi = %5.2f    |   pt = %5.2f",dxy,tp_dxy->at(it),TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y),tp_eta->at(it),tp_phi->at(it),tp_pt->at(it))<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
         htmp_tp_z.Fill(tmp_z0,tp_pt->at(it));
         htmp_tp_d.Fill(tmp_d0,tp_pt->at(it));
=======
         selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tp_dxy->at(it), tp_z0_prod->at(it), tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, tp_pdgid->at(it)));
>>>>>>> origin/main
      } // end tp loop

      // cout<<"End of tracking particle loop \n";

      // ----------------------------------------------------------------------------------------------------------------
      // Gen particle loop
      for (int it = 0; it < nPart; it++){
         if (partPt[it] < TP_minPt)
            continue;
         if (partPt[it] > TP_maxPt)
            continue;
         if (partStat[it] != 1)
            continue;
         if (abs(partId[it]) != 13)
            continue;

         // selectedTPs->push_back(new Track_Parameters(partPt[it], -99999, -99999, partEta[it], partPhi[it], partCh[it], it, partId[it]));
      }
      // cout << "End of Gen particle loop \n";

// --------------------------------------------------------------------------------------------
//                SELECTION CRITERIA -  Atleast 2 tracking particles 
// --------------------------------------------------------------------------------------------
      if (!(selectedTPs->size() >= 2)) continue;

   //* Primary Vertex Finding
      float zvtx_sliding = -999;
      float sigma_max = -999;
      int imax = -999;
      std::vector<int> found_z;
      std::vector<float> z_vertex;
      found_z.reserve(1);
      z_vertex.reserve(1);
      for (int ivtx = 0; ivtx < 1; ivtx++) {
         zvtx_sliding = -999;
         sigma_max = -999;
         imax = -999;
         for (int i = 2; i <= htmp_tp_z.GetNbinsX() - 1; i++) {
            float a0 = htmp_tp_z.GetBinContent(i - 1);
            float a1 = htmp_tp_z.GetBinContent(i);
            float a2 = htmp_tp_z.GetBinContent(i + 1);
            float sigma = a0 + a1 + a2;
            if ((sigma > sigma_max) && (find(found_z.begin(), found_z.end(), i) == found_z.end())) {
            sigma_max = sigma;
            imax = i;
            float z0 = htmp_tp_z.GetBinCenter(i - 1);
            float z1 = htmp_tp_z.GetBinCenter(i);
            float z2 = htmp_tp_z.GetBinCenter(i + 1);
            zvtx_sliding = (a0 * z0 + a1 * z1 + a2 * z2) / sigma;
            }
         }
         found_z.push_back(imax);
         z_vertex.push_back(zvtx_sliding);
      }
      // cout<<zvtx_sliding<<"\t";

      sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);
      Double_t x_dv = (tp_x->at((*selectedTPs)[0]->index));//+tp_x->at((*selectedTPs)[1]->index))/2.0;
      Double_t y_dv = (tp_y->at((*selectedTPs)[0]->index));//+tp_y->at((*selectedTPs)[1]->index))/2.0;
      Double_t z_dv = (tp_z->at((*selectedTPs)[0]->index));//+tp_z->at((*selectedTPs)[1]->index))/2.0;
      // float dist_12 = dist((*selectedTPs)[0]->x, (*selectedTPs)[0]->y,(*selectedTPs)[1]->x, (*selectedTPs)[1]->y);

      if(!(fabs(tp_x->at((*selectedTPs)[0]->index) - tp_x->at((*selectedTPs)[1]->index))<0.01)) continue;
      Vertex((*selectedTPs)[0],(*selectedTPs)[1],x_dv,y_dv); //(selectedTPs->size())>2?2:1
      // z_dv = ((*selectedTPs)[0]->z0 + (*selectedTPs)[1]->z0)/2.0;
      // (*selectedTPs)[0]->Propagate_Transverse();
      // (*selectedTPs)[1]->Propagate_Transverse();

      // if(dist(x_dv,y_dv,tp_x->at((*selectedTPs)[0]->index),tp_y->at((*selectedTPs)[0]->index))>0.5) continue;

      // if(fabs(tp_x->at((*selectedTPs)[0]->index) - x_dv)>0.1)
      //    cout<<Form("x_dv = %5.2f    |   y_dv = %5.2f    |  x[0] = %5.2f    |   y[0] = %5.2f    |   x[1] = %5.2f    |   y[1] = %5.2f    |   z[0] = -    |   z0[0] = -",x_dv,y_dv,tp_x->at((*selectedTPs)[0]->index),tp_y->at((*selectedTPs)[0]->index),tp_x->at((*selectedTPs)[1]->index),tp_y->at((*selectedTPs)[1]->index))<<endl;;//,(*selectedTPs)[0]->z((*selectedTPs)[0]->phi_T(tp_x->at((*selectedTPs)[0]->index))),(*selectedTPs)[0]->z0)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
      // cout <<x_dv<<"\t"<<y_dv<<"   Adding to vectors \n";
      for (int j = 0; j < selectedTPs->size(); j++){
         (*selectedTPs)[j]->dist_calc(x_dv,y_dv,tp_x->at((*selectedTPs)[j]->index),tp_y->at((*selectedTPs)[j]->index));
         // cout<<"Check 1\n";
         // if((fabs((*selectedTPs)[j]->z0 - z_dv)) < 1.0){
            selectedTPs_zmin->push_back(new Track_Parameters((*selectedTPs)[j]->pt, (*selectedTPs)[j]->d0, (*selectedTPs)[j]->z0, (*selectedTPs)[j]->eta, (*selectedTPs)[j]->phi, (*selectedTPs)[j]->charge, (*selectedTPs)[j]->index, (*selectedTPs)[j]->pdgid));
            (*selectedTPs_zmin)[selectedTPs_zmin->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedTPs_zmin)[selectedTPs_zmin->size()-1]->index),tp_y->at((*selectedTPs_zmin)[selectedTPs_zmin->size()-1]->index));
            (*selectedTPs_zmin)[selectedTPs_zmin->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedTPs_zmin)[selectedTPs_zmin->size()-1]);
         // }
         if(dist_Vertex(x_dv, y_dv,(*selectedTPs)[j]) < 0.5){
            // if((fabs(tp_z->at((*selectedTPs)[j]->index) - z_dv)) < 0.1){
            selectedTPs_dmin->push_back(new Track_Parameters((*selectedTPs)[j]->pt, (*selectedTPs)[j]->d0, (*selectedTPs)[j]->z0, (*selectedTPs)[j]->eta, (*selectedTPs)[j]->phi, (*selectedTPs)[j]->charge, (*selectedTPs)[j]->index, (*selectedTPs)[j]->pdgid));
            (*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->index),tp_y->at((*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->index));
            (*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedTPs_dmin)[selectedTPs_dmin->size()-1]);
         }
         // }

      }

      // Counters 
      // cout << i << "Start of Counters \n";

      int pt_check[5] = {0,0,0,0,0};
      int d0_check[5] = {0,0,0,0,0};
      int pt_d0_check[5][5];
      for (int k=0; k < 5; k++){
         for (int l=0; l < 5; l++){
            pt_d0_check[k][l] = 0;
         }
      }

      if ((selectedTracks->size() >= 2)){ nCount_2trk++;}

      if ((selectedTracks->size() >= 2)){ 
      for (int j = 0; j < selectedTracks->size()-1; j++){
         for (int k = j+1; k < selectedTracks->size(); k++){
            if(fabs((*selectedTracks)[j]->z0 - (*selectedTracks)[k]->z0) < 0.5 ){
               if ((*selectedTracks)[j]->pt   >   pt_cuts[0]){ pt_check[0]++; }
               if ((*selectedTracks)[j]->pt   >   pt_cuts[1]){ pt_check[1]++; }
               if ((*selectedTracks)[j]->pt   >   pt_cuts[2]){ pt_check[2]++; }
               if ((*selectedTracks)[j]->pt   >   pt_cuts[3]){ pt_check[3]++; }
               if ((*selectedTracks)[j]->pt   >   pt_cuts[4]){ pt_check[4]++; }

               if (fabs((*selectedTracks)[j]->d0)>d0_cuts[0]){ d0_check[0]++;}
               if (fabs((*selectedTracks)[j]->d0)>d0_cuts[1]){ d0_check[1]++;}
               if (fabs((*selectedTracks)[j]->d0)>d0_cuts[2]){ d0_check[2]++;}
               if (fabs((*selectedTracks)[j]->d0)>d0_cuts[3]){ d0_check[3]++;}
               if (fabs((*selectedTracks)[j]->d0)>d0_cuts[4]){ d0_check[4]++;}
               for (int l=0; l < 5; l++){
                  if ((*selectedTracks)[j]->pt   >   pt_cuts[0] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
                  if ((*selectedTracks)[j]->pt   >   pt_cuts[1] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
                  if ((*selectedTracks)[j]->pt   >   pt_cuts[2] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
                  if ((*selectedTracks)[j]->pt   >   pt_cuts[3] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
                  if ((*selectedTracks)[j]->pt   >   pt_cuts[4] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
               }
            }
         }
      }

      for(int k=0;k<5;k++){
         if(pt_check[k]>=1) nCount_trk_pt[k]++;
         if(d0_check[k]>=1) nCount_trk_d0[k]++;
         pt_check[k] = 0;
         d0_check[k] = 0;
         for (int l=0; l < 5; l++){
            if(pt_d0_check[k][l]>=1) nCount_trk_pt_d0[k][l]++;
            pt_d0_check[k][l] = 0;
         }
      }
      }

      // cout << i << "\t" << selectedTPs->size() << "\t" << (*selectedTPs)[0]->z0 << "\t" << (*selectedTPs)[0]->d0 << "\t" << (*selectedTPs)[0]->pt << "\t";
      
      if ((selectedTPs->size() == 2)){ nCount_2tp++;}

      if ((selectedTPs->size() >= 2)){ 
      for (int j = 0; j < selectedTPs->size()-1; j++){
         for (int k = j+1; k < selectedTPs->size(); k++){

            if(fabs((*selectedTPs)[j]->z0 - (*selectedTPs)[k]->z0) < 0.5 ){
               if ((*selectedTPs)[j]->pt   >   pt_cuts[0]){ pt_check[0]++; }
               if ((*selectedTPs)[j]->pt   >   pt_cuts[1]){ pt_check[1]++; }
               if ((*selectedTPs)[j]->pt   >   pt_cuts[2]){ pt_check[2]++; }
               if ((*selectedTPs)[j]->pt   >   pt_cuts[3]){ pt_check[3]++; }
               if ((*selectedTPs)[j]->pt   >   pt_cuts[4]){ pt_check[4]++; }

               if (fabs((*selectedTPs)[j]->d0)>d0_cuts[0]){ d0_check[0]++;}
               if (fabs((*selectedTPs)[j]->d0)>d0_cuts[1]){ d0_check[1]++;}
               if (fabs((*selectedTPs)[j]->d0)>d0_cuts[2]){ d0_check[2]++;}
               if (fabs((*selectedTPs)[j]->d0)>d0_cuts[3]){ d0_check[3]++;}
               if (fabs((*selectedTPs)[j]->d0)>d0_cuts[4]){ d0_check[4]++;}
               for (int l=0; l < 5; l++){
                  if ((*selectedTPs)[j]->pt   >   pt_cuts[0] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
                  if ((*selectedTPs)[j]->pt   >   pt_cuts[1] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
                  if ((*selectedTPs)[j]->pt   >   pt_cuts[2] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
                  if ((*selectedTPs)[j]->pt   >   pt_cuts[3] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
                  if ((*selectedTPs)[j]->pt   >   pt_cuts[4] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
               }
            }
         }
      }

      for(int k=0;k<5;k++){
         if(pt_check[k]>=1) nCount_tp_pt[k]++;
         if(d0_check[k]>=1) nCount_tp_d0[k]++;
         pt_check[k] = 0;
         d0_check[k] = 0;  
         for (int l=0; l < 5; l++){
            if(pt_d0_check[k][l]>=1) nCount_tp_pt_d0[k][l]++;
            pt_d0_check[k][l] = 0;
         }
      }
      }
      // cout << i << "End of Counters \n";


      sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);
      sort(selectedTPs_zmin->begin(), selectedTPs_zmin->end(), ComparePtTrack);
      sort(selectedTPs_dmin->begin(), selectedTPs_dmin->end(), ComparePtTrack);
      TP_minz0_nindex = -1;
      TP_mind0_nindex = -1;
      TP_minz0_pindex = -1;
      TP_mind0_pindex = -1;

      // cout << i << "\t" << selectedTPs->size() << "\t" << (*selectedTPs)[0]->z0 << "\t" << (*selectedTPs)[0]->d0 << "\t" << (*selectedTPs)[0]->pt << "\t";
      int check = 0;
      std::vector<int> mu_zmin_p_index;
      std::vector<int> mu_zmin_n_index;
      std::vector<int> mu_dmin_p_index;
      std::vector<int> mu_dmin_n_index;
      for (int j = 0; j < selectedTPs_zmin->size(); j++){
         if(abs((*selectedTPs_zmin)[j]->pdgid)==13){
            if((*selectedTPs_zmin)[j]->charge==1) {
               mu_zmin_p_index.push_back(selectedMu_zmin->size());
            }
            else
            {
               mu_zmin_n_index.push_back(selectedMu_zmin->size());
            }
            selectedMu->push_back(new Track_Parameters((*selectedTPs_zmin)[j]->pt, (*selectedTPs_zmin)[j]->d0, (*selectedTPs_zmin)[j]->z0, (*selectedTPs_zmin)[j]->eta, (*selectedTPs_zmin)[j]->phi, (*selectedTPs_zmin)[j]->charge, (*selectedTPs_zmin)[j]->index, (*selectedTPs_zmin)[j]->pdgid));
            (*selectedMu)[selectedMu->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedMu)[selectedMu->size()-1]->index),tp_y->at((*selectedMu)[selectedMu->size()-1]->index));
            (*selectedMu)[selectedMu->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedMu)[selectedMu->size()-1]);
            // selectedMu_dmin->push_back(new Track_Parameters((*selectedTPs_zmin)[j]->pt, (*selectedTPs_zmin)[j]->d0, (*selectedTPs_zmin)[j]->z0, (*selectedTPs_zmin)[j]->eta, (*selectedTPs_zmin)[j]->phi, (*selectedTPs_zmin)[j]->charge, j, (*selectedTPs_zmin)[j]->pdgid));
            selectedMu_zmin->push_back(new Track_Parameters((*selectedTPs_zmin)[j]->pt, (*selectedTPs_zmin)[j]->d0, (*selectedTPs_zmin)[j]->z0, (*selectedTPs_zmin)[j]->eta, (*selectedTPs_zmin)[j]->phi, (*selectedTPs_zmin)[j]->charge, (*selectedTPs_zmin)[j]->index, (*selectedTPs_zmin)[j]->pdgid));
            (*selectedMu_zmin)[selectedMu_zmin->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedMu_zmin)[selectedMu_zmin->size()-1]->index),tp_y->at((*selectedMu_zmin)[selectedMu_zmin->size()-1]->index));
            (*selectedMu_zmin)[selectedMu_zmin->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedMu_zmin)[selectedMu_zmin->size()-1]);
         }
      }
      for (int j = 0; j < selectedTPs_dmin->size(); j++){
         if(abs((*selectedTPs_dmin)[j]->pdgid)==13){
            if((*selectedTPs_dmin)[j]->charge==1) {
               mu_dmin_p_index.push_back(selectedMu_dmin->size());
            }
            else
            {
               mu_dmin_n_index.push_back(selectedMu_dmin->size());
            }
            selectedMu->push_back(new Track_Parameters((*selectedTPs_dmin)[j]->pt, (*selectedTPs_dmin)[j]->d0, (*selectedTPs_dmin)[j]->z0, (*selectedTPs_dmin)[j]->eta, (*selectedTPs_dmin)[j]->phi, (*selectedTPs_dmin)[j]->charge, (*selectedTPs_dmin)[j]->index, (*selectedTPs_dmin)[j]->pdgid));
            (*selectedMu)[selectedMu->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedMu)[selectedMu->size()-1]->index),tp_y->at((*selectedMu)[selectedMu->size()-1]->index));
            (*selectedMu)[selectedMu->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedMu)[selectedMu->size()-1]);
            selectedMu_dmin->push_back(new Track_Parameters((*selectedTPs_dmin)[j]->pt, (*selectedTPs_dmin)[j]->d0, (*selectedTPs_dmin)[j]->z0, (*selectedTPs_dmin)[j]->eta, (*selectedTPs_dmin)[j]->phi, (*selectedTPs_dmin)[j]->charge, (*selectedTPs_dmin)[j]->index, (*selectedTPs_dmin)[j]->pdgid));
            (*selectedMu_dmin)[selectedMu_dmin->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedMu_dmin)[selectedMu_dmin->size()-1]->index),tp_y->at((*selectedMu_dmin)[selectedMu_dmin->size()-1]->index));
            (*selectedMu_dmin)[selectedMu_dmin->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedMu_dmin)[selectedMu_dmin->size()-1]);
            // selectedMu_zmin->push_back(new Track_Parameters((*selectedTPs_dmin)[j]->pt, (*selectedTPs_dmin)[j]->d0, (*selectedTPs_dmin)[j]->z0, (*selectedTPs_dmin)[j]->eta, (*selectedTPs_dmin)[j]->phi, (*selectedTPs_dmin)[j]->charge, j, (*selectedTPs_dmin)[j]->pdgid));
         }
      }
      // cout<<endl<<i;
      // cout << "\tMuon selection \t" << selectedMu->size()<<"\t"<<mu_n_index.size()<<"\t"<<mu_p_index.size()<<"\t";
      if ((selectedMu_zmin->size()>=2) && (selectedMu_dmin->size()>=2)){
         nCount_2mu++;
      }
      if ((selectedMu_zmin->size() >= 2) && (mu_zmin_n_index.size() != 0) && (mu_zmin_p_index.size() != 0)){
         // nCount_2muOS++;
         
         float TP_minz0 = std::numeric_limits<float>::infinity();
         for (int j = 0; j < mu_zmin_p_index.size(); j++)
         {
            for (int k = 0; k < mu_zmin_n_index.size(); k++)
            {
               if (fabs((*selectedMu_zmin)[mu_zmin_p_index[j]]->z0 - (*selectedMu_zmin)[mu_zmin_n_index[k]]->z0) < TP_minz0)
               {
                  TP_minz0 = fabs((*selectedMu_zmin)[mu_zmin_p_index[j]]->z0 - (*selectedMu_zmin)[mu_zmin_n_index[k]]->z0);
                  TP_minz0_pindex = j;
                  TP_minz0_nindex = k;
               }
            }
         }
      }
      if ((selectedMu_dmin->size() >= 2) && (mu_dmin_n_index.size() != 0) && (mu_dmin_p_index.size() != 0)){
         nCount_2muOS++;
         float TP_mind0 = std::numeric_limits<float>::infinity();
         for (int j = 0; j < mu_dmin_p_index.size(); j++)
         {
            for (int k = 0; k < mu_dmin_n_index.size(); k++)
            {
               if (fabs((*selectedMu_dmin)[mu_dmin_p_index[j]]->d0 - (*selectedMu_dmin)[mu_dmin_n_index[k]]->d0) < TP_mind0)
               {
                  TP_mind0 = fabs((*selectedMu_dmin)[mu_dmin_p_index[j]]->d0 - (*selectedMu_dmin)[mu_dmin_n_index[k]]->d0);
                  TP_mind0_pindex = j;
                  TP_mind0_nindex = k;
               }
            }
         }
      }

         // cout << TP_minz0_nindex << "\t" << TP_mind0_nindex<<"\t";

         // if((*selectedMu_zmin)[mu_n_index[TP_minz0_nindex]]->charge==1)
            // cout<<" \t Charge = "<<(*selectedMu_zmin)[mu_p_index[TP_minz0_nindex]]->charge <<" \t PDGID = "<<(*selectedMu_zmin)[mu_p_index[TP_minz0_nindex]]->pdgid;

         // Swapping to get the Muon as the first element followed by the Anti-Muon
         /*
         std::iter_swap((*selectedMu_zmin)[mu_n_index[TP_minz0_nindex]], (*selectedMu_zmin)[0]);
         std::iter_swap((*selectedMu_zmin)[mu_p_index[TP_minz0_pindex]], (*selectedMu_zmin)[1]);

         std::iter_swap((*selectedMu_dmin)[mu_n_index[TP_mind0_nindex]], (*selectedMu_dmin)[0]);
         std::iter_swap((*selectedMu_dmin)[mu_p_index[TP_mind0_pindex]], (*selectedMu_dmin)[1]);
   */
      /*
      if (selectedMu->size()>=2){
      
         float TP_minz0 = std::numeric_limits<float>::infinity();
         float TP_mind0 = std::numeric_limits<float>::infinity();
         sort(selectedMu_zmin->begin(), selectedMu_zmin->end(), CompareZ0Track);
         sort(selectedMu_dmin->begin(), selectedMu_dmin->end(), CompareD0Track);
         for (int j = 0; j < selectedMu_zmin->size() - 1; j++)
         {
            if (fabs((*selectedMu_zmin)[j + 1]->z0 - (*selectedMu_zmin)[j]->z0) < TP_minz0)
            {
               TP_minz0 = fabs((*selectedMu_zmin)[j + 1]->z0 - (*selectedMu_zmin)[j]->z0);
               TP_minz0_nindex = j;
            }
            if (fabs((*selectedMu_dmin)[j + 1]->d0 - (*selectedMu_dmin)[j]->d0) < TP_mind0)
            {
               TP_mind0 = fabs((*selectedMu_dmin)[j + 1]->d0 - (*selectedMu_dmin)[j]->d0);
               TP_mind0_nindex = j;
            }
         }

         // cout << TP_minz0_index << "\t" << TP_mind0_index<<"\t";

         if ((*selectedMu_zmin)[TP_minz0_nindex + 1]->pt > (*selectedMu_zmin)[TP_minz0_nindex]->pt)
         {
            std::iter_swap((*selectedMu_zmin)[TP_minz0_nindex + 1], (*selectedMu_zmin)[0]);
            std::iter_swap((*selectedMu_zmin)[TP_minz0_nindex], (*selectedMu_zmin)[1]);
         }
         else
         {
            std::iter_swap((*selectedMu_zmin)[TP_minz0_nindex], (*selectedMu_zmin)[0]);
            std::iter_swap((*selectedMu_zmin)[TP_minz0_nindex + 1], (*selectedMu_zmin)[1]);
         }

         if ((*selectedMu_dmin)[TP_mind0_nindex + 1]->pt > (*selectedMu_dmin)[TP_mind0_nindex]->pt)
         {
            std::iter_swap((*selectedMu_dmin)[TP_mind0_nindex + 1], (*selectedMu_dmin)[0]);
            std::iter_swap((*selectedMu_dmin)[TP_mind0_nindex], (*selectedMu_dmin)[1]);
         }
         else
         {
            std::iter_swap((*selectedMu_dmin)[TP_mind0_nindex], (*selectedMu_dmin)[0]);
            std::iter_swap((*selectedMu_dmin)[TP_mind0_nindex + 1], (*selectedMu_dmin)[1]);
         }
      }
      // cout << "End of d0 and z0 sorting \t";

*/
      if (!(selectedTPs->size() >= 2)) continue;
      // cout<<selectedTPs_dmin->size()<<"\t"<<selectedTPs_dmin->size()<<endl;
      nAccept++;
         // ---------------------------------------------------------------------------------------------------------
         //Filling up Histograms

         // h_d0_Mu_n_Mu_p->Fill((*selectedMu_dmin)[mu_n_index[TP_mind0_nindex]]->d0,(*selectedMu_dmin)[mu_p_index[TP_mind0_pindex]]->d0);
         
         h_index=0;
         h_displaced_vertex_x->Fill( (tp_x->at((*selectedTPs)[0]->index)),x_dv);
         h_displaced_vertex_y->Fill( (tp_y->at((*selectedTPs)[0]->index)),y_dv);
         Hists[0][h_index++]->Fill(selectedTPs->size());
         Hists[0][h_index++]->Fill(selectedMu->size());
         for (int i = 0; i < selectedTPs->size(); i++)
         {
            Hists[0][h_index]->Fill((*selectedTPs)[i]->pt);
            Hists[0][h_index+1]->Fill((*selectedTPs)[i]->eta);
            Hists[0][h_index+2]->Fill((*selectedTPs)[i]->phi);
            Hists[0][h_index+3]->Fill((*selectedTPs)[i]->z0);
            Hists[0][h_index+4]->Fill((*selectedTPs)[i]->d0);
            Hists[0][h_index+5]->Fill((*selectedTPs)[i]->charge);
            Hists[0][h_index+6]->Fill((*selectedTPs)[i]->pdgid);
            Hists[0][h_index+7]->Fill((*selectedTPs)[i]->dxy);
         }
         h_index +=8;
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->dxy);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->dxy);
         Hists[0][h_index++]->Fill(fabs((*selectedTPs)[0]->eta - (*selectedTPs)[1]->eta));
         Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTPs)[0]->phi,(*selectedTPs)[1]->phi)));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
         Hists[0][h_index++]->Fill((deltaR((*selectedTPs)[0]->eta,(*selectedTPs)[0]->phi, (*selectedTPs)[1]->eta,(*selectedTPs)[1]->phi)));

         if ((selectedTPs_dmin->size() >= 2) && (selectedTPs_zmin->size() >= 2)) {
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->dxy);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->dxy);
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->eta - (*selectedTPs_zmin)[1]->eta));
         Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
         Hists[0][h_index++]->Fill((deltaR((*selectedTPs_zmin)[0]->eta, (*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->eta, (*selectedTPs_zmin)[1]->phi)));
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->dxy);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->dxy);
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->eta - (*selectedTPs_dmin)[1]->eta));
         Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
         Hists[0][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
         Hists[0][h_index++]->Fill((deltaR((*selectedTPs_dmin)[0]->eta,(*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->eta,(*selectedTPs_dmin)[1]->phi)));
         }
         // cout<<"TPs_dmin TPs_zmin \t"<<fabs((*selectedMu_dmin)[mu_n_index[TP_mind0_nindex]]->eta - (*selectedMu_dmin)[mu_p_index[TP_mind0_pindex]]->eta);

         if ((selectedMu_zmin->size() >= 2) && (selectedMu_dmin->size() >= 2) && (mu_zmin_n_index.size() != 0) && (mu_zmin_p_index.size() != 0) && (mu_dmin_n_index.size() != 0) && (mu_dmin_p_index.size() != 0)){
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->pt);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->z0);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->d0);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->charge);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->pdgid);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->dxy);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->pt);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->z0);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->d0);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->charge);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->pdgid);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->dxy);
         Hists[0][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta));
         Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi)));
         Hists[0][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->z0 - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->z0));
         Hists[0][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->d0 - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->d0));
         Hists[0][h_index++]->Fill((deltaR((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta, (*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi)));
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->pt);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->z0);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->d0);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->charge);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->pdgid);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_dmin_n_index[TP_mind0_nindex]]->dxy);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->pt);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->z0);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->d0);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->charge);
         Hists[0][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->pdgid);
         Hists[0][h_index++]->Fill((*selectedMu_zmin)[mu_dmin_p_index[TP_mind0_pindex]]->dxy);
         Hists[0][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta));
         Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi, (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi)));
         Hists[0][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->z0 - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->z0));
         Hists[0][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->d0 - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->d0));
         Hists[0][h_index++]->Fill((deltaR((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta,(*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi, (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta,(*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi)));
         }

         // cout<<"end of Histogram filling"<<"\n";

         // Atleast 1 pair of opposite sign muons
         if ((selectedMu_zmin->size() >= 2) && (selectedMu_dmin->size() >= 2) && (mu_zmin_n_index.size() != 0) && (mu_zmin_p_index.size() != 0) && (mu_dmin_n_index.size() != 0) && (mu_dmin_p_index.size() != 0) && (selectedTPs_dmin->size() >= 2) && (selectedTPs_zmin->size() >= 2))
         {
            // h_displaced_vertex->Fill(x_dv,y_dv);
            h_d0_Mu_n_Mu_p->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->d0,(*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->d0);
            h_z0_Mu_n_Mu_p->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->z0,(*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->z0);
            h_index = 0;
            Hists[1][h_index++]->Fill(selectedTPs->size());
            Hists[1][h_index++]->Fill(selectedMu->size());
            for (int i = 0; i < selectedTPs->size(); i++)
            {
               Hists[1][h_index]->Fill((*selectedTPs)[i]->pt);
               Hists[1][h_index + 1]->Fill((*selectedTPs)[i]->eta);
               Hists[1][h_index + 2]->Fill((*selectedTPs)[i]->phi);
               Hists[1][h_index + 3]->Fill((*selectedTPs)[i]->z0);
               Hists[1][h_index + 4]->Fill((*selectedTPs)[i]->d0);
               Hists[1][h_index + 5]->Fill((*selectedTPs)[i]->charge);
               Hists[1][h_index + 6]->Fill((*selectedTPs)[i]->pdgid);
               Hists[1][h_index + 7]->Fill((*selectedTPs)[i]->dxy);
            }
            h_index += 8;
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->dxy);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->dxy);
            Hists[1][h_index++]->Fill(fabs((*selectedTPs)[0]->eta - (*selectedTPs)[1]->eta));
            Hists[1][h_index++]->Fill(fabs(deltaPhi((*selectedTPs)[0]->phi,(*selectedTPs)[1]->phi)));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
            Hists[1][h_index++]->Fill((deltaR((*selectedTPs)[0]->eta,(*selectedTPs)[0]->phi, (*selectedTPs)[1]->eta,(*selectedTPs)[1]->phi)));

            
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->dxy);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->dxy);
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->eta - (*selectedTPs_zmin)[1]->eta));
            Hists[1][h_index++]->Fill(fabs(deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
            Hists[1][h_index++]->Fill((deltaR((*selectedTPs_zmin)[0]->eta, (*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->eta, (*selectedTPs_zmin)[1]->phi)));
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->dxy);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->dxy);
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->eta - (*selectedTPs_dmin)[1]->eta));
            Hists[1][h_index++]->Fill(fabs(deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
            Hists[1][h_index++]->Fill(fabs((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
            Hists[1][h_index++]->Fill((deltaR((*selectedTPs_dmin)[0]->eta,(*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->eta,(*selectedTPs_dmin)[1]->phi)));
            

            // if ((*selectedMu_zmin)[0]->charge == 1)
               // cout <<" \t Charge = " << (*selectedMu_zmin)[0]->charge << " \t PDGID = " << (*selectedMu_zmin)[0]->pdgid;

            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->pt);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->z0);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->d0);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->charge);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->pdgid);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->dxy);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->pt);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->z0);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->d0);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->charge);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->pdgid);
            Hists[1][h_index++]->Fill((*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->dxy);
            Hists[1][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta));
            Hists[1][h_index++]->Fill(fabs(deltaPhi((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi)));
            Hists[1][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->z0 - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->z0));
            Hists[1][h_index++]->Fill(fabs((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->d0 - (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->d0));
            Hists[1][h_index++]->Fill((deltaR((*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->eta, (*selectedMu_zmin)[mu_zmin_n_index[TP_minz0_nindex]]->phi, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->eta, (*selectedMu_zmin)[mu_zmin_p_index[TP_minz0_pindex]]->phi)));
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->pt);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->z0);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->d0);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->charge);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->pdgid);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->dxy);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->pt);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->z0);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->d0);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->charge);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->pdgid);
            Hists[1][h_index++]->Fill((*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->dxy);
            Hists[1][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta));
            Hists[1][h_index++]->Fill(fabs(deltaPhi((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi, (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi)));
            Hists[1][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->z0 - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->z0));
            Hists[1][h_index++]->Fill(fabs((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->d0 - (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->d0));
            Hists[1][h_index++]->Fill((deltaR((*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->eta, (*selectedMu_dmin)[mu_dmin_n_index[TP_mind0_nindex]]->phi, (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->eta, (*selectedMu_dmin)[mu_dmin_p_index[TP_mind0_pindex]]->phi)));
         }

      for (int l = 0; l < selectedTracks->size(); l++)
      {
         delete (*selectedTracks)[l];
         // delete (*selectedTracks_zmin)[l];
         // delete (*selectedTracks_dmin)[l];
      }

      for (int l = 0; l < selectedMu->size(); l++)
      {
         // delete (*selectedTPs)[l];
         delete (*selectedMu)[l];
      }
      for (int l = 0; l < selectedMu_zmin->size(); l++)
      {
         delete (*selectedMu_zmin)[l];
      }
      for (int l = 0; l < selectedMu_dmin->size(); l++)
      {
         delete (*selectedMu_dmin)[l];
      }
      for (int l = 0; l < selectedTPs->size(); l++)
      {
         delete (*selectedTPs)[l];
      }
      for (int l = 0; l < selectedTPs_dmin->size(); l++)
      {
         delete (*selectedTPs_dmin)[l];
      }
      for (int l = 0; l < selectedTPs_zmin->size(); l++)
      {
         delete (*selectedTPs_zmin)[l];
      }

      selectedTracks->clear();
      selectedTracks->shrink_to_fit();
      delete selectedTracks;
/*
      selectedTracks_zmin->clear();
      selectedTracks_zmin->shrink_to_fit();
      delete selectedTracks_zmin;
      selectedTracks_dmin->clear();
      selectedTracks_dmin->shrink_to_fit();
      delete selectedTracks_dmin;
*/
      selectedTPs->clear();
      selectedTPs->shrink_to_fit();
      delete selectedTPs;
      selectedMu->clear();
      selectedMu->shrink_to_fit();
      delete selectedMu;
      selectedMu_zmin->clear();
      selectedMu_zmin->shrink_to_fit();
      delete selectedMu_zmin;
      selectedMu_dmin->clear();
      selectedMu_dmin->shrink_to_fit();
      delete selectedMu_dmin;
      selectedTPs_dmin->clear();
      selectedTPs_dmin->shrink_to_fit();
      delete selectedTPs_dmin;
      selectedTPs_zmin->clear();
      selectedTPs_zmin->shrink_to_fit();
      delete selectedTPs_zmin;

      // cout<<endl;

   } // end event loop
   cout << endl
        << "from " << nevt << " events, " << nAccept << " events are accepted"<<endl;

   cout << "from " << nevt << " events, " << nCount_2mu << " events have 2 Muons" << endl;
   cout << "from " << nevt << " events, " << nCount_2muOS << " events have atleast 1 pair or opposite sign 2 Muons" << endl;
   // ---------------------------------------------------------------------------------------------------------
   //some Histograms

   char ctxt[500];
   TCanvas c;

   TString DIR = type_dir + "AnalyzerTrkPlots/";
   TString makedir = "mkdir -p " + DIR;
   const char *mkDIR = makedir.Data();
   gSystem->Exec(mkDIR);

   TFile *fout;
   fout = new TFile(type_dir + "output_" + type + ".root", "recreate");

   for (int k = 0; k < regions.size(); ++k)
   {
      makedir = "mkdir -p " + DIR + regions[k] + "/";
      gSystem->Exec(makedir.Data());
      for (int l = 0; l < vars.size(); ++l)
      {
         Hists[k][l]->GetXaxis()->SetRange(1, Hists[k][l]->GetNbinsX() + 2);
         Hists[k][l]->Draw();
         Hists[k][l]->GetYaxis()->SetTitle("Events");/*
         if (vars[l].Contains("pt"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("Pt");
         }
         else if (vars[l].Contains("eta"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("Eta");
         }
         else if (vars[l].Contains("phi"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("Phi");
         }
         else if (vars[l].Contains("z0"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("Z0");
         }
         else if (vars[l].Contains("d0"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("d0");
         }
         else if (vars[l].Contains("charge"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("charge");
         }*/
         Hists[k][l]->GetXaxis()->SetTitle(Hists[k][l]->GetName());
         Hists[k][l]->Write("", TObject::kOverwrite);
         c.SaveAs(DIR + regions[k] + "/" + type + "_" + Hists[k][l]->GetName() + ".pdf");
      }
   }
   for (int k = 0; k < regions.size(); ++k)
   {
      for (int l = 0; l < vars.size(); ++l)
      {
         delete Hists[k][l];
      }
   }
   // h_d0_Mu_n_Mu_p->GetXaxis()->SetRange(1, h_d0_Mu_n_Mu_p->GetNbinsX() + 2);
   h_d0_Mu_n_Mu_p->Draw("colz");
   h_d0_Mu_n_Mu_p->GetYaxis()->SetTitle("d0 Mu+");
   h_d0_Mu_n_Mu_p->GetXaxis()->SetTitle("d0 Mu-");
   h_d0_Mu_n_Mu_p->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_d0_Mu_n_Mu_p->GetName() + ".pdf");
   h_z0_Mu_n_Mu_p->Draw("colz");
   h_z0_Mu_n_Mu_p->GetYaxis()->SetTitle("z0 Mu+");
   h_z0_Mu_n_Mu_p->GetXaxis()->SetTitle("z0 Mu-");
   h_z0_Mu_n_Mu_p->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_z0_Mu_n_Mu_p->GetName() + ".pdf");
   delete h_d0_Mu_n_Mu_p;
   delete h_z0_Mu_n_Mu_p;
   h_displaced_vertex_x->Draw("colz");
   h_displaced_vertex_x->GetYaxis()->SetTitle("Leading p_t - x");
   h_displaced_vertex_x->GetXaxis()->SetTitle("Displaced vertex - x");
   h_displaced_vertex_x->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_displaced_vertex_x->GetName() + ".pdf");
   delete h_displaced_vertex_x;
   h_displaced_vertex_y->Draw("colz");
   h_displaced_vertex_y->GetYaxis()->SetTitle("Leading p_t - y");
   h_displaced_vertex_y->GetXaxis()->SetTitle("Displaced vertex - y");
   h_displaced_vertex_y->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_displaced_vertex_y->GetName() + ".pdf");
   delete h_displaced_vertex_y;
   fout->Close();

   ofstream myfile;
   myfile.open(type_dir+type+"_"+"Counters.txt");

   myfile<<endl<<"Total Number of Events = "<<nevt;
   myfile << endl;
   myfile << endl << "Tracking Particles"<<endl; 
   myfile<<endl<<"Number of Events with atleast 2 tracking particles = "<<nAccept;
   myfile<<endl<<"Number of Events with Exactly 2 tracking particles = "<<nCount_2tp;
   myfile<<endl<<"Number of Events with atleast 1 pair of muons = "<<nCount_2mu;
   myfile<<endl<<"Number of Events with atleast 1 pair of muons of opposite sign = "<<nCount_2muOS;

   myfile<<endl;
   std::stringstream txt;
   float rate;
   if(type.Contains("NeutrinoGun")){
      rate = 40000.0;
      txt << " khz ";
   }
   else{
      rate = 100.0;
      txt << " % ";
   }
   for(int k=0; k<5 ;k++){
      myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with "<<"Z0 < 0.5 and 1 pair of Tracking Particles with Pt greater than "<<pt_cuts[k]<<" = "<<nCount_tp_pt[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_pt[k]<<")";
   }
   for(int k=0; k<5 ;k++){
      myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with "<<"Z0 < 0.5 and 1 pair of Tracking Particles with D0 greater than "<<d0_cuts[k]<<" = "<<nCount_tp_d0[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_d0[k]<<")";
   }
   // myfile<<endl<<"Track Rate with "<< (txt.str()).c_str()<<endl;
   myfile<<"\n\n\n";
   int tablesize = 10;
   for (int k=-1; k<5; k++){
      if(k==-1){
         myfile << std::right << std::setw(tablesize*2)<<" d0 \\ Pt ";
         for (int l=0; l<5; l++){
            myfile << std::right << std::setw(tablesize*2)<< pt_cuts[l];
         }
      }
      else{
         for (int l=-1; l<5; l++){
            if(l==-1){ 
               myfile << std::right << std::setw(tablesize*2) << d0_cuts[k];
            }
            else{
            myfile <<std::fixed<<std::setprecision(3)<< std::right<< std::setw(tablesize) << nCount_tp_pt_d0[l][k]* rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_pt_d0[l][k]<<")";
            }
         }
      }
      myfile<<endl;
   }

   myfile << endl;
   myfile << endl << "Tracks"<<endl; 
   myfile<<endl<<"Number of Events with Atleast 2 tracks = "<<nCount_2trk;
   myfile<<endl;
   for(int k=0; k<5 ;k++){
      myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with chi2rz_dof<10, chi2rphi_dof<10, bendchi2<10, "<<"Z0 < 0.5 and 1 pair of Tracks with Pt greater than "<<pt_cuts[k]<<" = "<<nCount_trk_pt[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_pt[k]<<")";
   }
   for(int k=0; k<5 ;k++){
      myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with chi2rz_dof<10, chi2rphi_dof<10, bendchi2<10, "<<"Z0 < 0.5 and 1 pair of Tracks with D0 greater than "<<d0_cuts[k]<<" = "<<nCount_trk_d0[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_d0[k]<<")";
   }
   myfile<<"\n\n\n";
   for (int k=-1; k<5; k++){
      if(k==-1){
         myfile << std::right << std::setw(tablesize*2) <<" d0 \\ Pt ";
         for (int l=0; l<5; l++){
            myfile << std::right<< std::setw(tablesize*2) << pt_cuts[l];
         }
      }
      else{
         for (int l=-1; l<5; l++){
            if(l==-1){ 
               myfile <<std::fixed<<std::setprecision(3)<< std::right << std::setw(tablesize*2) << d0_cuts[k];
            }
            else{
            myfile <<std::fixed<<std::setprecision(3)<< std::right<< std::setw(tablesize) <<nCount_trk_pt_d0[l][k]* rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_pt_d0[l][k]<<")";
            }
         }
      }
      myfile<<endl;
   }

   myfile.close();
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
   gStyle->SetPadRightMargin(0.05);
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
   gStyle->SetOptStat(0);
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
