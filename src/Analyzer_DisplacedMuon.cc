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
   float TP_maxPt = 100.0;
   float TP_maxEta = 2.4;
   int TP_select_pdgid = 0;

   // ----------------------------------------------------------------------------------------------------------------
   // histograms
   // ----------------------------------------------------------------------------------------------------------------

   typedef vector<TH1F *> Dim1;
   typedef vector<Dim1> Dim2;
   typedef vector<Dim2> Dim3;

   std::vector<TString> regions{"alltrk", "trkeq2", "trkgeq2", "trkgeq2os", "tpeq2", "tpgeq2", "tpgeq2osMu"};
   std::vector<TString> vars{"trk_pt", "trk_eta", "trk_phi", "trk_z0", "trk_d0", "trk_rinv",
                             "t1pt", "t1eta", "t1phi", "t1z0", "t1d0", "t1rinv",
                             "t2pt", "t2eta", "t2phi", "t2z0", "t2d0", "t2rinv",
                             "ntrk", "ntp",
                             "tp1pt", "tp1eta", "tp1phi", "tp1z0", "tp1d0", "tp1rinv", "tp1pdgid",
                             "tp2pt", "tp2eta", "tp2phi", "tp2z0", "tp2d0", "tp2rinv", "tp2pdgid"};
   std::vector<int> nbins{50, 50, 50, 50, 50, 50,
                          50, 50, 50, 50, 50, 50,
                          50, 50, 50, 50, 50, 50,
                          15, 15,
                          50, 50, 50, 50, 50, 50, 350,
                          50, 50, 50, 50, 50, 50, 350};
   std::vector<float> lowEdge{0, -2.4, -4, -20, -10, -0.02,
                              0, -2.4, -4, -20, -10, -0.02,
                              0, -2.4, -4, -20, -10, -0.02,
                              0, 0,
                              0, -2.4, -4, -20, -10, -0.02, 0,
                              0, -2.4, -4, -20, -10, -0.02, 0};
   std::vector<float> highEdge{50, 2.4, 4, 20, 10, 0.02,
                               50, 2.4, 4, 20, 10, 0.02,
                               50, 2.4, 4, 20, 10, 0.02,
                               15, 15,
                               50, 2.4, 4, 20, 10, 0.02, 350,
                               50, 2.4, 4, 20, 10, 0.02, 350};

   typedef vector<TH1F *> Dim1;
   typedef vector<Dim1> Dim2;
   typedef vector<Dim2> Dim3;
   typedef vector<Dim3> Dim4;
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
   Long64_t nevt = fChain->GetEntries();

   for (Long64_t i=0; i<nevt; i++) {
      fChain->GetEntry(i);
      // if (Cut(ientry) < 0) continue;

      displayProgress(i, nevt);

      // ----------------------------------------------------------------------------------------------------------------
      // track particle loop

      int ntrk = 0, ntp = 0;
      selectedTracks = new std::vector<Track_Parameters *>();
      selectedTPs = new std::vector<Track_Parameters *>();

      for (int it = 0; it < (int)trk_pt->size(); it++)
      {
         if (abs(trk_eta->at(it)) > TP_maxEta)
            continue;
         if (trk_pt->at(it) < TP_minPt)
            continue;
         if (trk_pt->at(it) > TP_maxPt)
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
         ntrk++;
         selectedTracks->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), trk_rinv->at(it), it, -99999));
      }

      // ----------------------------------------------------------------------------------------------------------------
      // tracking particle loop

      for (int it = 0; it < (int)tp_pt->size(); it++)
      {
         if (std::abs(tp_d0->at(it)) > TP_maxD0)
            continue;
         if (std::abs(tp_d0->at(it)) < TP_minD0)
            continue;
         if (tp_pt->at(it) < 0.2)
            continue;
         if (tp_pt->at(it) > TP_maxPt)
            continue;
         if (std::abs(tp_eta->at(it)) > TP_maxEta)
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
            if (abs(tp_pdgid->at(it)) != abs(TP_select_pdgid))
               continue;
         }
         ntp++;
         selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tp_d0->at(it), tp_z0->at(it), tp_eta->at(it), tp_phi->at(it), matchtrk_rinv->at(it), it, tp_pdgid->at(it)));
      } // end tp loop

      sort(selectedTracks->begin(),selectedTracks->end(),ComparePtTrack);
      sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);

      // ---------------------------------------------------------------------------------------------------------
      //Filling up Histograms

      // All tracks
      for (int i = 0; i < ntrk; i++)
      {
         Hists[0][0]->Fill((*selectedTracks)[i]->pt);
         Hists[0][1]->Fill((*selectedTracks)[i]->eta);
         Hists[0][2]->Fill((*selectedTracks)[i]->phi);
         Hists[0][3]->Fill((*selectedTracks)[i]->z0);
         Hists[0][4]->Fill((*selectedTracks)[i]->d0);
         Hists[0][5]->Fill((*selectedTracks)[i]->rinv);
      }
      if (ntrk > 0)
      {
         nAccept++;
         Hists[0][6]->Fill((*selectedTracks)[0]->pt);
         Hists[0][7]->Fill((*selectedTracks)[0]->eta);
         Hists[0][8]->Fill((*selectedTracks)[0]->phi);
         Hists[0][9]->Fill((*selectedTracks)[0]->z0);
         Hists[0][10]->Fill((*selectedTracks)[0]->d0);
         Hists[0][11]->Fill((*selectedTracks)[0]->rinv);
      }
      if (ntrk > 1)
      {
         Hists[0][12]->Fill((*selectedTracks)[1]->pt);
         Hists[0][13]->Fill((*selectedTracks)[1]->eta);
         Hists[0][14]->Fill((*selectedTracks)[1]->phi);
         Hists[0][15]->Fill((*selectedTracks)[1]->z0);
         Hists[0][16]->Fill((*selectedTracks)[1]->d0);
         Hists[0][17]->Fill((*selectedTracks)[1]->rinv);
      }
      Hists[0][18]->Fill(ntrk);
      Hists[0][19]->Fill(ntp);
      if (ntp > 0)
      {
         Hists[0][20]->Fill((*selectedTPs)[0]->pt);
         Hists[0][21]->Fill((*selectedTPs)[0]->eta);
         Hists[0][22]->Fill((*selectedTPs)[0]->phi);
         Hists[0][23]->Fill((*selectedTPs)[0]->z0);
         Hists[0][24]->Fill((*selectedTPs)[0]->d0);
         Hists[0][25]->Fill((*selectedTPs)[0]->rinv);
         Hists[0][26]->Fill((*selectedTPs)[0]->pdgid);
      }
      if (ntp > 1)
      {
         Hists[0][27]->Fill((*selectedTPs)[1]->pt);
         Hists[0][28]->Fill((*selectedTPs)[1]->eta);
         Hists[0][29]->Fill((*selectedTPs)[1]->phi);
         Hists[0][30]->Fill((*selectedTPs)[1]->z0);
         Hists[0][31]->Fill((*selectedTPs)[1]->d0);
         Hists[0][32]->Fill((*selectedTPs)[1]->rinv);
         Hists[0][33]->Fill((*selectedTPs)[1]->pdgid);
      }

      // Number of Tracks == 2
      if (ntrk == 2)
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[1][0]->Fill((*selectedTracks)[i]->pt);
            Hists[1][1]->Fill((*selectedTracks)[i]->eta);
            Hists[1][2]->Fill((*selectedTracks)[i]->phi);
            Hists[1][3]->Fill((*selectedTracks)[i]->z0);
            Hists[1][4]->Fill((*selectedTracks)[i]->d0);
            Hists[1][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[1][6]->Fill((*selectedTracks)[0]->pt);
         Hists[1][7]->Fill((*selectedTracks)[0]->eta);
         Hists[1][8]->Fill((*selectedTracks)[0]->phi);
         Hists[1][9]->Fill((*selectedTracks)[0]->z0);
         Hists[1][10]->Fill((*selectedTracks)[0]->d0);
         Hists[1][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[1][12]->Fill((*selectedTracks)[1]->pt);
         Hists[1][13]->Fill((*selectedTracks)[1]->eta);
         Hists[1][14]->Fill((*selectedTracks)[1]->phi);
         Hists[1][15]->Fill((*selectedTracks)[1]->z0);
         Hists[1][16]->Fill((*selectedTracks)[1]->d0);
         Hists[1][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[1][18]->Fill(ntrk);
         Hists[1][19]->Fill(ntp);
         if (ntp >= 2)
         {
            Hists[1][20]->Fill((*selectedTPs)[0]->pt);
            Hists[1][21]->Fill((*selectedTPs)[0]->eta);
            Hists[1][22]->Fill((*selectedTPs)[0]->phi);
            Hists[1][23]->Fill((*selectedTPs)[0]->z0);
            Hists[1][24]->Fill((*selectedTPs)[0]->d0);
            Hists[1][25]->Fill((*selectedTPs)[0]->rinv);
            Hists[1][26]->Fill((*selectedTPs)[0]->pdgid);
            Hists[1][27]->Fill((*selectedTPs)[1]->pt);
            Hists[1][28]->Fill((*selectedTPs)[1]->eta);
            Hists[1][29]->Fill((*selectedTPs)[1]->phi);
            Hists[1][30]->Fill((*selectedTPs)[1]->z0);
            Hists[1][31]->Fill((*selectedTPs)[1]->d0);
            Hists[1][32]->Fill((*selectedTPs)[1]->rinv);
            Hists[1][33]->Fill((*selectedTPs)[1]->pdgid);
         }
      }
      // Number of Tracks >= 2
      if (ntrk >= 2)
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[2][0]->Fill((*selectedTracks)[i]->pt);
            Hists[2][1]->Fill((*selectedTracks)[i]->eta);
            Hists[2][2]->Fill((*selectedTracks)[i]->phi);
            Hists[2][3]->Fill((*selectedTracks)[i]->z0);
            Hists[2][4]->Fill((*selectedTracks)[i]->d0);
            Hists[2][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[2][6]->Fill((*selectedTracks)[0]->pt);
         Hists[2][7]->Fill((*selectedTracks)[0]->eta);
         Hists[2][8]->Fill((*selectedTracks)[0]->phi);
         Hists[2][9]->Fill((*selectedTracks)[0]->z0);
         Hists[2][10]->Fill((*selectedTracks)[0]->d0);
         Hists[2][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[2][12]->Fill((*selectedTracks)[1]->pt);
         Hists[2][13]->Fill((*selectedTracks)[1]->eta);
         Hists[2][14]->Fill((*selectedTracks)[1]->phi);
         Hists[2][15]->Fill((*selectedTracks)[1]->z0);
         Hists[2][16]->Fill((*selectedTracks)[1]->d0);
         Hists[2][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[2][18]->Fill(ntrk);
         Hists[2][19]->Fill(ntp);
         if (ntp > 0)
         {
            Hists[2][20]->Fill((*selectedTPs)[0]->pt);
            Hists[2][21]->Fill((*selectedTPs)[0]->eta);
            Hists[2][22]->Fill((*selectedTPs)[0]->phi);
            Hists[2][23]->Fill((*selectedTPs)[0]->z0);
            Hists[2][24]->Fill((*selectedTPs)[0]->d0);
            Hists[2][25]->Fill((*selectedTPs)[0]->rinv);
            Hists[2][26]->Fill((*selectedTPs)[0]->pdgid);
         }
         if (ntp > 1)
         {
            Hists[2][27]->Fill((*selectedTPs)[1]->pt);
            Hists[2][28]->Fill((*selectedTPs)[1]->eta);
            Hists[2][29]->Fill((*selectedTPs)[1]->phi);
            Hists[2][30]->Fill((*selectedTPs)[1]->z0);
            Hists[2][31]->Fill((*selectedTPs)[1]->d0);
            Hists[2][32]->Fill((*selectedTPs)[1]->rinv);
            Hists[2][33]->Fill((*selectedTPs)[1]->pdgid);
         }
      }
      // Number of Tracks >= 2 and opposite sign highest pt tracks
      if (ntrk >= 2 && ((*selectedTracks)[0]->rinv * (*selectedTracks)[1]->rinv < 0))
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[3][0]->Fill((*selectedTracks)[i]->pt);
            Hists[3][1]->Fill((*selectedTracks)[i]->eta);
            Hists[3][2]->Fill((*selectedTracks)[i]->phi);
            Hists[3][3]->Fill((*selectedTracks)[i]->z0);
            Hists[3][4]->Fill((*selectedTracks)[i]->d0);
            Hists[3][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[3][6]->Fill((*selectedTracks)[0]->pt);
         Hists[3][7]->Fill((*selectedTracks)[0]->eta);
         Hists[3][8]->Fill((*selectedTracks)[0]->phi);
         Hists[3][9]->Fill((*selectedTracks)[0]->z0);
         Hists[3][10]->Fill((*selectedTracks)[0]->d0);
         Hists[3][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[3][12]->Fill((*selectedTracks)[1]->pt);
         Hists[3][13]->Fill((*selectedTracks)[1]->eta);
         Hists[3][14]->Fill((*selectedTracks)[1]->phi);
         Hists[3][15]->Fill((*selectedTracks)[1]->z0);
         Hists[3][16]->Fill((*selectedTracks)[1]->d0);
         Hists[3][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[3][18]->Fill(ntrk);
         Hists[3][19]->Fill(ntp);
         if (ntp > 0)
         {
            Hists[3][20]->Fill((*selectedTPs)[0]->pt);
            Hists[3][21]->Fill((*selectedTPs)[0]->eta);
            Hists[3][22]->Fill((*selectedTPs)[0]->phi);
            Hists[3][23]->Fill((*selectedTPs)[0]->z0);
            Hists[3][24]->Fill((*selectedTPs)[0]->d0);
            Hists[3][25]->Fill((*selectedTPs)[0]->rinv);
            Hists[3][26]->Fill((*selectedTPs)[0]->pdgid);
         }
         if (ntp > 1)
         {
            Hists[3][27]->Fill((*selectedTPs)[1]->pt);
            Hists[3][28]->Fill((*selectedTPs)[1]->eta);
            Hists[3][29]->Fill((*selectedTPs)[1]->phi);
            Hists[3][30]->Fill((*selectedTPs)[1]->z0);
            Hists[3][31]->Fill((*selectedTPs)[1]->d0);
            Hists[3][32]->Fill((*selectedTPs)[1]->rinv);
            Hists[3][33]->Fill((*selectedTPs)[1]->pdgid);
         }
      }
      // True Number of Tracks = 2
      if (ntp == 2 && ntrk >= 2)
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[4][0]->Fill((*selectedTracks)[i]->pt);
            Hists[4][1]->Fill((*selectedTracks)[i]->eta);
            Hists[4][2]->Fill((*selectedTracks)[i]->phi);
            Hists[4][3]->Fill((*selectedTracks)[i]->z0);
            Hists[4][4]->Fill((*selectedTracks)[i]->d0);
            Hists[4][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[4][6]->Fill((*selectedTracks)[0]->pt);
         Hists[4][7]->Fill((*selectedTracks)[0]->eta);
         Hists[4][8]->Fill((*selectedTracks)[0]->phi);
         Hists[4][9]->Fill((*selectedTracks)[0]->z0);
         Hists[4][10]->Fill((*selectedTracks)[0]->d0);
         Hists[4][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[4][12]->Fill((*selectedTracks)[1]->pt);
         Hists[4][13]->Fill((*selectedTracks)[1]->eta);
         Hists[4][14]->Fill((*selectedTracks)[1]->phi);
         Hists[4][15]->Fill((*selectedTracks)[1]->z0);
         Hists[4][16]->Fill((*selectedTracks)[1]->d0);
         Hists[4][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[4][18]->Fill(ntrk);
         Hists[4][19]->Fill(ntp);
         Hists[4][20]->Fill((*selectedTPs)[0]->pt);
         Hists[4][21]->Fill((*selectedTPs)[0]->eta);
         Hists[4][22]->Fill((*selectedTPs)[0]->phi);
         Hists[4][23]->Fill((*selectedTPs)[0]->z0);
         Hists[4][24]->Fill((*selectedTPs)[0]->d0);
         Hists[4][25]->Fill((*selectedTPs)[0]->rinv);
         Hists[4][26]->Fill((*selectedTPs)[0]->pdgid);
         Hists[4][27]->Fill((*selectedTPs)[1]->pt);
         Hists[4][28]->Fill((*selectedTPs)[1]->eta);
         Hists[4][29]->Fill((*selectedTPs)[1]->phi);
         Hists[4][30]->Fill((*selectedTPs)[1]->z0);
         Hists[4][31]->Fill((*selectedTPs)[1]->d0);
         Hists[4][32]->Fill((*selectedTPs)[1]->rinv);
         Hists[4][33]->Fill((*selectedTPs)[1]->pdgid);
      }
      // True Number of Tracks >= 2
      if (ntp >= 2 && ntrk >= 2)
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[5][0]->Fill((*selectedTracks)[i]->pt);
            Hists[5][1]->Fill((*selectedTracks)[i]->eta);
            Hists[5][2]->Fill((*selectedTracks)[i]->phi);
            Hists[5][3]->Fill((*selectedTracks)[i]->z0);
            Hists[5][4]->Fill((*selectedTracks)[i]->d0);
            Hists[5][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[5][6]->Fill((*selectedTracks)[0]->pt);
         Hists[5][7]->Fill((*selectedTracks)[0]->eta);
         Hists[5][8]->Fill((*selectedTracks)[0]->phi);
         Hists[5][9]->Fill((*selectedTracks)[0]->z0);
         Hists[5][10]->Fill((*selectedTracks)[0]->d0);
         Hists[5][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[5][12]->Fill((*selectedTracks)[1]->pt);
         Hists[5][13]->Fill((*selectedTracks)[1]->eta);
         Hists[5][14]->Fill((*selectedTracks)[1]->phi);
         Hists[5][15]->Fill((*selectedTracks)[1]->z0);
         Hists[5][16]->Fill((*selectedTracks)[1]->d0);
         Hists[5][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[5][18]->Fill(ntrk);
         Hists[5][19]->Fill(ntp);
         Hists[5][20]->Fill((*selectedTPs)[0]->pt);
         Hists[5][21]->Fill((*selectedTPs)[0]->eta);
         Hists[5][22]->Fill((*selectedTPs)[0]->phi);
         Hists[5][23]->Fill((*selectedTPs)[0]->z0);
         Hists[5][24]->Fill((*selectedTPs)[0]->d0);
         Hists[5][25]->Fill((*selectedTPs)[0]->rinv);
         Hists[5][26]->Fill((*selectedTPs)[0]->pdgid);
         Hists[5][27]->Fill((*selectedTPs)[1]->pt);
         Hists[5][28]->Fill((*selectedTPs)[1]->eta);
         Hists[5][29]->Fill((*selectedTPs)[1]->phi);
         Hists[5][30]->Fill((*selectedTPs)[1]->z0);
         Hists[5][31]->Fill((*selectedTPs)[1]->d0);
         Hists[5][32]->Fill((*selectedTPs)[1]->rinv);
         Hists[5][33]->Fill((*selectedTPs)[1]->pdgid);
      }
      // True Number of Tracks >= 2 and highest pt tracks are opposite sign muons
      if (ntrk >= 2 && ntp >= 2 && ((*selectedTPs)[0]->rinv * (*selectedTPs)[1]->rinv < 0) && (abs((*selectedTPs)[0]->pdgid) == 13) && (abs((*selectedTPs)[1]->pdgid) == 13))
      {
         for (int i = 0; i < ntrk; i++)
         {
            Hists[6][0]->Fill((*selectedTracks)[i]->pt);
            Hists[6][1]->Fill((*selectedTracks)[i]->eta);
            Hists[6][2]->Fill((*selectedTracks)[i]->phi);
            Hists[6][3]->Fill((*selectedTracks)[i]->z0);
            Hists[6][4]->Fill((*selectedTracks)[i]->d0);
            Hists[6][5]->Fill((*selectedTracks)[i]->rinv);
         }
         Hists[6][6]->Fill((*selectedTracks)[0]->pt);
         Hists[6][7]->Fill((*selectedTracks)[0]->eta);
         Hists[6][8]->Fill((*selectedTracks)[0]->phi);
         Hists[6][9]->Fill((*selectedTracks)[0]->z0);
         Hists[6][10]->Fill((*selectedTracks)[0]->d0);
         Hists[6][11]->Fill((*selectedTracks)[0]->rinv);
         Hists[6][12]->Fill((*selectedTracks)[1]->pt);
         Hists[6][13]->Fill((*selectedTracks)[1]->eta);
         Hists[6][14]->Fill((*selectedTracks)[1]->phi);
         Hists[6][15]->Fill((*selectedTracks)[1]->z0);
         Hists[6][16]->Fill((*selectedTracks)[1]->d0);
         Hists[6][17]->Fill((*selectedTracks)[1]->rinv);
         Hists[6][18]->Fill(ntrk);
         Hists[6][19]->Fill(ntp);
         Hists[6][20]->Fill((*selectedTPs)[0]->pt);
         Hists[6][21]->Fill((*selectedTPs)[0]->eta);
         Hists[6][22]->Fill((*selectedTPs)[0]->phi);
         Hists[6][23]->Fill((*selectedTPs)[0]->z0);
         Hists[6][24]->Fill((*selectedTPs)[0]->d0);
         Hists[6][25]->Fill((*selectedTPs)[0]->rinv);
         Hists[6][26]->Fill((*selectedTPs)[0]->pdgid);
         Hists[6][27]->Fill((*selectedTPs)[1]->pt);
         Hists[6][28]->Fill((*selectedTPs)[1]->eta);
         Hists[6][29]->Fill((*selectedTPs)[1]->phi);
         Hists[6][30]->Fill((*selectedTPs)[1]->z0);
         Hists[6][31]->Fill((*selectedTPs)[1]->d0);
         Hists[6][32]->Fill((*selectedTPs)[1]->rinv);
         Hists[6][33]->Fill((*selectedTPs)[1]->pdgid);
      }

      for (int l = 0; l < selectedTracks->size(); l++)
      {
         delete (*selectedTracks)[l];
      }
      for (int l = 0; l < selectedTPs->size(); l++)
      {
         delete (*selectedTPs)[l];
      }
      selectedTracks->clear();
      selectedTracks->shrink_to_fit();
      delete selectedTracks;
      selectedTPs->clear();
      selectedTPs->shrink_to_fit();
      delete selectedTPs;
      
   } // end event loop
   cout << endl
        << "from " << nevt << " events, " << nAccept << " events are accepted" << endl;

   // ---------------------------------------------------------------------------------------------------------
   //some Histograms

   if (TP_select_pdgid != 0)
   {
      char pdgidtxt[500];
      sprintf(pdgidtxt, "_pdgid%i", TP_select_pdgid);
      type = type + pdgidtxt;
   }
   else if (TP_select_injet == 1)
      type = type + "_injet";
   else if (TP_select_injet == 2)
      type = type + "_injet_highpt";
   else if (TP_select_injet == 3)
      type = type + "_injet_vhighpt";

   if (TP_minPt > 2.0)
   {
      char pttxt[500];
      sprintf(pttxt, "_pt%.0f", TP_minPt);
      type = type + pttxt;
   }

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
         Hists[k][l]->Write("", TObject::kOverwrite);
         Hists[k][l]->Draw();
         Hists[k][l]->GetYaxis()->SetTitle("Events");
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
         else if (vars[l].Contains("rinv"))
         {
            Hists[k][l]->GetXaxis()->SetTitle("rInv");
         }
         // Hists[k][l]->GetXaxis()->SetTitle(+Hists[k][l]->GetName());
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
   fout->Close();
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
   gStyle->SetOptTitle(1);
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
