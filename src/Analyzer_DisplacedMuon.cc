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

Double_t deltaPhi(Double_t phi1, Double_t phi2)
{
   Double_t dPhi = phi1 - phi2;
   if (dPhi > TMath::Pi())
      dPhi -= 2. * TMath::Pi();
   if (dPhi < -TMath::Pi())
      dPhi += 2. * TMath::Pi();
   return dPhi;
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
   float TP_maxPt = 100.0;
   float TP_maxEta = 2.4;
   int TP_select_pdgid = 0;

   // ----------------------------------------------------------------------------------------------------------------
   // histograms
   // ----------------------------------------------------------------------------------------------------------------

   typedef vector<TH1F *> Dim1;
   typedef vector<Dim1> Dim2;
   typedef vector<Dim2> Dim3;

   std::vector<TString> regions{"tpgeq2", "tpgeq2osMu", "tpgeq2osMuzmin", "tpgeq2osMudmin"};
   std::vector<TString> vars{"ntp",
                             "tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_d0", "tp_charge", "tp_pdgid",
                             "tp1pt", "tp1eta", "tp1phi", "tp1z0", "tp1d0", "tp1charge", "tp1pdgid",
                             "tp2pt", "tp2eta", "tp2phi", "tp2z0", "tp2d0", "tp2charge", "tp2pdgid",
                             "tp_tp_z0","tp_tp_d0","tp_tp_phi",
                             "tpz1pt", "tpz1eta", "tpz1phi", "tpz1z0", "tpz1d0", "tpz1charge", "tpz1pdgid",
                             "tpz2pt", "tpz2eta", "tpz2phi", "tpz2z0", "tpz2d0", "tpz2charge", "tpz2pdgid",
                             "tpz_tpz_z0","tpz_tpz_d0","tpz_tpz_phi",
                             "tpd1pt", "tpd1eta", "tpd1phi", "tpd1z0", "tpd1d0", "tpd1charge", "tpd1pdgid",
                             "tpd2pt", "tpd2eta", "tpd2phi", "tpd2z0", "tpd2d0", "tpd2charge", "tpd2pdgid",
                             "tpd_tpd_z0","tpd_tpd_d0","tpd_tpd_phi"
                             };
   std::vector<int> nbins{50,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50, 50, 50, 50, 400,
                          50, 50, 50};
   std::vector<float> lowEdge{0,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              -20, -10, -4,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              -20, -10, -4,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              0, -2.4, -3.2, -20, -10, -2, 0,
                              -20, -10, -4};
   std::vector<float> highEdge{50,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               20, 10, 4,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               20, 10, 4,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               50, 2.4, 3.2, 20, 10, 2, 400,
                               20, 10, 4};

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
   // std::vector<Track_Parameters *> *selectedTracks;
   // std::vector<Track_Parameters *> *selectedTracks_zmin;
   // std::vector<Track_Parameters *> *selectedTracks_dmin;
   std::vector<Track_Parameters *> *selectedTPs;
   std::vector<Track_Parameters *> *selectedTPs_zmin;
   std::vector<Track_Parameters *> *selectedTPs_dmin;
   int Track_mind0_index = -1;
   int Track_minz0_index = -1;
   int TP_mind0_index = -1;
   int TP_minz0_index = -1;
   Long64_t nevt = fChain->GetEntries();

   for (Long64_t i=0; i<nevt; i++) {
      fChain->GetEntry(i);
      // if (Cut(ientry) < 0) continue;

      displayProgress(i, nevt);

      // ----------------------------------------------------------------------------------------------------------------
      // track particle loop

      int ntrk = 0, ntp = 0;
      int h_index;
      // selectedTracks = new std::vector<Track_Parameters *>();
      // selectedTracks_zmin = new std::vector<Track_Parameters *>();
      // selectedTracks_dmin = new std::vector<Track_Parameters *>();

      selectedTPs = new std::vector<Track_Parameters *>();
      selectedTPs_zmin = new std::vector<Track_Parameters *>();
      selectedTPs_dmin = new std::vector<Track_Parameters *>();
      /*

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
         ntrk++;
         selectedTracks->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
         selectedTracks_zmin->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
         selectedTracks_dmin->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
      }
*/
      // ----------------------------------------------------------------------------------------------------------------
      // tracking particle loop

      for (int it = 0; it < (int)tp_pt->size(); it++)
      {
         if (std::fabs(tp_d0->at(it)) > TP_maxD0)
            continue;
         if (std::fabs(tp_d0->at(it)) < TP_minD0)
            continue;
         if (tp_pt->at(it) < 0.2)
            continue;
         if (tp_pt->at(it) > TP_maxPt)
            continue;
         if (std::fabs(tp_eta->at(it)) > TP_maxEta)
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
         ntp++;
         selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tp_d0->at(it), tp_z0->at(it), tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, tp_pdgid->at(it)));
         selectedTPs_zmin->push_back(new Track_Parameters(tp_pt->at(it), tp_d0->at(it), tp_z0->at(it), tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, tp_pdgid->at(it)));
         selectedTPs_dmin->push_back(new Track_Parameters(tp_pt->at(it), tp_d0->at(it), tp_z0->at(it), tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, tp_pdgid->at(it)));
      } // end tp loop

      // sort(selectedTracks->begin(),selectedTracks->end(),ComparePtTrack);
      sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);
      Track_minz0_index = -1;
      Track_mind0_index = -1;
/*
      if (selectedTracks->size() >= 2)
      {
         float Track_minz0 = std::numeric_limits<float>::infinity();
         float Track_mind0 = std::numeric_limits<float>::infinity();
         sort(selectedTracks_zmin->begin(), selectedTracks_zmin->end(), CompareZ0Track);
         sort(selectedTracks_dmin->begin(), selectedTracks_dmin->end(), CompareD0Track);
         for (int j = 0; j < selectedTracks->size() - 1; j++)
         {
            if (fabs((*selectedTracks_zmin)[j+1]->z0 - (*selectedTracks_zmin)[j]->z0) < Track_minz0){
               Track_minz0 = fabs((*selectedTracks_zmin)[j + 1]->z0 - (*selectedTracks_zmin)[j]->z0);
               Track_minz0_index = j;
            }
            if (fabs((*selectedTracks_dmin)[j + 1]->d0 - (*selectedTracks_dmin)[j]->d0) < Track_mind0)
            {
               Track_mind0 = fabs((*selectedTracks_dmin)[j + 1]->d0 - (*selectedTracks_zmin)[j]->d0);
               Track_mind0_index = j;
            }
         }
      }
*/
      TP_minz0_index = -1;
      TP_mind0_index = -1;

      if (selectedTPs->size() >= 2)
      {
         float TP_minz0 = std::numeric_limits<float>::infinity();
         float TP_mind0 = std::numeric_limits<float>::infinity();
         sort(selectedTPs_zmin->begin(), selectedTPs_zmin->end(), CompareZ0Track);
         sort(selectedTPs_dmin->begin(), selectedTPs_dmin->end(), CompareD0Track);
         for (int j = 0; j < selectedTPs->size() - 1; j++)
         {
            if (fabs((*selectedTPs_zmin)[j + 1]->z0 - (*selectedTPs_zmin)[j]->z0) < TP_minz0)
            {
               TP_minz0 = fabs((*selectedTPs_zmin)[j + 1]->z0 - (*selectedTPs_zmin)[j]->z0);
               TP_minz0_index = j;
               if ((*selectedTPs_zmin)[j + 1]->pt > (*selectedTPs_zmin)[j]->pt)
                  std::iter_swap((*selectedTPs_zmin)[j + 1], (*selectedTPs_zmin)[j]);
            }
            if (fabs((*selectedTPs_dmin)[j + 1]->d0 - (*selectedTPs_dmin)[j]->d0) < TP_mind0)
            {
               TP_mind0 = fabs((*selectedTPs_dmin)[j + 1]->d0 - (*selectedTPs_zmin)[j]->d0);
               TP_mind0_index = j;
               if ((*selectedTPs_dmin)[j + 1]->pt > (*selectedTPs_dmin)[j]->pt)
                  std::iter_swap((*selectedTPs_dmin)[j + 1], (*selectedTPs_dmin)[j]);
            }
         }
      }

      if(ntp<2) continue;
      // cout<<ntrk<<"\t"<<ntp<<endl;
      nAccept++;
         // ---------------------------------------------------------------------------------------------------------
         //Filling up Histograms
         
         h_index=0;
         Hists[0][h_index++]->Fill(ntp);
         for (int i = 0; i < ntp; i++)
         {
            Hists[0][h_index]->Fill((*selectedTPs)[i]->pt);
            Hists[0][h_index+1]->Fill((*selectedTPs)[i]->eta);
            Hists[0][h_index+2]->Fill((*selectedTPs)[i]->phi);
            Hists[0][h_index+3]->Fill((*selectedTPs)[i]->z0);
            Hists[0][h_index+4]->Fill((*selectedTPs)[i]->d0);
            Hists[0][h_index+5]->Fill((*selectedTPs)[i]->charge);
            Hists[0][h_index+6]->Fill((*selectedTPs)[i]->pdgid);
         }
         h_index +=7;
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs)[1]->pdgid);
         Hists[0][h_index++]->Fill(((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
         Hists[0][h_index++]->Fill(((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
         Hists[0][h_index++]->Fill((deltaPhi((*selectedTPs)[0]->phi, (*selectedTPs)[1]->phi)));
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
         Hists[0][h_index++]->Fill(((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
         Hists[0][h_index++]->Fill(((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
         Hists[0][h_index++]->Fill((deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
         Hists[0][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
         Hists[0][h_index++]->Fill(((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
         Hists[0][h_index++]->Fill(((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
         Hists[0][h_index++]->Fill((deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));

         // True Number of Tracks >= 2 and highest pt tracks are opposite sign muons
         if (ntp >= 2 && ((*selectedTPs)[0]->charge * (*selectedTPs)[1]->charge < 0) && (abs((*selectedTPs)[0]->pdgid) == 13) && (abs((*selectedTPs)[1]->pdgid) == 13))
         {
            h_index = 0;
            Hists[1][h_index++]->Fill(ntp);
            for (int i = 0; i < ntp; i++)
            {
               Hists[1][h_index]->Fill((*selectedTPs)[i]->pt);
               Hists[1][h_index + 1]->Fill((*selectedTPs)[i]->eta);
               Hists[1][h_index + 2]->Fill((*selectedTPs)[i]->phi);
               Hists[1][h_index + 3]->Fill((*selectedTPs)[i]->z0);
               Hists[1][h_index + 4]->Fill((*selectedTPs)[i]->d0);
               Hists[1][h_index + 5]->Fill((*selectedTPs)[i]->charge);
               Hists[1][h_index + 6]->Fill((*selectedTPs)[i]->pdgid);
            }
            h_index += 7;
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs)[1]->pdgid);
            Hists[1][h_index++]->Fill(((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
            Hists[1][h_index++]->Fill(((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
            Hists[1][h_index++]->Fill((deltaPhi((*selectedTPs)[0]->phi, (*selectedTPs)[1]->phi)));
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
            Hists[1][h_index++]->Fill(((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
            Hists[1][h_index++]->Fill(((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
            Hists[1][h_index++]->Fill((deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
            Hists[1][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
            Hists[1][h_index++]->Fill(((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
            Hists[1][h_index++]->Fill(((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
            Hists[1][h_index++]->Fill((deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));
         }

      // True Number of Tracks >= 2 and minimum z0 tracks are opposite sign muons
      if (ntp >= 2 && ((*selectedTPs_zmin)[0]->charge * (*selectedTPs_zmin)[1]->charge < 0) && (abs((*selectedTPs_zmin)[0]->pdgid) == 13) && (abs((*selectedTPs_zmin)[1]->pdgid) == 13))
      {
         h_index = 0;
         Hists[2][h_index++]->Fill(ntp);
         for (int i = 0; i < ntp; i++)
         {
            Hists[2][h_index]->Fill((*selectedTPs)[i]->pt);
            Hists[2][h_index + 1]->Fill((*selectedTPs)[i]->eta);
            Hists[2][h_index + 2]->Fill((*selectedTPs)[i]->phi);
            Hists[2][h_index + 3]->Fill((*selectedTPs)[i]->z0);
            Hists[2][h_index + 4]->Fill((*selectedTPs)[i]->d0);
            Hists[2][h_index + 5]->Fill((*selectedTPs)[i]->charge);
            Hists[2][h_index + 6]->Fill((*selectedTPs)[i]->pdgid);
         }
         h_index += 7;
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs)[0]->pdgid);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs)[1]->pdgid);
         Hists[2][h_index++]->Fill(((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
         Hists[2][h_index++]->Fill(((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
         Hists[2][h_index++]->Fill((deltaPhi((*selectedTPs)[0]->phi, (*selectedTPs)[1]->phi)));
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
         Hists[2][h_index++]->Fill(((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
         Hists[2][h_index++]->Fill(((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
         Hists[2][h_index++]->Fill((deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
         Hists[2][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
         Hists[2][h_index++]->Fill(((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
         Hists[2][h_index++]->Fill(((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
         Hists[2][h_index++]->Fill((deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));
      }

      // True Number of Tracks >= 2 and minimum d0 tracks are opposite sign muons
      if (ntp >= 2 && ((*selectedTPs_dmin)[0]->charge * (*selectedTPs_dmin)[1]->charge < 0) && (abs((*selectedTPs_dmin)[0]->pdgid) == 13) && (abs((*selectedTPs_dmin)[1]->pdgid) == 13))
      {
         h_index = 0;
         Hists[3][h_index++]->Fill(ntp);
         for (int i = 0; i < ntp; i++)
         {
            Hists[3][h_index]->Fill((*selectedTPs)[i]->pt);
            Hists[3][h_index + 1]->Fill((*selectedTPs)[i]->eta);
            Hists[3][h_index + 2]->Fill((*selectedTPs)[i]->phi);
            Hists[3][h_index + 3]->Fill((*selectedTPs)[i]->z0);
            Hists[3][h_index + 4]->Fill((*selectedTPs)[i]->d0);
            Hists[3][h_index + 5]->Fill((*selectedTPs)[i]->charge);
            Hists[3][h_index + 6]->Fill((*selectedTPs)[i]->pdgid);
         }
         h_index += 7;
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs)[0]->pdgid);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs)[1]->pdgid);
         Hists[3][h_index++]->Fill(((*selectedTPs)[0]->z0 - (*selectedTPs)[1]->z0));
         Hists[3][h_index++]->Fill(((*selectedTPs)[0]->d0 - (*selectedTPs)[1]->d0));
         Hists[3][h_index++]->Fill((deltaPhi((*selectedTPs)[0]->phi, (*selectedTPs)[1]->phi)));
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[0]->pdgid);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs_zmin)[1]->pdgid);
         Hists[3][h_index++]->Fill(((*selectedTPs_zmin)[0]->z0 - (*selectedTPs_zmin)[1]->z0));
         Hists[3][h_index++]->Fill(((*selectedTPs_zmin)[0]->d0 - (*selectedTPs_zmin)[1]->d0));
         Hists[3][h_index++]->Fill((deltaPhi((*selectedTPs_zmin)[0]->phi, (*selectedTPs_zmin)[1]->phi)));
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[0]->pdgid);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->pt);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->eta);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->phi);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->z0);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->d0);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->charge);
         Hists[3][h_index++]->Fill((*selectedTPs_dmin)[1]->pdgid);
         Hists[3][h_index++]->Fill(((*selectedTPs_dmin)[0]->z0 - (*selectedTPs_dmin)[1]->z0));
         Hists[3][h_index++]->Fill(((*selectedTPs_dmin)[0]->d0 - (*selectedTPs_dmin)[1]->d0));
         Hists[3][h_index++]->Fill((deltaPhi((*selectedTPs_dmin)[0]->phi, (*selectedTPs_dmin)[1]->phi)));
      }
/*
      for (int l = 0; l < selectedTracks->size(); l++)
      {
         delete (*selectedTracks)[l];
         delete (*selectedTracks_zmin)[l];
         delete (*selectedTracks_dmin)[l];
      }
*/
      for (int l = 0; l < selectedTPs->size(); l++)
      {
         delete (*selectedTPs)[l];
         delete (*selectedTPs_zmin)[l];
         delete (*selectedTPs_dmin)[l];
      }
/*
      selectedTracks->clear();
      selectedTracks->shrink_to_fit();
      delete selectedTracks;
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
      selectedTPs_zmin->clear();
      selectedTPs_zmin->shrink_to_fit();
      delete selectedTPs_zmin;
      selectedTPs_dmin->clear();
      selectedTPs_dmin->shrink_to_fit();
      delete selectedTPs_dmin;

   } // end event loop
   cout << endl
        << "from " << nevt << " events, " << nAccept << " events are accepted\t" << Hists[1][0]->GetEntries() << "  " << Hists[2][0]->GetEntries() << "  " << Hists[3][0]->GetEntries() << endl;

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
         Hists[k][l]->GetXaxis()->SetTitle(+Hists[k][l]->GetName());
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
