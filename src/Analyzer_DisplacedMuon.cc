// ----------------------------------------------------------------------------------------------------------------
// Feasibility study of using L1 Tracks to identify Displaced Vertex
//
// By Bharadwaj Harikrishnan, May 2021
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
#include "TMath.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

bool VERBOSE[]={false,false,false,false};

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

Double_t dist_TPs(Track_Parameters *a, Track_Parameters *b); // Closest distance between 2 tracks

Int_t Vertex(Track_Parameters *a, Track_Parameters*b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx); 
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
                            float TP_minD0 = 0.0,
                            int TP_select_injet = 0)
{
   gROOT->SetBatch();
   gErrorIgnoreLevel = kWarning;

   SetPlotStyle();

   float TP_minPt = 2.0;
   float TP_maxPt = 10000.0;
   float TP_maxEta = 2.4;

   TH1F *h_trk_d0 = new TH1F("h_trk_d0","h_trk_d0; Track d_{0} Distribution (cm) ; Events",200,-10,10);
   
   // Chi2 plots
   TH1F *h_trk_chi2rphidof = new TH1F("h_trk_chi2rphidof","h_trk_chi2rphidof; Track #chi^{2}_{r#phi}/d.o.f ; Events",50,0,5);
   TH1F *h_trk_chi2rphidof_8 = new TH1F("h_trk_chi2rphidof_8","h_trk_chi2rphidof_8; Track #chi^{2}_{r#phi}/d.o.f (p_{T}>8 GeV); Events",50,0,5);   
   TH1F *h_trk_chi2rzdof = new TH1F("h_trk_chi2rzdof","h_trk_chi2rzdof; Track #chi^{2}_{rz}/d.o.f ; Events",50,0,5);
   TH1F *h_trk_chi2rzdof_8 = new TH1F("h_trk_chi2rzdof_8","h_trk_chi2rzdof_8; Track #chi^{2}_{rz}/d.o.f (p_{T}>8 GeV); Events",50,0,5);
   TH1F *h_trk_bendchi2 = new TH1F("h_trk_bendchi2","h_trk_bendchi2; Track bend #chi^{2} ; Events",50,0,5);
   TH1F *h_trk_bendchi2_8 = new TH1F("h_trk_bendchi2_8","h_trk_bendchi2_8; Track bend #chi^{2} (p_{T}>8 GeV); Events",50,0,5);

   // Efficiency of Identifying Tracks Plots
   TH1F* h_tp_pt = new TH1F("h_tp_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_match_tp_pt = new TH1F("h_match_tp_pt", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100, 0, 100.0);
   TH1F* h_tp_eta = new TH1F("h_tp_eta", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_match_tp_eta = new TH1F("h_match_tp_eta", ";Tracking particle #eta; Tracking particles / 0.1", 50, -2.5, 2.5);
   TH1F* h_tp_d0 = new TH1F("h_tp_d0", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);
   TH1F* h_match_tp_d0 = new TH1F("h_match_tp_d0", ";Tracking particle d_{0} [cm]; Tracking particles / 0.01 cm", 50, -2, 2);   

   // Displaced Vertex Plots
   TH1F *h_delta_dist_xy = new TH1F("h_delta_dist_xy","h_delta_dist_xy; (a) Distance between chosen TPs in x-y (cm) ; Events",50,0,1);
   TH1F *h_error_delta_x = new TH1F("h_error_delta_x","h_error_delta_x; (b) #Delta x_{displaced} with chosen TP true x (cm) ; Events",100,0,2);
   TH1F *h_delta_dist_z = new TH1F("h_delta_dist_z","h_delta_dist_z; (c) Distance between chosen TPs in z (cm) ; Events",40,0,2);
   TH1F *h_error_delta_z = new TH1F("h_error_delta_z","h_error_delta_z; (d) Chosen TP error in z (cm) ; Events",100,0,10);
   TH1F *h_delta_x = new TH1F("h_delta_x","h_delta_x; #Delta x between chosen TPs true x (cm) ; Events",40,0,2);

   TH1I *h_trk_Counter_TPcombination = new TH1I("h_trk_Counter_TPcombination","h_trk_Counter_TPcombination; Track combination chosen; Events",6,0,6);
   TH1F *h_trk_delta_dist_xy = new TH1F("h_trk_delta_dist_xy","h_trk_delta_dist_xy; Distance between chosen Tracks in x-y (cm) ; Events",50,0,1);
   TH1F *h_trk_delta_dist_z = new TH1F("h_trk_delta_dist_z","h_trk_delta_dist_z; Distance between chosen Tracks in z (cm) ; Events",40,0,2);

   TH1F *h_res_tp_trk_x = new TH1F("h_res_tp_trk_x","h_res_tp_trk_x; x of vertex (cm) ; Events",100,-1,1);
   TH1F *h_res_tp_trk_y = new TH1F("h_res_tp_trk_y","h_res_tp_trk_y; y of vertex (cm) ; Events",100,-1,1);
   TH1F *h_res_tp_trk_z = new TH1F("h_res_tp_trk_z","h_res_tp_trk_z; z of vertex (cm) ; Events",1000,-10,10);

   // Efficiency of Identifying Displaced Vertex Plots
   TH1F *h_all_trk_pt = new TH1F("h_all_trk_pt","h_all_trk_pt; p_{T} of Leading p_{T} trk (GeV); Events",50,0,50);
   TH1F *h_correct_trk_pt = new TH1F("h_correct_trk_pt","h_correct_trk_pt; p_{T} of Leading p_{T} trk (GeV); Events",50,0,50);
   TH1F *h_all_trk_eta = new TH1F("h_all_trk_eta","h_all_trk_eta; #eta of Leading p_{T} trk ; Events",50,-2.4,2.4);
   TH1F *h_correct_trk_eta = new TH1F("h_correct_trk_eta","h_correct_trk_eta; #eta of Leading p_{T} trk ; Events",50,-2.4,2.4);
   TH1F *h_all_trk_dxy = new TH1F("h_all_trk_dxy","h_all_trk_dxy; dxy of vertex (cm); Events",20,0,2);
   TH1F *h_correct_trk_dxy = new TH1F("h_correct_trk_dxy","h_correct_trk_dxy; dxy of vertex (cm) ; Events",20,0,2);

   // Trigger Rates Study 
   TH2F *h_Count_trk_pt_d0 = new TH2F("h_Count_trk_pt_d0","h_Count_trk_pt_d0; Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks
   TH2F *h_Count_trk_pt_d0_dv = new TH2F("h_Count_trk_pt_d0_dv","h_Count_trk_pt_d0_dv;(DV selection) Transverse Impact Parameter d_{0} (cm); Transverse Momentum p_{T}(GeV)",5,0.1,0.6,5,2,12); // Count of Selected tracks including the displaced vertex selection
   

   float pt_cuts[5] = {2.0,4.0,6.0,8.0,10.0};     // Cuts to control event rate
   float d0_cuts[5] = {0.1,0.2,0.3,0.4,0.5};

   if (fChain == 0) return;

   std::vector<Track_Parameters *> *selectedTracks;      // Tracks 
   std::vector<Track_Parameters *> *selectedTracks_dmin; // Tracks from DV
   std::vector<Track_Parameters *> *selectedTPs;         // Tracking particles
   std::vector<Track_Parameters *> *selectedTPs_dmin;    // Tracking particles from DV

   Long64_t nevt = fChain->GetEntries();

   for (Long64_t i=0; i<nevt; i++) {
       fChain->GetEntry(i);
       displayProgress(i, nevt);
       
       selectedTracks = new std::vector<Track_Parameters *>();
       selectedTracks_dmin = new std::vector<Track_Parameters *>();
       selectedTPs = new std::vector<Track_Parameters *>();
       selectedTPs_dmin = new std::vector<Track_Parameters *>();

       // ----------------------------------------------------------------------------------------------------------------
       // track loop

	   for (int it = 0; it < (int)trk_pt->size(); it++){
			if (fabs(trk_eta->at(it)) > TP_maxEta)
		   		continue;
			if (trk_pt->at(it) < TP_minPt)
        		continue;
        	if (trk_pt->at(it) > TP_maxPt)
        		continue;
			if (std::fabs(trk_d0->at(it)) > TP_maxD0)
         	   continue;
		
			int ndof = 2 * trk_nstub->at(it) - 5;
			float chi2rphidof = (float)trk_chi2rphi->at(it) / ndof;
			float chi2rzdof = (float)trk_chi2rz->at(it) / ndof;
         
         h_trk_d0->Fill(trk_d0->at(it));
         h_trk_chi2rphidof->Fill(chi2rphidof);
         h_trk_chi2rzdof->Fill(chi2rzdof);
         h_trk_bendchi2->Fill(trk_bendchi2 ->at(it));

         if(trk_pt->at(it) > 8){   
            h_trk_chi2rphidof_8->Fill(chi2rphidof);
            h_trk_chi2rzdof_8->Fill(chi2rzdof);
            h_trk_bendchi2_8->Fill(trk_bendchi2 ->at(it));
         }

			if(chi2rphidof > 2)
				continue;
			if(chi2rzdof   > 2)
				continue;
			if(trk_bendchi2 ->at(it) > 5)
				continue;

			selectedTracks->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
	   }

		bool check_d0 = true;

    	for (int j = 0; j < selectedTracks->size(); j++){
        	if(fabs((*selectedTracks)[j]->d0) > TP_minD0){
        		check_d0=false;
            	break;
			}
      }
		
		if(check_d0) continue;

      int pt_check[5] = {0,0,0,0,0};
      int d0_check[5] = {0,0,0,0,0};
      int pt_d0_check[5][5];
      for (int k=0; k < 5; k++){
         for (int l=0; l < 5; l++){
            pt_d0_check[k][l] = 0;
         }
      }

      if ((selectedTracks->size() >= 2)){ 
      for (int j = 0; j < selectedTracks->size(); j++){
         // for (int k = j+1; k < selectedTracks->size(); k++){
            // if((fabs((*selectedTracks)[j]->z(x_dv_trk,y_dv_trk,z_dv_trk) - z_dv_trk)) < 1.0){
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

		for (int it = 0; it < (int)tp_pt->size(); it++){
      
         float tmp_d0 = -tp_d0->at(it);	// Sign difference in the NTupleMaker
         float tmp_z0 = tp_z0->at(it);
         if (std::fabs(tmp_d0) > TP_maxD0)
            continue;
			if (tp_pt->at(it) < TP_minPt)
				continue;
			if (tp_pt->at(it) > TP_maxPt)
				continue;
			if (std::fabs(tp_eta->at(it)) > TP_maxEta)
				continue;

         h_tp_pt->Fill(tp_pt->at(it));
         h_tp_eta->Fill(tp_eta->at(it));
         h_tp_d0->Fill(tp_d0->at(it));

         if (tp_nmatch->at(it) >= 1){
            h_match_tp_pt->Fill(tp_pt->at(it));
            h_match_tp_eta->Fill(tp_eta->at(it));
            h_match_tp_d0->Fill(tp_d0->at(it));
         }
			
			selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tmp_d0, tmp_z0, tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, abs(tp_pdgid->at(it))));
		}

      // --------------------------------------------------------------------------------------------
      //         Vertex finding in Tracking Particles
      // --------------------------------------------------------------------------------------------
      if (!(selectedTracks->size() >= 2)) continue;
      bool true_DV = false;

      if(VERBOSE[0])
         std::cout<<"End of z-vertex Finding"<<endl;
      
      
      Double_t x_dv = -9999.0;// (tp_x->at((*selectedTPs)[0]->index));//+tp_x->at((*selectedTPs)[1]->index))/2.0;
      Double_t y_dv = -9999.0;// (tp_y->at((*selectedTPs)[0]->index));//+tp_y->at((*selectedTPs)[1]->index))/2.0;
      Double_t z_dv = -9999.0;// (tp_z->at((*selectedTPs)[0]->index));//+tp_z->at((*selectedTPs)[1]->index))/2.0;
            
      int Vertex_check = -1;
      int selected_track_j0 = -1;
      int selected_track_j1 = -1;

      Double_t x_tmp = x_dv;
      Double_t y_tmp = y_dv;
      Double_t z_tmp = z_dv;
      Int_t Vertex_check_tmp = -1;

      if(selectedTPs->size()>=2){

         sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);

         // Selection of leading 2 p_T tracks from true common vertex
         // if(!(!((fabs(tp_x->at((*selectedTPs)[0]->index) - tp_x->at((*selectedTPs)[1]->index))<0.01) && (fabs(tp_y->at((*selectedTPs)[0]->index) - tp_y->at((*selectedTPs)[1]->index))<0.01)))){nCount_leadingpt++;}// continue;

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
         
         h_delta_dist_xy->Fill(dist_TPs((*selectedTPs)[selected_track_j0],(*selectedTPs)[selected_track_j1]));
         h_error_delta_x->Fill(fabs((tp_x->at((*selectedTPs)[selected_track_j0]->index) - x_dv)));
         
         float z_dv_1 = (*selectedTPs)[selected_track_j0]->z(x_dv,y_dv);
         float z_dv_2 = (*selectedTPs)[selected_track_j1]->z(x_dv,y_dv);

         h_delta_dist_z->Fill(fabs(z_dv_1-z_dv_2));
         h_error_delta_z->Fill(fabs((tp_z->at((*selectedTPs)[selected_track_j0]->index) - z_dv_1)));
         if((fabs(x_dv-(tp_x->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF) && (fabs(z_dv-(tp_z->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF) && (fabs(z_dv-tp_z->at((*selectedTPs)[selected_track_j0]->index))<CUTOFF)){
            h_delta_x->Fill(fabs((tp_x->at((*selectedTPs)[selected_track_j0]->index) - (tp_x->at((*selectedTPs)[selected_track_j1]->index)))));
            true_DV = true;
         }
      }

      // --------------------------------------------------------------------------------------------
      //                Vertex finding in Tracks
      // --------------------------------------------------------------------------------------------

      sort(selectedTracks->begin(), selectedTracks->end(), ComparePtTrack);
      Double_t x_dv_trk = -9999.0;// (tp_x->at((*selectedTracks)[0]->index));//+tp_x->at((*selectedTracks)[1]->index))/2.0;
      Double_t y_dv_trk = -9999.0;// (tp_y->at((*selectedTracks)[0]->index));//+tp_y->at((*selectedTracks)[1]->index))/2.0;
      Double_t z_dv_trk = -9999.0;// (tp_z->at((*selectedTracks)[0]->index));//+tp_z->at((*selectedTracks)[1]->index))/2.0;
            
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
               h_trk_Counter_TPcombination->AddBinContent(3); 
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
      
      h_trk_delta_dist_xy->Fill(dist_TPs((*selectedTracks)[selected_track_j0],(*selectedTracks)[selected_track_j1]));
      
      float z_dv_trk_1 = (*selectedTracks)[selected_track_j0]->z(x_dv_trk,y_dv_trk);
      float z_dv_trk_2 = (*selectedTracks)[selected_track_j1]->z(x_dv_trk,y_dv_trk);

      h_trk_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2)); //* Use this condition
      
      if(true_DV){

      h_all_trk_pt->Fill((*selectedTracks)[selected_track_j0]->pt);
      h_all_trk_eta->Fill((*selectedTracks)[selected_track_j0]->eta);
      h_all_trk_dxy->Fill(dist(x_dv_trk,y_dv_trk));
      }

      // if((dist_Vertex(x_dv_trk,y_dv_trk,(*selectedTracks)[selected_track_j0])<CUTOFF) && (fabs(z_dv_trk_1-z_dv_trk)<CUTOFF)){
      if(true_DV && (fabs(x_dv-x_dv_trk)<CUTOFF) && (fabs(y_dv-y_dv_trk)<CUTOFF) && (fabs(z_dv-z_dv_trk)<CUTOFF)){
         h_correct_trk_pt->Fill((*selectedTracks)[selected_track_j0]->pt);
         h_correct_trk_eta->Fill((*selectedTracks)[selected_track_j0]->eta);
         h_correct_trk_dxy->Fill(dist(x_dv_trk,y_dv_trk));
            
         h_res_tp_trk_x->Fill(x_dv-x_dv_trk);
         h_res_tp_trk_y->Fill(y_dv-y_dv_trk);
         h_res_tp_trk_z->Fill(z_dv-z_dv_trk);

         for (int j = 0; j < selectedTracks->size(); j++){
            // for (int k = j+1; k < selectedTracks->size(); k++){
               // if((fabs((*selectedTracks)[j]->z(x_dv_trk,y_dv_trk,z_dv_trk) - z_dv_trk)) < 1.0){
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
               // }
         // }
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

   h_trk_chi2rphidof_8->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rphidof_8->GetXaxis()->SetRange(1, h_trk_chi2rphidof_8->GetNbinsX() + 2);
   h_trk_chi2rphidof_8->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rphidof_8->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rphidof_8->GetName() + ".jpeg");
   delete h_trk_chi2rphidof_8;

   h_trk_chi2rzdof_8->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_chi2rzdof_8->GetXaxis()->SetRange(1, h_trk_chi2rzdof_8->GetNbinsX() + 2);
   h_trk_chi2rzdof_8->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_chi2rzdof_8->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_chi2rzdof_8->GetName() + ".jpeg");
   delete h_trk_chi2rzdof_8;

   h_trk_bendchi2_8->GetYaxis()->SetNoExponent(kTRUE);
   h_trk_bendchi2_8->GetXaxis()->SetRange(1, h_trk_bendchi2_8->GetNbinsX() + 2);
   h_trk_bendchi2_8->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_trk_bendchi2_8->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_trk_bendchi2_8->GetName() + ".jpeg");
   delete h_trk_bendchi2_8;

   h_match_tp_pt->Sumw2();
   h_tp_pt->Sumw2();
   TH1F* h_eff_pt = (TH1F*)h_match_tp_pt->Clone();
   h_eff_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_pt->GetXaxis()->SetRange(1, h_eff_pt->GetNbinsX() + 2);
   h_eff_pt->SetName("eff_pt");
   h_eff_pt->GetXaxis()->SetTitle("Tracking Particle p_{T} (GeV)");
   h_eff_pt->GetYaxis()->SetTitle("Efficiency");
   h_eff_pt->Divide(h_match_tp_pt, h_tp_pt, 1.0, 1.0, "B");
   h_eff_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_pt->GetName() + ".jpeg");
   delete h_eff_pt;
   delete h_match_tp_pt;
   delete h_tp_pt;

   h_match_tp_eta->Sumw2();
   h_tp_eta->Sumw2();
   TH1F* h_eff_eta = (TH1F*)h_match_tp_eta->Clone();
   h_eff_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_eta->GetXaxis()->SetRange(1, h_eff_eta->GetNbinsX() + 2);
   h_eff_eta->SetName("eff_eta");
   h_eff_eta->GetXaxis()->SetTitle("Tracking Particle #eta");
   h_eff_eta->GetYaxis()->SetTitle("Efficiency");
   h_eff_eta->Divide(h_match_tp_eta, h_tp_eta, 1.0, 1.0, "B");
   h_eff_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_eta->GetName() + ".jpeg");
   delete h_eff_eta;
   delete h_match_tp_eta;
   delete h_tp_eta;

   h_match_tp_d0->Sumw2();
   h_tp_d0->Sumw2();
   TH1F* h_eff_d0 = (TH1F*)h_match_tp_d0->Clone();
   h_eff_d0->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_d0->GetXaxis()->SetRange(1, h_eff_d0->GetNbinsX() + 2);
   h_eff_d0->SetName("eff_d0");
   h_eff_d0->GetXaxis()->SetTitle("Tracking Particle d_{0} (cm)");
   h_eff_d0->GetYaxis()->SetTitle("Efficiency");
   h_eff_d0->Divide(h_match_tp_d0, h_tp_d0, 1.0, 1.0, "B");
   h_eff_d0->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_d0->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_d0->GetName() + ".jpeg");
   delete h_eff_d0;
   delete h_match_tp_d0;
   delete h_tp_d0;

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

   h_all_trk_pt->GetXaxis()->SetRange(1, h_all_trk_pt->GetNbinsX() + 2);
   h_all_trk_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trk_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trk_pt->GetName() + ".jpeg");

   h_correct_trk_pt->GetXaxis()->SetRange(1, h_correct_trk_pt->GetNbinsX() + 2);
   h_correct_trk_pt->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trk_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trk_pt->GetName() + ".jpeg");
   
   h_all_trk_pt->Sumw2();
   h_correct_trk_pt->Sumw2();
   TH1F* h_eff_trk_pt = (TH1F*)h_correct_trk_pt->Clone();
   h_eff_trk_pt->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trk_pt->GetXaxis()->SetRange(1, h_eff_trk_pt->GetNbinsX() + 2);
   h_eff_trk_pt->SetName("eff_trk_pt");
   h_eff_trk_pt->GetXaxis()->SetTitle("Track p_{T} (GeV)");
   h_eff_trk_pt->GetYaxis()->SetTitle("Efficiency");
   h_eff_trk_pt->Divide(h_correct_trk_pt,h_all_trk_pt, 1.0, 1.0, "B");
   h_eff_trk_pt->Draw();
   h_eff_trk_pt->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trk_pt->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trk_pt->GetName() + ".jpeg");
   delete h_eff_trk_pt;
   delete h_all_trk_pt;
   delete h_correct_trk_pt;

   h_all_trk_eta->GetXaxis()->SetRange(1, h_all_trk_eta->GetNbinsX() + 2);
   h_all_trk_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trk_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trk_eta->GetName() + ".jpeg");

   h_correct_trk_eta->GetXaxis()->SetRange(1, h_correct_trk_eta->GetNbinsX() + 2);
   h_correct_trk_eta->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trk_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trk_eta->GetName() + ".jpeg");

   h_all_trk_eta->Sumw2();
   h_correct_trk_eta->Sumw2();
   TH1F* h_eff_trk_eta = (TH1F*)h_correct_trk_eta->Clone();
   h_eff_trk_eta->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trk_eta->GetXaxis()->SetRange(1, h_eff_trk_eta->GetNbinsX() + 2);
   h_eff_trk_eta->SetName("eff_trk_eta");
   h_eff_trk_eta->GetXaxis()->SetTitle("Track #eta");
   h_eff_trk_eta->GetYaxis()->SetTitle("Efficiency");
   h_eff_trk_eta->Divide(h_correct_trk_eta,h_all_trk_eta, 1.0, 1.0, "B");
   h_eff_trk_eta->Draw();
   h_eff_trk_eta->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trk_eta->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trk_eta->GetName() + ".jpeg");
   delete h_eff_trk_eta;
   delete h_all_trk_eta;
   delete h_correct_trk_eta;

   h_all_trk_dxy->GetXaxis()->SetRange(1, h_all_trk_dxy->GetNbinsX() + 2);
   h_all_trk_dxy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_all_trk_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_all_trk_dxy->GetName() + ".jpeg");

   h_correct_trk_dxy->GetXaxis()->SetRange(1, h_correct_trk_dxy->GetNbinsX() + 2);
   h_correct_trk_dxy->Draw();
   mySmallText(0.4, 0.82, 1, ctxt);
   h_correct_trk_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_correct_trk_dxy->GetName() + ".jpeg");
   
   h_all_trk_dxy->Sumw2();
   h_correct_trk_dxy->Sumw2();
   TH1F* h_eff_trk_dxy = (TH1F*)h_correct_trk_dxy->Clone();
   h_eff_trk_dxy->GetYaxis()->SetNoExponent(kTRUE);
   h_eff_trk_dxy->GetXaxis()->SetRange(1, h_eff_trk_dxy->GetNbinsX() + 2);
   h_eff_trk_dxy->SetName("eff_trk_dxy");
   h_eff_trk_dxy->GetXaxis()->SetTitle("Track d_{xy} (cm)");
   h_eff_trk_dxy->GetYaxis()->SetTitle("Efficiency");
   h_eff_trk_dxy->Divide(h_correct_trk_dxy,h_all_trk_dxy, 1.0, 1.0, "B");
   h_eff_trk_dxy->Draw();
   h_eff_trk_dxy->SetAxisRange(0, 1.1, "Y");
   mySmallText(0.4, 0.82, 1, ctxt);
   h_eff_trk_dxy->Write("", TObject::kOverwrite);
   c.SaveAs(DIR + "/"+ h_eff_trk_dxy->GetName() + ".jpeg");
   delete h_eff_trk_dxy;
   delete h_all_trk_dxy;
   delete h_correct_trk_dxy;

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


Double_t dist_TPs(Track_Parameters *a, Track_Parameters *b){  
   float x1 = a->x0; //   Centers of the circles
   float y1 = a->y0; // 
   float x2 = b->x0; // 
   float y2 = b->y0; // 
   float R1 = a->rho;   // Radii of the circles
   float R2 = b->rho;
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

Int_t Vertex(Track_Parameters *a, Track_Parameters*b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx){
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
      
      
      float ap1 = a->phi_T(ix1,iy1);
      float bp1 = b->phi_T(ix1,iy1);
      float ap2 = a->phi_T(ix2,iy2);
      float bp2 = b->phi_T(ix2,iy2);/*
      float z11 = a->z(ap1);
      float z21 = b->z(bp1);
      float z12 = a->z(ap2);
      float z22 = b->z(bp2);
      */
      float z11 = a->z(ix1,iy1);
      float z21 = b->z(ix1,iy1);
      float z12 = a->z(ix2,iy2);
      float z22 = b->z(ix2,iy2);
      
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
   float az = a->z(a->phi_T(x_vtx,y_vtx));
   float bz = b->z(b->phi_T(x_vtx,y_vtx));
   z_vtx = az;//(az+bz)/2.0;
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
