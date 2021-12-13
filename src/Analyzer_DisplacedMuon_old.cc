// #define Analyzer_DisplacedMuon_cxx
// #include "Analyzer_DisplacedMuon.h"
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TSystem.h>

// #include "TLatex.h"
// #include "TFile.h"
// #include "TTree.h"
// #include "TChain.h"
// #include "TBranch.h"
// #include "TLeaf.h"
// #include "TCanvas.h"
// #include "TLegend.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TF1.h"
// #include "TMath.h"
// #include <sstream>
// #include <fstream>
// #include <iostream>
// #include <string>
// #include <vector>
// #include <algorithm>

// using namespace std;

// bool VERBOSE[]={false,false,false,false};

// float CUTOFF = 1.0; // Common z cut off in cm for Displaced Vertex and verification

// void SetPlotStyle();
// void mySmallText(Double_t x, Double_t y, Color_t color, char *text);

// bool ComparePtTrack(Track_Parameters *a, Track_Parameters *b) { return a->pt > b->pt; }
// bool CompareZ0Track(Track_Parameters *a, Track_Parameters *b) { return a->z0 > b->z0; }
// bool CompareD0Track(Track_Parameters *a, Track_Parameters *b) { return a->d0 > b->d0; }

// Double_t dist(Double_t x1, Double_t y1 , Double_t x2=0, Double_t y2=0){
//    return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)));
// }
// Double_t dist3(Double_t x1, Double_t y1 , Double_t z1, Double_t x2=0, Double_t y2=0 , Double_t z2=0){
//    return (TMath::Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)));
// }

// Double_t dist_TPs(Track_Parameters *a, Track_Parameters *b){
//    float x1 = a->x0; //   Centers of the circles
//    float y1 = a->y0; // 
//    float x2 = b->x0; // 
//    float y2 = b->y0; // 
//    float R1 = a->rho;   // Radii of the circles
//    float R2 = b->rho;
//    float R = dist(x1,y1,x2,y2); // Distance between centers
//    if((R>=(R1-R2)) && (R<=(R1+R2))){
//       return (0);
//    }
//    else if(R==0){
//       return (-99999.0);
//    }
//    else{

//       return(R-R1-R2);
//    }

// }

// //! Remove 3D functions
// Double_t dist_Vertex(Double_t x_vtx, Double_t y_vtx, Double_t z_vtx, Track_Parameters *a){
//    float R = fabs(dist(x_vtx,y_vtx,a->x0,a->y0) -(a->rho));
//    return TMath::Sqrt(R*R + (z_vtx - a->z0)*(z_vtx - a->z0));
// }

// Double_t dist_Vertex(Double_t x_vtx, Double_t y_vtx, Track_Parameters *a){
//    float R = dist(x_vtx,y_vtx,a->x0,a->y0);
//    return (fabs(R-(a->rho)));
// }

// Int_t Vertex(Track_Parameters *a, Track_Parameters*b, Double_t &x_vtx, Double_t &y_vtx, Double_t &z_vtx){
//    float x1 = a->x0; //   Centers of the circles
//    float y1 = a->y0; // 
//    float x2 = b->x0; // 
//    float y2 = b->y0; // 
//    float R1 = a->rho;   // Radii of the circles
//    float R2 = b->rho;
//    float del1, del2, x_11, x_12, x_21, x_22, y_11, y_12, y_21, y_22;
//    float R = dist(x1,y1,x2,y2); // Distance between centers
//    float centerdx = x1 - x2;
//    float centerdy = y1 - y2;
//    int retval = -1;

//    if((R>=(R1-R2)) && (R<=(R1+R2))){
//       // Circles Intersect
//       float R4 = R*R*R*R;
//       float A = (R1*R1 - R2*R2) / (2 * R*R);
//       float r2r2 = (R1*R1 - R2*R2);
//       float C = TMath::Sqrt(2 * (R1*R1 + R2*R2) / (R*R) - (r2r2 * r2r2) / R4 - 1);
      
//       float fx = (x1+x2) / 2 + A * (x2 - x1);
//       float gx = C * (y2 - y1) / 2;
//       float ix1 = fx + gx;
//       float ix2 = fx - gx;

//       float fy = (y1+y2) / 2 + A * (y2 - y1);
//       float gy = C * (x1 - x2) / 2;
//       float iy1 = fy + gy;
//       float iy2 = fy - gy;
      
      
//       float ap1 = a->phi_T(ix1,iy1);
//       float bp1 = b->phi_T(ix1,iy1);
//       float ap2 = a->phi_T(ix2,iy2);
//       float bp2 = b->phi_T(ix2,iy2);/*
//       float z11 = a->z(ap1);
//       float z21 = b->z(bp1);
//       float z12 = a->z(ap2);
//       float z22 = b->z(bp2);
//       */
//       float z11 = a->z(ix1,iy1);
//       float z21 = b->z(ix1,iy1);
//       float z12 = a->z(ix2,iy2);
//       float z22 = b->z(ix2,iy2);
      
//       float delz1 = fabs(z11 - z21);//fabs(fabs(ap1 - bp1)-TMath::Pi());// fabs(dist3(ix1,iy1,z11)-dist3(ix1,iy1,z21));
//       float delz2 = fabs(z12 - z22);//fabs(fabs(ap2 - bp2)-TMath::Pi());// fabs(dist3(ix2,iy2,z12)-dist3(ix2,iy2,z22));
      
//       if(VERBOSE[3]){// &&(fabs(z11-z_vtx)>0.1 && fabs(z12-z_vtx)>0.1 && fabs(z21-z_vtx)>0.1 && fabs(z22-z_vtx)>0.1)){
//          std::cout<<Form("ix1 = %5.2f    |   iy1 = %5.2f    |  ix2 = %5.2f    |   iy2 = %5.2f",ix1,iy1,ix2,iy2)<<endl;
//          // std::cout<<Form("ap1 = %5.2f    |   bp1 = %5.2f    |  ap2 = %5.2f    |   bp2 = %5.2f",ap1,bp1,ap2,bp2)<<endl;
//          std::cout<<Form("z11 = %5.2f    |   z21 = %5.2f    |  z12 = %5.2f    |   z22 = %5.2f    |   delz1 = %5.2f    |   delz2 = %5.2f",z11,z21,z12,z22,delz1,delz2)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
//       }
//       retval = 0;
//       if(gx==0 && gy==0){ // Only 1 intersection
//          x_vtx = ix1;
//          y_vtx = iy1;
//          if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
//             retval = 1;
//          }
//          z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
//          // cout<<"----1 intersection ----";
//          return retval;//true;
//       }
      
//       if(dist(ix1,iy1)>20){
//          x_vtx = ix2;
//          y_vtx = iy2;
//          if(delz2<CUTOFF){//(fabs(z12-tp_z)<1.0 || fabs(z22-tp_z)<1.0){//
//             retval = 2;
//          }
//          z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;//fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;// dist3(ix2,iy2,z12)<dist3(ix2,iy2,z22)?z12:z22;//fabs(a->z0-z12)<fabs(b->z0-z22)?z12:z22;//z12;
//          return retval;//true;
//       }
//       if(dist(ix2,iy2)>20){
//          x_vtx = ix1;
//          y_vtx = iy1;
//          if(delz1<CUTOFF){//(fabs(z11-tp_z)<1.0 || fabs(z21-tp_z)<1.0){//
//             retval = 2;
//          }
//          z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;//fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;//dist3(ix1,iy1,z11)<dist3(ix1,iy1,z21)?z11:z21;//fabs(a->z0-z11)<fabs(b->z0-z21)?z11:z21;//z11;
//          return retval;//true;
//       }
      
//       // Finding minimum z distance to decide between the 2 intersection vertex

//       // if((fabs(z11-tp_z)<1.0 || fabs(z12-tp_z)<1.0 || fabs(z21-tp_z)<1.0 || fabs(z22-tp_z)>1.0)){
//           retval = 3;
//       // }
 
//       if (delz1<=delz2){
//          x_vtx = ix1;
//          y_vtx = iy1;
//          z_vtx = z11;//fabs(z11)<fabs(z21)?z11:z21;// fabs(a->z0-z11)<fabs(a->z0-z21)?z11:z21;
//          // cout<<"----2 intersection ----";
//          return retval;//true;         
//       }
//       else{
//          x_vtx = ix2;
//          y_vtx = iy2;
//          z_vtx = z12;//fabs(z12)<fabs(z22)?z12:z22;// fabs(a->z0-z12)<fabs(a->z0-z22)?z12:z22;
//          return retval;//true;
//       }
      

//    }
//    else if(R==0){
//       // Circles with Same center. 
//       x_vtx = -999.0;//a->x();  // Considering First track coordinates as vertex
//       y_vtx = -999.0;//a->y();
//       return -2;
//    }/*
//    x_vtx = -99999.0;
//    y_vtx = -99999.0;
//    z_vtx = -99999.0;
//    return false;
//    else{
//       x_vtx = -9999.0;
//       y_vtx = -9999.0;
//       return;
//    }
//    */
//    //* Circles don't intersect. 

//    if(x1==x2){
//       if(y1==y2){
//          // Circles with Same center. 
//          //x_vtx = a->x();  // Considering First track coordinates as vertex
//          // y_vtx = a->y();
//          return -2;
//       }
//       x_11 = x1;
//       x_12 = x1;

//       x_21 = x2;
//       x_22 = x2;

//       y_11 = y1 + R1;
//       y_12 = y1 - R1;

//       y_21 = y2 + R2;
//       y_22 = y2 - R2;

//    }
//    else{

//       del1 = R1 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));
//       del2 = R2 / (TMath::Sqrt(1+((y1-y2)*(y1-y2)/((x1-x2)*(x1-x2)))));

//       x_11 = x1 + del1;
//       x_12 = x1 - del1;

//       x_21 = x2 + del2;
//       x_22 = x2 - del2;

//       y_11 = y1 + (x_11-x1) * ((y1-y2)/(x1-x2));
//       y_12 = y1 + (x_12-x1) * ((y1-y2)/(x1-x2));

//       y_21 = y2 + (x_21-x2) * ((y1-y2)/(x1-x2));
//       y_22 = y2 + (x_22-x2) * ((y1-y2)/(x1-x2));
//    }

//    if(dist(x_11,y_11,x2,y2) <= dist(x_12,y_12,x2,y2)){
//       x_vtx = x_11;
//       y_vtx = y_11;
//    }
//    else{
//       x_vtx = x_12;
//       y_vtx = y_12;
//    }

//    if(dist(x_21,y_21,x1,y1) <= dist(x_22,y_22,x1,y1)){
//       x_vtx = (x_vtx + x_21)/2;
//       y_vtx = (y_vtx + y_21)/2;
//    }
//    else{
//       x_vtx = (x_vtx + x_22)/2;
//       y_vtx = (y_vtx + y_22)/2;
//    }
//    //! Does this make sense for no intersection case???
//    float az = a->z(a->phi_T(x_vtx,y_vtx));
//    float bz = b->z(b->phi_T(x_vtx,y_vtx));
//    z_vtx = az;//(az+bz)/2.0;
//    if(VERBOSE[3]){
//       std::cout<<Form("z_vtx = %5.1f  |  az = %5.1f  |  bz = %5.1f",z_vtx,az,bz)<<endl;
//    }
//    return -1;
// }

// Double_t deltaPhi(Double_t phi1, Double_t phi2)
// {
//    Double_t dPhi = phi1 - phi2;
//    if (dPhi > TMath::Pi())
//       dPhi -= 2. * TMath::Pi();
//    if (dPhi < -TMath::Pi())
//       dPhi += 2. * TMath::Pi();
//    return dPhi;
// }

// Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
// {
//    Double_t dEta, dPhi;
//    dEta = eta1 - eta2;
//    dPhi = deltaPhi(phi1, phi2);
//    return sqrt(dEta * dEta + dPhi * dPhi);
// }

// void displayProgress(long current, long max)
// {
//    using std::cerr;
//    if (max < 2500)
//       return;
//    if (current % (max / 2500) != 0 && current < max - 1)
//       return;

//    int width = 52; // Hope the terminal is at least that wide.
//    int barWidth = width - 2;
//    cerr << "\x1B[2K";    // Clear line
//    cerr << "\x1B[2000D"; // Cursor left
//    cerr << '[';
//    for (int i = 0; i < barWidth; ++i)
//    {
//       if (i < barWidth * current / max)
//       {
//          cerr << '=';
//       }
//       else
//       {
//          cerr << ' ';
//       }
//    }
//    cerr << ']';
//    cerr << " " << Form("%8d/%8d (%5.2f%%)", (int)current, (int)max, 100.0 * current / max);
//    cerr.flush();
// }

// void Analyzer_DisplacedMuon::Loop(TString type,
//                             TString type_dir = "",
//                             float TP_maxD0 = 10.0,
//                             float TP_minD0 = 0.0,
//                             int TP_select_injet = 0)
// {
//    gROOT->SetBatch();
//    gErrorIgnoreLevel = kWarning;

//    SetPlotStyle();

//    float TP_minPt = 2.0;
//    float TP_maxPt = 10000.0;
//    float TP_maxEta = 2.4;
//    int TP_select_pdgid = 0;

//    // ----------------------------------------------------------------------------------------------------------------
//    // histograms
//    // ----------------------------------------------------------------------------------------------------------------

//    std::vector<TString> regions{"reg"};
//    std::vector<TString> vars{"ntrk",
//                              "trk_pt", "trk_eta", "trk_phi", "trk_z0", "trk_d0", "trk_charge", "trk_pdgid","trk_dxy",
//                              "trk1pt", "trk1eta", "trk1phi", "trk1z0", "trk1d0", "trk1charge", "trk1pdgid","trk1_dxy",
//                              "trk2pt", "trk2eta", "trk2phi", "trk2z0", "trk2d0", "trk2charge", "trk2pdgid","trk2_dxy",
//                              "trk_trk_eta", "trk_trk_phi", "trk_trk_z0", "trk_trk_d0", "trk_trk_dR","trk_trk_dist",
//                              "trkd1pt", "trkd1eta", "trkd1phi", "trkd1z0", "trkd1d0", "trkd1charge", "trkd1pdgid","trkd1_dxy",
//                              "trkd2pt", "trkd2eta", "trkd2phi", "trkd2z0", "trkd2d0", "trkd2charge", "trkd2pdgid","trkd2_dxy",
//                              "trk_d0_min_delta_eta", "trk_d0_min_delta_phi", "trk_d0_min_delta_z0", "trk_d0_min_delta_d0", "trk_d0_min_delta_dR","trk_d0_min_delta_dist",
//                              "matchmud1pt", "matchmud1eta", "matchmud1phi", "matchmud1z0", "matchmud1pdgid","matchmud1_dxy",
//                              "matchmud2pt", "matchmud2eta", "matchmud2phi", "matchmud2z0", "matchmud2pdgid","matchmud2_dxy",
//                              "matchmu_d0_min_delta_eta", "matchmu_d0_min_delta_phi", "matchmu_d0_min_delta_z0", "matchmu_d0_min_delta_dR","matchmu_d0_min_delta_dist",
//                              "mu1pt", "mu1eta", "mu1phi", "mu1z0", "mu1d0", "mu1charge", "mu1pdgid","mu1_dxy",
//                              "mu2pt", "mu2eta", "mu2phi", "mu2z0", "mu2d0", "mu2charge", "mu2pdgid","mu2_dxy",
//                              "mu_mu_eta", "mu_mu_phi", "mu_mu_z0", "mu_mu_d0", "mu_mu_dR","mu_mu_dist",};
//    std::vector<int> nbins{100,
//                           20, 100, 100, 100, 100, 100, 400, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100,
//                           100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100, 50, 100,
//                           100, 100, 100, 100, 100, 100};
//    std::vector<float> lowEdge{0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, 0, 0, 0, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, 0, 0, 0, 0, 0,
//                               0, -2.4, -3.2, -20, 0, 0,
//                               0, -2.4, -3.2, -20, 0, 0,
//                               0, 0, 0, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, -2.4, -3.2, -20, -1, -2, 0, 0,
//                               0, 0, 0, 0, 0, 0};
//    std::vector<float> highEdge{100,
//                                20,  2.4, 3.2, 20, 1, 2, 400, 5,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                5, 4, 0.5, 0.5, 5, 1,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                5, 4, 0.5, 0.5, 5, 1,
//                                100, 2.4, 3.2, 20, 50, 5,
//                                100, 2.4, 3.2, 20, 50, 5,
//                                5, 4, 0.5, 5, 1,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                100, 2.4, 3.2, 20, 1, 2, 50, 5,
//                                5, 4, 0.5, 0.5, 5, 1};

//    typedef vector<TH1F *> Dim1;
//    typedef vector<Dim1> Dim2;
//    typedef vector<Dim2> Dim3;
//    typedef vector<Dim3> Dim4;
//    TH1I *h_Counter_TPcombination = new TH1I("h_Counter_TPcombination","h_Counter_TPcombination; TP combination chosen; Events",6,0,6);
//    TH1F *h_delta_dist_xy = new TH1F("h_delta_dist_xy","h_delta_dist_xy; (a) Distance between chosen TPs in x-y (cm) ; Events",50,0,1);
//    TH1F *h_error_delta_x = new TH1F("h_error_delta_x","h_error_delta_x; (b) #Delta x_{displaced} with chosen TP true x (cm) ; Events",100,0,2);
//    TH1F *h_delta_x = new TH1F("h_delta_x","h_delta_x; #Delta x between chosen TPs true x (cm) ; Events",40,0,2);
//    TH1F *h_delta_dist_z = new TH1F("h_delta_dist_z","h_delta_dist_z; (c) Distance between chosen TPs in z (cm) ; Events",40,0,2);
//    TH1F *h_error_delta_z = new TH1F("h_error_delta_z","h_error_delta_z; (d) Chosen TP error in z (cm) ; Events",100,0,10);
//    TH2F *h_d0_Mu_n_Mu_p = new TH2F("h_d0_Mu_n_Mu_p","h_d0_Mu_n_Mu_p",50,-5,5,50,-5,5);
//    TH2F *h_displaced_vertex_x = new TH2F("h_displaced_vertex_x","(b) displaced_vertex_x",40,-1,1,40,-1,1);
//    TH2F *h_displaced_vertex_y = new TH2F("h_displaced_vertex_y","(b) displaced_vertex_y",40,-1,1,40,-1,1);
//    TH2F *h_displaced_vertex_z = new TH2F("h_displaced_vertex_z","(b) displaced_vertex_z",40,-10,10,40,-10,10);
   
//    TH1I *h_trk_Counter_TPcombination = new TH1I("h_trk_Counter_TPcombination","h_trk_Counter_TPcombination; Track combination chosen; Events",6,0,6);
//    TH1F *h_trk_delta_dist_xy = new TH1F("h_trk_delta_dist_xy","h_trk_delta_dist_xy; (a) Distance between chosen Tracks in x-y (cm) ; Events",50,0,1);
//    TH1F *h_trk_error_delta_x = new TH1F("h_trk_error_delta_x","h_trk_error_delta_x; (b) #Delta x_{displaced} with chosen Track true x (cm) ; Events",100,0,2);
//    TH1F *h_trk_delta_x = new TH1F("h_trk_delta_x","h_trk_delta_x; #Delta x between chosen Track true x (cm) ; Events",40,0,2);
//    TH1F *h_trk_delta_dist_z = new TH1F("h_trk_delta_dist_z","h_trk_delta_dist_z; (c) Distance between chosen Tracks in z (cm) ; Events",40,0,2);
//    TH1F *h_trk_error_delta_z = new TH1F("h_trk_error_delta_z","h_trk_error_delta_z; (d) Chosen Track error in z (cm) ; Events",100,0,10);
//    TH2F *h_trk_d0_Mu_n_Mu_p = new TH2F("h_trk_d0_Mu_n_Mu_p","h_trk_d0_Mu_n_Mu_p",50,-5,5,50,-5,5);
//    TH2F *h_trk_displaced_vertex_x = new TH2F("h_trk_displaced_vertex_x","(b) displaced_vertex_x",40,-20,20,40,-20,20);
//    TH2F *h_trk_displaced_vertex_y = new TH2F("h_trk_displaced_vertex_y","(b) displaced_vertex_y",40,-20,20,40,-20,20);
//    TH2F *h_trk_displaced_vertex_z = new TH2F("h_trk_displaced_vertex_z","(b) displaced_vertex_z",60,-30,30,60,-30,30);

//    TH1F *h_eff_trk_pt = new TH1F("h_eff_trk_pt","h_eff_trk_pt; p_{T} of Leading p_{T} trk (GeV); Efficiency",50,0,50);
//    TH1F *h_all_trk_pt = new TH1F("h_all_trk_pt","h_all_trk_pt; p_{T} of Leading p_{T} trk (GeV); Events",50,0,50);
//    TH1F *h_correct_trk_pt = new TH1F("h_correct_trk_pt","h_correct_trk_pt; p_{T} of Leading p_{T} trk (GeV); Events",50,0,50);
//    TH1F *h_eff_trk_eta = new TH1F("h_eff_trk_eta","h_eff_trk_eta; #eta of Leading p_{T} trk ; Efficiency",50,-2.4,2.4);
//    TH1F *h_all_trk_eta = new TH1F("h_all_trk_eta","h_all_trk_eta; #eta of Leading p_{T} trk ; Events",50,-2.4,2.4);
//    TH1F *h_correct_trk_eta = new TH1F("h_correct_trk_eta","h_correct_trk_eta; #eta of Leading p_{T} trk ; Events",50,-2.4,2.4);
//    TH1F *h_eff_trk_dxy = new TH1F("h_eff_trk_dxy","h_eff_trk_dxy; dxy of vertex (cm); Efficiency",20,0,2);
//    TH1F *h_all_trk_dxy = new TH1F("h_all_trk_dxy","h_all_trk_dxy; dxy of vertex (cm); Events",20,0,2);
//    TH1F *h_correct_trk_dxy = new TH1F("h_correct_trk_dxy","h_correct_trk_dxy; dxy of vertex (cm) ; Events",20,0,2);

//    TH1F *h_res_tp_trk_x = new TH1F("h_res_tp_trk_x","h_res_tp_trk_x; x of vertex (cm) ; Events",100,-1,1);
//    TH1F *h_res_tp_trk_y = new TH1F("h_res_tp_trk_y","h_res_tp_trk_y; y of vertex (cm) ; Events",100,-1,1);
//    TH1F *h_res_tp_trk_z = new TH1F("h_res_tp_trk_z","h_res_tp_trk_z; z of vertex (cm) ; Events",1000,-10,10);

//    TH2F *h_Count_trk_pt_d0 = new TH2F("h_Count_trk_pt_d0","h_Count_trk_pt_d0",5,0.1,0.6,5,2,12);
//    // TH1F *h_tp_xy_cond = new TH1F("h_tp_xy_cond","h_tp_xy_cond; x-y distance; Events",50,0,5);
//    // TH1F *h_tp_z_cond = new TH1F("h_tp_z_cond","h_tp_z_cond; z distance; Events",50,0,5);
//    // TH1F *h_trk_xy_cond = new TH1F("h_trk_xy_cond","h_trk_xy_cond; Track x-y distance; Events",50,0,5);
//    // TH1F *h_trk_z_cond = new TH1F("h_trk_z_cond","h_trk_z_cond; Track z distance; Events",50,0,5);

//    TH1F *h_tp_xy_cond1 = new TH1F("h_tp_xy_cond1","h_tp_xy_cond1; x-y distance; Events",50,0,5);
//    TH1F *h_tp_z_cond1 = new TH1F("h_tp_z_cond1","h_tp_z_cond1; z distance; Events",50,0,5);
//    TH1F *h_trk_xy_cond1 = new TH1F("h_trk_xy_cond1","h_trk_xy_cond1; Track x-y distance; Events",50,0,5);
//    TH1F *h_trk_z_cond1 = new TH1F("h_trk_z_cond1","h_trk_z_cond1; Track z distance; Events",50,0,5);

//    TH1F *h_tp_xy_cond2 = new TH1F("h_tp_xy_cond2","h_tp_xy_cond2; x-y distance; Events",50,0,5);
//    TH1F *h_tp_z_cond2 = new TH1F("h_tp_z_cond2","h_tp_z_cond2; z distance; Events",50,0,5);
//    TH1F *h_trk_xy_cond2 = new TH1F("h_trk_xy_cond2","h_trk_xy_cond2; Track x-y distance; Events",50,0,5);
//    TH1F *h_trk_z_cond2 = new TH1F("h_trk_z_cond2","h_trk_z_cond2; Track z distance; Events",50,0,5);
   

//    Dim2 Hists(regions.size(), Dim1(vars.size()));
//    std::stringstream name;
//    TH1F *h_test;
//    Float_t **Output_vars = new Float_t *[regions.size()];
//    for (int k = 0; k < regions.size(); ++k)
//    {
//       Output_vars[k] = new Float_t[vars.size()];
//       for (int l = 0; l < vars.size(); ++l)
//       {
//          name << regions[k] << "_" << vars[l];
//          h_test = new TH1F((name.str()).c_str(), (name.str()).c_str(), nbins[l], lowEdge[l], highEdge[l]);
//          h_test->StatOverflows(kTRUE);
//          h_test->Sumw2(kTRUE);
//          Hists[k][l] = h_test;
//          name.str("");
//          Output_vars[k][l] = 0;
//       }
//    }

//    if (fChain == 0) return;

//    // Counters 
//    int nAccept = 0;           // Acceptance 
//    int nCount_2trk = 0;       // atleast 2 tracks
//    int nCount_leadingpt = 0;  // Leading pt TPs come from same true vertex
//    int nCount_intersect = 0;  // Tracks intersect in x-y plane
//    int nCount_no_intersect_z = 0;  // Tracks intersect in x-y plane but not in z
//    int nCount_1intersect = 0;        // only one intersection
//    int nCount_only1accept = 0;        // only one x,y within 20cm
//    int nCount_accept_both = 0;        // Both x,y are within 20cm
//    int nCount_no_intersect = 0; // Tracks don't intersect in x-y 
//    int nCount_z = 0;          // Identified z is true z
//    int nCount_matchtrk = 0; // Both Tracking Particles are matched to tracks
//    int nCount_xyz = 0; // Both chosen TPs within 1.0cm of (x_dv,y_dv,z_dv)
//    int nCount_tp_x = 0;   // TP and true have same x within 0.5
//    int nCount_tp_y = 0;   // TP and true have same y within 0.5
//    int nCount_tp_z = 0;   // TP and true have same z within 0.5
//    // Tracks
//    int nCount_trk_intersect = 0;  // Tracks intersect in x-y plane
//    int nCount_trk_no_intersect_z = 0;  // Tracks intersect in x-y plane but not in z
//    int nCount_trk_1intersect = 0;        // only one intersection
//    int nCount_trk_only1accept = 0;        // only one x,y within 20cm
//    int nCount_trk_accept_both = 0;        // Both x,y are within 20cm
//    int nCount_trk_no_intersect = 0; // Tracks don't intersect in x-y 
//    int nCount_trk_xyz = 0; // Both chosen tracks within 1.0cm of (x_dv,y_dv,z_dv)

//    int nCount_tp_trk_x = 0;   // TP and trk have same x within 0.5
//    int nCount_tp_trk_y = 0;   // TP and trk have same y within 0.5
//    int nCount_tp_trk_z = 0;   // TP and trk have same z within 0.5

//    float pt_cuts[5] = {2.0,4.0,6.0,8.0,10.0};     // Cuts to control event rate
//    float d0_cuts[5] = {0.1,0.2,0.3,0.4,0.5};
//    // float z0_cuts[5] = {0.2,0.4,0.6,0.8,1.0};
//    int nCount_trk_pt[5] = {0,0,0,0,0};
//    int nCount_trk_d0[5] = {0,0,0,0,0};
//    int nCount_trk_pt_d0[5][5] = {0};
//    // int nCount_trk_z0[5] = {0,0,0,0,0};
//    int nCount_2tp = 0;        // atleast 2 tracking particles
//    int nCount_tp_pt[5] = {0,0,0,0,0};
//    int nCount_tp_d0[5] = {0,0,0,0,0};
//    int nCount_tp_pt_d0[5][5] = {0};
//    for (int k=0; k < 5; k++){
//       for (int l=0; l < 5; l++){
//          nCount_tp_pt_d0[k][l] = 0;
//          nCount_trk_pt_d0[k][l] = 0;
//       }
//    }
//    // int nCount_tp_z0[5] = {0,0,0,0,0};
//    int nCount_OS=0;       // atleast 2 muons
//    int nCount_2muOS = 0;   // atleast 2 muons, with one pair of opposite sign

//    std::vector<Track_Parameters *> *selectedTracks;      // Tracks 
//    std::vector<Track_Parameters *> *selectedTracks_dmin; // Tracks from DV
//    std::vector<Track_Parameters *> *selectedTPs;         // Tracking particles
//    std::vector<Track_Parameters *> *selectedTPs_dmin;    // Tracking particles from DV
//    int nindex_mu_trigpass = -1;
//    int pindex_mu_trigpass = -1;
//    int nindex_os_trigpass = -1;
//    int pindex_os_trigpass = -1;

//    Long64_t nevt = fChain->GetEntries();
   
//    Long64_t i_in = 0;

//    for (Long64_t i=i_in; i<nevt; i++) { // 44,36,32,30,24,22,16,13,8,5
//       fChain->GetEntry(i);
//       int h_index = 0;

//       if(VERBOSE[0]||VERBOSE[2]){
//          std::cout<<"------------------------------------------------------------------------------------------------------"<<endl;
//          std::cout<<"Event Number  = "<<i<<endl;
//          if(i>i_in+50) break;
//       }
//       else if(VERBOSE[1]){
//          std::cout<<"------------------------------------------------------------------------------------------------------"<<endl;
//          std::cout<<"Event Number  = "<<i<<endl;
//          if(i>i_in+50) break;
//       }
//       else{
//          displayProgress(i, nevt);
//       }

//       selectedTracks = new std::vector<Track_Parameters *>();
//       selectedTracks_dmin = new std::vector<Track_Parameters *>();

//       selectedTPs = new std::vector<Track_Parameters *>();
//       selectedTPs_dmin = new std::vector<Track_Parameters *>();

//       // ----------------------------------------------------------------------------------------------------------------
//       // track loop

//       for (int it = 0; it < (int)trk_pt->size(); it++)
//       {
         
//          if (fabs(trk_eta->at(it)) > TP_maxEta)
//             continue;
//          if (trk_pt->at(it) < TP_minPt)
//             continue;
//          if (trk_pt->at(it) > TP_maxPt)
//             continue;
//          // if (std::fabs(trk_d0->at(it)) > TP_maxD0)
//          //    continue;
//          // if (std::fabs(trk_d0->at(it)) < TP_minD0)
//          //    continue;
// /*
//          int ndof = 2 * trk_nstub->at(it) - 5;
//          float chi2rphidof = (float)trk_chi2rphi->at(it) / ndof;
//          float chi2rzdof = (float)trk_chi2rz->at(it) / ndof;

//          if(chi2rphidof > 2)
//             continue;
//          if(chi2rzdof   > 2)
//             continue;
//          if(trk_bendchi2 ->at(it) > 5)
//             continue;
// */
//          selectedTracks->push_back(new Track_Parameters(trk_pt->at(it), trk_d0->at(it), trk_z0->at(it), trk_eta->at(it), trk_phi->at(it), -trk_rinv->at(it), it, -99999));
//       }

//       bool check_d0 = true;

//       for (int j = 0; j < selectedTracks->size(); j++){
//          if(fabs((*selectedTracks)[j]->d0) > TP_minD0){
//             check_d0=false;
//             break;
//          }
//       }
//       // if(check_d0) continue;

//       if(VERBOSE[0])
//          std::cout<<"End of Track loop"<<endl;
   
//       // ----------------------------------------------------------------------------------------------------------------
//       // tracking particle loop

//       for (int it = 0; it < (int)tp_pt->size(); it++)
//       {
//          float tmp_d0 = -tp_d0->at(it);
//          float tmp_z0 = tp_z0->at(it);
//          if (std::fabs(tmp_d0) > TP_maxD0)
//             continue;
//          // if (std::fabs(tmp_d0) < TP_minD0)
//          //    continue;
//          if (tp_pt->at(it) < TP_minPt)
//             continue;
//          if (tp_pt->at(it) > TP_maxPt)
//             continue;
//          if (std::fabs(tp_eta->at(it)) > TP_maxEta)
//             continue;
//           //! Issue with NaN. Verify NTuple.
//          // if (TMath::IsNaN(tmp_d0))
//          //    continue;
//          // if (TMath::IsNaN(tmp_z0))
//          //    continue;
         
         
//          // Only keep tracks of particular particles using PDGID
//          if (TP_select_pdgid != 0)
//          {
//             if (fabs(tp_pdgid->at(it)) != fabs(TP_select_pdgid))
//                continue;
//          }
//          selectedTPs->push_back(new Track_Parameters(tp_pt->at(it), tmp_d0, tmp_z0, tp_eta->at(it), tp_phi->at(it), tp_charge->at(it), it, abs(tp_pdgid->at(it))));

//          //----------------------------------------------------------------------------------------------
//          // matchtrk 
//          /*
//          float tmp_d0 = matchtrk_d0->at(it);
//          float tmp_z0 = matchtrk_z0->at(it);
//          if (std::fabs(tmp_d0) > TP_maxD0)
//             continue;
//          if (std::fabs(tmp_d0) < TP_minD0)
//             continue;
//          if (matchtrk_pt->at(it) < TP_minPt)
//             continue;
//          if (matchtrk_pt->at(it) > TP_maxPt)
//             continue;
//          if (std::fabs(matchtrk_eta->at(it)) > TP_maxEta)
//             continue;
//          //  //! Issue with NaN. Verify NTuple.
//          // if (TMath::IsNaN(tmp_d0))
//          //    continue;
//          // if (TMath::IsNaN(tmp_z0))
//          //    continue;
         
 
//          selectedTPs->push_back(new Track_Parameters(matchtrk_pt->at(it), tmp_d0, tmp_z0, matchtrk_eta->at(it), matchtrk_phi->at(it), -matchtrk_rinv->at(it), it, -99999));
//          */
//       } // end tp loop

//       if(VERBOSE[0]){
//          std::cout<<"End of TP loop"<<endl;
//          for (int j = 0; j < selectedTPs->size(); j++){
//             std::cout<<"Tracking Particle = "<<j<<endl;
//             std::cout<<Form("pt = %5.1f  |  d0 = %5.1f  |  z0 = %5.1f  |  pdgid = %5d  |  rho = %5.1f  |  x0 = %5.1f  |  y0 = %5.1f  |  atan = %5.1f  |  phi = %5.1f  |  Dphi = %5.1f",(*selectedTPs)[j]->pt,(*selectedTPs)[j]->d0,(*selectedTPs)[j]->z0,(*selectedTPs)[j]->pdgid,(*selectedTPs)[j]->rho,(*selectedTPs)[j]->x0,(*selectedTPs)[j]->y0,TMath::ATan2(-(*selectedTPs)[j]->x0,(*selectedTPs)[j]->y0),(*selectedTPs)[j]->phi,deltaPhi((*selectedTPs)[j]->phi,TMath::ATan2(-(*selectedTPs)[j]->x0,(*selectedTPs)[j]->y0)))<<endl;
//          }
//       }

//       // ----------------------------------------------------------------------------------------------------------------
//       // Gen particle loop
//       for (int it = 0; it < nPart; it++){
//          if (partPt[it] < TP_minPt)
//             continue;
//          if (partPt[it] > TP_maxPt)
//             continue;
//          if (partStat[it] != 1)
//             continue;
//          if (abs(partId[it]) != 13)
//             continue;

//          // selectedTPs->push_back(new Track_Parameters(partPt[it], -99999, -99999, partEta[it], partPhi[it], partCh[it], it, partId[it]));
//       } // end gen loop

//       if(VERBOSE[0])
//          std::cout<<"End of Gen Particle loop"<<endl;

// // --------------------------------------------------------------------------------------------
// //                Vertex finding in Tracking Particles
// // --------------------------------------------------------------------------------------------
//       if (!(selectedTracks->size() >= 2)) continue;
//       bool true_DV = false;

//       if(VERBOSE[0])
//          std::cout<<"End of z-vertex Finding"<<endl;
      
      
//       Double_t x_dv = -9999.0;// (tp_x->at((*selectedTPs)[0]->index));//+tp_x->at((*selectedTPs)[1]->index))/2.0;
//       Double_t y_dv = -9999.0;// (tp_y->at((*selectedTPs)[0]->index));//+tp_y->at((*selectedTPs)[1]->index))/2.0;
//       Double_t z_dv = -9999.0;// (tp_z->at((*selectedTPs)[0]->index));//+tp_z->at((*selectedTPs)[1]->index))/2.0;
            
//       int Vertex_check = -1;
//       int selected_track_j0 = -1;
//       int selected_track_j1 = -1;

//       Double_t x_tmp = x_dv;
//       Double_t y_tmp = y_dv;
//       Double_t z_tmp = z_dv;
//       Int_t Vertex_check_tmp = -1;

//       if(selectedTPs->size()>=2){

//       sort(selectedTPs->begin(), selectedTPs->end(), ComparePtTrack);

//       // Temporary selection of leading 2 p_T tracks from common vertex
//       if(!(!((fabs(tp_x->at((*selectedTPs)[0]->index) - tp_x->at((*selectedTPs)[1]->index))<0.01) && (fabs(tp_y->at((*selectedTPs)[0]->index) - tp_y->at((*selectedTPs)[1]->index))<0.01)))){nCount_leadingpt++;}// continue;

//       Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[1],x_dv,y_dv,z_dv); // 0+1
//       Vertex_check_tmp = Vertex_check;

//       if (Vertex_check>0){
//          selected_track_j0 = 0;
//          selected_track_j1 = 1;
//          h_Counter_TPcombination->AddBinContent(1); 
//       } 
//       else if(selectedTPs->size()>=3){
//          Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[2],x_dv,y_dv,z_dv);
//          if(Vertex_check>0){ //0+2
//             selected_track_j0 = 0;
//             selected_track_j1 = 2;
//             h_Counter_TPcombination->AddBinContent(2); 
//          }
//          else if((Vertex_check=Vertex((*selectedTPs)[1],(*selectedTPs)[2],x_dv,y_dv,z_dv))>0){ //1+2
//             selected_track_j0 = 1;
//             selected_track_j1 = 2;
//             h_Counter_TPcombination->AddBinContent(3); 
//          }
//          else if(selectedTPs->size()>=4){
//             Vertex_check = Vertex((*selectedTPs)[0],(*selectedTPs)[3],x_dv,y_dv,z_dv);
//             if(Vertex_check>0){ //0+3
//                selected_track_j0 = 0;
//                selected_track_j1 = 3;
//                h_Counter_TPcombination->AddBinContent(4); 
//             }
//             else if((Vertex_check=Vertex((*selectedTPs)[1],(*selectedTPs)[3],x_dv,y_dv,z_dv))>0){ //1+3
//                selected_track_j0 = 1;
//                selected_track_j1 = 3;
//                h_Counter_TPcombination->AddBinContent(5); 
//             }
//             else if((Vertex_check=Vertex((*selectedTPs)[2],(*selectedTPs)[3],x_dv,y_dv,z_dv))>0){ //2+3
//                selected_track_j0 = 2;
//                selected_track_j1 = 3;
//                h_Counter_TPcombination->AddBinContent(6); 
//             }
//             else{
               
//                x_dv = x_tmp;
//                y_dv = y_tmp;
//                z_dv = z_tmp;//(*selectedTPs)[0]->z0;
//                Vertex_check = Vertex_check_tmp;
               
//                selected_track_j0 = 0;
//                selected_track_j1 = 1;
//                h_Counter_TPcombination->AddBinContent(0); 
//             }
//          }
//          else{
            
//             x_dv = x_tmp;
//             y_dv = y_tmp;
//             z_dv = z_tmp;//(*selectedTPs)[0]->z0;
//             Vertex_check = Vertex_check_tmp;
            
//             selected_track_j0 = 0;
//             selected_track_j1 = 1;
//             h_Counter_TPcombination->AddBinContent(0); 
//          }
//       } 
//       else{
//          x_dv = x_tmp;
//          y_dv = y_tmp;
//          z_dv = z_tmp;//(*selectedTPs)[0]->z0;
//          Vertex_check = Vertex_check_tmp;
         
//          selected_track_j0 = 0;
//          selected_track_j1 = 1;
//          h_Counter_TPcombination->AddBinContent(0); 
//       }
      

//       if(tp_nmatch->at((*selectedTPs)[selected_track_j0]->index)>=1 && (*selectedTPs)[selected_track_j1]->index >=1){
//          nCount_matchtrk++;
//       }
//       if(fabs(z_dv-tp_z->at((*selectedTPs)[selected_track_j0]->index))<CUTOFF){
//          nCount_z++;
//       }
//       if(Vertex_check>=0){
//          nCount_intersect++;
//          if(Vertex_check==0){
//             nCount_no_intersect_z++;
//          }
//          if(Vertex_check==1){
//             nCount_1intersect++;
//          }
//          else if(Vertex_check==2){
//             nCount_only1accept++;
//          }
//          else if(Vertex_check==3){
//             nCount_accept_both++;
//          }
//       }
//       else if(Vertex_check==-1){
//          nCount_no_intersect++;
//       }
//       if((fabs(x_dv-(tp_x->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF)) nCount_tp_x++;
//       if((fabs(y_dv-(tp_y->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF)) nCount_tp_y++;
//       if((fabs(z_dv-(tp_z->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF)) nCount_tp_z++;

//       if(VERBOSE[0])
//          std::cout<<"Selected Tracks = "<<selected_track_j0<<"\t and  "<<selected_track_j1<<endl;      
      
//       // h_delta_dist_xy->Fill(fabs(dist_Vertex(x_dv, y_dv,(*selectedTPs)[0]) - dist_Vertex(x_dv, y_dv,(*selectedTPs)[1])));
//       // if(Vertex_check>0){
//       h_delta_dist_xy->Fill(dist_TPs((*selectedTPs)[selected_track_j0],(*selectedTPs)[selected_track_j1]));
//       h_error_delta_x->Fill(fabs((tp_x->at((*selectedTPs)[selected_track_j0]->index) - x_dv)));

//       float z_dv_1 = (*selectedTPs)[selected_track_j0]->z(x_dv,y_dv);
//       float z_dv_2 = (*selectedTPs)[selected_track_j1]->z(x_dv,y_dv);

//       h_delta_dist_z->Fill(fabs(z_dv_1-z_dv_2));

//       h_error_delta_z->Fill(fabs((tp_z->at((*selectedTPs)[selected_track_j0]->index) - z_dv_1)));

//       // h_displaced_vertex_x->Fill( x_dv,(tp_x->at((*selectedTPs)[selected_track_j0]->index)));
//       // h_displaced_vertex_y->Fill( y_dv,(tp_y->at((*selectedTPs)[selected_track_j0]->index)));
//       // h_displaced_vertex_z->Fill( z_dv,(tp_z->at((*selectedTPs)[selected_track_j0]->index)));
//       // }
//       h_tp_xy_cond1->Fill(fabs(x_dv-(tp_x->at((*selectedTPs)[selected_track_j0]->index))));
//       h_tp_z_cond1->Fill(fabs(z_dv-(tp_z->at((*selectedTPs)[selected_track_j0]->index))));      
//       h_tp_xy_cond2->Fill((dist_Vertex(x_dv,y_dv,(*selectedTPs)[selected_track_j1])));
//       h_tp_z_cond2->Fill((fabs(z_dv_2-z_dv)));
      
//       //h_tp_xy_cond->Fill((dist_Vertex(x_dv,y_dv,(*selectedTPs)[selected_track_j0]))+ (dist_Vertex(x_dv,y_dv,(*selectedTPs)[selected_track_j0])));
//       //h_tp_z_cond->Fill((fabs(z_dv_1-z_dv)));
      
//       // if((dist_Vertex(x_dv,y_dv,(*selectedTPs)[selected_track_j0])<CUTOFF) && (fabs(z_dv_1-z_dv)<CUTOFF) && (dist_Vertex(x_dv,y_dv,(*selectedTPs)[selected_track_j1])<CUTOFF) && (fabs(z_dv_2-z_dv)<CUTOFF) ){
//       if((fabs(x_dv-(tp_x->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF) && (fabs(z_dv-(tp_z->at((*selectedTPs)[selected_track_j0]->index)))<CUTOFF) && (fabs(z_dv-tp_z->at((*selectedTPs)[selected_track_j0]->index))<CUTOFF)){
//          nCount_xyz++;
//          h_delta_x->Fill(fabs((tp_x->at((*selectedTPs)[selected_track_j0]->index) - (tp_x->at((*selectedTPs)[selected_track_j1]->index)))));
//          true_DV = true;
//          // h_correct_tp_pt->Fill((*selectedTPs)[0]->pt);
//          // h_correct_tp_eta->Fill((*selectedTPs)[0]->eta);
//          // h_correct_tp_dxy->Fill(dist(x_dv,y_dv));
//       }
//       }

// // --------------------------------------------------------------------------------------------
// //                Vertex finding in Tracks
// // --------------------------------------------------------------------------------------------


//       sort(selectedTracks->begin(), selectedTracks->end(), ComparePtTrack);
//       Double_t x_dv_trk = -9999.0;// (tp_x->at((*selectedTracks)[0]->index));//+tp_x->at((*selectedTracks)[1]->index))/2.0;
//       Double_t y_dv_trk = -9999.0;// (tp_y->at((*selectedTracks)[0]->index));//+tp_y->at((*selectedTracks)[1]->index))/2.0;
//       Double_t z_dv_trk = -9999.0;// (tp_z->at((*selectedTracks)[0]->index));//+tp_z->at((*selectedTracks)[1]->index))/2.0;
            
//       Vertex_check = -1;
//       selected_track_j0 = -1;
//       selected_track_j1 = -1;
//       Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[1],x_dv_trk,y_dv_trk,z_dv_trk); // 0+1
//       x_tmp = x_dv_trk;
//       y_tmp = y_dv_trk;
//       z_tmp = z_dv_trk;
//       Vertex_check_tmp = Vertex_check;
//       if (Vertex_check>0){
//          selected_track_j0 = 0;
//          selected_track_j1 = 1;
//          h_trk_Counter_TPcombination->AddBinContent(1); 
//       } 
//       else if(selectedTracks->size()>=3){
//          Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[2],x_dv_trk,y_dv_trk,z_dv_trk);
//          if(Vertex_check>0){ //0+2
//             selected_track_j0 = 0;
//             selected_track_j1 = 2;
//             h_trk_Counter_TPcombination->AddBinContent(2); 
//          }
//          else if((Vertex_check=Vertex((*selectedTracks)[1],(*selectedTracks)[2],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //1+2
//             selected_track_j0 = 1;
//             selected_track_j1 = 2;
//             h_trk_Counter_TPcombination->AddBinContent(3); 
//          }
//          else if(selectedTracks->size()>=4){
//             Vertex_check = Vertex((*selectedTracks)[0],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk);
//             if(Vertex_check>0){ //0+3
//                selected_track_j0 = 0;
//                selected_track_j1 = 3;
//                h_trk_Counter_TPcombination->AddBinContent(4); 
//             }
//             else if((Vertex_check=Vertex((*selectedTracks)[1],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //1+3
//                selected_track_j0 = 1;
//                selected_track_j1 = 3;
//                h_trk_Counter_TPcombination->AddBinContent(5); 
//             }
//             else if((Vertex_check=Vertex((*selectedTracks)[2],(*selectedTracks)[3],x_dv_trk,y_dv_trk,z_dv_trk))>0){ //2+3
//                selected_track_j0 = 2;
//                selected_track_j1 = 3;
//                h_trk_Counter_TPcombination->AddBinContent(6); 
//             }
//             else{
               
//                x_dv_trk = x_tmp;
//                y_dv_trk = y_tmp;
//                z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
//                Vertex_check = Vertex_check_tmp;
               
//                selected_track_j0 = 0;
//                selected_track_j1 = 1;
//                h_trk_Counter_TPcombination->AddBinContent(0); 
//             }
//          }
//          else{
            
//             x_dv_trk = x_tmp;
//             y_dv_trk = y_tmp;
//             z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
//             Vertex_check = Vertex_check_tmp;
            
//             selected_track_j0 = 0;
//             selected_track_j1 = 1;
//             h_trk_Counter_TPcombination->AddBinContent(0); 
//          }
//       } 
//       else{
//          x_dv_trk = x_tmp;
//          y_dv_trk = y_tmp;
//          z_dv_trk = z_tmp;//(*selectedTracks)[0]->z0;
//          Vertex_check = Vertex_check_tmp;
         
//          selected_track_j0 = 0;
//          selected_track_j1 = 1;
//          h_trk_Counter_TPcombination->AddBinContent(0); 
//       }

//       //! Due to difference in d0 signs for trks and tps - fix the d0 sign directly
//       // x_dv_trk *= -1;
//       // y_dv_trk *= -1;

//       h_trk_delta_dist_xy->Fill(dist_TPs((*selectedTracks)[selected_track_j0],(*selectedTracks)[selected_track_j1]));
//       // h_trk_error_delta_x->Fill(fabs((tp_x->at((*selectedTracks)[selected_track_j0]->index) - x_dv_trk)));
//       // h_trk_delta_x->Fill(fabs((tp_x->at((*selectedTracks)[selected_track_j0]->index) - (tp_x->at((*selectedTracks)[selected_track_j1]->index)))));

//       float z_dv_trk_1 = (*selectedTracks)[selected_track_j0]->z(x_dv_trk,y_dv_trk);
//       float z_dv_trk_2 = (*selectedTracks)[selected_track_j1]->z(x_dv_trk,y_dv_trk);

//       h_trk_delta_dist_z->Fill(fabs(z_dv_trk_1-z_dv_trk_2)); //* Use this condition
//       // h_trk_delta_dist_z->Fill(fabs(x_dv - x_dv_trk));

//       // h_trk_error_delta_z->Fill(fabs((tp_z->at((*selectedTracks)[selected_track_j0]->index) - z_dv_trk_1)));

      
//       h_trk_xy_cond1->Fill(fabs(x_dv-x_dv_trk));
//       h_trk_z_cond1->Fill(fabs(z_dv-z_dv_trk));      
//       h_trk_xy_cond2->Fill((dist_Vertex(x_dv_trk,y_dv_trk,(*selectedTracks)[selected_track_j1])));
//       h_trk_z_cond2->Fill((fabs(z_dv_trk_2-z_dv_trk)));

//       h_displaced_vertex_x->Fill( x_dv,x_dv_trk);
//       h_displaced_vertex_y->Fill( y_dv,y_dv_trk);
//       h_displaced_vertex_z->Fill( z_dv,z_dv_trk);
      
//       if((dist_Vertex(x_dv_trk,y_dv_trk,(*selectedTracks)[selected_track_j0])<CUTOFF) && (fabs(z_dv_trk_1-z_dv_trk)<CUTOFF) && (dist_Vertex(x_dv_trk,y_dv_trk,(*selectedTracks)[selected_track_j1])<CUTOFF) && (fabs(z_dv_trk_2-z_dv_trk)<CUTOFF) ){
//          nCount_trk_xyz++;
//       }
//       if(Vertex_check>=0){
//          nCount_trk_intersect++;
//          if(Vertex_check==0){
//             nCount_trk_no_intersect_z++;
//          }
//          if(Vertex_check==1){
//             nCount_trk_1intersect++;
//          }
//          else if(Vertex_check==2){
//             nCount_trk_only1accept++;
//          }
//          else if(Vertex_check==3){
//             nCount_trk_accept_both++;
//          }
//       }
//       else if(Vertex_check==-1){
//          nCount_trk_no_intersect++;
//       }
//       if((fabs(x_dv-x_dv_trk)<CUTOFF)) nCount_tp_trk_x++;
//       if((fabs(y_dv-y_dv_trk)<CUTOFF)) nCount_tp_trk_y++;
//       if((fabs(z_dv-z_dv_trk)<CUTOFF)) nCount_tp_trk_z++;
      
//       if(true_DV){

//       h_all_trk_pt->Fill((*selectedTracks)[selected_track_j0]->pt);
//       h_all_trk_eta->Fill((*selectedTracks)[selected_track_j0]->eta);
//       h_all_trk_dxy->Fill(dist(x_dv_trk,y_dv_trk));
//       }

//       // if((dist_Vertex(x_dv_trk,y_dv_trk,(*selectedTracks)[selected_track_j0])<CUTOFF) && (fabs(z_dv_trk_1-z_dv_trk)<CUTOFF)){
//       if(true_DV && (fabs(x_dv-x_dv_trk)<CUTOFF) && (fabs(y_dv-y_dv_trk)<CUTOFF) && (fabs(z_dv-z_dv_trk)<CUTOFF)){
//          h_correct_trk_pt->Fill((*selectedTracks)[selected_track_j0]->pt);
//          h_correct_trk_eta->Fill((*selectedTracks)[selected_track_j0]->eta);
//          h_correct_trk_dxy->Fill(dist(x_dv_trk,y_dv_trk));
            
//          h_res_tp_trk_x->Fill(x_dv-x_dv_trk);
//          h_res_tp_trk_y->Fill(y_dv-y_dv_trk);
//          h_res_tp_trk_z->Fill(z_dv-z_dv_trk);
//       }
      


//    //---------------------------------------------------------------------------------------------------------
      
//       if(VERBOSE[0])
//          std::cout<<"End of Displaced Vertex Finding"<<endl;
//       // if(fabs(tp_x->at((*selectedTPs)[0]->index) - x_dv)>0.1)
//       if(VERBOSE[1]){
//          std::cout<<Form("x_dv = %5.2f  |  y_dv = %5.2f  |  z_dv = %5.2f  |  Vertex_check = %2d",x_dv,y_dv,z_dv,Vertex_check)<<endl;;//,(*selectedTPs)[0]->z((*selectedTPs)[0]->phi_T(tp_x->at((*selectedTPs)[0]->index))),(*selectedTPs)[0]->z0)<<endl;//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
//       }
//       if(VERBOSE[2]){
//          for (int j = 0; j < selectedTPs->size() && j<4; j++){
//             std::cout<<"Tracking Particle = "<<j<<endl;
//             std::cout<<Form("pt = %5.1f  |  d0 = %5.1f  |  z0 = %5.1f  |  pdgid = %5d  |  rho = %5.1f  |  x = %5.1f  |  y = %5.1f  |  z = %5.1f",(*selectedTPs)[j]->pt,(*selectedTPs)[j]->d0,(*selectedTPs)[j]->z0,(*selectedTPs)[j]->pdgid,(*selectedTPs)[j]->rho,tp_x->at((*selectedTPs)[j]->index),tp_y->at((*selectedTPs)[j]->index),tp_z->at((*selectedTPs)[j]->index))<<endl;
//          }
//       }

//       // Counters 
//       if(VERBOSE[0])
//          std::cout << "Start of Counters"<<endl;

//       int pt_check[5] = {0,0,0,0,0};
//       int d0_check[5] = {0,0,0,0,0};
//       int pt_d0_check[5][5];
//       for (int k=0; k < 5; k++){
//          for (int l=0; l < 5; l++){
//             pt_d0_check[k][l] = 0;
//          }
//       }

//       if ((selectedTracks->size() == 2)){ nCount_2trk++;}

//       if ((selectedTracks->size() >= 2)){ 
//       for (int j = 0; j < selectedTracks->size(); j++){
//          // for (int k = j+1; k < selectedTracks->size(); k++){
//             // if((fabs((*selectedTracks)[j]->z(x_dv_trk,y_dv_trk,z_dv_trk) - z_dv_trk)) < 1.0){
//                if ((*selectedTracks)[j]->pt   >   pt_cuts[0]){ pt_check[0]++; }
//                if ((*selectedTracks)[j]->pt   >   pt_cuts[1]){ pt_check[1]++; }
//                if ((*selectedTracks)[j]->pt   >   pt_cuts[2]){ pt_check[2]++; }
//                if ((*selectedTracks)[j]->pt   >   pt_cuts[3]){ pt_check[3]++; }
//                if ((*selectedTracks)[j]->pt   >   pt_cuts[4]){ pt_check[4]++; }

//                if (fabs((*selectedTracks)[j]->d0)>d0_cuts[0]){ d0_check[0]++;}
//                if (fabs((*selectedTracks)[j]->d0)>d0_cuts[1]){ d0_check[1]++;}
//                if (fabs((*selectedTracks)[j]->d0)>d0_cuts[2]){ d0_check[2]++;}
//                if (fabs((*selectedTracks)[j]->d0)>d0_cuts[3]){ d0_check[3]++;}
//                if (fabs((*selectedTracks)[j]->d0)>d0_cuts[4]){ d0_check[4]++;}
//                for (int l=0; l < 5; l++){
//                   if ((*selectedTracks)[j]->pt   >   pt_cuts[0] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
//                   if ((*selectedTracks)[j]->pt   >   pt_cuts[1] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
//                   if ((*selectedTracks)[j]->pt   >   pt_cuts[2] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
//                   if ((*selectedTracks)[j]->pt   >   pt_cuts[3] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
//                   if ((*selectedTracks)[j]->pt   >   pt_cuts[4] && fabs((*selectedTracks)[j]->d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
//                }
//             // }
//         // }
//       }

//       for(int k=0;k<5;k++){
//          if(pt_check[k]>=2) nCount_trk_pt[k]++;
//          if(d0_check[k]>=2) nCount_trk_d0[k]++;
//          pt_check[k] = 0;
//          d0_check[k] = 0;
//          for (int l=0; l < 5; l++){
//             if(pt_d0_check[k][l]>=2){
//                nCount_trk_pt_d0[k][l]++;
//                h_Count_trk_pt_d0->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0->GetBinContent(l+1,k+1) + 1));
//             }
//             pt_d0_check[k][l] = 0;
//          }
//       }
//       }
      
//       if ((selectedTPs->size() >= 2)){ nCount_2tp++;}

//       if ((selectedTPs->size() >= 2)){ 
//       for (int j = 0; j < selectedTPs->size(); j++){
//         // for (int k = j+1; k < selectedTPs->size(); k++){
//             // if((fabs((*selectedTPs)[j]->z(x_dv,y_dv,z_dv) - z_dv)) < 1.0){
//                if ((*selectedTPs)[j]->pt   >   pt_cuts[0]){ pt_check[0]++; }
//                if ((*selectedTPs)[j]->pt   >   pt_cuts[1]){ pt_check[1]++; }
//                if ((*selectedTPs)[j]->pt   >   pt_cuts[2]){ pt_check[2]++; }
//                if ((*selectedTPs)[j]->pt   >   pt_cuts[3]){ pt_check[3]++; }
//                if ((*selectedTPs)[j]->pt   >   pt_cuts[4]){ pt_check[4]++; }

//                if (fabs((*selectedTPs)[j]->d0)>d0_cuts[0]){ d0_check[0]++;}
//                if (fabs((*selectedTPs)[j]->d0)>d0_cuts[1]){ d0_check[1]++;}
//                if (fabs((*selectedTPs)[j]->d0)>d0_cuts[2]){ d0_check[2]++;}
//                if (fabs((*selectedTPs)[j]->d0)>d0_cuts[3]){ d0_check[3]++;}
//                if (fabs((*selectedTPs)[j]->d0)>d0_cuts[4]){ d0_check[4]++;}
//                for (int l=0; l < 5; l++){
//                   if ((*selectedTPs)[j]->pt   >   pt_cuts[0] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[0][l]++;}
//                   if ((*selectedTPs)[j]->pt   >   pt_cuts[1] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[1][l]++;}
//                   if ((*selectedTPs)[j]->pt   >   pt_cuts[2] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[2][l]++;}
//                   if ((*selectedTPs)[j]->pt   >   pt_cuts[3] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[3][l]++;}
//                   if ((*selectedTPs)[j]->pt   >   pt_cuts[4] && fabs((*selectedTPs)[j]->d0)>d0_cuts[l]) { pt_d0_check[4][l]++;}
//                }
//             // }
//         // }
//       }

//       for(int k=0;k<5;k++){
//          if(pt_check[k]>=2) nCount_tp_pt[k]++;
//          if(d0_check[k]>=2) nCount_tp_d0[k]++;
//          pt_check[k] = 0;
//          d0_check[k] = 0;  
//          for (int l=0; l < 5; l++){
//             if(pt_d0_check[k][l]>=2) nCount_tp_pt_d0[k][l]++;
//             pt_d0_check[k][l] = 0;
//          }
//       }
//       }

//       if(VERBOSE[0])
//          std::cout << "End of Counters"<<endl;

//       for (int j = 0; j < selectedTPs->size(); j++){
//          (*selectedTPs)[j]->dist_calc(x_dv,y_dv,tp_x->at((*selectedTPs)[j]->index),tp_y->at((*selectedTPs)[j]->index));
//          // (*selectedTPs)[selectedTPs->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedTPs)[selectedTPs->size()-1]);
//          if(dist_Vertex(x_dv, y_dv,(*selectedTPs)[j]) < 1.0){
//           if((fabs((*selectedTPs)[j]->z(x_dv,y_dv,z_dv) - z_dv)) < 1.0){
//             selectedTPs_dmin->push_back(new Track_Parameters((*selectedTPs)[j]->pt, (*selectedTPs)[j]->d0, (*selectedTPs)[j]->z0, (*selectedTPs)[j]->eta, (*selectedTPs)[j]->phi, (*selectedTPs)[j]->charge, (*selectedTPs)[j]->index, (*selectedTPs)[j]->pdgid));
//             (*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->dist_calc(x_dv,y_dv,tp_x->at((*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->index),tp_y->at((*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->index));
//             // (*selectedTPs_dmin)[selectedTPs_dmin->size()-1]->dxy = dist_Vertex(x_dv, y_dv,(*selectedTPs_dmin)[selectedTPs_dmin->size()-1]);
//           }
//          }

//       }

//       for (int j = 0; j < selectedTracks->size(); j++){
//          // (*selectedTracks)[j]->dist_calc(x_dv_trk,y_dv_trk,tp_x->at((*selectedTracks)[j]->index),tp_y->at((*selectedTracks)[j]->index));
//          // (*selectedTracks)[selectedTracks->size()-1]->dxy = dist_Vertex(x_dv_trk, y_dv_trk,(*selectedTracks)[selectedTracks->size()-1]);
//          if(dist_Vertex(x_dv_trk, y_dv_trk,(*selectedTracks)[j]) < 1.0){
//           if((fabs((*selectedTracks)[j]->z(x_dv_trk,y_dv_trk,z_dv_trk) - z_dv_trk)) < 1.0){
//             selectedTracks_dmin->push_back(new Track_Parameters((*selectedTracks)[j]->pt, (*selectedTracks)[j]->d0, (*selectedTracks)[j]->z0, (*selectedTracks)[j]->eta, (*selectedTracks)[j]->phi, (*selectedTracks)[j]->charge, (*selectedTracks)[j]->index, (*selectedTracks)[j]->pdgid));
//             // (*selectedTracks_dmin)[selectedTracks_dmin->size()-1]->dist_calc(x_dv_trk,y_dv_trk,tp_x->at((*selectedTracks_dmin)[selectedTracks_dmin->size()-1]->index),tp_y->at((*selectedTracks_dmin)[selectedTracks_dmin->size()-1]->index));
//             // (*selectedTracks_dmin)[selectedTracks_dmin->size()-1]->dxy = dist_Vertex(x_dv_trk, y_dv_trk,(*selectedTracks_dmin)[selectedTracks_dmin->size()-1]);
//           }
//          }
//       }

//       pindex_os_trigpass = -1;
//       nindex_os_trigpass = -1;

//       bool triggerpass = false;
//       for (int j = 0; j < selectedTracks_dmin->size() && !triggerpass; j++){
//          for (int k = j+1; k < selectedTracks_dmin->size() && !triggerpass; k++){
//             if((*selectedTracks_dmin)[j]->charge*(*selectedTracks_dmin)[k]->charge == -1){
//                triggerpass = true;
//                if((*selectedTracks_dmin)[j]->charge == 1){
//                   pindex_os_trigpass = j;
//                   nindex_os_trigpass = k;
//                }
//                else{
//                   pindex_os_trigpass = k;
//                   nindex_os_trigpass = j;
//                }
               
//             }
//          }
//       }
      
//       if (triggerpass){
//          nCount_OS++;
//       }

//       pindex_mu_trigpass = -1;
//       nindex_mu_trigpass = -1;

//       bool triggerpass_mu = false;
//       for (int j = 0; j < selectedTracks_dmin->size() && !triggerpass_mu; j++){
//          for (int k = j+1; k < selectedTracks_dmin->size() && !triggerpass_mu; k++){
//             if((*selectedTracks_dmin)[j]->charge*(*selectedTracks_dmin)[k]->charge == -1 && fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[j]->index)) == 13 && fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[k]->index)) == 13){
//                triggerpass_mu = true;
//                if((*selectedTracks_dmin)[j]->charge == 1){
//                   pindex_mu_trigpass = j;
//                   nindex_mu_trigpass = k;
//                }
//                else{
//                   pindex_mu_trigpass = k;
//                   nindex_mu_trigpass = j;
//                }
               
//             }
//          }
//       }

//       if (triggerpass_mu){
//          nCount_2muOS++;
//       }

//       if(VERBOSE[0])
//          std::cout << "Selected the Muons"<<endl;

//       if (!(selectedTracks->size() >= 2)) continue;
//       nAccept++;
//          // ---------------------------------------------------------------------------------------------------------
//          //Filling up Histograms
         
//          h_index=0;
//          Hists[0][h_index++]->Fill(selectedTracks->size());
//          for (int i = 0; i < selectedTracks->size(); i++)
//          {
//             Hists[0][h_index]->Fill((*selectedTracks)[i]->pt);
//             Hists[0][h_index+1]->Fill((*selectedTracks)[i]->eta);
//             Hists[0][h_index+2]->Fill((*selectedTracks)[i]->phi);
//             Hists[0][h_index+3]->Fill((*selectedTracks)[i]->z0);
//             Hists[0][h_index+4]->Fill((*selectedTracks)[i]->d0);
//             Hists[0][h_index+5]->Fill((*selectedTracks)[i]->charge);
//             Hists[0][h_index+6]->Fill((*selectedTracks)[i]->pdgid);
//             Hists[0][h_index+7]->Fill((*selectedTracks)[i]->dxy);
//          }
//          h_index +=8;
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->charge);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->pdgid);
//          Hists[0][h_index++]->Fill((*selectedTracks)[0]->dxy);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->charge);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->pdgid);
//          Hists[0][h_index++]->Fill((*selectedTracks)[1]->dxy);
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks)[0]->eta - (*selectedTracks)[1]->eta));
//          Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTracks)[0]->phi,(*selectedTracks)[1]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks)[0]->z0 - (*selectedTracks)[1]->z0));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks)[0]->d0 - (*selectedTracks)[1]->d0));
//          Hists[0][h_index++]->Fill((deltaR((*selectedTracks)[0]->eta,(*selectedTracks)[0]->phi, (*selectedTracks)[1]->eta,(*selectedTracks)[1]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks)[0]->dxy - (*selectedTracks)[1]->dxy));

//          if(VERBOSE[0])
//             std::cout << "All tracks Histograms filled"<<endl;

//          if (triggerpass) {
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->charge);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->pdgid);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_os_trigpass]->dxy);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->charge);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->pdgid);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_os_trigpass]->dxy);
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_os_trigpass]->eta - (*selectedTracks_dmin)[pindex_os_trigpass]->eta));
//          Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTracks_dmin)[nindex_os_trigpass]->phi, (*selectedTracks_dmin)[pindex_os_trigpass]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_os_trigpass]->z0 - (*selectedTracks_dmin)[pindex_os_trigpass]->z0));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_os_trigpass]->d0 - (*selectedTracks_dmin)[pindex_os_trigpass]->d0));
//          Hists[0][h_index++]->Fill((deltaR((*selectedTracks_dmin)[nindex_os_trigpass]->eta,(*selectedTracks_dmin)[nindex_os_trigpass]->phi, (*selectedTracks_dmin)[pindex_os_trigpass]->eta,(*selectedTracks_dmin)[pindex_os_trigpass]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_os_trigpass]->dxy - (*selectedTracks_dmin)[pindex_os_trigpass]->dxy));

//          if(VERBOSE[0])
//             std::cout << "OS tracks Histograms filled"<<endl;
         
//          Hists[0][h_index++]->Fill(trk_matchtp_pt->at((*selectedTracks_dmin)[nindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_eta->at((*selectedTracks_dmin)[nindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_phi->at((*selectedTracks_dmin)[nindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_z0->at((*selectedTracks_dmin)[nindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[nindex_os_trigpass]->index)));
//          Hists[0][h_index++]->Fill(trk_matchtp_dxy->at((*selectedTracks_dmin)[nindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_pt->at((*selectedTracks_dmin)[pindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_eta->at((*selectedTracks_dmin)[pindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_phi->at((*selectedTracks_dmin)[pindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(trk_matchtp_z0->at((*selectedTracks_dmin)[pindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[pindex_os_trigpass]->index)));
//          Hists[0][h_index++]->Fill(trk_matchtp_dxy->at((*selectedTracks_dmin)[pindex_os_trigpass]->index));
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_eta->at((*selectedTracks_dmin)[nindex_os_trigpass]->index) - trk_matchtp_eta->at((*selectedTracks_dmin)[pindex_os_trigpass]->index)));
//          Hists[0][h_index++]->Fill(fabs(deltaPhi(trk_matchtp_phi->at((*selectedTracks_dmin)[nindex_os_trigpass]->index), trk_matchtp_phi->at((*selectedTracks_dmin)[pindex_os_trigpass]->index))));
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_z0->at((*selectedTracks_dmin)[nindex_os_trigpass]->index) - trk_matchtp_z0->at((*selectedTracks_dmin)[pindex_os_trigpass]->index)));
//          Hists[0][h_index++]->Fill((deltaR(trk_matchtp_eta->at((*selectedTracks_dmin)[nindex_os_trigpass]->index),trk_matchtp_phi->at((*selectedTracks_dmin)[nindex_os_trigpass]->index), trk_matchtp_eta->at((*selectedTracks_dmin)[pindex_os_trigpass]->index),trk_matchtp_phi->at((*selectedTracks_dmin)[pindex_os_trigpass]->index))));
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_dxy->at((*selectedTracks_dmin)[nindex_os_trigpass]->index) - trk_matchtp_dxy->at((*selectedTracks_dmin)[pindex_os_trigpass]->index)));
//          }

//          if(VERBOSE[0])
//             std::cout << "Matched OS tracks Histograms filled"<<endl;

//          if (triggerpass_mu){
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->charge);
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[nindex_mu_trigpass]->index)));
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[nindex_mu_trigpass]->dxy);
//          //std::cout << "Negative charge"<<endl;
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->pt);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->eta);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->phi);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->z0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->d0);
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->charge);
//          Hists[0][h_index++]->Fill(fabs(trk_matchtp_pdgid->at((*selectedTracks_dmin)[pindex_mu_trigpass]->index)));
//          Hists[0][h_index++]->Fill((*selectedTracks_dmin)[pindex_mu_trigpass]->dxy);
//          // std::cout << "Positive charge"<<endl;
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_mu_trigpass]->eta - (*selectedTracks_dmin)[pindex_mu_trigpass]->eta));
//          Hists[0][h_index++]->Fill(fabs(deltaPhi((*selectedTracks_dmin)[nindex_mu_trigpass]->phi, (*selectedTracks_dmin)[pindex_mu_trigpass]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_mu_trigpass]->z0 - (*selectedTracks_dmin)[pindex_mu_trigpass]->z0));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_mu_trigpass]->d0 - (*selectedTracks_dmin)[pindex_mu_trigpass]->d0));
//          Hists[0][h_index++]->Fill((deltaR((*selectedTracks_dmin)[nindex_mu_trigpass]->eta,(*selectedTracks_dmin)[nindex_mu_trigpass]->phi, (*selectedTracks_dmin)[pindex_mu_trigpass]->eta,(*selectedTracks_dmin)[pindex_mu_trigpass]->phi)));
//          Hists[0][h_index++]->Fill(fabs((*selectedTracks_dmin)[nindex_mu_trigpass]->dxy - (*selectedTracks_dmin)[pindex_mu_trigpass]->dxy));
//          }

//          if(VERBOSE[0])
//             std::cout << "First Region Histogram Filled"<<endl;

//       for (int l = 0; l < selectedTracks->size(); l++)
//       {
//          delete (*selectedTracks)[l];
//       }
//       for (int l = 0; l < selectedTracks_dmin->size(); l++)
//       {
//          delete (*selectedTracks_dmin)[l];
//       }
//       for (int l = 0; l < selectedTPs->size(); l++)
//       {
//          delete (*selectedTPs)[l];
//       }
//       for (int l = 0; l < selectedTPs_dmin->size(); l++)
//       {
//          delete (*selectedTPs_dmin)[l];
//       }

//       selectedTracks->clear();
//       selectedTracks->shrink_to_fit();
//       delete selectedTracks;
//       selectedTracks_dmin->clear();
//       selectedTracks_dmin->shrink_to_fit();
//       delete selectedTracks_dmin;
//       selectedTPs->clear();
//       selectedTPs->shrink_to_fit();
//       delete selectedTPs;
//       selectedTPs_dmin->clear();
//       selectedTPs_dmin->shrink_to_fit();
//       delete selectedTPs_dmin;
//       // cout<<endl;

//    } // end event loop
//    cout << endl
//         << "from " << nevt << " events, " << nAccept << " events are accepted"<<endl;

//    cout << "from " << nevt << " events, " << nCount_OS << " events have 1 pair of opposite sign tracks" << endl;
//    cout << "from " << nevt << " events, " << nCount_2muOS << " events have atleast 1 pair or opposite sign 2 Muons" << endl;
//    // ---------------------------------------------------------------------------------------------------------
//    //some Histograms

//    char ctxt[500];
//    if(type.Contains("cT0")){
//       sprintf(ctxt, "Dark Photon, PU=0, #tau=0mm");
//    }
//    else if(type.Contains("cT10000")){
//       sprintf(ctxt, "Dark Photon, PU=0, #tau=10000mm");
//    }
//    else if(type.Contains("cT5000")){
//       sprintf(ctxt, "Dark Photon, PU=0, #tau=5000mm");
//    }   
//    else if(type.Contains("cT100")){
//       sprintf(ctxt, "Dark Photon, PU=0, #tau=100mm");
//    }
//    else if(type.Contains("cT10")){
//       if(type.Contains("PU200")){
//          sprintf(ctxt, "Dark Photon, PU=200, #tau=10mm");
//       }
//       else{
//          sprintf(ctxt, "Dark Photon, PU=0, #tau=10mm");
//       }
//    }
//    else if(type.Contains("NeutrinoGun")){
//       sprintf(ctxt, "Neutrino Gun, PU=200");
//    }
//    else if(type.Contains("DisplacedMu")){
//       if(type.Contains("PU200")){
//          sprintf(ctxt, "Displaced Mu, PU=200");
//       }
//       else{
//          sprintf(ctxt, "Displaced Mu, PU=0");
//       }
//    }
//    else{
//       sprintf(ctxt, "");
//    }
//    TCanvas c;

//    TString DIR = type_dir + "AnalyzerTrkPlots/";
//    TString makedir = "mkdir -p " + DIR;
//    const char *mkDIR = makedir.Data();
//    gSystem->Exec(mkDIR);

//    TFile *fout;
//    fout = new TFile(type_dir + "output_" + type + ".root", "recreate");

//    for (int k = 0; k < regions.size(); ++k)
//    {
//       makedir = "mkdir -p " + DIR + regions[k] + "/";
//       gSystem->Exec(makedir.Data());
//       for (int l = 0; l < vars.size(); ++l)
//       {
//          Hists[k][l]->GetXaxis()->SetRange(1, Hists[k][l]->GetNbinsX() + 2);
//          Hists[k][l]->Draw();
//          mySmallText(0.4, 0.82, 1, ctxt);
//          Hists[k][l]->GetYaxis()->SetTitle("Events");/*
//          if (vars[l].Contains("pt"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("Pt");
//          }
//          else if (vars[l].Contains("eta"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("Eta");
//          }
//          else if (vars[l].Contains("phi"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("Phi");
//          }
//          else if (vars[l].Contains("z0"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("Z0");
//          }
//          else if (vars[l].Contains("d0"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("d0");
//          }
//          else if (vars[l].Contains("charge"))
//          {
//             Hists[k][l]->GetXaxis()->SetTitle("charge");
//          }*/
//          Hists[k][l]->GetXaxis()->SetTitle(Hists[k][l]->GetName());
//          Hists[k][l]->Write("", TObject::kOverwrite);
//          c.SaveAs(DIR + regions[k] + "/" + type + "_" + Hists[k][l]->GetName() + ".jpeg");
//       }
//    }
//    for (int k = 0; k < regions.size(); ++k)
//    {
//       for (int l = 0; l < vars.size(); ++l)
//       {
//          delete Hists[k][l];
//       }
//    }

//    char res[1000];
//    float rms = 0;
//    TF1* fit;
//    fit = new TF1("fit", "gaus", -1, 1);
//    h_res_tp_trk_x->Fit("fit", "0");
//    h_res_tp_trk_x->Draw();
//    rms = fit->GetParameter(2);
//    sprintf(res, "RMS = %.4f", rms);
//    mySmallText(0.22, 0.82, 1, res);
//    mySmallText(0.4, 0.42, 1, ctxt);
//    h_res_tp_trk_x->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_res_tp_trk_x->GetName() + ".jpeg");
//    delete h_res_tp_trk_x;
//    delete fit;

//    fit = new TF1("fit", "gaus", -1, 1);
//    h_res_tp_trk_y->Fit("fit", "0");
//    h_res_tp_trk_y->Draw();
//    rms = fit->GetParameter(2);
//    sprintf(res, "RMS = %.4f", rms);
//    mySmallText(0.22, 0.82, 1, res);
//    mySmallText(0.4, 0.42, 1, ctxt);
//    h_res_tp_trk_y->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_res_tp_trk_y->GetName() + ".jpeg");
//    delete h_res_tp_trk_y;
//    delete fit;

//    fit = new TF1("fit", "gaus", -10, 10);
//    h_res_tp_trk_z->Fit("fit", "0");
//    h_res_tp_trk_z->Draw();
//    rms = fit->GetParameter(2);
//    sprintf(res, "RMS = %.4f", rms);
//    mySmallText(0.22, 0.82, 1, res);
//    mySmallText(0.4, 0.42, 1, ctxt);
//    h_res_tp_trk_z->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_res_tp_trk_z->GetName() + ".jpeg");
//    delete h_res_tp_trk_z;
//    delete fit;

//    char binlabel[1000];
//    sprintf(binlabel, "0 + 1,   0 + 2,   1 + 2,   0 + 3,   1 + 3,   2 + 3");

//    h_all_trk_pt->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_all_trk_pt->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_all_trk_pt->GetName() + ".jpeg");

//    h_correct_trk_pt->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_correct_trk_pt->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_correct_trk_pt->GetName() + ".jpeg");
   
//    h_all_trk_pt->Sumw2();
//    h_correct_trk_pt->Sumw2();
//    h_eff_trk_pt->Divide(h_correct_trk_pt,h_all_trk_pt, 1.0, 1.0, "B");
//    h_eff_trk_pt->Draw();
//    h_eff_trk_pt->SetAxisRange(0, 1.1, "Y");
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_eff_trk_pt->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_eff_trk_pt->GetName() + ".jpeg");
//    delete h_eff_trk_pt;
//    delete h_all_trk_pt;
//    delete h_correct_trk_pt;

//    h_all_trk_eta->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_all_trk_eta->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_all_trk_eta->GetName() + ".jpeg");

//    h_correct_trk_eta->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_correct_trk_eta->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_correct_trk_eta->GetName() + ".jpeg");

//    h_all_trk_eta->Sumw2();
//    h_correct_trk_eta->Sumw2();
//    h_eff_trk_eta->Divide(h_correct_trk_eta,h_all_trk_eta, 1.0, 1.0, "B");
//    h_eff_trk_eta->Draw();
//    h_eff_trk_eta->SetAxisRange(0, 1.1, "Y");
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_eff_trk_eta->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_eff_trk_eta->GetName() + ".jpeg");
//    delete h_eff_trk_eta;
//    delete h_all_trk_eta;
//    delete h_correct_trk_eta;

//    h_all_trk_dxy->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_all_trk_dxy->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_all_trk_dxy->GetName() + ".jpeg");

//    h_correct_trk_dxy->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_correct_trk_dxy->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_correct_trk_dxy->GetName() + ".jpeg");
   
//    h_all_trk_dxy->Sumw2();
//    h_correct_trk_dxy->Sumw2();
//    h_eff_trk_dxy->Divide(h_correct_trk_dxy,h_all_trk_dxy, 1.0, 1.0, "B");
//    h_eff_trk_dxy->Draw();
//    h_eff_trk_dxy->SetAxisRange(0, 1.1, "Y");
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_eff_trk_dxy->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_eff_trk_dxy->GetName() + ".jpeg");
//    delete h_eff_trk_dxy;
//    delete h_all_trk_dxy;
//    delete h_correct_trk_dxy;


//    h_Counter_TPcombination->GetYaxis()->SetNoExponent(kTRUE);
//    h_Counter_TPcombination->GetXaxis()->SetRange(0, h_Counter_TPcombination->GetNbinsX() + 1);
//    h_Counter_TPcombination->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    mySmallText(0.25, 0.95, 1, binlabel);
//    h_Counter_TPcombination->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_Counter_TPcombination->GetName() + ".jpeg");
//    delete h_Counter_TPcombination;

//    // h_d0_Mu_n_Mu_p->GetXaxis()->SetRange(1, h_d0_Mu_n_Mu_p->GetNbinsX() + 2);
//    h_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
//    h_delta_dist_xy->GetXaxis()->SetRange(1, h_delta_dist_xy->GetNbinsX() + 2);
//    h_delta_dist_xy->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_delta_dist_xy->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_delta_dist_xy->GetName() + ".jpeg");
//    delete h_delta_dist_xy;
//    h_error_delta_x->GetXaxis()->SetRange(1, h_error_delta_x->GetNbinsX() + 2);
//    h_error_delta_x->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_error_delta_x->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_error_delta_x->GetName() + ".jpeg");
//    delete h_error_delta_x;
//    h_delta_x->GetXaxis()->SetRange(1, h_delta_x->GetNbinsX() + 2);
//    h_delta_x->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_delta_x->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_delta_x->GetName() + ".jpeg");
//    delete h_delta_x;
//    h_delta_dist_z->GetXaxis()->SetRange(1, h_delta_dist_z->GetNbinsX() + 2);
//    h_delta_dist_z->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_delta_dist_z->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_delta_dist_z->GetName() + ".jpeg");
//    delete h_delta_dist_z;
//    h_error_delta_z->GetXaxis()->SetRange(1, h_error_delta_z->GetNbinsX() + 2);
//    h_error_delta_z->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_error_delta_z->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_error_delta_z->GetName() + ".jpeg");
//    delete h_error_delta_z;
//    // gStyle->SetPadRightMargin(0.2);
//    h_d0_Mu_n_Mu_p->Draw("colz");
//    h_d0_Mu_n_Mu_p->GetYaxis()->SetTitle("d0 Mu+");
//    h_d0_Mu_n_Mu_p->GetXaxis()->SetTitle("d0 Mu-");
//    h_d0_Mu_n_Mu_p->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_d0_Mu_n_Mu_p->GetName() + ".jpeg");
//    delete h_d0_Mu_n_Mu_p;
//    h_displaced_vertex_x->GetXaxis()->SetRange(0, h_displaced_vertex_x->GetNbinsX() + 2);
//    h_displaced_vertex_x->GetYaxis()->SetRange(0, h_displaced_vertex_x->GetNbinsX() + 2);
//    h_displaced_vertex_x->Draw("colz");
//    h_displaced_vertex_x->GetYaxis()->SetTitle("Track Displaced vertex - x");
//    h_displaced_vertex_x->GetXaxis()->SetTitle("Tracking Particle Displaced vertex - x");
//    h_displaced_vertex_x->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_displaced_vertex_x->GetName() + ".jpeg");
//    delete h_displaced_vertex_x;
//    h_displaced_vertex_y->GetXaxis()->SetRange(0, h_displaced_vertex_y->GetNbinsX() + 2);
//    h_displaced_vertex_y->GetYaxis()->SetRange(0, h_displaced_vertex_y->GetNbinsX() + 2);
//    h_displaced_vertex_y->Draw("colz");
//    h_displaced_vertex_y->GetYaxis()->SetTitle("Track Displaced vertex - y");
//    h_displaced_vertex_y->GetXaxis()->SetTitle("Tracking Particle Displaced vertex - y");
//    h_displaced_vertex_y->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_displaced_vertex_y->GetName() + ".jpeg");
//    delete h_displaced_vertex_y;
//    h_displaced_vertex_z->GetXaxis()->SetRange(0, h_displaced_vertex_z->GetNbinsX() + 2);
//    h_displaced_vertex_z->GetYaxis()->SetRange(0, h_displaced_vertex_z->GetNbinsX() + 2);
//    h_displaced_vertex_z->Draw("colz");
//    h_displaced_vertex_z->GetYaxis()->SetTitle("Track Displaced vertex - z");
//    h_displaced_vertex_z->GetXaxis()->SetTitle("Tracking Particle Displaced vertex - z");
//    h_displaced_vertex_z->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_displaced_vertex_z->GetName() + ".jpeg");
//    delete h_displaced_vertex_z;

//    h_trk_Counter_TPcombination->GetYaxis()->SetNoExponent(kTRUE);
//    h_trk_Counter_TPcombination->GetXaxis()->SetRange(0, h_trk_Counter_TPcombination->GetNbinsX() + 1);
//    h_trk_Counter_TPcombination->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    mySmallText(0.25, 0.95, 1, binlabel);
//    h_trk_Counter_TPcombination->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_Counter_TPcombination->GetName() + ".jpeg");
//    delete h_trk_Counter_TPcombination;

//    // h_trk_d0_Mu_n_Mu_p->GetXaxis()->SetRange(1, h_trk_d0_Mu_n_Mu_p->GetNbinsX() + 2);
//    h_trk_delta_dist_xy->GetYaxis()->SetNoExponent(kTRUE);
//    h_trk_delta_dist_xy->GetXaxis()->SetRange(1, h_trk_delta_dist_xy->GetNbinsX() + 2);
//    h_trk_delta_dist_xy->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_delta_dist_xy->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_delta_dist_xy->GetName() + ".jpeg");
//    delete h_trk_delta_dist_xy;
//    // h_trk_error_delta_x->GetXaxis()->SetRange(1, h_trk_error_delta_x->GetNbinsX() + 2);
//    // h_trk_error_delta_x->Draw();
//    // mySmallText(0.4, 0.82, 1, ctxt);
//    // h_trk_error_delta_x->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_error_delta_x->GetName() + ".jpeg");
//    delete h_trk_error_delta_x;
//    // h_trk_delta_x->GetXaxis()->SetRange(1, h_trk_delta_x->GetNbinsX() + 2);
//    // h_trk_delta_x->Draw();
//    // mySmallText(0.4, 0.82, 1, ctxt);
//    // h_trk_delta_x->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_delta_x->GetName() + ".jpeg");
//    delete h_trk_delta_x;
//    h_trk_delta_dist_z->GetXaxis()->SetRange(1, h_trk_delta_dist_z->GetNbinsX() + 2);
//    h_trk_delta_dist_z->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_delta_dist_z->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_delta_dist_z->GetName() + ".jpeg");
//    delete h_trk_delta_dist_z;
//    // h_trk_error_delta_z->GetXaxis()->SetRange(1, h_trk_error_delta_z->GetNbinsX() + 2);
//    // h_trk_error_delta_z->Draw();
//    // mySmallText(0.4, 0.82, 1, ctxt);
//    // h_trk_error_delta_z->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_error_delta_z->GetName() + ".jpeg");
//    delete h_trk_error_delta_z;
//    // // gStyle->SetPadRightMargin(0.2);
//    // h_trk_d0_Mu_n_Mu_p->Draw("colz");
//    // h_trk_d0_Mu_n_Mu_p->GetYaxis()->SetTitle("d0 Mu+");
//    // h_trk_d0_Mu_n_Mu_p->GetXaxis()->SetTitle("d0 Mu-");
//    // h_trk_d0_Mu_n_Mu_p->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_d0_Mu_n_Mu_p->GetName() + ".jpeg");
//    delete h_trk_d0_Mu_n_Mu_p;
//    // h_trk_displaced_vertex_x->GetXaxis()->SetRange(0, h_trk_displaced_vertex_x->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_x->GetYaxis()->SetRange(0, h_trk_displaced_vertex_x->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_x->Draw("colz");
//    // h_trk_displaced_vertex_x->GetYaxis()->SetTitle("Leading p_t - x");
//    // h_trk_displaced_vertex_x->GetXaxis()->SetTitle("Displaced vertex - x");
//    // h_trk_displaced_vertex_x->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_displaced_vertex_x->GetName() + ".jpeg");
//    delete h_trk_displaced_vertex_x;
//    // h_trk_displaced_vertex_y->GetXaxis()->SetRange(0, h_trk_displaced_vertex_y->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_y->GetYaxis()->SetRange(0, h_trk_displaced_vertex_y->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_y->Draw("colz");
//    // h_trk_displaced_vertex_y->GetYaxis()->SetTitle("Leading p_t - y");
//    // h_trk_displaced_vertex_y->GetXaxis()->SetTitle("Displaced vertex - y");
//    // h_trk_displaced_vertex_y->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_displaced_vertex_y->GetName() + ".jpeg");
//    delete h_trk_displaced_vertex_y;
//    // h_trk_displaced_vertex_z->GetXaxis()->SetRange(0, h_trk_displaced_vertex_z->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_z->GetYaxis()->SetRange(0, h_trk_displaced_vertex_z->GetNbinsX() + 2);
//    // h_trk_displaced_vertex_z->Draw("colz");
//    // h_trk_displaced_vertex_z->GetYaxis()->SetTitle("Leading p_t - z");
//    // h_trk_displaced_vertex_z->GetXaxis()->SetTitle("Displaced vertex - z");
//    // h_trk_displaced_vertex_z->Write("", TObject::kOverwrite);
//    // c.SaveAs(DIR + "/"+ h_trk_displaced_vertex_z->GetName() + ".jpeg");
//    delete h_trk_displaced_vertex_z;

//    h_tp_xy_cond1->GetXaxis()->SetRange(1, h_tp_xy_cond1->GetNbinsX() + 2);
//    h_tp_xy_cond1->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_tp_xy_cond1->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_tp_xy_cond1->GetName() + ".jpeg");
//    delete h_tp_xy_cond1;
//    h_tp_xy_cond2->GetXaxis()->SetRange(1, h_tp_xy_cond2->GetNbinsX() + 2);
//    h_tp_xy_cond2->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_tp_xy_cond2->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_tp_xy_cond2->GetName() + ".jpeg");
//    delete h_tp_xy_cond2;
//    h_tp_z_cond1->GetXaxis()->SetRange(1, h_tp_z_cond1->GetNbinsX() + 2);
//    h_tp_z_cond1->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_tp_z_cond1->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_tp_z_cond1->GetName() + ".jpeg");
//    delete h_tp_z_cond1;
//    h_tp_z_cond2->GetXaxis()->SetRange(1, h_tp_z_cond2->GetNbinsX() + 2);
//    h_tp_z_cond2->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_tp_z_cond2->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_tp_z_cond2->GetName() + ".jpeg");
//    delete h_tp_z_cond2;
   
//    h_trk_xy_cond1->GetXaxis()->SetRange(1, h_trk_xy_cond1->GetNbinsX() + 2);
//    h_trk_xy_cond1->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_xy_cond1->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_xy_cond1->GetName() + ".jpeg");
//    delete h_trk_xy_cond1;
//    h_trk_xy_cond2->GetXaxis()->SetRange(1, h_trk_xy_cond2->GetNbinsX() + 2);
//    h_trk_xy_cond2->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_xy_cond2->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_xy_cond2->GetName() + ".jpeg");
//    delete h_trk_xy_cond2;
//    h_trk_z_cond1->GetXaxis()->SetRange(1, h_trk_z_cond1->GetNbinsX() + 2);
//    h_trk_z_cond1->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_z_cond1->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_z_cond1->GetName() + ".jpeg");
//    delete h_trk_z_cond1;
//    h_trk_z_cond2->GetXaxis()->SetRange(1, h_trk_z_cond2->GetNbinsX() + 2);
//    h_trk_z_cond2->Draw();
//    mySmallText(0.4, 0.82, 1, ctxt);
//    h_trk_z_cond2->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_trk_z_cond2->GetName() + ".jpeg");
//    delete h_trk_z_cond2;

//    // h_Count_trk_pt_d0->GetXaxis()->SetRange(0, h_Count_trk_pt_d0->GetNbinsX() + 2);
//    // h_Count_trk_pt_d0->GetYaxis()->SetRange(0, h_Count_trk_pt_d0->GetNbinsX() + 2);
   
//    std::stringstream txt;
//    float rate;
//    if(type.Contains("NeutrinoGun")){
//       rate = 40000.0;
//       txt << " kHz ";
//    }
//    else{
//       rate = 100.0;
//       txt << " % ";
//    }

//    for(int k=0;k<5;k++){
//       for (int l=0; l < 5; l++){
//          h_Count_trk_pt_d0->SetBinContent(l+1,k+1,(h_Count_trk_pt_d0->GetBinContent(l+1,k+1) * rate/nevt));
//       }
//    }
   
//    h_Count_trk_pt_d0->Draw("colz");
//    h_Count_trk_pt_d0->SetMarkerSize(2);
//    h_Count_trk_pt_d0->Draw("textsame");
//    h_Count_trk_pt_d0->GetXaxis()->SetNdivisions(5);
//    h_Count_trk_pt_d0->GetYaxis()->SetNdivisions(5);
//    h_Count_trk_pt_d0->GetXaxis()->CenterLabels();
//    h_Count_trk_pt_d0->GetYaxis()->CenterLabels();
//    h_Count_trk_pt_d0->GetYaxis()->SetTitle("Transverse Momentum p_T (GeV)");
//    h_Count_trk_pt_d0->GetXaxis()->SetTitle("Transverse Impact Parameter d_0 (cm)");
//    h_Count_trk_pt_d0->Write("", TObject::kOverwrite);
//    c.SaveAs(DIR + "/"+ h_Count_trk_pt_d0->GetName() + ".pdf");
//    delete h_Count_trk_pt_d0;


//    fout->Close();

//    ofstream myfile;
//    myfile.open(type_dir+type+"_"+"Counters.txt");

//    myfile<<endl<<"Total Number of Events = "<<nevt;
//    myfile << endl;
//    myfile << endl << "Tracking Particles"<<endl; 
//    myfile<<endl<<"Number of Events with atleast 2 tracking particles = "<<nAccept;
//    myfile<<endl<<"Number of Events with Exactly 2 tracking particles = "<<nCount_2tp;
//    myfile<<endl<<"Number of Events with atleast 1 pair of opposite sign tracks = "<<nCount_OS;
//    myfile<<endl<<"Number of Events with atleast 1 pair of muons of opposite sign = "<<nCount_2muOS;
//    myfile<<endl;
//    myfile<<endl<<"Number of Events with leading pT TPs from same true vertex = "<<nCount_leadingpt;
//    myfile<<endl<<"Number of Events with selected z as true z = "<<nCount_z;
//    myfile<<endl<<"Number of Events with both TPs within 1.0cm of DV = "<<nCount_xyz;
//    myfile<<endl<<"Number of Events with no intersection in x-y = "<<nCount_no_intersect;
//    myfile<<endl<<"Number of Events with intersection in x-y = "<<nCount_intersect;
//    myfile<<endl<<"Number of Events with intersection in x-y but not in z = "<<nCount_no_intersect_z;
//    myfile<<endl<<"Number of Events with only 1 intersection = "<<nCount_1intersect;
//    myfile<<endl<<"Number of Events with 2 intersections but only 1 detectable = "<<nCount_only1accept;
//    myfile<<endl<<"Number of Events with 2 intersections but both detectable = "<<nCount_accept_both;
//    myfile<<endl<<"Number of Events with both identified TPs matched to tracks = "<<nCount_matchtrk;
//    myfile<<endl;
   
//    myfile<<endl<<"Number of Events with difference in x_dv for TP and true value < 1.0 = "<<nCount_tp_x;
//    myfile<<endl<<"Number of Events with difference in y_dv for TP and true value < 1.0 = "<<nCount_tp_y;
//    myfile<<endl<<"Number of Events with difference in z_dv for TP and true value < 1.0 = "<<nCount_tp_z;

//    myfile<<endl;

//    for(int k=0; k<5 ;k++){
//       myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with "<<"Z0 < 0.5 and 1 pair of Tracking Particles with Pt greater than "<<pt_cuts[k]<<" = "<<nCount_tp_pt[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_pt[k]<<")";
//    }
//    for(int k=0; k<5 ;k++){
//       myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with "<<"Z0 < 0.5 and 1 pair of Tracking Particles with D0 greater than "<<d0_cuts[k]<<" = "<<nCount_tp_d0[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_d0[k]<<")";
//    }
//    // myfile<<endl<<"Track Rate with "<< (txt.str()).c_str()<<endl;
//    myfile<<"\n\n\n";
//    int tablesize = 10;
//    for (int k=-1; k<5; k++){
//       if(k==-1){
//          myfile << std::right << std::setw(tablesize*2)<<" d0 \\ Pt ";
//          for (int l=0; l<5; l++){
//             myfile << std::right << std::setw(tablesize*2)<< pt_cuts[l];
//          }
//       }
//       else{
//          for (int l=-1; l<5; l++){
//             if(l==-1){ 
//                myfile << std::right << std::setw(tablesize*2) << d0_cuts[k];
//             }
//             else{
//             myfile <<std::fixed<<std::setprecision(3)<< std::right<< std::setw(tablesize) << nCount_tp_pt_d0[l][k]* rate/nevt<<(txt.str()).c_str()<<"("<<nCount_tp_pt_d0[l][k]<<")";
//             }
//          }
//       }
//       myfile<<endl;
//    }

//    myfile << endl;
//    myfile << endl << "Tracks"<<endl; 
//    myfile<<endl<<"Number of Events with Atleast 2 tracks = "<<nCount_2trk;
//    myfile<<endl;
//    myfile<<endl<<"Number of Events with both tracks within 1.0cm of DV = "<<nCount_trk_xyz;
//    myfile<<endl<<"Number of Events with no intersection in x-y = "<<nCount_trk_no_intersect;
//    myfile<<endl<<"Number of Events with intersection in x-y = "<<nCount_trk_intersect;
//    myfile<<endl<<"Number of Events with intersection in x-y but not in z = "<<nCount_trk_no_intersect_z;
//    myfile<<endl<<"Number of Events with only 1 intersection = "<<nCount_trk_1intersect;
//    myfile<<endl<<"Number of Events with 2 intersections but only 1 detectable = "<<nCount_trk_only1accept;
//    myfile<<endl<<"Number of Events with 2 intersections but both detectable = "<<nCount_trk_accept_both;
//    myfile<<endl;

//    myfile<<endl<<"Number of Events with difference in x_dv for TP and tracks < 1.0 = "<<nCount_tp_trk_x<< "("<<(float)nCount_tp_trk_x*100.0/nevt<<")";
//    myfile<<endl<<"Number of Events with difference in y_dv for TP and tracks < 1.0 = "<<nCount_tp_trk_y<< "("<<(float)nCount_tp_trk_y*100.0/nevt<<")";
//    myfile<<endl<<"Number of Events with difference in z_dv for TP and tracks < 1.0 = "<<nCount_tp_trk_z<< "("<<(float)nCount_tp_trk_z*100.0/nevt<<")";

//    myfile<<endl;

//    for(int k=0; k<5 ;k++){
//       myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with chi2rz_dof<2, chi2rphi_dof<2, bendchi2<5, "<<"Z0 < 0.5 and 1 pair of Tracks with Pt greater than "<<pt_cuts[k]<<" = "<<nCount_trk_pt[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_pt[k]<<")";
//    }
//    for(int k=0; k<5 ;k++){
//       myfile<<endl<<std::fixed<<std::setprecision(3)<<"Track Rate with chi2rz_dof<2, chi2rphi_dof<2, bendchi2<5, "<<"Z0 < 0.5 and 1 pair of Tracks with D0 greater than "<<d0_cuts[k]<<" = "<<nCount_trk_d0[k] * rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_d0[k]<<")";
//    }
//    myfile<<"\n\n\n";
//    for (int k=-1; k<5; k++){
//       if(k==-1){
//          myfile << std::right << std::setw(tablesize*2) <<" d0 \\ Pt ";
//          for (int l=0; l<5; l++){
//             myfile << std::right<< std::setw(tablesize*2) << pt_cuts[l];
//          }
//       }
//       else{
//          for (int l=-1; l<5; l++){
//             if(l==-1){ 
//                myfile <<std::fixed<<std::setprecision(3)<< std::right << std::setw(tablesize*2) << d0_cuts[k];
//             }
//             else{
//             myfile <<std::fixed<<std::setprecision(3)<< std::right<< std::setw(tablesize) <<nCount_trk_pt_d0[l][k]* rate/nevt<<(txt.str()).c_str()<<"("<<nCount_trk_pt_d0[l][k]<<")";
//             }
//          }
//       }
//       myfile<<endl;
//    }

//    myfile.close();
// }

// void SetPlotStyle()
// {
//    // from ATLAS plot style macro

//    // use plain black on white colors
//    gStyle->SetFrameBorderMode(0);
//    gStyle->SetFrameFillColor(0);
//    gStyle->SetCanvasBorderMode(0);
//    gStyle->SetCanvasColor(0);
//    gStyle->SetPadBorderMode(0);
//    gStyle->SetPadColor(0);
//    gStyle->SetStatColor(0);
//    gStyle->SetHistLineColor(1);

//    gStyle->SetPalette(1);

//    // set the paper & margin sizes
//    gStyle->SetPaperSize(20, 26);
//    gStyle->SetPadTopMargin(0.05);
//    gStyle->SetPadRightMargin(0.15);
//    gStyle->SetPadBottomMargin(0.16);
//    gStyle->SetPadLeftMargin(0.16);

//    // set title offsets (for axis label)
//    gStyle->SetTitleXOffset(1.4);
//    gStyle->SetTitleYOffset(1.4);

//    // use large fonts
//    gStyle->SetTextFont(42);
//    gStyle->SetTextSize(0.05);
//    gStyle->SetLabelFont(42, "x");
//    gStyle->SetTitleFont(42, "x");
//    gStyle->SetLabelFont(42, "y");
//    gStyle->SetTitleFont(42, "y");
//    gStyle->SetLabelFont(42, "z");
//    gStyle->SetTitleFont(42, "z");
//    gStyle->SetLabelSize(0.05, "x");
//    gStyle->SetTitleSize(0.05, "x");
//    gStyle->SetLabelSize(0.05, "y");
//    gStyle->SetTitleSize(0.05, "y");
//    gStyle->SetLabelSize(0.05, "z");
//    gStyle->SetTitleSize(0.05, "z");

//    // use bold lines and markers
//    gStyle->SetMarkerStyle(20);
//    gStyle->SetMarkerSize(1.2);
//    gStyle->SetHistLineWidth(2.);
//    gStyle->SetLineStyleString(2, "[12 12]");

//    // get rid of error bar caps
//    gStyle->SetEndErrorSize(0.);

//    // do not display any of the standard histogram decorations
//    gStyle->SetOptTitle(0);
//    gStyle->SetOptStat(0);
//    gStyle->SetOptFit(0);

//    // put tick marks on top and RHS of plots
//    gStyle->SetPadTickX(1);
//    gStyle->SetPadTickY(1);
// }

// void mySmallText(Double_t x, Double_t y, Color_t color, char *text)
// {
//    Double_t tsize = 0.044;
//    TLatex l;
//    l.SetTextSize(tsize);
//    l.SetNDC();
//    l.SetTextColor(color);
//    l.DrawLatex(x, y, text);
// }
