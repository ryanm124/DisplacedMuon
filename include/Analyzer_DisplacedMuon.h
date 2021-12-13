//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 21 23:43:52 2020 by ROOT version 6.22/00
// from TTree eventTree/Event tree
// found on file: events_Dark_Photon_cT100.root
//////////////////////////////////////////////////////////

#ifndef Analyzer_DisplacedMuon_h
#define Analyzer_DisplacedMuon_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TVector3.h>

// Header file for the classes stored in the TTree if any.
#include "TMath.h"
#include "vector"
#include "vector"
#include "vector"

using std::vector;

class Track_Parameters
{
public:
   float pt;
   float d0;
   float dxy = -99999;
   float z0;
   float eta;
   float phi;
   float charge;
   float rho;
   int index;
   int pdgid = -99999;
   float x0;
   float y0;
   float dist_calc(float x_dv, float y_dv, float x, float y){
      dxy = TMath::Sqrt((x_dv-x)*(x_dv-x) + (y_dv-y)*(y_dv-y));
   }
   float x(float phi_T=0){
      return (-charge * rho * TMath::Sin(phi - charge*phi_T) + (d0 + charge * rho) * TMath::Sin(phi));
   }
   float y(float phi_T=0){
      return ( charge * rho * TMath::Cos(phi - charge*phi_T) - (d0 + charge * rho) * TMath::Cos(phi));
   }
   float z(float phi_T=0){
      float theta = 2 * TMath::ATan(TMath::Exp(-eta));
      return (z0 + rho*phi_T/TMath::Tan(theta));
   }
   
   float deltaPhi_T(Double_t phi1, Double_t phi2)
   {
      Double_t dPhi = phi1 - phi2;
      if (dPhi >= TMath::Pi())
         dPhi -= 2. * TMath::Pi();
      if (dPhi <= -TMath::Pi())
         dPhi += 2. * TMath::Pi();
      return dPhi;
   }
   float phi_T(float x, float y){
      // if(fabs(x)<0.01 && fabs(y)<0.01) return (0); //! Understand why I get Pi instead of 0 for this case!
      float num = x - (d0 + charge * rho) * TMath::Sin(phi);
      float den = y + (d0 + charge * rho) * TMath::Cos(phi);
      // std::cout<<Form("x = %5.1f  |  y = %5.1f  |  num = %5.1f  |  den = %5.1f  |  atan = %5.1f  |  phi = %5.1f  |  dphi = %5.1f  |  phi-atan = %5.1f",x,y,num,den,TMath::ATan2(num,-den)/charge,phi/charge,deltaPhi_T(phi/charge,TMath::ATan2(num,-den)/charge),(phi - TMath::ATan2(num,-den))/charge)<<std::endl;
      return ((phi-TMath::ATan2(num,-den))/charge);   

/*
      if(fabs(1-num/den)>1.1){
         std::cout<<Form("num = %5.2f    |   den = %5.2f    |   check = %10.8f",num,den,fabs(1-num/den))<<"\n";//<<" dxy = "<<tp_dxy->at(it)<<"\t \t dist = "<<TMath::Sqrt((*selectedTPs)[j]->x*(*selectedTPs)[j]->x + (*selectedTPs)[j]->y*(*selectedTPs)[j]->y)<<" \t eta = "<<tp_eta->at(it)<<" \t phi = "<<tp_phi->at(it)<<" \t pt = "<<tp_pt->at(it)<<endl;
         return (-99999.0);
      }
      return(TMath::ACos(1-num/den)/charge);

      

      if(x==-99999){
         if(fabs((y+((d0 + charge * rho) * TMath::Cos(phi)))/(charge*rho))>1)
            return(-99999.0);
         return((phi-TMath::ACos((y+((d0 + charge * rho) * TMath::Cos(phi)))/(charge*rho)))/charge);
      }
      else if(y==-99999){
         if(fabs((x-((d0 + charge * rho) * TMath::Sin(phi)))/(-charge*rho))>1)
            return (-99999.0);
         return((phi-TMath::ASin((x-((d0 + charge * rho) * TMath::Sin(phi)))/(-charge*rho)))/charge);
      }
      else{
         return(-99999);
      }
      */
   }
   float z(float x, float y, float compare_z=9999.0){
      float orig_phiT = phi_T(x,y); 
      float temp_z = z(orig_phiT);// (z_val==9999.0)?z(phiT):z_val;
      if(fabs(temp_z)<50){
         if(compare_z==9999.0)    return (temp_z);
         else if(fabs(temp_z-compare_z)<1.0) return (temp_z);
      }
      std::vector<float> z_val;
      z_val.push_back(temp_z);
      float phiT = (TMath::Pi()/charge) + orig_phiT; 
      z_val.push_back(z(phiT));
      if(fabs(z_val.back())<50 ){
         if(compare_z==9999.0)    return (z_val.back());
         else if(fabs(z_val.back()-compare_z)<1.0) return (z_val.back());
      }
      phiT = (-TMath::Pi()/charge) + orig_phiT; 
      z_val.push_back(z(phiT));
      if(fabs(z_val.back())<50){
         if(compare_z==9999.0)    return (z_val.back());
         else if(fabs(z_val.back()-compare_z)<1.0) return (z_val.back());
      }
      /*
      for(int n=1;abs(n)<2;n>0?n++:n--){
         phiT = (n*TMath::Pi()/charge) + phi_T(x,y); 
         z_val.push_back(z(phiT));
         if(fabs(z_val.back())>fabs(temp_z) && n==1){
            n*=-1;
            n++;
         }
         else if(fabs(z_val.back())<50){
            if(compare_z==9999.0){
               return (z_val.back());
            }
            else{
               for(int j=0;j<z_val.size();j++){
                  if(fabs(z_val[j]-compare_z)<1.0){
                     return (z_val[j]);
                  }
               }
            }
         }
      }
      */
      return (temp_z);
   }
   Track_Parameters(float pt_in, float d0_in, float z0_in, float eta_in, float phi_in, float charge_in, int index_in, int pdgid_in)
   {
      pt = pt_in;
      d0 = d0_in;
      z0 = z0_in;
      eta = eta_in;
      phi = phi_in;
      if(charge_in > 0){
         charge = 1;
      }
      else if (charge_in < 0){
         charge = -1;
      }
      else{
         charge = 0;
      }
      index = index_in;
      pdgid = pdgid_in;
      rho = 100*pt / (0.3 * 3.8);
      x0 =  (d0 + charge * rho)*TMath::Sin(phi);
      y0 = -(d0 + charge * rho)*TMath::Cos(phi);
   }
   void Propagate_Transverse(){
      /*
      float theta = 2 * TMath::ATan(TMath::Exp(-eta));
      float phi_T = -z0/(rho/TMath::Tan(theta));
      x0 = -charge * rho * TMath::Sin(phi - charge*phi_T) + (d0 + charge * rho) * TMath::Sin(phi);
      y0 =  charge * rho * TMath::Cos(phi - charge*phi_T) - (d0 + charge * rho) * TMath::Cos(phi);
      
      
      rho = pt / (0.3 * 3.8);
      float x0 = -d0 * TMath::Sin(phi);
      float y0 =  d0 * TMath::Cos(phi);
      float lambda = TMath::PiOver2() - 2 * TMath::ATan(TMath::Exp(-eta));
      float s = -z0/TMath::Sin(lambda);
      
      x = x0 + rho *(TMath::Cos(phi + charge * s * TMath::Cos(lambda)/rho) - TMath::Cos(phi));
      y = y0 + rho *(TMath::Sin(phi + charge * s * TMath::Cos(lambda)/rho) - TMath::Sin(phi));
*/
      // return (phi_T);
   }
   ~Track_Parameters(){};
};
/*
void Track_Parameters::Propagate_Transverse(){
   float rho = pt / (0.3 * 3.8);
   if(charge>0){
      charge = 1;
   }
   else if (charge < 0){
      charge = -1;
   }
   float theta = 2 * TMath::ATan(TMath::Exp(-eta));
   float phi_T = -z0/(rho/TMath::Tan(theta));
   x = -charge * rho * TMath::Sin(phi - charge*phi_T) + (d0 + charge * rho) * TMath::Sin(phi);
   y =  charge * rho * TMath::Cos(phi - charge*phi_T) - (d0 + charge * rho) * TMath::Cos(phi);
}
*/
class Analyzer_DisplacedMuon {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *trk_pt;
   vector<float>   *trk_eta;
   vector<float>   *trk_phi;
   vector<float>   *trk_d0;
   vector<float>   *trk_rinv;
   vector<float>   *trk_z0;
   vector<float>   *trk_chi2;
   vector<float>   *trk_chi2rphi;
   vector<float>   *trk_chi2rz;
   vector<float>   *trk_bendchi2;
   vector<int>     *trk_nstub;
   vector<int>     *trk_lhits;
   vector<int>     *trk_dhits;
   vector<int>     *trk_seed;
   vector<int>     *trk_hitpattern;
   vector<unsigned int> *trk_phiSector;
   vector<int>     *trk_genuine;
   vector<int>     *trk_loose;
   vector<int>     *trk_unknown;
   vector<int>     *trk_combinatoric;
   vector<int>     *trk_fake;
   vector<int>     *trk_matchtp_pdgid;
   vector<float>   *trk_matchtp_pt;
   vector<float>   *trk_matchtp_eta;
   vector<float>   *trk_matchtp_phi;
   vector<float>   *trk_matchtp_z0;
   vector<float>   *trk_matchtp_dxy;
   vector<int>     *trk_injet;
   vector<int>     *trk_injet_highpt;
   vector<int>     *trk_injet_vhighpt;
   vector<float>   *tp_pt;
   vector<float>   *tp_eta;
   vector<float>   *tp_phi;
   vector<float>   *tp_dxy;
   vector<float>   *tp_d0;
   vector<float>   *tp_z0;
   vector<float>   *tp_x;
   vector<float>   *tp_y;
   vector<float>   *tp_z;
   vector<float>   *tp_d0_prod;
   vector<float>   *tp_z0_prod;
   vector<int>     *tp_pdgid;
   vector<int>     *tp_nmatch;
   vector<int>     *tp_nstub;
   vector<int>     *tp_eventid;
   vector<int>     *tp_charge;
   vector<int>     *tp_injet;
   vector<int>     *tp_injet_highpt;
   vector<int>     *tp_injet_vhighpt;
   vector<float>   *matchtrk_pt;
   vector<float>   *matchtrk_eta;
   vector<float>   *matchtrk_phi;
   vector<float>   *matchtrk_z0;
   vector<float>   *matchtrk_d0;
   vector<float>   *matchtrk_rinv;
   vector<float>   *matchtrk_chi2;
   vector<float>   *matchtrk_chi2rphi;
   vector<float>   *matchtrk_chi2rz;
   vector<float>   *matchtrk_bendchi2;
   vector<int>     *matchtrk_nstub;
   vector<int>     *matchtrk_lhits;
   vector<int>     *matchtrk_dhits;
   vector<int>     *matchtrk_seed;
   vector<int>     *matchtrk_hitpattern;
   vector<int>     *matchtrk_injet;
   vector<int>     *matchtrk_injet_highpt;
   vector<int>     *matchtrk_injet_vhighpt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_pt;
   vector<float>   *jet_tp_sumpt;
   vector<float>   *jet_trk_sumpt;
   vector<float>   *jet_matchtrk_sumpt;
 //L1Analysis::L1AnalysisGeneratorDataFormat *Generator;
   Float_t         weight;
   Float_t         pthat;
   Int_t           nVtx;
   Int_t           nMeanPU;
   Int_t           nPart;
   vector<int>     partId;
   vector<int>     partStat;
   vector<int>     partParent;
   vector<float>   partPt;
   vector<float>   partEta;
   vector<float>   partPhi;
   vector<float>   partE;
   vector<int>     partCh;
   Int_t           nJet;
   vector<float>   jetPt;
   vector<float>   jetEta;
   vector<float>   jetPhi;
   vector<float>   jetM;

   // List of branches
   TBranch        *b_trk_pt;   //!
   TBranch        *b_trk_eta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_d0;   //!
   TBranch        *b_trk_rinv;   //!
   TBranch        *b_trk_z0;   //!
   TBranch        *b_trk_chi2;   //!
   TBranch        *b_trk_chi2rphi;   //!
   TBranch        *b_trk_chi2rz;   //!
   TBranch        *b_trk_bendchi2;   //!
   TBranch        *b_trk_nstub;   //!
   TBranch        *b_trk_lhits;   //!
   TBranch        *b_trk_dhits;   //!
   TBranch        *b_trk_seed;   //!
   TBranch        *b_trk_hitpattern;   //!
   TBranch        *b_trk_phiSector;   //!
   TBranch        *b_trk_genuine;   //!
   TBranch        *b_trk_loose;   //!
   TBranch        *b_trk_unknown;   //!
   TBranch        *b_trk_combinatoric;   //!
   TBranch        *b_trk_fake;   //!
   TBranch        *b_trk_matchtp_pdgid;   //!
   TBranch        *b_trk_matchtp_pt;   //!
   TBranch        *b_trk_matchtp_eta;   //!
   TBranch        *b_trk_matchtp_phi;   //!
   TBranch        *b_trk_matchtp_z0;   //!
   TBranch        *b_trk_matchtp_dxy;   //!
   TBranch        *b_trk_injet;   //!
   TBranch        *b_trk_injet_highpt;   //!
   TBranch        *b_trk_injet_vhighpt;   //!
   TBranch        *b_tp_pt;   //!
   TBranch        *b_tp_eta;   //!
   TBranch        *b_tp_phi;   //!
   TBranch        *b_tp_dxy;   //!
   TBranch        *b_tp_d0;   //!
   TBranch        *b_tp_z0;   //!
   TBranch        *b_tp_x;   //!
   TBranch        *b_tp_y;   //!
   TBranch        *b_tp_z;   //!
   TBranch        *b_tp_d0_prod;   //!
   TBranch        *b_tp_z0_prod;   //!
   TBranch        *b_tp_pdgid;   //!
   TBranch        *b_tp_nmatch;   //!
   TBranch        *b_tp_nstub;   //!
   TBranch        *b_tp_eventid;   //!
   TBranch        *b_tp_charge;   //!
   TBranch        *b_tp_injet;   //!
   TBranch        *b_tp_injet_highpt;   //!
   TBranch        *b_tp_injet_vhighpt;   //!
   TBranch        *b_matchtrk_pt;   //!
   TBranch        *b_matchtrk_eta;   //!
   TBranch        *b_matchtrk_phi;   //!
   TBranch        *b_matchtrk_z0;   //!
   TBranch        *b_matchtrk_d0;   //!
   TBranch        *b_matchtrk_rinv;   //!
   TBranch        *b_matchtrk_chi2;   //!
   TBranch        *b_matchtrk_chi2rphi;   //!
   TBranch        *b_matchtrk_chi2rz;   //!
   TBranch        *b_matchtrk_bendchi2;   //!
   TBranch        *b_matchtrk_nstub;   //!
   TBranch        *b_matchtrk_lhits;   //!
   TBranch        *b_matchtrk_dhits;   //!
   TBranch        *b_matchtrk_seed;   //!
   TBranch        *b_matchtrk_hitpattern;   //!
   TBranch        *b_matchtrk_injet;   //!
   TBranch        *b_matchtrk_injet_highpt;   //!
   TBranch        *b_matchtrk_injet_vhighpt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_tp_sumpt;   //!
   TBranch        *b_jet_trk_sumpt;   //!
   TBranch        *b_jet_matchtrk_sumpt;   //!
   TBranch        *b_Generator_weight;   //!
   TBranch        *b_Generator_pthat;   //!
   TBranch        *b_Generator_nVtx;   //!
   TBranch        *b_Generator_nMeanPU;   //!
   TBranch        *b_Generator_nPart;   //!
   TBranch        *b_Generator_partId;   //!
   TBranch        *b_Generator_partStat;   //!
   TBranch        *b_Generator_partParent;   //!
   TBranch        *b_Generator_partPt;   //!
   TBranch        *b_Generator_partEta;   //!
   TBranch        *b_Generator_partPhi;   //!
   TBranch        *b_Generator_partE;   //!
   TBranch        *b_Generator_partCh;   //!
   TBranch        *b_Generator_nJet;   //!
   TBranch        *b_Generator_jetPt;   //!
   TBranch        *b_Generator_jetEta;   //!
   TBranch        *b_Generator_jetPhi;   //!
   TBranch        *b_Generator_jetM;   //!

   Analyzer_DisplacedMuon(TTree *tree=0);
   virtual ~Analyzer_DisplacedMuon();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString, TString, float, float, int);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyzer_DisplacedMuon_cxx
Analyzer_DisplacedMuon::Analyzer_DisplacedMuon(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("/eos/user/b/bharikri/NTuples/L1Tracking/Dark_Photon_cT10/events_Dark_Photon_cT10.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("events_Dark_Photon_cT10.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("events_Dark_Photon_cT10.root:/L1TrackNtuple");
      dir->GetObject("eventTree",tree);

   }
   Init(tree);
}

Analyzer_DisplacedMuon::~Analyzer_DisplacedMuon()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyzer_DisplacedMuon::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyzer_DisplacedMuon::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyzer_DisplacedMuon::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trk_pt = 0;
   trk_eta = 0;
   trk_phi = 0;
   trk_d0 = 0;
   trk_rinv = 0;
   trk_z0 = 0;
   trk_chi2 = 0;
   trk_chi2rphi = 0;
   trk_chi2rz = 0;
   trk_bendchi2 = 0;
   trk_nstub = 0;
   trk_lhits = 0;
   trk_dhits = 0;
   trk_seed = 0;
   trk_hitpattern = 0;
   trk_phiSector = 0;
   trk_genuine = 0;
   trk_loose = 0;
   trk_unknown = 0;
   trk_combinatoric = 0;
   trk_fake = 0;
   trk_matchtp_pdgid = 0;
   trk_matchtp_pt = 0;
   trk_matchtp_eta = 0;
   trk_matchtp_phi = 0;
   trk_matchtp_z0 = 0;
   trk_matchtp_dxy = 0;
   trk_injet = 0;
   trk_injet_highpt = 0;
   trk_injet_vhighpt = 0;
   tp_pt = 0;
   tp_eta = 0;
   tp_phi = 0;
   tp_dxy = 0;
   tp_d0 = 0;
   tp_z0 = 0;
   tp_x = 0;
   tp_y = 0;
   tp_z = 0;
   tp_d0_prod = 0;
   tp_z0_prod = 0;
   tp_pdgid = 0;
   tp_nmatch = 0;
   tp_nstub = 0;
   tp_eventid = 0;
   tp_charge = 0;
   tp_injet = 0;
   tp_injet_highpt = 0;
   tp_injet_vhighpt = 0;
   matchtrk_pt = 0;
   matchtrk_eta = 0;
   matchtrk_phi = 0;
   matchtrk_z0 = 0;
   matchtrk_d0 = 0;
   matchtrk_rinv = 0;
   matchtrk_chi2 = 0;
   matchtrk_chi2rphi = 0;
   matchtrk_chi2rz = 0;
   matchtrk_bendchi2 = 0;
   matchtrk_nstub = 0;
   matchtrk_lhits = 0;
   matchtrk_dhits = 0;
   matchtrk_seed = 0;
   matchtrk_hitpattern = 0;
   matchtrk_injet = 0;
   matchtrk_injet_highpt = 0;
   matchtrk_injet_vhighpt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_pt = 0;
   jet_tp_sumpt = 0;
   jet_trk_sumpt = 0;
   jet_matchtrk_sumpt = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
   fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_d0", &trk_d0, &b_trk_d0);
   fChain->SetBranchAddress("trk_rinv", &trk_rinv, &b_trk_rinv);
   fChain->SetBranchAddress("trk_z0", &trk_z0, &b_trk_z0);
   fChain->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
   fChain->SetBranchAddress("trk_chi2rphi", &trk_chi2rphi, &b_trk_chi2rphi);
   fChain->SetBranchAddress("trk_chi2rz", &trk_chi2rz, &b_trk_chi2rz);
   fChain->SetBranchAddress("trk_bendchi2", &trk_bendchi2, &b_trk_bendchi2);
   fChain->SetBranchAddress("trk_nstub", &trk_nstub, &b_trk_nstub);
   fChain->SetBranchAddress("trk_lhits", &trk_lhits, &b_trk_lhits);
   fChain->SetBranchAddress("trk_dhits", &trk_dhits, &b_trk_dhits);
   fChain->SetBranchAddress("trk_seed", &trk_seed, &b_trk_seed);
   fChain->SetBranchAddress("trk_hitpattern", &trk_hitpattern, &b_trk_hitpattern);
   fChain->SetBranchAddress("trk_phiSector", &trk_phiSector, &b_trk_phiSector);
   fChain->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
   fChain->SetBranchAddress("trk_loose", &trk_loose, &b_trk_loose);
   fChain->SetBranchAddress("trk_unknown", &trk_unknown, &b_trk_unknown);
   fChain->SetBranchAddress("trk_combinatoric", &trk_combinatoric, &b_trk_combinatoric);
   fChain->SetBranchAddress("trk_fake", &trk_fake, &b_trk_fake);
   fChain->SetBranchAddress("trk_matchtp_pdgid", &trk_matchtp_pdgid, &b_trk_matchtp_pdgid);
   fChain->SetBranchAddress("trk_matchtp_pt", &trk_matchtp_pt, &b_trk_matchtp_pt);
   fChain->SetBranchAddress("trk_matchtp_eta", &trk_matchtp_eta, &b_trk_matchtp_eta);
   fChain->SetBranchAddress("trk_matchtp_phi", &trk_matchtp_phi, &b_trk_matchtp_phi);
   fChain->SetBranchAddress("trk_matchtp_z0", &trk_matchtp_z0, &b_trk_matchtp_z0);
   fChain->SetBranchAddress("trk_matchtp_dxy", &trk_matchtp_dxy, &b_trk_matchtp_dxy);
   fChain->SetBranchAddress("trk_injet", &trk_injet, &b_trk_injet);
   fChain->SetBranchAddress("trk_injet_highpt", &trk_injet_highpt, &b_trk_injet_highpt);
   fChain->SetBranchAddress("trk_injet_vhighpt", &trk_injet_vhighpt, &b_trk_injet_vhighpt);
   fChain->SetBranchAddress("tp_pt", &tp_pt, &b_tp_pt);
   fChain->SetBranchAddress("tp_eta", &tp_eta, &b_tp_eta);
   fChain->SetBranchAddress("tp_phi", &tp_phi, &b_tp_phi);
   fChain->SetBranchAddress("tp_dxy", &tp_dxy, &b_tp_dxy);
   fChain->SetBranchAddress("tp_d0", &tp_d0, &b_tp_d0);
   fChain->SetBranchAddress("tp_z0", &tp_z0, &b_tp_z0);
   fChain->SetBranchAddress("tp_x", &tp_x, &b_tp_x);
   fChain->SetBranchAddress("tp_y", &tp_y, &b_tp_y);
   fChain->SetBranchAddress("tp_z", &tp_z, &b_tp_z);
   fChain->SetBranchAddress("tp_d0_prod", &tp_d0_prod, &b_tp_d0_prod);
   fChain->SetBranchAddress("tp_z0_prod", &tp_z0_prod, &b_tp_z0_prod);
   fChain->SetBranchAddress("tp_pdgid", &tp_pdgid, &b_tp_pdgid);
   fChain->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
   fChain->SetBranchAddress("tp_nstub", &tp_nstub, &b_tp_nstub);
   fChain->SetBranchAddress("tp_eventid", &tp_eventid, &b_tp_eventid);
   fChain->SetBranchAddress("tp_charge", &tp_charge, &b_tp_charge);
   fChain->SetBranchAddress("tp_injet", &tp_injet, &b_tp_injet);
   fChain->SetBranchAddress("tp_injet_highpt", &tp_injet_highpt, &b_tp_injet_highpt);
   fChain->SetBranchAddress("tp_injet_vhighpt", &tp_injet_vhighpt, &b_tp_injet_vhighpt);
   fChain->SetBranchAddress("matchtrk_pt", &matchtrk_pt, &b_matchtrk_pt);
   fChain->SetBranchAddress("matchtrk_eta", &matchtrk_eta, &b_matchtrk_eta);
   fChain->SetBranchAddress("matchtrk_phi", &matchtrk_phi, &b_matchtrk_phi);
   fChain->SetBranchAddress("matchtrk_z0", &matchtrk_z0, &b_matchtrk_z0);
   fChain->SetBranchAddress("matchtrk_d0", &matchtrk_d0, &b_matchtrk_d0);
   fChain->SetBranchAddress("matchtrk_rinv", &matchtrk_rinv, &b_matchtrk_rinv);
   fChain->SetBranchAddress("matchtrk_chi2", &matchtrk_chi2, &b_matchtrk_chi2);
   fChain->SetBranchAddress("matchtrk_chi2rphi", &matchtrk_chi2rphi, &b_matchtrk_chi2rphi);
   fChain->SetBranchAddress("matchtrk_chi2rz", &matchtrk_chi2rz, &b_matchtrk_chi2rz);
   fChain->SetBranchAddress("matchtrk_bendchi2", &matchtrk_bendchi2, &b_matchtrk_bendchi2);
   fChain->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
   fChain->SetBranchAddress("matchtrk_lhits", &matchtrk_lhits, &b_matchtrk_lhits);
   fChain->SetBranchAddress("matchtrk_dhits", &matchtrk_dhits, &b_matchtrk_dhits);
   fChain->SetBranchAddress("matchtrk_seed", &matchtrk_seed, &b_matchtrk_seed);
   fChain->SetBranchAddress("matchtrk_hitpattern", &matchtrk_hitpattern, &b_matchtrk_hitpattern);
   fChain->SetBranchAddress("matchtrk_injet", &matchtrk_injet, &b_matchtrk_injet);
   fChain->SetBranchAddress("matchtrk_injet_highpt", &matchtrk_injet_highpt, &b_matchtrk_injet_highpt);
   fChain->SetBranchAddress("matchtrk_injet_vhighpt", &matchtrk_injet_vhighpt, &b_matchtrk_injet_vhighpt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_tp_sumpt", &jet_tp_sumpt, &b_jet_tp_sumpt);
   fChain->SetBranchAddress("jet_trk_sumpt", &jet_trk_sumpt, &b_jet_trk_sumpt);
   fChain->SetBranchAddress("jet_matchtrk_sumpt", &jet_matchtrk_sumpt, &b_jet_matchtrk_sumpt);
   fChain->SetBranchAddress("weight", &weight, &b_Generator_weight);
   fChain->SetBranchAddress("pthat", &pthat, &b_Generator_pthat);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_Generator_nVtx);
   fChain->SetBranchAddress("nMeanPU", &nMeanPU, &b_Generator_nMeanPU);
   fChain->SetBranchAddress("nPart", &nPart, &b_Generator_nPart);
   fChain->SetBranchAddress("partId", &partId, &b_Generator_partId);
   fChain->SetBranchAddress("partStat", &partStat, &b_Generator_partStat);
   fChain->SetBranchAddress("partParent", &partParent, &b_Generator_partParent);
   fChain->SetBranchAddress("partPt", &partPt, &b_Generator_partPt);
   fChain->SetBranchAddress("partEta", &partEta, &b_Generator_partEta);
   fChain->SetBranchAddress("partPhi", &partPhi, &b_Generator_partPhi);
   fChain->SetBranchAddress("partE", &partE, &b_Generator_partE);
   fChain->SetBranchAddress("partCh", &partCh, &b_Generator_partCh);
   fChain->SetBranchAddress("nJet", &nJet, &b_Generator_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_Generator_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_Generator_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_Generator_jetPhi);
   fChain->SetBranchAddress("jetM", &jetM, &b_Generator_jetM);
   Notify();
}

Bool_t Analyzer_DisplacedMuon::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyzer_DisplacedMuon::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyzer_DisplacedMuon::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyzer_DisplacedMuon_cxx
