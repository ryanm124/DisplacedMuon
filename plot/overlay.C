#include "TROOT.h"
#include "TStyle.h"
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
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text); 


// ----------------------------------------------------------------------------------------------------------------
// Main script
void overlay(int numfiles = 2) {
  TString what;

  TString makedir = "mkdir -p Overlay/";
  const char *mkDIR = makedir.Data();
  gSystem->Exec(mkDIR);

  std::vector<TString> regions{"reg"};
  for (int k = 0; k < regions.size(); ++k)
  {
    makedir = "mkdir -p Overlay/"+ regions[k] + "/";
    gSystem->Exec(makedir.Data());
  }
  TString path;
    gROOT->SetBatch();
    gErrorIgnoreLevel = kWarning;

    SetPlotStyle();
    cout<<"Numfiles = "<<numfiles<<endl;

    TFile *tree1, *tree2, *tree3, *tree4, *tree5;

    switch(numfiles){
      default:
      case 5:tree5 = new TFile("cT10000/output_events_Dark_Photon_cT10000.root");
      case 4:tree4 = new TFile("cT5000/output_events_Dark_Photon_cT5000.root");
      case 3:tree3 = new TFile("cT100/output_events_Dark_Photon_cT100.root");
      case 2:tree2 = new TFile("cT10/output_events_Dark_Photon_cT10.root");
      case 1:tree1 = new TFile("cT10_PU200/output_events_Dark_Photon_cT10_PU200.root");
    }
    /*
    TFile *tree1 = new TFile("cT0/output_events_Dark_Photon_cT0.root");
    TFile *tree2 = new TFile("cT10/output_events_Dark_Photon_cT10.root");
    TFile *tree3 = new TFile("cT100/output_events_Dark_Photon_cT100.root");
    TFile *tree4 = new TFile("cT5000/output_events_Dark_Photon_cT5000.root");
    TFile *tree5 = new TFile("cT10000/output_events_Dark_Photon_cT10000.root");*/
    // TFile* tree4 = new TFile("output_TTbar_PU200_hybrid_displaced_500evt.root");

    TList *list = tree1->GetListOfKeys();
    TObject *obj = list->First();

    while ((obj = list->After(obj)))
    {
      what = obj->GetName();

      TH1F *h1,*h2,*h3,*h4,*h5;

      switch(numfiles){
        default:
        case 5:h5 = (TH1F *)tree5->Get(what);
        case 4:h4 = (TH1F *)tree4->Get(what);
        case 3:h3 = (TH1F *)tree3->Get(what);
        case 2:h2 = (TH1F *)tree2->Get(what);
        case 1:h1 = (TH1F *)tree1->Get(what);
      }
      /*
      TH1F *h1 = (TH1F *)tree1->Get(what);
      TH1F *h2 = (TH1F *)tree2->Get(what);
      TH1F *h3 = (TH1F *)tree3->Get(what);
      TH1F *h4 = (TH1F *)tree4->Get(what);
      TH1F *h5 = (TH1F *)tree5->Get(what);*/

      // h1->Sumw2();
      // h2->Sumw2();
      // h3->Sumw2();
      // h4->Sumw2();
      // h5->Sumw2();

      if (!what.Contains("Mu_n_Mu_p")){

      /*
      h1->GetXaxis()->SetRange(1, h1->GetNbinsX() + 2);
      h2->GetXaxis()->SetRange(1, h2->GetNbinsX() + 2);
      h3->GetXaxis()->SetRange(1, h3->GetNbinsX() + 2);
      h4->GetXaxis()->SetRange(1, h4->GetNbinsX() + 2);
      h5->GetXaxis()->SetRange(1, h5->GetNbinsX() + 2);*/
      }

      if(!what.Contains("eff")){
      switch(numfiles){
        default:
        case 5:if (h5->GetEntries()) h5->Scale(1.0 / h5->GetEntries());
        case 4:if (h4->GetEntries()) h4->Scale(1.0 / h4->GetEntries());
        case 3:if (h3->GetEntries()) h3->Scale(1.0 / h3->GetEntries());
        case 2:if (h2->GetEntries()) h2->Scale(1.0 / h2->GetEntries());
        case 1:if (h1->GetEntries()) h1->Scale(1.0 / h1->GetEntries());
      }
      }
      /*
      if (h1->GetEntries())
        h1->Scale(1.0 / h1->GetEntries());
      if (h2->GetEntries())
        h2->Scale(1.0 / h2->GetEntries());
      if (h3->GetEntries())
        h3->Scale(1.0 / h3->GetEntries());
      if (h4->GetEntries())
        h4->Scale(1.0 / h4->GetEntries());
      if (h5->GetEntries())
        h5->Scale(1.0 / h5->GetEntries());
*/
      /*
      switch(numfiles){
        default:
        case 5:h5->Rebin();
        case 4:h4->Rebin();
        case 3:h3->Rebin();
        case 2:h2->Rebin();
        case 1:h1->Rebin();
      }*/
      
      switch(numfiles){
        default:
        case 5:h5->GetXaxis()->SetRange(0, h5->GetNbinsX() + 2);
        case 4:h4->GetXaxis()->SetRange(0, h4->GetNbinsX() + 2);
        case 3:h3->GetXaxis()->SetRange(0, h3->GetNbinsX() + 2);
        case 2:h2->GetXaxis()->SetRange(0, h2->GetNbinsX() + 2);
        case 1:h1->GetXaxis()->SetRange(0, h1->GetNbinsX() + 2);
      }
  

      THStack hs("hs", what);

      TCanvas c;
      switch(numfiles){
        default:
        case 5:
                h5->SetLineColor(28);
                h5->SetMarkerColor(28);
                h5->SetMarkerStyle(29);
                hs.Add(h5);
        case 4:
                h4->SetLineColor(8);
                h4->SetMarkerColor(8);
                h4->SetMarkerStyle(23);
                hs.Add(h4);
        case 3:
                h3->SetLineColor(4);
                h3->SetMarkerColor(4);
                h3->SetMarkerStyle(22);
                hs.Add(h3);
        case 2:
                h2->SetLineColor(2);
                h2->SetMarkerColor(2);
                h2->SetMarkerStyle(24);
                hs.Add(h2);
        case 1:
                h1->SetLineColor(1);
                h1->SetMarkerColor(1);
                h1->SetMarkerStyle(8);
                hs.Add(h1);
      }
      /*
      h1->SetLineColor(1);
      h1->SetMarkerColor(1);
      h1->SetMarkerStyle(8);

      h2->SetLineColor(2);
      h2->SetMarkerColor(2);
      h2->SetMarkerStyle(24);

      h3->SetLineColor(4);
      h3->SetMarkerColor(4);
      h3->SetMarkerStyle(22);

      h4->SetLineColor(8);
      h4->SetMarkerColor(8);
      h4->SetMarkerStyle(23);

      h5->SetLineColor(28);
      h5->SetMarkerColor(28);
      h5->SetMarkerStyle(29);
      

      hs.Add(h1);
      hs.Add(h2);
      hs.Add(h3);
      //hs.Add(h4);
      //hs.Add(h5);*/
      hs.Draw("nostack");
      // hs.GetXaxis()->SetLimits(, 3);
      hs.GetXaxis()->SetRange(0, h1->GetNbinsX() + 2);

      hs.GetYaxis()->SetTitle("Normalized Events");
      
      /*  if (what.Contains("pt"))
  {
    hs.GetXaxis()->SetTitle("Pt");
  }
  else if (what.Contains("eta"))
  {
    hs.GetXaxis()->SetTitle("Eta");
  }
  else if (what.Contains("phi"))
  {
    hs.GetXaxis()->SetTitle("Phi");
  }
  else if (what.Contains("z0"))
  {
    hs.GetXaxis()->SetTitle("Z0");
  }
  else if (what.Contains("d0"))
  {
    hs.GetXaxis()->SetTitle("d0");
  }
  else if (what.Contains("rinv"))
  {
    hs.GetXaxis()->SetTitle("rInv");
  }
  else{*/
      hs.GetXaxis()->SetTitle(what);
      //}
      if(what.Contains("d0"))
        hs.GetXaxis()->SetTitle("d_{0} (cm)");
      
      if (what.Contains("eff")){
        hs.GetYaxis()->SetTitle("Efficiency");
        if (what.Contains("eta")){
          hs.GetXaxis()->SetTitle("Pseudorapidity #eta");
        }
        else if (what.Contains("pt")){
          hs.GetXaxis()->SetTitle("Transverse Momentum p_{T} (GeV)");
        }
        else if (what.Contains("dxy")){
          hs.GetXaxis()->SetTitle("Transverse displacement, dxy (cm)");
        }
      }
      else if (what.Contains("res")){
        // hs.GetYaxis()->SetTitle("Events");
        if (what.Contains("x")){
          hs.GetXaxis()->SetTitle("x resolution (cm)");
        }
        else if (what.Contains("z")){
          hs.GetXaxis()->SetTitle("z resolution (cm)");
        }
        else if (what.Contains("y")){
          hs.GetXaxis()->SetTitle("y resolution (cm)");
        }
      }

      TLegend *l;
      char sample[500];
      sprintf(sample, "Dark Photon");

      if (what.Contains("eff"))
      {
        //l = new TLegend(0.55,0.22,0.85,0.40);
        mySmallText(0.25, 0.22, 1, sample);
      }
      else
      {
        //l = new TLegend(0.2,0.72,0.5,0.9);
        mySmallText(0.4, 0.82, 1, sample);
      }
      //l = new TLegend(0.15, 0.22, "", "brNDC");
      l = new TLegend(0.45, 0.40, 0.65, 0.55,"","brNDC");
      l->SetFillColor(0);
      l->SetLineColor(0);
      l->SetTextSize(0.05);
      switch(numfiles){
        default:
        case 5:l->AddEntry(h5, "#tau 10000", "lep");
        case 4:l->AddEntry(h4, "#tau 5000", "lep");
        case 3:l->AddEntry(h1, "#tau = 100mm, PU = 0", "lep");
        case 2:l->AddEntry(h2, "#tau = 10mm, PU = 0", "lep");
        case 1:l->AddEntry(h1, "#tau = 10mm, PU = 200", "lep");
      }
      // l->AddEntry(h1, "cT0", "lep");
      // l->AddEntry(h2, "cT10", "lep");
      // l->AddEntry(h3, "cT100", "lep");
      // l->AddEntry(h4, "cT5000", "lep");
      // l->AddEntry(h5, "cT10000", "lep");
      // l->AddEntry(h1, "Seed L3L4L2", "lep");
      // l->AddEntry(h2, "Seed L5L6L4", "lep");
      // l->AddEntry(h3, "Seed L2L3D1", "lep");
      // l->AddEntry(h4, "Seed D1D2L2", "lep");

      l->SetTextFont(42);
      l->Draw();

      // gPad->SetGridy();
      path="";
      for (int k = 0; k < regions.size(); ++k)
      {

        if(what.Contains(regions[k]+"_")){
          path=regions[k];
          break;
        }
      }
      gPad->Update();
      cout<<"Saving "<<what<<".jpeg"<<endl;
       c.SaveAs("Overlay/" + path + "/overlay_" + what + ".jpeg");

        delete l;
      }
  }

void SetPlotStyle() {

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
  gStyle->SetPaperSize(20,26);
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
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

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


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}


