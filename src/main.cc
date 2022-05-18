
#include "Analyzer_DisplacedMuon.h"

// The Ntuples are located in /eos/user/b/bharikri/NTuples/L1Tracking/ (https://cernbox.cern.ch/index.php/s/jIGxA1m0a5BvWXE)

int main()
{
    // TString type_dir = "~/prj/Thesis/L1Tracking/DisplacedMuon/plot/";
    TString type_dir = "/eos/user/r/rmccarth/dispVert/";
    // TString file_dir = "~/prj/Thesis/Results_L1/2020_10_22_Dark_Photon_Files_GenParticles/";
    // TString file_dir = "/afs/cern.ch/user/b/bharikri/Projects/Thesis/Results_L1/2021_01_11_Dark_Photon_Files_TP_xyz/";
    //TString file_dir = "/eos/home-b/bharikri/NTuples/L1Tracking/";
    // TString file_dir = "~/prj/Thesis/Results_L1/2020_11_19_Neutrino_Files/";

    TString file_dir_nu = "/afs/cern.ch/work/r/rmccarth/private/DarkPhotonT10/";

    
/*
    TChain *ch1 = new TChain("L1TrackNtuple/eventTree");
    ch1->Add(file_dir + "Dark_Photon_cT0/events_Dark_Photon_cT0.root");
    Analyzer_DisplacedMuon t1(ch1);
    t1.Loop("events_Dark_Photon_cT0",type_dir+"cT0/",10.0,0.0,0);

    TChain *ch2 = new TChain("L1TrackNtuple/eventTree");
    ch2->Add(file_dir + "Dark_Photon_cT10/events_Dark_Photon_cT10.root");
    Analyzer_DisplacedMuon t2(ch2);
    t2.Loop("events_Dark_Photon_cT10", type_dir + "cT10/",10.0,0.0,0);
/*>
    TChain *ch3 = new TChain("L1TrackNtuple/eventTree");
    ch3->Add(file_dir + "Dark_Photon_cT100/events_Dark_Photon_cT100.root");
    Analyzer_DisplacedMuon t3(ch3);
    t3.Loop("events_Dark_Photon_cT100", type_dir + "cT100/",10.0,0.0,0);

    TChain *ch7 = new TChain("L1TrackNtuple/eventTree");
    ch7->Add(file_dir + "Dark_Photon_cT10_PU200/events_Dark_Photon_cT10_PU200.root");
    Analyzer_DisplacedMuon t7(ch7);
    t7.Loop("events_Dark_Photon_cT10_PU200", type_dir + "cT10_PU200/",10.0,0.0,0);
/*
    TChain *ch4 = new TChain("L1TrackNtuple/eventTree");
    ch4->Add(file_dir + "Dark_Photon_cT5000/events_Dark_Photon_cT5000.root");
    Analyzer_DisplacedMuon t4(ch4);
    t4.Loop("events_Dark_Photon_cT5000", type_dir + "cT5000/",10.0,0.0,0);

    TChain *ch5 = new TChain("L1TrackNtuple/eventTree");
    ch5->Add(file_dir + "Dark_Photon_cT10000/events_Dark_Photon_cT10000.root");
    Analyzer_DisplacedMuon t5(ch5);
    t5.Loop("events_Dark_Photon_cT10000", type_dir + "cT10000/",10.0,0.0,0);
    
    
    TChain *ch8 = new TChain("L1TrackNtuple/eventTree");
    ch8->Add(file_dir + "DisplacedMu_PU0/events_DisplacedMu_PU0.root");
    Analyzer_DisplacedMuon t8(ch8);
    t8.Loop("events_DisplacedMu_PU0", type_dir + "DisplacedMu_PU0/",30.0,0.0,0);
    
    TChain *ch9 = new TChain("L1TrackNtuple/eventTree");
    ch9->Add(file_dir + "DisplacedMu_PU200/events_DisplacedMu_PU200.root");
    Analyzer_DisplacedMuon t9(ch9);
    t9.Loop("events_DisplacedMu_PU200", type_dir + "DisplacedMu_PU200/",30.0,0.0,0);
*/
    TChain *ch6 = new TChain("L1TrackNtuple/eventTree");
    ch6->Add(file_dir_nu + "events_Dark_Photon_cT10.root");
    Analyzer_DisplacedMuon t6(ch6);
    t6.Loop("events_DarkPhoton_cT10", type_dir + "DarkPhoton_cT10_CMSPresentation/", 10.0, 0.0004196, 0);
    return 0;
}
