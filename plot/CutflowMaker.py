import ROOT
from ROOT import THStack
from ROOT import TGaxis
from ROOT import TColor
from array import array
import math
import gc
import sys
import numpy as np
import copy
import os
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1;")
ROOT.TH1.AddDirectory(ROOT.kFALSE)
ROOT.gStyle.SetOptStat(0)
import argparse
TGaxis.SetMaxDigits(1)

def cutFlowTable(hists, samples, regions,totals):
#    table = '\\begin{sidewaystable*}' + "\n"
    table = '\\begin{table}' + "\n"
    table += '\\centering' + "\n"
    table += '\\caption{' +"CutFlow"+ "}\n"
    table += '\\resizebox{\\textwidth}{!}{ \n'
    table += '\\begin{tabular}{|l|l|' + ''.join([('' if False else 'l|') for c in regions]).strip() + '}' + "\n"
    table += '\\hline' + "\n"
    table += 'Samples ' + ''.join([('' if False else ' & '+c) for c in regions]).strip() + '\\\\' + "\n"
    table += '\\hline' + "\n"

    for ids, s in enumerate(samples):
        table += s
        table += (' & '+ str(totals[ids]))
        for idr, r in enumerate(regions):
            table += (' & ' + str(round(hists[ids][idr][0].GetEntries(), 0)))
        table += '\\\\' + "\n"
    table += '\\hline' + "\n"
    table += '\\end{tabular}}' + "\n"
    table += '\\end{table}' + "\n"
#    table += '\\end{sidewaystable*}' + "\n"
    return table

regions = ["tpgeq2", "tpgeq2osMu"]
variables = ["ntp", "nMu",
             "tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_d0", "tp_charge", "tp_pdgid",\
             "tp1pt", "tp1eta", "tp1phi", "tp1z0", "tp1d0", "tp1charge", "tp1pdgid",\
             "tp2pt", "tp2eta", "tp2phi", "tp2z0", "tp2d0", "tp2charge", "tp2pdgid",\
             "tp_tp_eta", "tp_tp_phi", "tp_tp_z0", "tp_tp_d0", "tp_tp_dR",\
             "tpz1pt", "tpz1eta", "tpz1phi", "tpz1z0", "tpz1d0", "tpz1charge", "tpz1pdgid",\
             "tpz2pt", "tpz2eta", "tpz2phi", "tpz2z0", "tpz2d0", "tpz2charge", "tpz2pdgid",\
             "tpz_tpz_eta", "tpz_tpz_phi", "tpz_tpz_z0", "tpz_tpz_d0", "tpz_tpz_dR",\
             "tpd1pt", "tpd1eta", "tpd1phi", "tpd1z0", "tpd1d0", "tpd1charge", "tpd1pdgid",\
             "tpd2pt", "tpd2eta", "tpd2phi", "tpd2z0", "tpd2d0", "tpd2charge", "tpd2pdgid",\
             "tpd_tpd_eta", "tpd_tpd_phi", "tpd_tpd_z0", "tpd_tpd_d0", "tpd_tpd_dR"]
# variablesName = []

variablesName = variables

# set up an argument parser
parser = argparse.ArgumentParser()

parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--l', dest='LOCATION',default='/afs/cern.ch/user/b/bharikri/Projects/Thesis/L1Tracking/DisplacedMuon/plot/')

ARGS = parser.parse_args()

verbose = ARGS.VERBOSE
HistAddress = ARGS.LOCATION


Samples = ['events_Dark_Photon_cT0.root', 'events_Dark_Photon_cT10.root','events_Dark_Photon_cT100.root','events_Dark_Photon_cT5000.root','events_Dark_Photon_cT10000.root']#,'events_NeutrinoGun_PU200.root']
SamplesName = ['cT0','cT10','cT100','cT5000','cT10000','NeutrinoGun_PU200']
SamplesNameLatex = ['cT0', 'cT10', 'cT100', 'cT5000', 'cT10000','NeutrinoGun\_PU200']
Totals = [99998, 99998, 99998, 77883, 77883, 872538]

Hists = []
Files = []
for f in range(len(Samples)):
    l0=[]
    Files.append(ROOT.TFile.Open(HistAddress +SamplesName[f]+"/output_"+Samples[f]))
    print(HistAddress + SamplesName[f]+"/output_"+Samples[f])
    for numreg, namereg in enumerate(regions):
        l1=[]
        for numvar, namevar in enumerate(variables):
            h= Files[f].Get(namereg + '_' + namevar)
            l1.append(h)
        l0.append(l1)
    Hists.append(l0)       

le = '\\documentclass{article}' + "\n"
le += '\\usepackage{rotating}' + "\n"
le += '\\begin{document}' + "\n"

le += cutFlowTable(Hists, SamplesNameLatex, regions, Totals)

le += '\\end{document}' + "\n"

print le

with open('Cutflow.tex', 'w') as fout:
    for i in range(len(le)):
        fout.write(le[i])

