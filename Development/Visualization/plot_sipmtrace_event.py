import ROOT
import os

filepath = os.path.join("..", "tmpsipm72SiPMs_40Yield_part1.root")
f = ROOT.TFile.Open(filepath, "READ")
fTree = f.Get("fTree")
nSiPM = f.Get("NSiPM").GetVal()

theReader = ROOT.TTreeReader("fTree", f)
# Branches
sipm_branches = []
for i in range(nSiPM):
    branch = ROOT.TTreeReaderValue(int)(theReader, "SiPM{}".format(i));
    sipm_branches.append(branch)
