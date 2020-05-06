# Script:   map_event_opdetected
# Status:   deprecated
# Info:     Given an input simulation file as csv export,
#           this script adds a field `OPs` (detection efficiency from OP map).
#           This feature is now integrated in `Preprocessing.C`
import ROOT
import os, time
import pandas as pd
import numpy as np

# IO params
map_dir = os.path.join("..", "Data", "root")
map_filebase = "OpticalMapL200XeD.14String.5mm.root"
sim_dir = os.path.join("..", "Data", "cutROI")
sim_filebase= "output2eROI_all"

# Get OP map
print("[Info] Loading optical map: {}".format(map_filebase))
map_file = ROOT.TFile.Open(os.path.join(map_dir, map_filebase))
h_map = ROOT.TH3D()
map_file.GetObject("ProbMapInterior", h_map)

# Get event
print("[Info] Loading input file: {}.csv".format(sim_filebase))
df = pd.read_csv(os.path.join(sim_dir, "{}.csv".format(sim_filebase)), index_col=False)
ops = np.zeros((df.shape[0], 1))
time0 = time.time()
for i, entry in df.iterrows():
    if i % 1000000 == 0:
        print("[Info] Time: {}, entry: {}".format(time.time()-time0, i+1))
    if entry.energydeposition <= 0:
        continue
    bin = h_map.FindBin(entry.x, entry.y, entry.z)
    ops[i] = h_map.GetBinContent(bin)
df["OPs"] = ops
import ipdb
ipdb.set_trace()
df.to_csv("{}_wt_ops.csv".format(sim_filebase))
print("[Info] Done.")
