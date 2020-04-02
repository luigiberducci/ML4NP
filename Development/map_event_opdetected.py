import ROOT
import os
import pandas as pd
import numpy as np

# IO params
map_dir = os.path.join("..", "Data", "root")
map_filebase = "OpticalMapL200XeD.14String.5mm.root"
sim_dir = os.path.join("..", "Data", "cutROI")
sim_filebase= "output2eROI_all.csv"

# Get OP map
map_file = ROOT.TFile.Open(os.path.join(map_dir, map_filebase))
h_map = ROOT.TH3D()
map_file.GetObject("ProbMapInterior", h_map)

# Get event
df = pd.read_csv(os.path.join(sim_dir, sim_filebase), index_col=False)
df = df[df.energydeposition>0]
eventnumbers = [555, 19343]
for eventnumber in eventnumbers:
    print("[Info] Processing event {}".format(eventnumber))
    event = df[(df.eventnumber==eventnumber)].sort_values("time")
    OPs = []
    for entry in event.iterrows():
        bin = h_map.FindBin(entry.x, entry.y, entry.z)
        ops = h_map.GetBinContent(bin)
        OPs.append(ops)
    event["OPs"] = OPs
    event.to_csv("event{}_OPs.csv".format(eventnumber))
    print("[Info] Saved file: event{}_OPs.csv".format(eventnumber))
print("[Info] Done.")
