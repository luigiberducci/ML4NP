import ROOT
import os
import numpy as np
import pandas as pd
import matplotlib as matlib
import matplotlib.pyplot as plt

# Input params
filepath = os.path.join("..", "tmpsipm72SiPMs_40Yield_part1.root")
eventnumber = 100
delta_interval = 100000  # ns
discrete_interval = '1us'

# Open tree
f = ROOT.TFile.Open(filepath, "READ")
fTree = f.Get("fTree")
nSiPM = f.Get("NSiPM").GetVal()
sipm_cols = ["SiPM{}".format(i) for i in range(nSiPM)]  # Assume: SiPM_, where _ is int
# Select event
eventTree = fTree.CopyTree("(eventnumber=={})".format(eventnumber))

# Use Dataframe
data, columns = eventTree.AsMatrix(return_labels=True)
df = pd.DataFrame(data=data, columns=columns)
print(columns)

min_time = df.time.min() - delta_interval
max_time = min_time + 2 * delta_interval

delimiter_rows = np.zeros((2, len(columns)))
delimiter_rows[:, 0] = [eventnumber] * 2
delimiter_rows[:, 1] = [min_time, max_time]
delimiter_df = pd.DataFrame(delimiter_rows, columns=columns)

df = df[df.time <= max_time]
df = df.append(delimiter_df)
df["datetime"] = pd.to_timedelta(df.time - min_time, unit='ns')
df = df.set_index("datetime")

print(df)
df = df.resample(discrete_interval).sum()
print(df)

df = df.reset_index()
time = df.datetime
for sipm_col in sipm_cols:
    if df[sipm_col].max() > 0:
        plt.step(time, df[sipm_col], label=sipm_col)
    else:
        plt.step(time, df[sipm_col])
plt.legend()
plt.show()
