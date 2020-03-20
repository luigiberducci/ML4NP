import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import math
import os

particles = {'m': 13, 'n':2112, 'g':22}  # dict particle names -> ids

parser = argparse.ArgumentParser()
parser.add_argument("eventnumber", type=int, default=0, nargs='?', help="Event number (id)")
args = parser.parse_args()

event=args.eventnumber
plot_interval = 10000
file_path = os.path.join("..", "Data", "output123456789.csv")
file_path = os.path.join("..", "Data", "output10000.csv")
df = pd.read_csv(file_path, index_col=False)
df = df[df['eventnumber'] == event]
print("[Info] Event {} loaded. All entries: {}\nPlot only muons, electrons, neutrons, photons".format(event, len(df)))

threedee = plt.figure().gca(projection='3d')
threedee.set_xlabel('X')
threedee.set_ylabel('Y')
threedee.set_zlabel('Z')

# define sipm topology
n_sipm = 30     # num sipm for top/bottom
alphas = np.linspace(0, 2*math.pi, n_sipm)
radius = 1500
x = 0 + radius * np.cos(alphas)  # x-coord
z = 0 + radius * np.sin(alphas)  # z-coord, y-coord is fixed for top and bottom
y_top = + 1500
y_bottom = - 1500
# threedee.scatter(x, y_top, z, c="orange")
# threedee.scatter(x, y_bottom, z, c="orange")

for i in range(len(df)):
    x = df.iloc[i,:]['x']
    y = df.iloc[i,:]['y']
    z = df.iloc[i,:]['z']
    pid = df.iloc[i,:]['PID']
   # parent = df.iloc[i,:]['PID_parent']
    color = "black"    # all particles
    if pid == 11:
        color = "red"  # electrons
    elif pid == 13:
        color = "blue"  # muon
    elif pid == 22:
        color = "green"   # gamma
    elif pid == 2112:
        color = "black" # neutron
    # if pd.isna(parent) and pid!=13: #not muon
    #     color = "purple"
    if color=="black":
        continue
    threedee.scatter(x, y, z, c=color)
    if i>0 and i % plot_interval == 0:
        plt.pause(0.001)
threedee.set_xlim(-1500,1500)
threedee.set_ylim(-1500,1500)
threedee.set_zlim(-1500,1500)
plt.show()
