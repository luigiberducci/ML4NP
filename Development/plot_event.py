import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import argparse
import os
import warnings

particles = {'e-': 11, 'e+': -11, 'm': 13, 'n': 2112, 'g': 22}  # dict particle names -> ids
labels = {'e-': 'electron', 'e+': 'positron', 'm': 'muon', 'n': 'neutron', 'g': 'photon', 'ot': 'others'}  # dict particle id->name
colors = {'e-': 'r', 'e+': 'b', 'm': 'orange', 'n': 'y', 'g': 'g', 'ot': 'black'}   # dict particle colors
edep_scale = 0.1    # scaling factor to avoid "huge" markers

default_infile = os.path.join("..", "Data", "output_Muon_10000.csv")
full_df = pd.DataFrame()

def plot_event(file_path, event):
    df = full_df[full_df['eventnumber'] == event]
    available_pid = sorted(list(df.PID.unique()))
    selected_pid = sorted(set(particles.values()).intersection(available_pid))

    print("[Info]\tEvent {} loaded.\tNum entries: {}".format(event, len(df)))
    print("[Info]\tFound particles: {}\n\tPlot particles: {}".format(available_pid, selected_pid))

    threedee = plt.figure().gca(projection='3d')
    threedee.set_title('Event {} from {}'.format(event, file_path.name))
    threedee.set_xlabel('X')
    threedee.set_ylabel('Y')
    threedee.set_zlabel('Z')
    handles = []
    for id in labels.keys():
        patch = mpatches.Patch(color=colors[id], label=labels[id])
        handles.append(patch)
    threedee.legend(handles=handles)

    df_gb = df.sort_values("time").groupby("PID")
    muons = df_gb.get_group(particles["m"]) if particles['m'] in selected_pid else pd.DataFrame()
    posit = df_gb.get_group(particles["e+"]) if particles['e+'] in selected_pid else pd.DataFrame()
    elect = df_gb.get_group(particles["e-"]) if particles['e-'] in selected_pid else pd.DataFrame()
    neutr = df_gb.get_group(particles["n"]) if particles['n'] in selected_pid else pd.DataFrame()
    gamma = df_gb.get_group(particles["g"]) if particles['g'] in selected_pid else pd.DataFrame()

    # plot all
    others = df[~(df.PID.isin(particles.values()))]
    threedee.scatter(others.x, others.y, others.z, c=colors['ot'], marker='*', s=edep_scale*others.energydeposition)
    # plot muon trajectory (main trajectory)
    if not muons.empty:
        threedee.plot(muons.x, muons.y, muons.z, c=colors['m'])
    # the other particles, plot each single track
    if not elect.empty:
        elect_gb = elect.groupby("tracknumber")
        elect_gb.apply(lambda track: threedee.plot(track.x, track.y, track.z, c=colors['e-']))
        elect_gb.apply(lambda track: threedee.scatter(track.x, track.y, track.z, c=colors['e-'],
                                                      marker='*', s=edep_scale*track.energydeposition))
    if not posit.empty:
        posit_gb = posit.groupby("tracknumber")
        posit_gb.apply(lambda track: threedee.plot(track.x, track.y, track.z, c=colors['e+']))
        posit_gb.apply(lambda track: threedee.scatter(track.x, track.y, track.z, c=colors['e+'],
                                                      marker='*', s=edep_scale*track.energydeposition))
    if not neutr.empty:
        neutr_gb = neutr.groupby("tracknumber")
        neutr_gb.apply(lambda track: threedee.plot(track.x, track.y, track.z, c=colors['n']))
        neutr_gb.apply(lambda track: threedee.scatter(track.x, track.y, track.z, c=colors['n'],
                                                      marker='*', s=edep_scale*track.energydeposition))
    if not gamma.empty:
        gamma_gb = gamma.groupby("tracknumber")
        gamma_gb.apply(lambda track: threedee.plot(track.x, track.y, track.z, c=colors['g']))
        gamma_gb.apply(lambda track: threedee.scatter(track.x, track.y, track.z, c=colors['g'],
                                                      marker='*', s=edep_scale*track.energydeposition))
    threedee.set_xlim(-1950, +1950)
    threedee.set_ylim(-1950, +1950)
    threedee.set_zlim(-1950, +1950)
    plt.show()

def main():
    global full_df
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=argparse.FileType('r'), default=default_infile, help="Input filepath")
    parser.add_argument("eventnumber", type=int, default=None, nargs='?', help="Event number (id)")
    args = parser.parse_args()

    file_path = args.infile
    event = args.eventnumber
    print("[Info]\tLoad file {}".format(file_path.name))
    full_df = pd.read_csv(file_path, index_col=False)

    if event is not None:
        plot_event(file_path, event)
    else:
        eventnumbers = full_df.eventnumber.unique()
        print("[Info]\tLoaded. Num events: {}\nAvailable event numbers:\n{}".format(len(eventnumbers), eventnumbers))

        while True:
            try:
                event = int(input("[Input]\tInsert event number or 'q'..."))
                if event not in eventnumbers:
                    warnings.warn("no event {} in file {}".format(event, file_path.name))
            except ValueError:
                break
            plot_event(file_path, event)
    print("[Info] End.")

if __name__=="__main__":
    main()
