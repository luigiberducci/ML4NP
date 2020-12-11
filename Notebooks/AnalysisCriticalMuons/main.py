import os
import random
import time

import ROOT
import argparse
from joblib import load
import sklearn
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

in_file = None
tmp_file = None
fTree = None
eventTree = None


def get_config():
  config, types = {}, {}
  config['input_dir'] = os.path.join('/', 'home', 'luigi', 'Development', 'ML4NP', 'Data', 'OutputProcessing',
                                     'UnseenTestData_11_10_2020', 'Muons', 'Muons_Sliced/')
  config['sim_file'] = ''
  config['model'] = ''
  config['n_pileup'] = 1
  for key, value in config.items():
    types[key] = type(value)
  return config, types


def load_simulation_file(sim_file):
  global in_file, fTree
  try:
    in_file = ROOT.TFile.Open(sim_file, "READ")
    fTree = in_file.Get("fTree")
  except Exception:
    print(f'[Error] Fail to read simulation file: {sim_file}')
    exit(-1)


def get_event_set(keep_zero_detection_events=True):
  events = set()
  for i in range(0, fTree.GetEntries()):
    fTree.GetEntry(i)
    event = getattr(fTree, "eventnumber")
    npe = getattr(fTree, "pedetected")
    if keep_zero_detection_events:
      events.add(event)
    elif npe > 0:
      events.add(event)
  return events


def print_summary_events(max_numbers_per_line=20):
  events = get_event_set()
  print('[Info] List of event numbers:')
  for i, event in enumerate(events):
    end = ',\n' if i > 0 and i % max_numbers_per_line == 0 else ', '
    print(event, end=end)
  print("\n")


def load_prediction_model(model_file):
  try:
    model = load(model_file)
  except Exception:
    print(f'[Error] Fail to read saved model: {model_file}')
    exit(-1)
  return model


def get_input_event(events):
  event_nr = None
  try:
    inp = input("[Input] Insert an event number (-1 exit, 'events' for summary, 'r' for random):")
    if inp == "events":
      print_summary_events()
    elif inp == "r" or inp == "":
      event_nr = random.sample(events, 1)[0]
    else:
      event_nr = int(inp)
    if event_nr < 0:
        raise Exception
  except Exception:
    exit(0)
  return event_nr


def filter_event_entries(eventnumber):
  global eventTree, tmp_file
  tmp_file = ROOT.TFile.Open("tmp", "RECREATE")
  eventTree = fTree.CopyTree(f'eventnumber=={eventnumber}')


def plot_roi_deposition(coords, energies, npes, title="", energy_scale=500):
  f = plt.figure(figsize=(7, 22))
  gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
  ax1 = f.add_subplot(gs[0])
  ax1.set_title(title)
  ax1.scatter(coords[:, 0], coords[:, 2], s=energies / energy_scale, c=np.where(npes>0, 'g', 'r'))
  ax1.set_xlim(-710, +710)
  ax1.set_ylim(-850, +850)
  ax1.set_xlabel('X')
  ax1.set_ylabel('Z')

  ax2 = f.add_subplot(gs[1])
  ax2.scatter(coords[:, 0], coords[:, 1], s=energies / energy_scale, c=np.where(npes>0, 'g', 'r'))
  ax2.set_xlim(-710, +710)
  ax2.set_ylim(-710, +710)
  ax2.set_xlabel('X')
  ax2.set_ylabel('Y')
  plt.show()

def get_detections_on_time_window(eventTree, t0, time_window, n_inner_slices=12, n_outer_slices=20):
  inner_detections, outer_detections = np.zeros(n_inner_slices, dtype=int), np.zeros(n_outer_slices, dtype=int)
  for i in range(eventTree.GetEntries()):
    eventTree.GetEntry(i)
    t = getattr(eventTree, 'time')
    if t0 <= t <= t0 + time_window:
      for j in range(12):
        inner_detections[j] += getattr(eventTree, f'InnerSlice{j}')
      for j in range(20):
        outer_detections[j] += getattr(eventTree, f'OuterSlice{j}')
  return inner_detections, outer_detections

def plot_detections(inner_detections, outer_detections, title="", inner_r=175, outer_r=295):
  n_inner_slices, n_outer_slices = inner_detections.shape[0], outer_detections.shape[0]
  f = plt.figure(figsize=(10, 10))
  ax = f.gca()
  ax.set_title(title)
  inner_theta, outer_theta = 2 * np.pi / n_inner_slices, 2 * np.pi / n_outer_slices
  for n_slices, theta, radius, detections in zip([n_inner_slices, n_outer_slices], [inner_theta, outer_theta],
                                                 [inner_r, outer_r], [inner_detections, outer_detections]):
    shroud = plt.Circle((0, 0), radius, color='k', fill=False)
    ax.add_artist(shroud)
    for i in range(n_slices):
      ang = i * theta
      x, y = np.cos(ang) * radius, np.sin(ang) * radius
      ax.scatter(x, y, c='k')
      ang = i * theta + .5 * theta
      x, y = np.cos(ang) * radius, np.sin(ang) * radius
      ax.text(x, y, f'{int(detections[i])}', fontsize=25)
  ax.set_xlim(-310, +310)
  ax.set_ylim(-310, +310)
  plt.show()

def compute_rf_features(inner_detections, outer_detections):
  import feature_extraction as feat
  shroud_features = ["PEDetected", "NActiveSlices", "MeanNPEActive", "SpatialVar", "SpatialRange"]
  global_features = ["PEDetected", "NActiveSlices", "MeanNPEActive"]
  shroud_lambdas = [feat.pe_detected, feat.nr_active_slices, feat.mean_npe_active, feat.spatial_var, feat.range_detections]
  global_lambdas = [feat.pe_detected, feat.nr_active_slices, feat.mean_npe_active, feat.spatial_var]
  data = {}
  for shroud_name, detections in zip(["outer", "inner"], [outer_detections, inner_detections]):
    for feature_name, feature_fun in zip(shroud_features, shroud_lambdas):
      data[f'{feature_name}_{shroud_name}'] = feature_fun(detections)
  for feature_name, feature_fun in zip(global_features, global_lambdas):
    data[f'{feature_name}_tot'] = feature_fun(np.concatenate([inner_detections, outer_detections]))
  return data

def plot_model_classification(model, inner_detections, outer_detections, title="", flag_plot_decisions=False):
  feature_data = compute_rf_features(inner_detections, outer_detections)
  muon_proba = np.zeros(2+model.n_estimators) + 0.1
  colors = ['grey' for i in range(2+model.n_estimators)]
  ticks = ["<4NAS"] + [f"DT{i+1}" for i in range(model.n_estimators)] + [">60PE"]
  if feature_data["NActiveSlices_outer"] < 4:
    pred = np.array([[1, 0]])
    muon_proba[0] += 0
    colors[0] = 'r'
  elif feature_data["PEDetected_tot"] > 60:
    pred = np.array([[0, 1]])
    muon_proba[-1] += 1
    colors[-1] = 'g'
  else:
    X = np.array(list(feature_data.values())).reshape(1, -1)
    pred = model.predict_proba(X)
    estimators_proba = [model.estimators_[i].predict_proba(X)[0][1] for i in range(model.n_estimators)]
    colors[1: -1] = ['g' if estimators_proba[i]>.5 else 'r' for i in range(model.n_estimators)]
    muon_proba[1:-1] += estimators_proba
  if flag_plot_decisions:
    plt.title(title + ", (P[mu]={:.3f})".format(pred[0][1]))
    plt.bar(range(2+model.n_estimators), muon_proba, color=colors)
    plt.xticks(range(2+model.n_estimators), ticks)
    plt.yticks([0.1, .6, 1.1], ["Ar39", "boundary", "Muon"])
    plt.hlines(0.6, 0, 2+model.n_estimators, colors='r', linestyles='dashed')
    plt.ylim(0, 1.1)
    plt.show()
  return pred


def get_first_detection_time(eventTree):
  first_detection_t, time_window = np.Inf, 10000
  for i in range(0, eventTree.GetEntries()):
    eventTree.GetEntry(i)
    npe = getattr(eventTree, 'pedetected')
    t = getattr(eventTree, 'time')
    if npe > 0 and t < first_detection_t:
      first_detection_t = t
  return first_detection_t

def get_energy_depositions(eventTree):
  npes, energies, coords, times = [], [], [], []
  for i in range(0, eventTree.GetEntries()):
    eventTree.GetEntry(i)
    edep = getattr(eventTree, 'energydeposition')
    npe = getattr(eventTree, 'pedetected')
    x = getattr(eventTree, 'x')
    y = getattr(eventTree, 'y')
    z = getattr(eventTree, 'z')
    t = getattr(eventTree, 'time')
    npes.append(npe)
    energies.append(edep)
    coords.append([x, y, z])
    times.append(t)
  return np.array(npes), np.array(energies), np.array(coords), np.array(times)

def classify_single_file(input_filepath, model, n_pileup=1, n_inner_slices=12, n_outer_slices=20,
                         interactive=False, flag_plot_edep=False, flag_plot_detections=False, flag_plot_decisions=False):
  classes = ["Ar39", "Muon"]
  load_simulation_file(input_filepath)
  if interactive:
    print_summary_events(20)
  events = get_event_set()
  preds = []
  while True:
    events_to_show = []
    while len(events_to_show) < n_pileup:
      if interactive:
        event_nr = get_input_event(events)
        if event_nr is None:
          continue
        if event_nr not in events:
          print("[Info] Event not found.")
          continue
      elif len(events)>0:
        event_nr = events.pop()
      else:
        break
      events_to_show.append(event_nr)
    all_inner_detections, all_outer_detections = np.zeros(n_inner_slices, dtype=int), np.zeros(n_outer_slices, dtype=int)
    all_npes, all_energies, all_coords, all_times = [], [], np.array([]).reshape((0,3)), []
    for event_nr in events_to_show:
      filter_event_entries(event_nr)
      t0 = get_first_detection_time(eventTree)
      inner_detections, outer_detections = get_detections_on_time_window(eventTree, t0, 10000)
      all_inner_detections += inner_detections
      all_outer_detections += outer_detections
      if flag_plot_edep:
        npes, energies, coords, time = get_energy_depositions(eventTree)
        all_npes = np.concatenate([all_npes, npes])
        all_energies = np.concatenate([all_energies, energies])
        all_coords = np.concatenate([all_coords, coords])
        all_times = np.concatenate([all_times, time])
    if flag_plot_edep:
      plot_roi_deposition(all_coords, all_energies, all_npes, f'Events: {events_to_show}, Tot Energy Dep: {np.sum(all_energies)} KeV')
    if flag_plot_detections:
      plot_detections(all_inner_detections, all_outer_detections, f'Events: {events_to_show}, SiPM Detections ({np.sum(all_npes)} PE)')
    pred = plot_model_classification(model, all_inner_detections, all_outer_detections,
                                     f'Events {events_to_show}', flag_plot_decisions)
    preds.append(np.argmax(pred))
    if interactive:
      print('[Info] Events {}: ({:d} nzero entries): {:s}\n'.format(events_to_show, len(all_npes), classes[np.argmax(pred)]))
  return preds

def main(args):
  model = load_prediction_model(args.model)
  import glob
  all_preds = []
  if args.sim_file:
    preds = classify_single_file(args.sim_file, model, args.n_pileup, interactive=True,
                                 flag_plot_edep=True, flag_plot_detections=True, flag_plot_decisions=True)
    all_preds += preds
  elif args.input_dir:
    for i, file in enumerate(glob.glob(f'{args.input_dir}/*root')):
      print(f"[Info] File {i+1}: {file}")
      preds = classify_single_file(file, model, args.n_pileup)
      all_preds += preds
  else:
    raise NotImplemented("either input_dir or sim_file must be valued")
    exit(-1)
  print()
  print(f"[Result] Accuracy: {sum(all_preds)/len(all_preds)}")
  print(f"[Result] Tot Time: {sum(times)} sec, Avg Time File: {sum(times)/len(times)} sec")

if __name__ == '__main__':
  config, types = get_config()
  parser = argparse.ArgumentParser()
  for key, value in config.items():
    parser.add_argument(f'--{key}', type=types[key], default=value, nargs='?')
  args = parser.parse_args()
  main(args)
