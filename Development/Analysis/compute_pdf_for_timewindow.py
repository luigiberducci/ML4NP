import numpy as np
import argparse as parse

def fact(k):
    if k <= 1:
        return 1
    return k * fact(k-1)

def pois(l, k):
    return l**k * np.exp(-l) / fact(k)

def convert_ns_to_years(ns):
    yy = ns / 1000000000 / 60 / 60 / 24 / 365
    return yy

parser = parse.ArgumentParser()
parser.add_argument("rate", help="event rate in 1 second", type=float)
parser.add_argument("time", help="time window in nanoseconds (ns)", type=float)
parser.add_argument("--maxk", help="max number of events to consider", type=int, default=10)
parser.add_argument("--multxp", help="probability factor of coincident (independent) event", type=float, default=1)
args = parser.parse_args()

print("[Info] Computing the event's PDF over time interval {} ns, assuming rate {}/sec...".format(args.time, args.rate))

rate_t = args.rate / 10**9 * args.time
timewindow = args.time
multxp = args.multxp    # used to see prob of events in the same time (muon || ar39)
print("[Debug] Rate in T: {}".format(rate_t))

print("K\tp~Pois(k,rate)\tE[X=1|X~Gm(p)]\tExp.Time (Years)\tExp.Time (Days)")
for k in range(args.maxk + 1):
    prob_kevents_in_window = multxp * pois(rate_t, k)
    exp_snapshot_to_first = np.ceil(1 / prob_kevents_in_window)
    exp_time_to_first = convert_ns_to_years(exp_snapshot_to_first * timewindow)
    exp_time_to_first_days = exp_time_to_first * 365
    print("{}\t{:.8e}\t{:e}\t{}\t{}".format(k, prob_kevents_in_window, exp_snapshot_to_first, exp_time_to_first, exp_time_to_first_days))
