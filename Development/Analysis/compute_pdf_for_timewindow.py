import numpy as np
import argparse as parse

parser = parse.ArgumentParser()
parser.add_argument("rate", help="event rate in 1 second", type=int)
parser.add_argument("time", help="time window in nanoseconds (ns)", type=float)
args = parser.parse_args()

print("[Info] Computing the event's PDF over time interval {} ns, assuming rate {}/sec...".format(args.time, args.rate))

rate_t = args.rate / 10**9 * args.time
print("[Debug] Rate in T: {}".format(rate_t))

def fact(k):
    if k == 1:
        return 1
    return k * fact(k-1)

def pois(l, k):
    return l**k * np.exp(-l) / fact(k)

for k in range(10):
    print(pois(rate_t, k))
