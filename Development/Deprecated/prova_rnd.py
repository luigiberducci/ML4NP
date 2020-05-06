import numpy as np
import time

init = time.time()
N = 10000000
A = np.random.normal(0.5, 0.5, N)
print("Time: {}".format(time.time() - init))
