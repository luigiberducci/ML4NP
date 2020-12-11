import numpy as np

def pe_detected(row):
    return np.sum(row)

def nr_active_slices(row):
    return np.nonzero(row)[0].shape[0]

def mean_npe(row):
    return np.mean(row)

def mean_npe_active(row):
    rownonzero = np.nonzero(row)[0]
    return np.mean(row[rownonzero]) if rownonzero.shape[0]>0 else -1

def std_npe(row):
    return np.std(row)

def std_npe_active(row):
    rownonzero = np.nonzero(row)[0]
    return np.std(row[rownonzero]) if rownonzero.shape[0]>0 else -1

def range_detections(row):
    rownonzero = np.nonzero(row)[0]
    return rownonzero[-1] - rownonzero[0] + 1 if rownonzero.shape[0]>0 else -1

def spatial_var(row):
    ids = np.repeat(np.argwhere(row>0), row[row>0])
    return np.var(ids) if ids.shape[0]>0 else -1

def spatial_std(row):
    ids = np.repeat(np.argwhere(row>0), row[row>0])
    return np.std(ids) if ids.shape[0]>0 else -1