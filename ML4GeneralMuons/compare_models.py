import numpy as np, pandas as pd, os
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score, f1_score
from sklearn.metrics import confusion_matrix as cv
import datetime
import tensorflow as tf
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False    # no deprecation warning print
import logging
logger = tf.get_logger()
logger.setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
plt.rcParams['figure.figsize'] = [25, 15]
plt.rcParams.update({'font.size': 18})

def load_dataset():
    ar39filename = "Ar39_{}Pileup_cut265PE.csv"
    muonfilename = "MarginalMuons_wt_{}ar39_cut265PE.csv"
    # load ar39, mark various groups
    dfar39 = pd.read_csv(os.path.join(DATADIR, "Ar39", ar39filename.format(1)), index_col=False)
    dfar39['group'] = -1
    for i in range(2, 7+1):
        df = pd.read_csv(os.path.join(DATADIR, "Ar39", ar39filename.format(i)), index_col=False)
        df['group'] = -(i)
        dfar39 = pd.concat([dfar39, df])
    # load muons, mark various groups
    dfmuons = pd.read_csv(os.path.join(DATADIR, "Muons", muonfilename.format(0)), index_col=False)
    dfmuons['group'] = 0
    for i in range(1, 3+1):
        df = pd.read_csv(os.path.join(DATADIR, "Muons", muonfilename.format(i)), index_col=False)
        df['group'] = i
        dfmuons = pd.concat([dfmuons, df])
    print("[Info] Loaded {} Ar39 events, {} Muon events.".format(len(dfar39), len(dfmuons)))
    # set boolean label
    dfar39['y'] = 0
    dfmuons['y'] = 1
    dfall = pd.concat([dfmuons, dfar39]).to_numpy()
    return dfall.reshape((dfall.shape[0], dfall.shape[1], 1))

def get_twoconv_model_sigmoid01(n_delta_t=1, n_slices=72):
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Conv1D(filters=32, kernel_size=6, activation='relu', input_shape=(n_slices, n_delta_t)))
    model.add(tf.keras.layers.Conv1D(filters=32, kernel_size=3, activation='relu'))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=3))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(16, activation='relu'))
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
    return model

def get_twoconv_model_tanh_11(n_delta_t=1, n_slices=72):
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Conv1D(filters=32, kernel_size=6, activation='relu', input_shape=(n_slices, n_delta_t)))
    model.add(tf.keras.layers.Conv1D(filters=32, kernel_size=3, activation='relu'))
    model.add(tf.keras.layers.MaxPooling1D(pool_size=3))
    model.add(tf.keras.layers.Flatten())
    model.add(tf.keras.layers.Dense(16, activation='relu'))
    model.add(tf.keras.layers.Dense(1, activation='tanh'))
    return model

def get_my_confusion_matrix_as_string(g_test, y_pred):
    # Note: group <0 are ar39 (-1, -2, ..., -7), others are mu (0, 1, 2, 3)
    cv_string = "Test Set Size: {} ({} Ar39, {} Mu)\n".format(g_test.shape[0], g_test[g_test<0].shape[0], g_test[g_test>=0].shape[0])
    cv_string += "\t\tPred\tPred\t| Nr\n"
    cv_string += "\t\tAr39\tMu\t| samples\n"
    cv_string += "--------------------------------------------\n"

    ar39_pred_ids = np.nonzero(y_pred.round()==0)[0]
    for i_pileup in range(1, 7 + 1):
        ar39_pileup_ids = np.nonzero(g_test==-i_pileup)[0]
        tp = len(np.intersect1d(ar39_pred_ids, ar39_pileup_ids))
        cv_string += "{} Ar39\t\t{}\t{}\t| {}\n".format(i_pileup, tp, len(ar39_pileup_ids)-tp, len(ar39_pileup_ids))
    cv_string += "--------------------------------------------\n"
    for i_mu in range(0, 3 + 1):
        mu_pileup_ids = np.nonzero(g_test == i_mu)[0]
        fn = len(np.intersect1d(ar39_pred_ids, mu_pileup_ids))
        cv_string += "Mu + {} Ar39\t{}\t{}\t| {}\n".format(i_mu, fn, len(mu_pileup_ids)-fn, len(mu_pileup_ids))
    cv_string += "--------------------------------------------\n"
    cv_string += "Total\t\t{}\t{}\t| {}\n".format(len(ar39_pred_ids), len(y_pred)-len(ar39_pred_ids), len(y_pred))
    return cv_string

def plot_roc_curve_and_get_auc(y_test, y_pred, decision_threshold):
    fpr, tpr, thresholds_roc = roc_curve(y_test, y_pred)
    precision, recall, thresholds_prec = precision_recall_curve(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred)
#    plt.subplot(1, 2, 1)
#    plt.plot(fpr, tpr)
#    threshold_id = np.argmin(abs(thresholds_roc-decision_threshold))
#    plt.scatter(fpr[threshold_id], tpr[threshold_id], color='r', label="Threshold: {:.2f}".format(thresholds_roc[threshold_id]))
#    plt.text(.5, .5, "Training Set Size: {},\nTest Set Size: {}\nAUC: {:.3f}".format(len(X_train), len(X_test), auc))
#    plt.xlabel("False Positive Rate (FPR)")
#    plt.ylabel("True Positive Rate (TPR)")
#    plt.title("ROC Curve")
#    plt.legend()

 #   plt.subplot(1, 2, 2)
    plt.plot(recall, precision)
    threshold_id = np.argmin(abs(thresholds_prec-decision_threshold))
    plt.scatter(recall[threshold_id], precision[threshold_id], color='r', label="Threshold: {:.2f}".format(thresholds_prec[threshold_id]))
    plt.text(.25, .5, "Training Set Size: {},\nTest Set Size: {}".format(len(X_train), len(X_test)))
    plt.xlabel("Efficiency (or Recall)")
    plt.ylabel("Purity (or Precision)")
    plt.title("Purity vs Efficiency Curve")
    plt.legend()
    plt.show()
    return auc

def get_compiled_model(model_identifier, optimizer, loss, epochs, batch):
    model_generator = get_twoconv_model_sigmoid01
    if model_identifier=="model_2conv_1maxpool_1dense":
        if "hinge" in loss:
            model_generator = get_twoconv_model_tanh_11
        decision_threshold = 0.0 if "hinge" in loss else 0.5
        model = model_generator()
        binary_accuracy = tf.keras.metrics.BinaryAccuracy(threshold=decision_threshold)
        model.compile(loss=loss, optimizer=optimizer, metrics=[binary_accuracy])
        model_name_template = "model_2conv_1maxpool_1dense_{}epochs_{}batch_{}_{}_{}decthreshold_{}"
        model_name = model_name_template.format(epochs, batch, optimizer, loss, decision_threshold,
                                                datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    else:
        raise ValueError("model {} not supported".format(model_identifier))
    print(model_name)
    return model, model_name, decision_threshold

# Constants
LOGDIR = "logs"
MODELDIR = "models"
DATADIR = "data"
OUTDIR = "out"

# Parameters
test_set_perc = 0.1
dev_set_perc = 0.1
# Load dataset and label
all_data = load_dataset()
X, g, y = all_data[:, 2:-2], all_data[:, -2], all_data[:, -1]    # exclude edep, pedetected
# Prepare data
X_train, X_test, id_train, id_test = train_test_split(X, np.arange(0, len(y)), test_size=test_set_perc, random_state=1234)
y_train, y_test = y[id_train], y[id_test]
g_train, g_test = g[id_train], g[id_test]
X_train, X_dev, y_train, y_dev = train_test_split(X_train, y_train, test_size=dev_set_perc, random_state=1)

# Configure model (parametrize for multi testing)
model_identifier = "model_2conv_1maxpool_1dense"
optimizer = "adam"
n_epochs = 25
batch_size = 32
losses = ["hinge", "squared_hinge", "binary_crossentropy"]
models, model_names, decision_thresholds = [], [], []
for loss in losses:
    model, model_name, decision_threshold = get_compiled_model(model_identifier, optimizer, loss, n_epochs, batch_size)
    models.append(model)
    model_names.append(model_name)
    decision_thresholds.append(decision_threshold)

# loop on models
metrics_on_testset = []
for model, model_name, decision_threshold in zip(models, model_names, decision_thresholds):
    # train
    print("[Info] Start training {}".format(model_name))
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=os.path.join(LOGDIR, model_name), histogram_freq=1)
    model.fit(X_train, y_train, epochs=n_epochs, batch_size=batch_size,
              verbose=2, callbacks=[tensorboard_callback], validation_data=(X_dev, y_dev))
    # test
    y_pred = model.predict(X_test)
    y_binary_pred = [1 if y >= decision_threshold else 0 for y in y_pred]  #0,1 predictions
    # evaluate on test set
    tn, fp, fn, tp = cv(y_test, y_binary_pred ).ravel()
    accuracy = (tp + tn) / (tp + fp + fn + tn)
    recall = tp / (tp + fn)       # true positive rate
    precision = tp / (tp + fp)    # positive predictive power
    f1score = f1_score(y_test, y_binary_pred)    # armonic mean precision, recall
    auc = plot_roc_curve_and_get_auc(y_test, y_pred, decision_threshold)
    metrics_on_testset.append({"cv": [tn, fp, fn, tp],
                               "accuracy": accuracy,
                               "precision": precision,
                               "recall": recall,
                               "f1score": f1score,
                               "auc": auc})

# print out results
for model_name, metrics in zip(model_names, metrics_on_testset):
    output = "[Info] Model: {}.\n".format(model_name)
    output += get_my_confusion_matrix_as_string(g_test, y_pred)
    output += "[Test Result] tn, fp, fn, tp: {}\n".format(metrics["cv"])
    output += "[Test Result] Accuracy: {:.3f}, Precision: {:.3f}, Recall: {:.3f}, " \
              "F1: {:.3f}, AUC: {:.3f}\n".format(metrics["accuracy"], metrics["precision"], metrics["recall"],
                                               metrics["f1score"], metrics["auc"])
    outfile = os.path.join(OUTDIR, "{}.out".format(model_name))
    with open(outfile, "w+") as file:
        file.write(output)
    print("[Info] Output report in {}.".format(outfile))
    print(output)
