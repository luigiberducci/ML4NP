![ML4NP Logo](https://github.com/luigiberducci/ML4NP/blob/master/ml4np_logo_muoni.png)
# ML4NP
ML4NP is a project founded by the GARR Grant 2020 for research and study activities. The project aims to develop a fast-decision trigger to tag cosmogenic events (muon passages) for the experiment LEGEND-200.

# Motivation
The experiment **LEGEND-200** aims to study a rare neutrino-physics process that occurs at very low evergy, using Germanium detectors in liquid Argon. The experiment results extremely sensitive to background events such as radioactive decays, cosmogenic background and eletromagnetic noise, that could change the signals in the low-energy spectrum in which we are interested in. For this reason, we are interested in the rejection of cosmogenic events but recognize these events only based on the amount of energy deposited is not trivial because of the radioactivity of natural liquid Argon. We propose the adoption the machine learning to learn a model based on spatio-temporal information about energy depositions.

# Proposed approach
In this project, we aim to investigate the application of machine learning for the automatic synthesis of trigger logic starting from simulated data. We want to derive a classification model that receive the ADC readouts of SiPMs as input and trigger when a cosmogenic event is detected. The resulting trigger system has to be implemented in a real-time environment, where the decision must be taken in the order of O(microsec). To this aim, we would like to consider models as simple as possible, with a reasonable tradeoff between performance and accuracy.

The project is held by **Luigi Berducci**, in collaboration with the *Sezione INFN di Roma 3*.
