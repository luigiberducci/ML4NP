![ML4NP Logo](https://github.com/luigiberducci/ML4NP/blob/master/ml4np_logo.png)
# ML4NP
ML4NP is a project founded by the GARR Grant 2020 for research and study activities. The project aims to develop a fast-decision trigger to tag cosmogenic events (muon passages) for the experiment LEGEND-200.

# Motivation
The experiment LEGEND-200 aims to study a very rare physics process that occurs at very low evergy. For this reason, the experiment results extremely sensitive to background events such as radioactive decays and eletromagnetic noise. In particular, we are interested in the rejection of cosmogenic events that generates signals in the main detector and lead to write a large amount of meaningless data.

# Proposed approach
In this project, we aim to investigate the application of machine learning for the automatic synthesis of trigger logic. We want to derive a classification model starting from ADC readouts of SiPMs. The resulting model has to be implemented in a real-time environment, where the decision must be taken in the order of O(microsec). To this aim, we would like to consider models as simple as possible, with a reasonable tradeoff between performance and accuracy.

The project is held by Luigi Berducci, in collaboration with the Sezione INFN di Roma 3.
