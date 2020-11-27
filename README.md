# Deep-loop
Identification of chromatin loops using deep learning

File Description
1. data
The ‘data’ file contains all datasets used in this work, including training set and independent set in K562 and MCF-7 cell lines. In each sub file, ‘>+’ indicates that the following sequence is positive sample, whereas ‘>-’ indicates that the following sequence is negative sample.
2. featureCode
The ‘featureCode’ file contains all codes for feature extraction. For running these codes, python version 3 is required.
3. model
The ‘model’ file contains codes for training model (train_CNN_model.ipynb), testing model (test_model.ipynb). And all model files named with the .h5 suffix. ‘FF’, ‘FR’, ‘RF’, and ‘RR’ represent forward-forward orientation pairs, forward-reverse orientation pairs, reverse-forward orientation pairs, and reverse-reverse orientation pairs, respectively.
4. obtainNegativeSamples
The ‘obtainNegativeSamples’ file contains the code for getting negative samples.
5. dataOfIMR90
The ‘dataOfIMR90’ file contains data for evaluation the prediction ability of cell line/CBS pair model on IMR90 cell line.

