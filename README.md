# SeaLevel_RF_LSTM_EGU_OS_2025
Scripts to reproduce Heuzé et al. (2025) "Drivers of high frequency extreme sea level around Northern Europe – Synergies between recurrent neural networks and Random Forest", Ocean Science

1. Run prep_sealevel.m to generate levels_homogenised.mat
2. Run extract_fromERA5.m to generate the de-trended timeseries of all the tested predictors, at each location
   
Now it depends on whether you want to run the Random Forest (option a) of the LSTM (b)
3a. Run dataprep_RF.m to generate for each city a csv with the timeseries of all the predictors, on all the variables' common period, with all the delays tested, that has the right format for the RF
3b. Run dataprep_LSTM to generate for each city a csv with the timeseries of all the predictors min-max normalised, on all the variables' common period, that has the right format for the LSTM
Both scripts also feature a short second section to verify whether the predictor series are correlated to each other

Move to Jupyter

4a. Run RF_all.ipynb for the entire Random Forest part. By default, generates an ensemble of 100 models. Choose within the script whether to run with all predictors or to activate the feature importance
4b. Run LSTM_all.ipynb for the entire LSTM part. By default, generates an ensemble of 100 models. Third cell does the permutation feature importance.

Et voilà! The rest is just a bit of plotting that you can do on your own :)
