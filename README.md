# Previ-sais-MF
script to analyse seasonal prediction for Triatlas (done for May startdate might need adaptation for november)

Map-comp-FOR-NDG_clean.ipynb: script to compare the coupled "home reanalysis" (dcppA-assim) with the forced reanalysis which is used to nudged the coupled simulation (JRA55doV1.5_Acycl: NEMO PISCES forced by the JRA55 reanalysis).
Corr-maps-verif_compare_s8_clean.ipynb: script to compare the hindcasts (May) with the system 8 for SST
Corr-maps-verif_obs_sat_clean.ipynb: script to compare (maps of climatology, biases and correlations) the hindcasts (May) at seasonal scale (JJA) with satelital obervations for chlorophyll, intpp or SST Might need some adaptation for november startdate.
corr-maps_Biochem-NDG_clean.ipynb: script to compare the hindcasts (May) at seasonal scale (6months) with the "home reanalysis" (dcppA-assim) for global map of correlations and biases. Might need some adaptation for november startdate.
Forecast_veri.py: fonctions for analysing seasonal forecasts (anomalies, climatologies, anomalies in cross_validations)
Function_read.py: fonctions to read files, extract region, compute area averaging
skill-reg-Biochem_clean.ipynb: script to compare the hindcasts (May) at seasonal scale (6months) with the "home reanalysis" (dcppA-assim) and satelital obervations for variables averaged in a given region (raw timeserie, anomalies, climatology and correlation). Might need some adaptation for november startdate.
Skill-reg_compare_s8_clean.ipynb: script to compare in several region the skill of the SST from MF operational system and the skill of the triatlas hindcast in a given region for the May startdate, Might need some adaptation for november startdate.
valid-Biochem-HIST-NDG-obs_clean.ipynb: script to compare the "home reanalysis" (dcppA-assim) with the obervation and historical simulation in order to validate (for physics) and estimate added value (for biogeochemistry) of nudging for initialisation for a given region.
