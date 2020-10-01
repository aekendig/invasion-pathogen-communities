# invasion-pathogen-communities
This repository holds the data and code for the manuscript "Native perennial and non-native annual grasses shape pathogen community composition and disease severity in a California grassland" by Amy E. Kendig, Erin R. Spear, S. Caroline Daws, S. Luke Flory, and Erin A. Mordecai published in Journal of Ecology.

### Contents
- code: descriptions below 
- original-data: descriptions below  
- intermediate-data: datasets produced by processing original data  
- output: outputs of R scripts that are used as inputs to other R scripts or as figures or tables in manuscript
- metadata: metadata for each data file in original-data, see manuscript for details on data collection
- invasion-pathogen-communities.Rproj: RStudio project for running R scripts

|code                                     |desription |
|:----------------------------------------|:-------------------------------------------------------------------------------------------------|
|fungal_data_processing.R                 |R script to format foliar fungal pathogen isolate data from JRBP compilation |
|damage_data_processing.R                 |R script to format disease severity data from JRBP compilation |
|background_plant_data_procssing.R        |R script to format grass density data from the observational study and manipulated experiment |
|fungal_community_differences.R           |R script to assess differences between foliar fungal pathogen communities associated with native perennial and non-native annual grasses (Question 1a) |
|damage_differences.R                     |R script to assess differences in disease severity between native perennial and non-native annual grasses (Question 1b) |
|fungal_community_host_density.R          |R script to assess the effects of native perennial and non-native annual grass densities on foliar fungal pathogen communities (Question 2a) |
|damage_host_density.R                    |R script to assess the effects of native perennial and non-native annual grass densities on disease severity (Question 2b) |
|density_figures.R                        |R script to create figures for Questions 2a and 2b |
|manipulated_experiment_map.R             |R script to create Fig. S2 |
                                          

|data                                                   |description |
|:------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|2016ExperimentPlotList_26Apr16_PlotAssignments.csv     |treatments applied to plots in the manipulated experiment at JRBP |
|BgPlantSizeStatus_JepsonCalFloraJRBP_021120.csv        |species-level data for grasses at JRBP |
|competition_plot_seeds_composition_damage              |plant density data collected from the manipulated experiment at JRBP in 2016 |
|Data File 3 - Full Dataset.csv                         |foliar fungal isolates, sequences, OTUs, and estimated species identities collected from grasses at Jasper Ridge Biological Preserve (JRBP) in 2015-2017 |
|JRBP_PathogenDamage_21July2016_Full.csv                |disease severity data collected from grasses in the manipulated experiment and observational study at JRBP in 2016 |
|JRBP_plant_count_combined_2015-04-10.csv               |plant density data collected from the observational study at JRBP in 2015 |
|JRBPmidAprPathPrcntDamCulturesData_25Apr2015ERS.csv    |disease severity data collected from grasses in the observational study at JRBP in April 2015 |
|JRBPmidMarPathPrcntDamData_5Apr2015ERS.csv             |disease severity data collected from grasses in the observational study at JRBP in March 2015 |
|JRBPtransects11to14PathPrcntDam28Apr_7May2015ERS.csv   |disease severity data collected from grasses in transects at JRBP in April 2015 |
|transect_seeds_composition_damage.csv                  |plant density data collected from the observational study at JRBP in 2016 |


## Comments
Please contact me with any questions about raw data if you are interested in using it for your own project.


## Citation
Kendig, A. E., E. R. Spear, S. C. Daws, S. L. Flory, and E. A. Mordecai. (2020). Data from: Native perennial and non-native annual grasses shape pathogen community composition and disease severity in a California grassland. (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.4062434
