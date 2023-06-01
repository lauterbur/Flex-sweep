# Flex-sweep
Versatile tool for detecting selective sweeps with a variety of ages, strengths, starting allele frequencies, and completeness.

Creates the following directory structure:
~~~
outputDir
      |
      -- training_data
              |
              -- neutral (includes a file for each simulation)
                      |
                      -- stats
                              |
                              -- bins
                              -- ...
              -- sweep (includes a file for each simulation)
                      |
                      -- stats
                              |
                              -- ...
      -- fvs (includes neutral.fv and sweep.fv)
      -- model (includes model files, history, model test predictions)
              |
              -- predictions
      -- data_windows (includes a separate file and subdirectory for each window, and a fv for each window)
              |
              -- window_subdirectories
                      |
                      -- stats
                              |
                              -- ...
      -- classification (will include a predictions file)
~~~
### Singularity container to run Flex-sweep and two pre-trained models:
https://zenodo.org/record/7860595

### Configuration file for making simulation array ###
It is recommended to train Flex-sweep using simulations generated with a wide range of mutation rates, recombination rates,
sweep strengths, sweep ages, swept allele starting frequencies, and swept allele ending frequencies. These should be chosen
based on reasonable estimates that reflect your species.

Choose a demographic model that represents reasonable expectations for your population.

### Changelog 24 April 2023
- changed order of statistics in feature vector (was: [iHS, nsl, iSAFE, DIND, hapDAF-o, hapDAF-s, highfreq, lowfreq, Sratio, HAF, H12], now: [DIND, HAF, hapDAF-o, iSAFE, highfreq, hapDAF-s, nsl, Sratio, lowfreq, iHS, H12]) reflecting new analyses in revision

- upload new singularity (apptainer) image to zenodo with new statistic order

- upload new pre-trained models with new statistic order

- simulation config file can now take upper and lower bounds for normal distributions of parameters


### See wiki for full documentation
