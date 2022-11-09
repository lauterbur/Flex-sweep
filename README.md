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
### Configuration file for making simulation array ###
It is recommended to train Flex-sweep using simulations generated with a wide range of mutation rates, recombination rates,
sweep strengths, sweep ages, swept allele starting frequencies, and swept allele ending frequencies. These should be chosen
based on reasonable estimates that reflect your species.

Choose a demographic model that represents reasonable expectations for your population.

See wiki for full documentation
