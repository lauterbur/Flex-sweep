# Flex-sweep
Versatile tool for detecting selective sweeps with a variety of ages, strengths, starting allele frequencies, and completeness.

Creates the following directory structure:
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

###### Configuration file for making simulation array ######
It is recommended to train Flex-sweep using simulations generated with a wide range of mutation rates, recombination rates,
sweep strengths, sweep ages, swept allele starting frequencies, and swept allele ending frequencies. These should be chosen
based on reasonable estimates that reflect your species.

Choose a demographic model that represents reasonable expectations for your population.

###### Sequence to SIMULATE, TRAIN, and CLASSIFY ######
# singularity TODO FlexSweep.py simulate outputDir NE numberSims numberChroms configFile nodeFile (--continue) (--num_task)
# singularity TODO FlexSweep.py fv outputDir nodeFile (--continue) (--num_task) (--rmap </path/to/rmap>) (--keep)
# singularity TODO FlexSweep.py train outputDir nodeFile (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir nodeFile (--threshold) (--num_task)

###### Sequence to make FVs using existing data, TRAIN, and CLASSIFY ######
# singularity TODO FlexSweep.py fv outputDir nodeFile --data_loc </path/to/training_data_dir> (--continue) (--num_task) (--rmap </path/to/rmap>) (--keep)
# singularity TODO FlexSweep.py train outputDir nodeFile --fv_loc </path/to/feature_vector_dir> (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir nodeFile --model_loc </path/to/model_dir> (--threshold) (--num_task)

###### Sequence to TRAIN using existing fvs and CLASSIFY ######
# singularity TODO FlexSweep.py train outputDir nodeFile --fv_loc </path/to/feature_vector_dir> (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir nodeFile --model_loc </path/to/model_dir> --norm_loc </path/to/normalization_data> (--threshold) (--num_task)

###### Sequence to CLASSIFY with existing model (must include the neutral bins for normalizing data to classify, has to be same normalization data as used for training)
# singularity TODO FlexSweep.py classify outputDir nodeFile --model_loc </path/to/model_dir> --norm_loc </path/to/normalization_data> (--threshold) (--num_task)
