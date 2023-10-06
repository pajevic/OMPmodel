# OMPmodel

Python code implementing the OMP model described in Pajevic et al, "Oligodendrocyte-mediated Myelin Plasticity", eLife, 2023 [1]. The code in ompsimulator.py, omputils.py, generate_runs.py was originally used to simulate the data and obtain the results reported in the manuscript. These files are provided mainly for replicating the results in [1]. In the future use ompmodel.py for the latest implementation of the OMP model and that file will be available in February of 2024.

------------------------------------------------------------------------------------------------------------------------------

## Installation

The current code can be used simply as a script. Make sure you have numpy, scipy, and matplotlib installed, and that ompsimulator.py and omputils.py are placed in the same directory. Run test_runs.sh, which should produce the matching results to the ones already stored in results (with orig-results suffix).

## Generating the scripts and simulation data/results

Since a large number of runs have been saved we are only providing the scripts for running simulations, i.e., the tools to generate those scripts. To generate scripts for running all simulations use generate_runs.py and auxilary text file specifying the grid of values to be explored. In scripts/ directory we provide the template file (DefaultValues.txt) which can then be changed based on the Appendix Tables in the manuscript. We provide an example .txt file for generating runs in Appendix table 2 (left and right). For the table on the left only the runs for nol=2 are specified, since nepochs depend on it.

After running :> generate_runs.p AppendixTable2left_nol2    (or, generate_runs.py scripts/AppendixTable2left_nol2.txt ) the following files are created:
 AppendixTable2left_nol2-subminfo.npy  # saves a dictionary of parameter values used for each run
 AppendixTable2left_nol2-params.npy    # saves the original code and the parameter grid used for this simulation set
 AppendixTable2left_nol2-1.swarm       # individual commands to be executed, either in parallel (using swarm utility or sequentually

 For most tables you will get more than 1000 runs, e.g., for the one on the right in Table A2, a grid of 6912 different runs will be created and they are chunked into individual submission files with maximum of 1000 lines (due to a limit of the number of parallel jobs to be submitted on our NIH HPC Biowulf cluster). Running "generate_runs.py AppendixTable2right" will generate scripts divided into seven files:  AppendixTable2right-1.swarm ...  AppendixTable2right-7.swarm.

 For each of the tables in the Appendix modify the values of DefaultValues.txt for the listed parameters, to reproduce the results. Use the provided table that maps the variable names in the paper and in the code.

The synchronization profiles are saved in the numpy file results/"runname"-results.npy. When using saver=5 only the basic
information will be saved. Use np.load(``filename'').tolist() to load the results file as a dictionary. The key 'tstdarr' containes the
synchronization measure, $\sigma_\tau$, saved as a numpy array with shape (nreps, nepochs, nol, ngroups), where ngroups
indicates the number of different mixed signals. Instructions on using saver to generate model high temporal resolution histories of the OMP variables will be provided later (yielding enormous files), and will change in the new ompmodel.py distribution, which will use OMPmodel class.

## References:

1. Sinisa PajevicDietmar PlenzPeter J BasserR Douglas Fields (2023) Oligodendrocyte-mediated myelin plasticity and its role in neural synchronization eLife 12:e81982.  https://doi.org/10.7554/eLife.81982
