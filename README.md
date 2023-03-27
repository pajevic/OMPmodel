# OMPmodel

Python code implementing the OMP model described in Pajevic et al, "Oligodendrocyte-mediated Myelin Plasticity", eLife, 2023. The full code and scripts will become available after the acceptance and publication of the manuscript.
------------------------------------------------------------------------------------------------------------------------------
The code in ompsimulator.py, omputils.py, generate_runs.py  was originally used to simulate the data and obgtain the results reported in the manuscript. These files are provided mainly for replicating the results in [1]. In the future use ompmodel.py for the latest implementation of the OMP model, that just provides an interface to running the model with user specified parameters. The ompmodel.py file will be available in May, 2023.

To generate scripts for running all simulations use generate_runs.py and auxilary text file specifying the grid of values to be explored.
An example for generating runs with the Appendix table 1 parameters  is provided in the directory scripts in file xxx.

References:

1. Pajevic S, Plenz D, Basser PJ, Fields, RD, "Oligodendrocyte-mediated Myelin Plasticity", eLife, xxx,  2023
