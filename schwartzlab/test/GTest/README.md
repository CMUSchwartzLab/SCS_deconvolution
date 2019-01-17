This directory is storing the scripts for submitting the job to queue.

The scripts are used to solve the deconvolution problem using phylogeny-based method with the implementation of ILP, we used Gurobi to solve it.

There are two example script runDeTumorGurG07T3A00N00.sh and runDeTumorGurG33T3A00N00.sh to call DecomposeSolver.py. It contains some internal documentation indicating what arguments the user has to change for the user's installation of this entire repository.

The way we named the script is as follows:
- Gur = using gurobi
- G07 = work on GBM07  
- T3 = deconvolve 3 tumor samples
- A02 = regularization parameter $\alpha$ for the penalty equals to 0.2
- N00 = noise level is 0.0

The user might change the arguments inside the script and name it accordingly, the principle is that the user should know which script runs which job.
