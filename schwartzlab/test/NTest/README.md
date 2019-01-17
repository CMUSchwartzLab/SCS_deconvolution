This directory is storing the scripts for submitting the job to queue.

The scripts are used to solve the deconvolution problem using phylogeny-free method with the implementation of modified NMF.

There are two example scripts runDeTumorNMFG07T3A02.sh and runDeTumorNMFG33T3A02.sh to call DecomposeSolver.py. It contains some internal documentation indicating what arguments the user has to change for the user's installation of this entire repository.

The way we named the script is as follows:
- NMF = using modified Non-Negative Matrix Factorization algorithm
- G07 = working on GBM07
- T3 = deconvolve 3 tumor samples
- A02 = regularization parameter $\alpha$ for the penalty equals to 0.2
- N00 = Noise level is 0.0

The user might change the arguments inside the script and name it accordingly, the principle is that the user should know which script runs which job.
