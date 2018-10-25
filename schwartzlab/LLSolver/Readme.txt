This directory is part of the repository
github.com/CMUSchwartzLab/SCS_deconvolution

that contains software to solve several formulations of the problem of deconvolving bulk tumor data into subpopulations,
as described in (Lei et al., in preparation)

The manuscript describes two methods called "phylogeny-free' and "phylogeny-based". Much of the documentation below assumes that the
user has at least skimmed the manuscript and is familiar with the terminology therein.

The present structure of the repository is that, except for README.md, all files are in the subdirectory
schwartzlab/LLSolver

The setup assumes that the user will create another subdirectory
schwartzlab/test

so that LLSolver and test are parallel subdirectories.
Various other subdirectories get created to store the simulated data for deconvolution tests

This LLsolver subdirectory contains the main prorgams to 
1) simulate data
2) solve the deconvolution problem in two different ways
3) generate the figures for the manuscript

Programs were written initially by Haoyun Lei and Bchuan Lyu. Programs were modified by Haoyun Lei and E. Michael Gertz. Programs were tested by
E. Michael Gertz, Haoyun Lei, Bchuan Lyu, and Alejandro Schaffer.

All programs are written in python3 (not python2). Some programs assume the avilability of the Gurobi package to solve optimization problems.
All programs assume the available of the numpy package.

###################################################################################################################################

DataSimulation.py 
  #This program simulates bulk tumor data with desired number of samples from the single cell sequencing data
  #Simulation is based on the Geometric Model described in the manuscript and uses a Dirichlet distribution

  The arguments to DataSimulation.py are as follows:

  ParentDirectory: specify a directory that contains the LLSolver folder
  DateFolder: allows for different simulations run on different dates to be stored n different subdirectories (a.k.a. folders)
  TumorName: pick a tumor from which you choose the single cell data; the available data are from two tumors named GBM07 or GBM33 (GBM is short for glioblastoma multiforme)
  tumor_number: choose how many tumor samples you want to simulate, in the main experiments in the manuscript, we chose 3, 6, or 9 samples
  alpha: the alpha parameter for the Dirichlet distribution, describe in the manuscript
  N: the number of simulated replicates (the manuscript uses N=40)
  Cap: True of False, whether or not to cap the copy numbers larger than 10 at 10 (this was set to True for all runs in the manuscript and must be set to True for the example data provided)
  

After running DataSimulation.py the folder/subdirectory structure should be abstractly as follows:
/some/path/to/the/ParentDirectory:
                                  /data/single cell sequencing data
                                  /LLSolver/DataSimulation.py
                                           /DecomposeSolver.py
                                           /....
                                           
                                  /simulation/DateFolder/GBM07/3/simulateData1
                                                                /simulateData2
                                                                ...
                                                                /simulateDataN
                                                                
                                  /test/GTest/TestCase1.sh
                                             /TestCase2.sh
                                             ...
                                             /TestCaseN.sh
                                       /NTest/TestCase1.sh
                                             /TestCase2.sh
                                             ...
                                             /TestCaseN.sh

Here, GBM07 is the tumor name selected 9and could be GBM33 instead). Here, 3 is number of samples.
GTest refers to the phylogeny-based method, which uses the Gurobi (hence the G) python library, while Ntest refers to the phylogeny-free method, which is also called non-negative
matrix factorization (NMF, hence the N).


  ImportSCData(ImportSCData(ParentDirectory, TumorName, IntCNV=True, Cap=False)
    #this is to import the data from single cell sequencing data
    IntCNV: if or not choose the integer Copy number
  
  SimulateData(N=40, Cap=Cap)
  #This is to simulate the desired number of simulated data and save them in an independent folder as described above

####################################################################################################################################
DecomposeSolver.py
  #This is to decompose the bulk tumor data to resolve the Copy number in loci in each fundamental cell type. 
  #The code will retrieve the simulation data if you set up the input arguments correctly (see below)
  #The code will solve all simulateData when you specify the number of tumor samples, e.g. if the number of tumor samples is 3, it will       solve all the simulateData saved in the folder and save the results for each of the simulateData
    
  ParentDirectory: specify a directory that contains the LLSolver folder
  DateFolder: you may do different simulation at different time, this is just for you to record different date, consistent with the one               you set up with DataSimulation.py
  TumorName: pick a tumor from which you choose the single cell data, now the available single cell data are from GBM07 or GBM33
  TumorNUmber: specify a tumor sample number so it can get data from that subfolder of the simulation folder
  reg1: regularization parameter of the penalty in nmf, only will be effecitve if the solver is nmf
  alpha: regularization parameter of the penalty in ILP, only will be effective if the solver is gurobi
  solver: choose nmf or gurobi to solve the problem
  
  The folder structure will be as following:
  /some/path/to/the/ParentDirectory:
                                  /data/single cell sequencing data
                                  /LLSolver/DataSimulation.py
                                           /DecomposeSolver.py
                                           /....
                                           
                                  /simulation/DateFolder/GBM07/3/simulateData1
                                                                /simulateData2
                                                                ...
                                                                /simulateDataN
                                                                
                                  /test/GTest/TestCase1.sh
                                             /TestCase2.sh
                                             ...
                                             /TestCaseN.sh
                                       /NTest/TestCase1.sh
                                             /TestCase2.sh
                                             ...
                                             /TestCaseN.sh
                                             
                                  /results/DataFolder/GBM07/3/nmf/result_for_simulateData1
                                                                 /result_for_simulateData2
                                                                 ...
                                                                 /result_for_simulateDataN
                                                              
                                                             /gurobi/result_for_simulateData1
                                                                    /result_for_simulateData2
                                                                    ...
                                                                    /result_for_simulateDataN
  
  extractValue(directory)
    #This is to retrieve the data in the replicate in the simulated data generated by DataSimulation.py
    directory: the path to each simulateData
    
  SolveDecomposition(AllDataPaths, solver):
    #This is to solve the problem with desired solver
    AllDataPaths: a list contains the path to each simulateData : [/path/to/simulateData1, /path/to/simulateData2, ...,      
                                                                    /path/to/simulateDataN ]
    solver: choose nmf or gurobi to solve the problem

###################################################################################################################################
  
  
Reference:
Haoyun Lei, Bochuan Lyu, E. Michael Gertz, Alejandro A. Schaffer, Xulian Shi, Kui Wu, Guibo Li, Liqin Xu, Yong Hou, Michael Dean, Russell Schwartz,
Tumor Copy Number Deconvolution Integrating Bulk and Single-Cell Sequencing Data, in preparation.
