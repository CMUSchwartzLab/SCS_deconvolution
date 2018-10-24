This directory is part of the repository
github.com/CMUSchwartzLab/SCS_deconvolution

that contains software to solve several formulations of the problem of deconvolving bulk tumor data into subpopulations,
as described in (Lei et al., in preparaion)

This directory contains the main prorgams to 
1) simulate data
2) solve the deconvolution problem in two different ways
3) generate the figures for the manuscript

Programs were written initially by Haoyun Lei and Bchuan Lyu. Programs were modified by Haoyun Lei and E. Michael Gertz. Programs were tested by
E. Michael Gertz, Haoyun Lei, Bchuan Lyu, and Alejandro Schaffer.

###################################################################################################################################
### Before you call any function, please move the test folder into the parent folder of LLSolver, 
### so that test/ and LLSolver/ will be parallel in the same parent folder (See ReadMe.txt in test folder)

DataSimulation.py 
  #This is to simulate bulk tumor data with desired number of samples from the single cell sequencing data
  #Simulation is based Geometric Mix Model
  
  ParentDirectory: specify a directory that contains the LLSolver folder
  DateFolder: you may do different simulation at different time, this is just for you to record different date
  TumorName: pick a tumor from which you choose the single cell data, now the available single cell data are from GBM07 or GBM33
  tumor_number: choose how many tumor samples you want to simulate, in the current experiment, we do 3, 6, 9 samples
  alpha: the alpha parameter for the dirichlet distribution, this is F matrix in the objective, please refer to the paper
  N: the number of simulated data you want to get
  Cap: True of False, if or not cap the CN larger than 10 to 10
  

The folder structure now will be as following:
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
